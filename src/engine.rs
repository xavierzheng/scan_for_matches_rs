//! Port of ggpunit.c: the pattern parser and the backtracking matcher.
//!
//! The code follows the C original statement by statement.
//!
//! The C program keeps its pattern in fixed-size static arrays (100 pattern
//! units, 600 code bytes, 10500 weights) and its state in globals.  When a
//! pattern is too big for these arrays the C code silently writes past them
//! into the neighbouring arrays, which changes its results.  To give the same
//! results in those cases too, all of that state lives in one memory block
//! laid out exactly like the static data of the reference build
//! (`gcc -std=gnu89 -O2`, Apple clang, macOS arm64).  Pointers are stored
//! in it as 8-byte addresses, so overlapping writes behave as they do in C.
//!
//! Memory used by the emulated program:
//!   * the static block (globals, pu_s, cv, iv; reference addresses
//!     0x1_0000_c000..0x1_0002_8000) followed by a read-only __LINKEDIT
//!     page.  It is mapped 16 KiB aligned, so pointer bits that C exposes
//!     stay the same (higher bits are random in C too).
//!   * the coded sequence buffer (C: `cdata`), with the slack before it and
//!     the low 24 address bits of the reference runs.
//!   * malloc(16) blocks for rule sets.
//!
//! Emulated addresses are host addresses, so a checked pointer is used
//! directly.  Every access is range checked unless the pointer was checked
//! before (see `Pol`); an access outside these blocks ends the program the
//! way the C program ends (SIGSEGV / SIGBUS).

use crate::sys::{segv, sigbus};

pub const PEPTIDE: i32 = 1;
pub const DNA: i32 = 2;

/// Longest sequence: lengths and positions are C `int`s.
pub const MAX_SEQ_LEN: usize = i32::MAX as usize - 2;
/// Address space reserved for the coded sequence buffer (memory is only
/// used for the pages that are touched).
const CD_CAP: i64 = 1 << 32;
/// Zero bytes reserved before the coded sequence buffer.
const LOW_SLACK: i64 = 1 << 20;
/// Names p0..p50 and rule sets r0..r50 (the C program accepts 50 but has
/// only 50 entries).
const MAX_NAMES: i32 = 50;
const N_NAMES: i64 = MAX_NAMES as i64 + 1;

// ---- layout of the static block -----------------------------------------
// The C program's globals and pattern arrays, with room for every pattern
// that fits on a pattern line (31 999 characters): one unit needs at least
// one character, a code byte or a weight needs at least one character.
const REF_S: i64 = 0x1_0000_c000;
const A_KNOWN_CHAR: i64 = REF_S;
const A_KNOWN_CHAR_INDEX: i64 = REF_S + 0x10;
const A_INITIALIZED: i64 = REF_S + 0x100;
const A_AD_PU_S: i64 = REF_S + 0x108;
const A_CODE_TO_PUNIT: i64 = REF_S + 0x200;
const A_NAMES: i64 = REF_S + 0x300; // 51 pointers
const A_PAST_LAST: i64 = A_NAMES + 8 * N_NAMES + 8;
const A_SEQ_TYPE: i64 = A_PAST_LAST + 8;
const A_P2C: i64 = REF_S + 0x600;
const A_RULE_SETS: i64 = REF_S + 0x700; // 51 pointers
const A_START_SRCH: i64 = A_RULE_SETS + 8 * N_NAMES + 8;
const A_PU_S: i64 = REF_S + 0x1000;
/// pattern unit slots (C: 100)
const N_SLOTS: i64 = 32_768;
const A_CV: i64 = A_PU_S + N_SLOTS * 264;
/// pattern code bytes (C: 600)
const CV_LEN: i64 = 1 << 17;
const A_IV: i64 = A_CV + CV_LEN;
/// weights (C: 10 500)
const IV_LEN: i64 = 1 << 16;
const REF_S_END: i64 = A_IV + 4 * IV_LEN;
const S_LEN: usize = (REF_S_END - REF_S) as usize;
/// read-only page after the block
const LINKEDIT_LEN: usize = 0x4000;
const PAGE: usize = 0x4000;

/// at most one 16-byte block per rule set r0..r50
const HEAP_CAP: usize = 16 * 64;

// ---- struct punit layout (sizeof = 264) -----------------------------------
const SZ: i64 = 264;
const O_TYPE: i64 = 0;
const O_NXT: i64 = 8;
const O_PREV: i64 = 16;
const O_BR: i64 = 24;
const O_ANCH: i64 = 32;
const O_HIT: i64 = 40;
const O_MLEN: i64 = 48;
// union
const O_U0: i64 = 56; // or1, exact.len, code_matrix, range.min, ins, llim_vec
const O_U4: i64 = 60; // range.width, del
const O_U8: i64 = 64; // or2, exact.code, range.nxt, mis, wvec.vec
const O_U12: i64 = 68; // of, sim.len
const O_U16: i64 = 72; // or.SR, rule_set, sim.code, wvec.cutoff
const O_U20: i64 = 76; // wvec.tupsz
const O_U24: i64 = 80; // or.alt, wvec.maxwt
const O_LLIM_BOUND: i64 = 260;

const EXACT_PUNIT: i32 = 0;
const RANGE_PUNIT: i32 = 1;
const COMPL_PUNIT: i32 = 2;
const REPEAT_PUNIT: i32 = 3;
const SIM_PUNIT: i32 = 4;
const WEIGHT_PUNIT: i32 = 5;
const OR_PUNIT: i32 = 6;
const ANY_PUNIT: i32 = 8;
const LLIM_PUNIT: i32 = 9;
const INV_REP_PUNIT: i32 = 10;
const MATCH_START: i32 = 11;
const MATCH_END: i32 = 12;

const A_BIT: u8 = 0x01;
const C_BIT: u8 = 0x02;
const G_BIT: u8 = 0x04;
const T_BIT: u8 = 0x08;

const KNOWN_CHAR: [u8; 16] = [0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];
const KNOWN_CHAR_INDEX: [i8; 16] = [-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1];

/// Field of a range unit (free in the C union for ranges): the unit after
/// it, when gap lengths may be skipped (see `next_start`); else 0.
const O_SKIP: i64 = 72;
/// k-mer length of the per-record position index
const KIDX_K: usize = 5;
/// below this many candidate positions a plain scan is used
const SCAN_LIMIT: i64 = 48;

/// Positions of every k-mer of plain bases (A, C, G, T) in the coded
/// sequence: `pos[starts[key]..starts[key + 1]]`, in increasing order.
pub struct KIndex {
    generation: u64,
    starts: Vec<u32>,
    pos: Vec<u32>,
}

/// The coded sequence of one strand of a long record, shared by all the
/// engines that search its pieces (instead of a copy in each engine).  It
/// is read only: the matcher only reads the sequence (a write would be a
/// SIGSEGV, not a change seen by the other threads).
pub struct SharedSeq {
    /// start of the mapping: LOW_SLACK zero bytes, then the sequence
    map: *mut u8,
}

// The mapping is read only after it is made, and unmapped on drop.
unsafe impl Send for SharedSeq {}
unsafe impl Sync for SharedSeq {}

impl SharedSeq {
    fn cd(&self) -> *mut u8 {
        unsafe { self.map.add(LOW_SLACK as usize) }
    }
}

impl Drop for SharedSeq {
    fn drop(&mut self) {
        crate::sys::unmap(self.map, (LOW_SLACK + CD_CAP) as usize);
    }
}

/// A k-mer index shared by the engines that search pieces of one record
/// strand (built by the first that needs it).
pub type SharedIndex = std::sync::Arc<std::sync::OnceLock<KIndex>>;

/// 2-bit value of a plain base code, or None.
#[inline(always)]
fn base2(c: u8) -> Option<u32> {
    match c & 15 {
        1 => Some(0),
        2 => Some(1),
        4 => Some(2),
        8 => Some(3),
        _ => None,
    }
}

/// As `build_index`, with `threads` threads: each counts, then fills, the
/// k-mers that start in its part of `seq`.  Parts are in sequence order,
/// so every bucket is sorted as with one thread (same index).
fn build_index_par(seq: &[u8], threads: usize) -> KIndex {
    const MIN_PART: usize = 1 << 20;
    let parts = threads.min(seq.len() / MIN_PART).max(1);
    if parts == 1 {
        return build_index(seq, None);
    }
    let nkeys = 1usize << (2 * KIDX_K);
    let mask = (nkeys - 1) as u32;
    let bounds: Vec<usize> = (0..=parts).map(|p| p * seq.len() / parts).collect();
    // call `f(start, key)` for each k-mer of plain bases that starts in
    // lo..hi (the k-1 bases before lo are read too)
    let each = |lo: usize, hi: usize, f: &mut dyn FnMut(usize, u32)| {
        let from = lo.saturating_sub(KIDX_K - 1);
        let to = (hi + KIDX_K - 1).min(seq.len());
        let (mut key, mut run) = (0u32, 0usize);
        for (i, &c) in seq[from..to].iter().enumerate() {
            match base2(c) {
                Some(b) => {
                    key = (key << 2 | b) & mask;
                    run += 1;
                }
                None => run = 0,
            }
            if run >= KIDX_K {
                let at = from + i + 1 - KIDX_K;
                if at >= lo {
                    f(at, key);
                }
            }
        }
    };
    // pass 1: counts of each part
    let counts: Vec<Vec<u32>> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..parts)
            .map(|p| {
                let (lo, hi) = (bounds[p], bounds[p + 1]);
                let each = &each;
                sc.spawn(move || {
                    let mut c = vec![0u32; nkeys];
                    each(lo, hi, &mut |_, k| c[k as usize] += 1);
                    c
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    // bucket starts, and where each part writes in each bucket
    let mut starts = vec![0u32; nkeys + 1];
    let mut slots = vec![vec![0u32; nkeys]; parts];
    let mut at = 0u32;
    for k in 0..nkeys {
        starts[k] = at;
        for p in 0..parts {
            slots[p][k] = at;
            at += counts[p][k];
        }
    }
    starts[nkeys] = at;
    let mut pos = vec![0u32; at as usize];
    // pass 2: fill (parts write to separate slots)
    struct Ptr(*mut u32);
    unsafe impl Send for Ptr {}
    unsafe impl Sync for Ptr {}
    let out = Ptr(pos.as_mut_ptr());
    std::thread::scope(|sc| {
        for (p, mut slot) in slots.into_iter().enumerate() {
            let (lo, hi) = (bounds[p], bounds[p + 1]);
            let each = &each;
            let out = &out;
            sc.spawn(move || {
                each(lo, hi, &mut |x, k| {
                    let s = &mut slot[k as usize];
                    unsafe { *out.0.add(*s as usize) = x as u32 };
                    *s += 1;
                });
            });
        }
    });
    KIndex {
        generation: 0,
        starts,
        pos,
    }
}

/// Positions of every k-mer of plain bases in the coded sequence `seq`
/// (the buffers of `reuse` are used again).
fn build_index(seq: &[u8], reuse: Option<KIndex>) -> KIndex {
    let nkeys = 1usize << (2 * KIDX_K);
    let mask = (nkeys - 1) as u32;
    let mut idx = reuse.unwrap_or(KIndex {
        generation: 0,
        starts: Vec::new(),
        pos: Vec::new(),
    });
    idx.starts.clear();
    idx.starts.resize(nkeys + 1, 0);
    // pass 1: count, pass 2: fill
    for pass in 0..2 {
        let mut key = 0u32;
        let mut run = 0usize;
        for (i, &c) in seq.iter().enumerate() {
            match base2(c) {
                Some(b) => {
                    key = (key << 2 | b) & mask;
                    run += 1;
                }
                None => run = 0,
            }
            if run >= KIDX_K {
                let at = i + 1 - KIDX_K;
                if pass == 0 {
                    idx.starts[key as usize + 1] += 1;
                } else {
                    let slot = &mut idx.starts[key as usize];
                    idx.pos[*slot as usize] = at as u32;
                    *slot += 1;
                }
            }
        }
        if pass == 0 {
            for k in 0..nkeys {
                idx.starts[k + 1] += idx.starts[k];
            }
            idx.pos.clear();
            idx.pos.resize(idx.starts[nkeys] as usize, 0);
        } else {
            // starts[k] now holds the end of bucket k: shift back
            for k in (0..nkeys).rev() {
                idx.starts[k + 1] = idx.starts[k];
            }
            idx.starts[0] = 0;
        }
    }
    idx
}

/// rev_compl_data(): one code of the reverse complement, or with a rule
/// set the index of the base in the rule table.
#[inline(always)]
fn compl_code(i: u32, rule_set: i32) -> u8 {
    if rule_set == -1 {
        ((i >> 4) & 15) as u8
    } else {
        ((KNOWN_CHAR_INDEX[(i & 15) as usize] as i32) << 2) as u8
    }
}

/// Size of loose_match's on-stack buffer for the complemented operand
/// (longer operands use the heap).
const MAX_CODES: i32 = 100;

/// The sequence buffer (C: `data`).  It is reused for every sequence;
/// bytes past the current sequence keep what an earlier, longer sequence
/// left there, as in C.
pub struct Buf {
    pub v: Vec<u8>,
}

impl Buf {
    pub fn new() -> Buf {
        Buf { v: Vec::new() }
    }

    /// Store `body` followed by a NUL byte.
    pub fn store(&mut self, body: &[u8]) {
        if self.v.len() < body.len() + 1 {
            self.v.resize(body.len() + 1, 0);
        }
        self.v[..body.len()].copy_from_slice(body);
        self.v[body.len()] = 0;
    }

    #[inline(always)]
    pub fn get(&self, off: i64) -> u8 {
        if off >= 0 && (off as usize) < self.v.len() {
            self.v[off as usize]
        } else {
            0
        }
    }
}

/// Access policy for the matcher.  `Fast` reads the fields of a pattern
/// unit without range checks; it is only used while the unit pointer has
/// been checked to lie inside the static block.  `Checked` treats every
/// access like any other C memory access.
trait Pol {
    const CHECKED: bool;
}
struct Fast;
struct Checked;
impl Pol for Fast {
    const CHECKED: bool = false;
}
impl Pol for Checked {
    const CHECKED: bool = true;
}

#[inline(always)]
fn hr32(p: *mut u8, off: i64) -> i32 {
    unsafe { (p.add(off as usize) as *const i32).read_unaligned() }
}

#[inline(always)]
fn hr64(p: *mut u8, off: i64) -> i64 {
    unsafe { (p.add(off as usize) as *const i64).read_unaligned() }
}

#[inline(always)]
fn hw32(p: *mut u8, off: i64, v: i32) {
    unsafe { (p.add(off as usize) as *mut i32).write_unaligned(v) }
}

#[inline(always)]
fn hw64(p: *mut u8, off: i64, v: i64) {
    unsafe { (p.add(off as usize) as *mut i64).write_unaligned(v) }
}

/// Does a whole `struct punit` at `p` lie inside the static block?
#[inline(always)]
fn pu_ok_at(sb: i64, p: i64) -> bool {
    (p.wrapping_sub(sb) as u64) <= (S_LEN as i64 - SZ) as u64
}

/// How the fast matcher stopped.
enum Out {
    Done(i32),
    /// a unit pointer outside the static block: continue in checked mode
    Bail {
        back: bool,
        cr: i64,
        sr: i64,
        br: i64,
    },
}

/// Everything the rare (out-of-line) memory paths need.  It sits in a Box
/// so the hot code only carries one pointer to it.
struct SlowMem {
    s: *mut u8,
    cd: *mut u8,
    heap: *mut u8,
    heap_len: usize,
}

#[inline(always)]
fn rb_full(sm: &SlowMem, a: i64) -> u8 {
    let o = a.wrapping_sub(sm.cd as i64);
    if (-LOW_SLACK..CD_CAP).contains(&o) {
        // the bytes before the buffer are zero in our allocation
        return unsafe { *(a as *const u8) };
    }
    let o = a.wrapping_sub(sm.s as i64) as u64;
    if o < (S_LEN + LINKEDIT_LEN) as u64 {
        // __LINKEDIT reads as zero here
        return unsafe { *(a as *const u8) };
    }
    let o = a.wrapping_sub(sm.heap as i64);
    if o >= 0 && (o as usize) < sm.heap_len {
        return unsafe { *(a as *const u8) };
    }
    segv()
}

#[inline(always)]
fn wb_full(sm: &SlowMem, a: i64, v: u8) {
    let o = a.wrapping_sub(sm.s as i64) as u64;
    if o < S_LEN as u64 {
        unsafe { *(a as *mut u8) = v };
        return;
    }
    if o < (S_LEN + LINKEDIT_LEN) as u64 {
        sigbus(); // read-only segment
    }
    let o = a.wrapping_sub(sm.heap as i64);
    if o >= 0 && (o as usize) < sm.heap_len {
        unsafe { *(a as *mut u8) = v };
        return;
    }
    let o = a.wrapping_sub(sm.cd as i64) as u64;
    if o < CD_CAP as u64 {
        unsafe { *(a as *mut u8) = v };
        return;
    }
    segv()
}

#[cold]
#[inline(never)]
fn rb_slow(sm: *const SlowMem, a: i64) -> u8 {
    rb_full(unsafe { &*sm }, a)
}

#[cold]
#[inline(never)]
fn wb_slow(sm: *const SlowMem, a: i64, v: u8) {
    wb_full(unsafe { &*sm }, a, v)
}

#[cold]
#[inline(never)]
fn r32_slow(sm: *const SlowMem, a: i64) -> i32 {
    let sm = unsafe { &*sm };
    i32::from_le_bytes([
        rb_full(sm, a),
        rb_full(sm, a + 1),
        rb_full(sm, a + 2),
        rb_full(sm, a + 3),
    ])
}

#[cold]
#[inline(never)]
fn r64_slow(sm: *const SlowMem, a: i64) -> i64 {
    let sm = unsafe { &*sm };
    let mut b = [0u8; 8];
    for (k, x) in b.iter_mut().enumerate() {
        *x = rb_full(sm, a + k as i64);
    }
    i64::from_le_bytes(b)
}

#[cold]
#[inline(never)]
fn wn_slow(sm: *const SlowMem, a: i64, v: u64, n: usize) {
    let sm = unsafe { &*sm };
    for (k, x) in v.to_le_bytes().iter().take(n).enumerate() {
        wb_full(sm, a + k as i64, *x);
    }
}

/// Raw view of the C program's memory: the static block and the coded
/// sequence buffer, plus a pointer to the rest.  It is small and `Copy`
/// so the hot matching code keeps it in registers.
#[derive(Clone, Copy)]
struct M {
    s: *mut u8,
    cd: *mut u8,
    sm: *const SlowMem,
    pep: bool,
}

impl M {
    /// Address of a C static object (given by its reference address) in
    /// this engine's static block.
    #[inline(always)]
    fn sa(self, a: i64) -> i64 {
        self.s as i64 + (a - REF_S)
    }

    #[inline(always)]
    fn rb(self, a: i64) -> u8 {
        let o = a.wrapping_sub(self.cd as i64) as u64;
        if o < CD_CAP as u64 {
            return unsafe { *(a as *const u8) };
        }
        let o = a.wrapping_sub(self.s as i64) as u64;
        if o < S_LEN as u64 {
            return unsafe { *(a as *const u8) };
        }
        rb_slow(self.sm, a)
    }

    #[inline(always)]
    fn wb(self, a: i64, v: u8) {
        let o = a.wrapping_sub(self.s as i64) as u64;
        if o < S_LEN as u64 {
            unsafe { *(a as *mut u8) = v };
            return;
        }
        wb_slow(self.sm, a, v)
    }

    #[inline(always)]
    fn r32(self, a: i64) -> i32 {
        let o = a.wrapping_sub(self.s as i64) as u64;
        if o + 4 <= S_LEN as u64 {
            return unsafe { (a as *const i32).read_unaligned() };
        }
        r32_slow(self.sm, a)
    }

    #[inline(always)]
    fn r64(self, a: i64) -> i64 {
        let o = a.wrapping_sub(self.s as i64) as u64;
        if o + 8 <= S_LEN as u64 {
            return unsafe { (a as *const i64).read_unaligned() };
        }
        r64_slow(self.sm, a)
    }

    #[inline(always)]
    fn w32(self, a: i64, v: i32) {
        let o = a.wrapping_sub(self.s as i64) as u64;
        if o + 4 <= S_LEN as u64 {
            unsafe { (a as *mut i32).write_unaligned(v) };
            return;
        }
        wn_slow(self.sm, a, v as u32 as u64, 4);
    }

    #[inline(always)]
    fn w64(self, a: i64, v: i64) {
        let o = a.wrapping_sub(self.s as i64) as u64;
        if o + 8 <= S_LEN as u64 {
            unsafe { (a as *mut i64).write_unaligned(v) };
            return;
        }
        wn_slow(self.sm, a, v as u64, 8);
    }

    #[inline(always)]
    fn pr32<P: Pol>(self, a: i64) -> i32 {
        if P::CHECKED {
            self.r32(a)
        } else {
            unsafe { (a as *const i32).read_unaligned() }
        }
    }

    #[inline(always)]
    fn pr64<P: Pol>(self, a: i64) -> i64 {
        if P::CHECKED {
            self.r64(a)
        } else {
            unsafe { (a as *const i64).read_unaligned() }
        }
    }

    /// Host pointer to the unit at `pu` (only valid when m.pu_ok(pu)).
    #[inline(always)]
    fn up(self, pu: i64) -> *mut u8 {
        pu as *mut u8
    }

    #[inline(always)]
    fn pu_ok(self, p: i64) -> bool {
        pu_ok_at(self.s as i64, p)
    }

    /// names[i] (static slot unless `i` is out of range)
    #[inline(always)]
    fn names_p<P: Pol>(self, i: i32) -> i64 {
        if !P::CHECKED && (0..=MAX_NAMES).contains(&i) {
            self.pr64::<Fast>(self.sa(A_NAMES) + 8 * i as i64)
        } else {
            self.names(i)
        }
    }

    /// a field of a unit reached through a pointer that has not been checked
    #[inline(always)]
    fn pu_r64<P: Pol>(self, pu: i64, off: i64) -> i64 {
        if !P::CHECKED && self.pu_ok(pu) {
            self.pr64::<Fast>(pu + off)
        } else {
            self.r64(pu + off)
        }
    }

    #[inline(always)]
    fn pu_r32<P: Pol>(self, pu: i64, off: i64) -> i32 {
        if !P::CHECKED && self.pu_ok(pu) {
            self.pr32::<Fast>(pu + off)
        } else {
            self.r32(pu + off)
        }
    }

    /// Are the `n` bytes from `a` inside the coded sequence buffer?
    #[inline(always)]
    fn cd_ok(self, a: i64, n: i64) -> bool {
        let o = a.wrapping_sub(self.cd as i64);
        o >= 0 && n >= 0 && o.wrapping_add(n) <= CD_CAP
    }

    #[inline(always)]
    fn names(self, i: i32) -> i64 {
        self.r64(self.sa(A_NAMES) + 8 * i as i64)
    }

    /// KnownChar(C)
    #[inline(always)]
    fn known_char(self, c: u8) -> bool {
        if self.pep {
            true
        } else {
            KNOWN_CHAR[c as usize] != 0
        }
    }

    /// Matches(C1, C2)
    #[inline(always)]
    fn matches(self, c1: u8, c2: u8) -> bool {
        if self.pep {
            c1 == c2 || c2 == b'X'
        } else {
            let a = c1 & 15;
            KNOWN_CHAR[a as usize] != 0 && (a & (c2 & 15)) == a
        }
    }

    /// ExMatches(RuleSet, C1, C2)
    #[inline(always)]
    fn ex_matches(self, rule_set: i32, c1: u8, c2: u8) -> bool {
        if rule_set == -1 {
            self.matches(c1, c2)
        } else {
            let base = self.r64(self.sa(A_RULE_SETS) + 8 * rule_set as i64);
            let idx = c2 as i64 + KNOWN_CHAR_INDEX[(c1 & 15) as usize] as i64;
            self.rb(base + idx) != 0
        }
    }

    #[inline(always)]
    fn one_get(self, one: &One, i: i64) -> u8 {
        match one {
            One::Mem(base) => self.rb(base + i),
            One::Bytes(b) => {
                if i >= 0 && (i as usize) < b.len() {
                    b[i as usize]
                } else {
                    0
                }
            }
        }
    }
}

#[inline(always)]
fn isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

#[inline(always)]
fn isalpha(c: u8) -> bool {
    c.is_ascii_alphabetic()
}

#[inline(always)]
fn at(l: &[u8], p: usize) -> u8 {
    if p < l.len() { l[p] } else { 0 }
}

fn ws(l: &[u8], mut p: usize) -> usize {
    while isspace(at(l, p)) {
        p += 1;
    }
    p
}

/// num(): returns the value and the new position.
fn num(l: &[u8], mut p: usize) -> Option<(i32, usize)> {
    let sign: i32 = if at(l, p) == b'-' {
        p += 1;
        -1
    } else {
        1
    };
    let c = at(l, p);
    if !c.is_ascii_digit() {
        return None;
    }
    let mut i: i32 = (c - b'0') as i32;
    p += 1;
    while at(l, p).is_ascii_digit() {
        i = i.wrapping_mul(10).wrapping_add((at(l, p) - b'0') as i32);
        p += 1;
    }
    Some((i.wrapping_mul(sign), p))
}

fn name_id(l: &[u8], p: usize) -> Option<(i32, usize)> {
    if at(l, p) != b'p' {
        return None;
    }
    match num(l, p + 1) {
        Some((n, p)) if (0..=MAX_NAMES).contains(&n) => Some((n, p)),
        _ => None,
    }
}

fn name_assgn(l: &[u8], p: usize) -> Option<(i32, usize)> {
    let (n, p) = name_id(l, p)?;
    let p = ws(l, p);
    if at(l, p) == b'=' {
        Some((n, p + 1))
    } else {
        None
    }
}

fn rule_id(l: &[u8], p: usize) -> Option<(i32, usize)> {
    if at(l, p) != b'r' {
        return None;
    }
    match num(l, p + 1) {
        Some((n, p)) if (0..=MAX_NAMES).contains(&n) => Some((n, p)),
        _ => None,
    }
}

/// returns (mis, ins, del, p)
fn misinsdel(l: &[u8], p: usize) -> Option<(i32, i32, i32, usize)> {
    if at(l, p) != b'[' {
        return None;
    }
    let (mis, p) = num(l, p + 1)?;
    if at(l, p) != b',' {
        return None;
    }
    let (ins, p) = num(l, p + 1)?;
    if at(l, p) != b',' {
        return None;
    }
    let (del, p) = num(l, p + 1)?;
    if at(l, p) != b']' {
        return None;
    }
    // negative counts are not allowed (C treats them as "no limit")
    if mis < 0 || ins < 0 || del < 0 {
        return None;
    }
    Some((mis, ins, del, p + 1))
}

fn elipses(l: &[u8], p: usize) -> Option<usize> {
    if at(l, p) == b'.' && at(l, p + 1) == b'.' && at(l, p + 2) == b'.' {
        Some(p + 3)
    } else {
        None
    }
}

fn word(s: &[u8], l: &[u8], mut p: usize) -> Option<usize> {
    for &c in s {
        if c != at(l, p) {
            return None;
        }
        p += 1;
    }
    Some(p)
}

fn range_pat(l: &[u8], p: usize) -> Option<(i32, i32, usize)> {
    let (min, p) = num(l, p)?;
    let p = ws(l, p);
    let p = elipses(l, p)?;
    let p = ws(l, p);
    let (max, p) = num(l, p)?;
    // negative lengths are not allowed (C reads outside the sequence)
    if min < 0 || max < 0 {
        return None;
    }
    Some((min, max, p))
}

fn dna_char(l: &[u8], p: usize) -> Option<(u8, usize)> {
    match at(l, p) {
        b'a' | b'A' => Some((0, p + 1)),
        b'c' | b'C' => Some((1, p + 1)),
        b'g' | b'G' => Some((2, p + 1)),
        b't' | b'T' | b'u' | b'U' => Some((3, p + 1)),
        _ => None,
    }
}

fn bond(l: &[u8], p: usize, rs: &mut [u8; 16]) -> Option<usize> {
    let (c1, p) = dna_char(l, p)?;
    let (c2, p) = dna_char(l, p)?;
    rs[((c1 << 2) + c2) as usize] = 1;
    Some(p)
}

fn parse_bonds(l: &[u8], p: usize, rs: &mut [u8; 16]) -> Option<usize> {
    let mut p = bond(l, p, rs)?;
    loop {
        if at(l, p) != b',' {
            return Some(p);
        }
        {
            let q = bond(l, p + 1, rs)?;
            p = q
        }
    }
}

/// C `compl()` from scan_for_matches.c (note: 's' maps to 'S').
#[inline(always)]
pub fn compl(c: u8) -> u8 {
    COMPL[c as usize]
}

static COMPL: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = compl_slow(i as u8);
        i += 1;
    }
    t
};

const fn compl_slow(c: u8) -> u8 {
    match c {
        b'a' => b't',
        b'A' => b'T',
        b'c' => b'g',
        b'C' => b'G',
        b'g' => b'c',
        b'G' => b'C',
        b't' | b'u' => b'a',
        b'T' | b'U' => b'A',
        b'm' => b'k',
        b'M' => b'K',
        b'r' => b'y',
        b'R' => b'Y',
        b'w' => b'w',
        b'W' => b'W',
        b's' => b's',
        b'S' => b'S',
        b'y' => b'r',
        b'Y' => b'R',
        b'k' => b'm',
        b'K' => b'M',
        b'b' => b'v',
        b'B' => b'V',
        b'd' => b'h',
        b'D' => b'H',
        b'h' => b'd',
        b'H' => b'D',
        b'v' => b'b',
        b'V' => b'B',
        b'n' => b'n',
        b'N' => b'N',
        _ => c,
    }
}

enum Lm {
    Top,
    Ins,
    Del,
}

#[derive(Clone, Copy, Default)]
struct StackEnt {
    p1: i64,
    p2: i64,
    n1: i32,
    n2: i32,
    mis: i32,
    ins: i32,
    del: i32,
    next_choice: i32,
}

/// Where loose_match's first operand lives.
enum One<'a> {
    /// C memory (sequence codes or pattern codes)
    Mem(i64),
    /// a private copy (complemented / reversed operand)
    Bytes(&'a [u8]),
}

pub struct Engine {
    /// malloc'd rule-set blocks (accessed through `mem`)
    _heap: Vec<u8>,
    slow: Box<SlowMem>,
    mem: M,
    /// punit_sequence_type (also stored in `s`)
    seq_type: i32,
    /// parse cursors (C: *pu_s, *cv, *iv)
    pup: i64,
    cvp: i64,
    ivp: i64,
    br1: i64,
    end_srch: i64,
    /// scratch for collect_hits (C: revhits[MAX_PUNITS] on the stack)
    revhits: Vec<i64>,
    /// counts changes of the coded sequence (for the k-mer index)
    cd_gen: u64,
    /// bytes of the coded sequence buffer used since `release_data`
    cd_hi: usize,
    /// the engine's own coded sequence buffer (`mem.cd` points to a
    /// `SharedSeq` while one is set)
    own_cd: *mut u8,
    /// k-mer index of the coded sequence, built when first needed
    kidx: Option<KIndex>,
    /// used instead of `kidx` while searching a piece of a long record
    shared_kidx: Option<SharedIndex>,
    /// threads that build a shared index
    index_threads: usize,
    /// highest start position of a hit (the first unit, when it is not
    /// anchored); i64::MAX except while searching a piece of a record
    start_lim: i64,
    /// scratch: masks of the unit after a range
    masks: Vec<u8>,
    /// `masks[k..k + 8]` packed in one word, as (k, word), when the masks
    /// checked are all plain bases; empty: check byte by byte.
    mask_words: Vec<(usize, u64)>,
    /// `masks[..mask_tol]` (a reverse complement with mismatches only)
    /// may have `mask_mis` mismatches; the other masks must match
    mask_tol: usize,
    mask_mis: i32,
    /// masks checked one by one in `masks_match_at` (the ones not covered
    /// by the index key, the tolerant part or `mask_words`)
    mask_bytes: [(usize, usize); 2],
    /// (range, units after its skip target that `next_start` checks too)
    skip_chains: Vec<(i64, Vec<i64>)>,
    /// punit_to_code for indexes 0..127 (never overwritten)
    p2c_lo: [u8; 128],
    /// for each unit slot: byte offset in the pattern line where the unit
    /// starts (output formats only; the parse result does not use it)
    unit_offs: Vec<u32>,
    /// record `hit_slots` for each hit (output formats only)
    track_units: bool,
    /// the unit slot of each entry of the last hit
    hit_slots: Vec<u32>,
}

// An engine owns its memory blocks (the raw pointers point only into
// them) and is used by one thread at a time.
unsafe impl Send for Engine {}

impl Engine {
    pub fn new() -> Engine {
        // static block followed by a read-only page
        let s_ptr = crate::sys::map_zeroed(0, S_LEN + LINKEDIT_LEN + 2 * PAGE);
        let s_ptr = {
            let a = s_ptr as usize;
            unsafe { s_ptr.add((PAGE - a % PAGE) % PAGE) }
        };
        let mut heap = vec![0u8; HEAP_CAP];
        // coded sequence buffer, with zero bytes before it
        let cd_ptr = unsafe {
            crate::sys::map_zeroed(0, (LOW_SLACK + CD_CAP) as usize).add(LOW_SLACK as usize)
        };
        let slow = Box::new(SlowMem {
            s: s_ptr,
            cd: cd_ptr,
            heap: heap.as_mut_ptr(),
            heap_len: 0,
        });
        let mem = M {
            s: slow.s,
            cd: slow.cd,
            sm: &*slow as *const SlowMem,
            pep: false,
        };
        let mut e = Engine {
            _heap: heap,
            slow,
            mem,
            seq_type: 0,
            pup: mem.sa(A_PU_S),
            cvp: mem.sa(A_CV),
            ivp: mem.sa(A_IV),
            br1: 0,
            end_srch: 0,
            revhits: Vec::with_capacity(256),
            cd_gen: 0,
            cd_hi: 0,
            own_cd: cd_ptr,
            kidx: None,
            shared_kidx: None,
            index_threads: 1,
            start_lim: i64::MAX,
            masks: Vec::new(),
            mask_words: Vec::new(),
            mask_tol: 0,
            mask_mis: 0,
            mask_bytes: [(0, 0); 2],
            skip_chains: Vec::new(),
            p2c_lo: [0; 128],
            unit_offs: Vec::new(),
            track_units: false,
            hit_slots: Vec::new(),
        };
        for (i, v) in KNOWN_CHAR.iter().enumerate() {
            e.wb(e.mem.sa(A_KNOWN_CHAR) + i as i64, *v);
        }
        for (i, v) in KNOWN_CHAR_INDEX.iter().enumerate() {
            e.wb(e.mem.sa(A_KNOWN_CHAR_INDEX) + i as i64, *v as u8);
        }
        e
    }

    /// Address of cdata[0] (C: the `cdata` pointer).
    pub fn cdata_base(&self) -> i64 {
        self.mem.cd as i64
    }

    // ------------------------------------------------------------------
    // memory

    #[inline(always)]
    fn rb(&self, a: i64) -> u8 {
        self.mem.rb(a)
    }

    #[inline(always)]
    fn wb(&mut self, a: i64, v: u8) {
        self.mem.wb(a, v)
    }

    #[inline(always)]
    fn r32(&self, a: i64) -> i32 {
        self.mem.r32(a)
    }

    #[inline(always)]
    fn r64(&self, a: i64) -> i64 {
        self.mem.r64(a)
    }

    #[inline(always)]
    fn w32(&mut self, a: i64, v: i32) {
        self.mem.w32(a, v)
    }

    #[inline(always)]
    fn w64(&mut self, a: i64, v: i64) {
        self.mem.w64(a, v)
    }

    fn malloc16(&mut self) -> i64 {
        let a = self.slow.heap as i64 + self.slow.heap_len as i64;
        self.slow.heap_len += 16;
        a
    }

    // globals that other C objects alias
    #[inline(always)]
    fn names(&self, i: i32) -> i64 {
        self.mem.names(i)
    }

    #[inline(always)]
    fn past_last(&self) -> i64 {
        self.r64(self.mem.sa(A_PAST_LAST))
    }

    #[inline(always)]
    fn start_srch(&self) -> i64 {
        self.r64(self.mem.sa(A_START_SRCH))
    }

    /// punit_to_code[c] where C subscripts with a signed char.
    #[inline(always)]
    fn p2c(&self, c: u8) -> u8 {
        // bytes >= 0x80 are not nucleotides (C indexes the table with a
        // signed char and reads other variables)
        if c < 0x80 { self.p2c_lo[c as usize] } else { 0 }
    }

    fn build_conversion_tables(&mut self) {
        for the_char in 0..256i64 {
            let lc = if the_char < 128 {
                (the_char as u8).to_ascii_lowercase()
            } else {
                the_char as u8
            };
            let mut v: u8 = match lc {
                b'a' => A_BIT,
                b'c' => C_BIT,
                b'g' => G_BIT,
                b't' => T_BIT,
                b'u' => T_BIT,
                b'm' => A_BIT | C_BIT,
                b'r' => A_BIT | G_BIT,
                b'w' => A_BIT | T_BIT,
                b's' => C_BIT | G_BIT,
                b'y' => C_BIT | T_BIT,
                b'k' => G_BIT | T_BIT,
                b'b' => C_BIT | G_BIT | T_BIT,
                b'd' => A_BIT | G_BIT | T_BIT,
                b'h' => A_BIT | C_BIT | T_BIT,
                b'v' => A_BIT | C_BIT | G_BIT,
                b'n' => A_BIT | C_BIT | G_BIT | T_BIT,
                _ => 0,
            };
            if v & A_BIT != 0 {
                v |= T_BIT << 4;
            }
            if v & C_BIT != 0 {
                v |= G_BIT << 4;
            }
            if v & G_BIT != 0 {
                v |= C_BIT << 4;
            }
            if v & T_BIT != 0 {
                v |= A_BIT << 4;
            }
            self.wb(self.mem.sa(A_P2C) + the_char, v);
            if the_char < 128 {
                self.p2c_lo[the_char as usize] = v;
            }
        }
        for the_char in 0..256i64 {
            let c = match the_char & 15 {
                0x1 => Some(b'A'),
                0x2 => Some(b'C'),
                0x4 => Some(b'G'),
                0x8 => Some(b'T'),
                0x3 => Some(b'M'),
                0x5 => Some(b'R'),
                0x9 => Some(b'W'),
                0x6 => Some(b'S'),
                0xA => Some(b'Y'),
                0xC => Some(b'K'),
                0xE => Some(b'B'),
                0xD => Some(b'D'),
                0xB => Some(b'H'),
                0x7 => Some(b'V'),
                0xF => Some(b'N'),
                _ => None,
            };
            if let Some(c) = c {
                self.wb(self.mem.sa(A_CODE_TO_PUNIT) + the_char, c);
            }
        }
        self.w32(self.mem.sa(A_INITIALIZED), 1);
    }

    /// comp_data(data, cdata): translate characters to codes, stop at NUL.
    pub fn comp_data(&mut self, data: &Buf) {
        self.set_shared_seq(None);
        self.cd_gen += 1;
        let n = data.v.iter().position(|&c| c == 0).expect("NUL at the end");
        self.cd_hi = self.cd_hi.max(n + 1);
        let cd = unsafe { std::slice::from_raw_parts_mut(self.mem.cd, n + 1) };
        let mut t = [0u8; 256];
        t[..128].copy_from_slice(&self.p2c_lo);
        for (d, &c) in cd[..n].iter_mut().zip(&data.v[..n]) {
            *d = t[c as usize];
        }
        cd[n] = 0;
    }

    /// strcpy(cdata, data)
    /// Pattern letters are upper case, so the data is matched in upper
    /// case too (the output shows the data as it is).
    pub fn copy_data(&mut self, data: &Buf) {
        self.set_shared_seq(None);
        self.cd_gen += 1;
        let mut k = 0usize;
        loop {
            let c = data.v[k].to_ascii_uppercase();
            unsafe { *self.mem.cd.add(k) = c };
            if c == 0 {
                break;
            }
            k += 1;
        }
        self.cd_hi = self.cd_hi.max(k + 1);
    }

    /// Code `data` (up to its NUL) as `comp_data` / `copy_data` do, into a
    /// new read-only buffer that engines can share (`set_shared_seq`).
    pub fn shared_seq(&self, data: &Buf) -> SharedSeq {
        let map = crate::sys::map_zeroed(0, (LOW_SLACK + CD_CAP) as usize);
        let s = SharedSeq { map };
        let n = data.v.iter().position(|&c| c == 0).expect("NUL at the end");
        let cd = unsafe { std::slice::from_raw_parts_mut(s.cd(), n + 1) };
        if self.seq_type == PEPTIDE {
            for (d, &c) in cd[..n].iter_mut().zip(&data.v[..n]) {
                *d = c.to_ascii_uppercase();
            }
        } else {
            let mut t = [0u8; 256];
            t[..128].copy_from_slice(&self.p2c_lo);
            for (d, &c) in cd[..n].iter_mut().zip(&data.v[..n]) {
                *d = t[c as usize];
            }
        }
        cd[n] = 0;
        crate::sys::protect_read_only(s.map, (LOW_SLACK + CD_CAP) as usize);
        s
    }

    /// Search `seq` (None: the engine's own buffer again) from now on.
    /// The caller keeps `seq` alive while it is set.
    pub fn set_shared_seq(&mut self, seq: Option<&SharedSeq>) {
        let cd = seq.map_or(self.own_cd, |s| s.cd());
        if cd != self.mem.cd {
            self.mem.cd = cd;
            self.slow.cd = cd;
            self.cd_gen += 1;
        }
    }

    /// Give the memory of the coded sequence and of the k-mer index back
    /// to the system (several patterns: an engine that is not used for a
    /// while must not hold a long record).  The next sequence is loaded
    /// again with `comp_data` / `copy_data`.
    pub fn release_data(&mut self) {
        self.set_shared_seq(None);
        self.kidx = None;
        self.shared_kidx = None;
        self.cd_gen += 1;
        crate::sys::release_pages(self.own_cd, self.cd_hi.next_multiple_of(PAGE));
        self.cd_hi = 0;
    }

    // ------------------------------------------------------------------
    // parser

    fn dna_pat(&mut self, l: &[u8], mut p: usize) -> Option<usize> {
        if self.seq_type != DNA {
            return None;
        }
        let start = self.cvp;
        let mut p1 = start;
        loop {
            let c = at(l, p);
            if c == 0 || c == b' ' || c == b'\t' || c == b'[' || c == b'\n' || c == b')' {
                break;
            }
            let code = self.p2c(c);
            self.wb(p1, code);
            if self.rb(p1) != 0 {
                p1 += 1;
            } else {
                return None;
            }
            p += 1;
        }
        if p1 > start {
            self.wb(p1, 0);
            self.cvp = p1 + 1;
            Some(p)
        } else {
            None
        }
    }

    fn char_pat(&mut self, l: &[u8], mut p: usize) -> Option<usize> {
        if self.seq_type == PEPTIDE && isalpha(at(l, p)) {
            let mut p1 = self.cvp;
            while isalpha(at(l, p)) {
                self.wb(p1, at(l, p).to_ascii_uppercase());
                p1 += 1;
                p += 1;
            }
            self.wb(p1, 0);
            self.cvp = p1 + 1;
            Some(p)
        } else {
            None
        }
    }

    /// returns (mis, ins, del, p)
    fn sim_pat(&mut self, l: &[u8], p: usize) -> Option<(i32, i32, i32, usize)> {
        let p1 = if self.seq_type == DNA {
            self.dna_pat(l, p)
        } else if self.seq_type == PEPTIDE {
            self.char_pat(l, p)
        } else {
            None
        }?;
        match misinsdel(l, p1) {
            Some(r) => Some(r),
            None => Some((0, 0, 0, p1)),
        }
    }

    /// returns (rs, n, p)
    fn compl_id(&self, l: &[u8], p: usize) -> Option<(i32, i32, usize)> {
        if self.seq_type != DNA {
            return None;
        }
        if at(l, p) == b'~' {
            let (n, p) = name_id(l, p + 1)?;
            Some((-1, n, p))
        } else {
            let (rs, p) = rule_id(l, p)?;
            if at(l, p) == b'~' {
                let (n, p) = name_id(l, p + 1)?;
                Some((rs, n, p))
            } else {
                None
            }
        }
    }

    /// returns (rs, n, mis, ins, del, p)
    fn compl_pat(&self, l: &[u8], p: usize) -> Option<(i32, i32, i32, i32, i32, usize)> {
        if self.seq_type != DNA {
            return None;
        }
        let (rs, n, p) = self.compl_id(l, p)?;
        match misinsdel(l, p) {
            Some((mis, ins, del, p1)) => Some((rs, n, mis, ins, del, p1)),
            None => Some((rs, n, 0, 0, 0, p)),
        }
    }

    /// returns (n, mis, ins, del, p)
    fn repeat_pat(&self, l: &[u8], p: usize) -> Option<(i32, i32, i32, i32, usize)> {
        let (n, p) = name_id(l, p)?;
        match misinsdel(l, p) {
            Some((mis, ins, del, p1)) => Some((n, mis, ins, del, p1)),
            None => Some((n, 0, 0, 0, p)),
        }
    }

    fn inv_rep_pat(&self, l: &[u8], p: usize) -> Option<(i32, i32, i32, i32, usize)> {
        if at(l, p) != b'<' {
            return None;
        }
        let (n, p) = name_id(l, p + 1)?;
        if self.seq_type == DNA
            && let Some((mis, ins, del, p1)) = misinsdel(l, p)
        {
            return Some((n, mis, ins, del, p1));
        }
        Some((n, 0, 0, 0, p))
    }

    /// n_tuple(n, t, p): stores the numbers at the int pointer t.
    fn n_tuple(&mut self, l: &[u8], p: usize, t: i64) -> Option<(i32, usize)> {
        if at(l, p) != b'(' {
            return None;
        }
        let p = ws(l, p + 1);
        let (v, mut p) = num(l, p)?;
        self.w32(t, v);
        let mut n: i32 = 1;
        loop {
            if at(l, p) == b',' {
                p += 1;
                p = ws(l, p);
                match num(l, p) {
                    Some((v, q)) => {
                        self.w32(t + 4 * n as i64, v);
                        n += 1;
                        p = ws(l, q);
                    }
                    None => return None,
                }
            } else {
                let c = at(l, p);
                p += 1;
                return if c == b')' { Some((n, p)) } else { None };
            }
        }
    }

    /// wt_template(): returns (tupsz, n, p).  The iv cursor moves even if
    /// the template fails, as in C.
    fn wt_template(&mut self, l: &[u8], p: usize) -> Option<(i32, i32, usize)> {
        let mut n: i32 = 0;
        if at(l, p) != b'{' {
            return None;
        }
        let p = ws(l, p + 1);
        let (tupsz, mut p) = self.n_tuple(l, p, self.ivp)?;
        self.ivp += 4 * tupsz as i64;
        n += 1;
        loop {
            p = ws(l, p);
            if at(l, p) != b',' {
                break;
            }
            p += 1;
            p = ws(l, p);
            match self.n_tuple(l, p, self.ivp) {
                Some((sz, q)) if sz == tupsz => {
                    n += 1;
                    self.ivp += 4 * sz as i64;
                    p = q;
                }
                _ => return None,
            }
        }
        if at(l, p) == b'}' && (tupsz == 4 || tupsz == 20 || tupsz == 21) {
            Some((tupsz, n, p + 1))
        } else {
            None
        }
    }

    /// returns (maxwt, cutoff, tupsz, n, p)
    fn wt_pat(&mut self, l: &[u8], p: usize) -> Option<(i32, i32, i32, i32, usize)> {
        if let Some((maxwt, p1)) = num(l, p) {
            let p1 = ws(l, p1);
            if at(l, p1) == b'>' {
                let p1 = ws(l, p1 + 1);
                if let Some((tupsz, n, p1)) = self.wt_template(l, p1) {
                    let p1 = ws(l, p1);
                    if at(l, p1) == b'>' {
                        let p1 = ws(l, p1 + 1);
                        if let Some((cutoff, p1)) = num(l, p1) {
                            return Some((maxwt, cutoff, tupsz, n, p1));
                        }
                    }
                }
            }
        }
        let maxwt = 1_000_000;
        let (tupsz, n, p) = self.wt_template(l, p)?;
        let p = ws(l, p);
        if at(l, p) != b'>' {
            return None;
        }
        let p = ws(l, p + 1);
        let (cutoff, p) = num(l, p)?;
        Some((maxwt, cutoff, tupsz, n, p))
    }

    /// name_id() storing through a pointer: num() writes *n before the range
    /// check, so the store happens even when the name is out of range.
    fn name_id_store(&mut self, l: &[u8], p: usize, dst: i64) -> Option<usize> {
        if at(l, p) != b'p' {
            return None;
        }
        let (n, q) = num(l, p + 1)?;
        self.w32(dst, n);
        if (0..=MAX_NAMES).contains(&n) {
            Some(q)
        } else {
            None
        }
    }

    /// llim_pat(v, max, p) with v = &pu->info.llim.llim_vec,
    /// max = &pu->info.llim.bound
    fn llim_pat(&mut self, l: &[u8], p: usize, pu: i64) -> Option<usize> {
        let v = pu + O_U0;
        self.w32(v, 1);
        let p = word(b"length(", l, p)?;
        let mut p = self.name_id_store(l, p, v + 4)?;
        loop {
            let p1 = ws(l, p);
            if at(l, p1) != b'+' {
                break;
            }
            let p1 = ws(l, p1 + 1);
            let cnt = self.r32(v);
            self.w32(v, cnt.wrapping_add(1));
            let dst = v + 4 + 4 * cnt as i64;
            match self.name_id_store(l, p1, dst) {
                Some(q) => p = q,
                None => break,
            }
        }
        let p = ws(l, p);
        if at(l, p) != b')' {
            return None;
        }
        let p = ws(l, p + 1);
        if at(l, p) != b'<' {
            return None;
        }
        let p = ws(l, p + 1);
        let (max, p) = num(l, p)?;
        self.w32(pu + O_LLIM_BOUND, max);
        Some(p)
    }

    /// any_pat(): returns (code matrix, p)
    fn any_pat(&self, l: &[u8], mut p: usize) -> Option<(i64, usize)> {
        if self.seq_type == PEPTIDE && word(b"any(", l, p).is_some() {
            p += 4;
            let mut cm: i64 = 0;
            while isalpha(at(l, p)) {
                cm |= (1i32 << (at(l, p).to_ascii_uppercase() - b'A')) as i64;
                p += 1;
            }
            let c = at(l, p);
            p += 1;
            if c == b')' { Some((cm, p)) } else { None }
        } else if self.seq_type == PEPTIDE && word(b"notany(", l, p).is_some() {
            p += 7;
            let mut cm: i64 = 0x3ffffff;
            while isalpha(at(l, p)) {
                cm &= !((1i32 << (at(l, p).to_ascii_uppercase() - b'A')) as i64);
                p += 1;
            }
            let c = at(l, p);
            p += 1;
            if c == b')' { Some((cm, p)) } else { None }
        } else {
            None
        }
    }

    fn or_pat(&mut self, l: &[u8], p: usize) -> Option<usize> {
        let or_pu = self.pup;
        self.pup += SZ;
        let r = (|| {
            if at(l, p) != b'(' {
                return None;
            }
            let p = ws(l, p + 1);
            let p1 = self.pup;
            let p = self.punit_cons_list(l, p)?;
            let p = ws(l, p);
            if at(l, p) != b'|' {
                return None;
            }
            let p = ws(l, p + 1);
            let p2 = self.pup;
            let p = self.punit_cons_list(l, p)?;
            let p = ws(l, p);
            if at(l, p) != b')' {
                return None;
            }
            Some((p1, p2, p + 1))
        })();
        match r {
            Some((p1, p2, p)) => {
                self.w32(or_pu + O_TYPE, OR_PUNIT);
                self.w64(or_pu + O_U0, p1);
                self.w64(or_pu + O_U8, p2);
                self.w64(p2 + O_PREV, or_pu);
                self.w64(p1 + O_PREV, or_pu);
                self.w64(or_pu + O_PREV, 0);
                self.w64(or_pu + O_NXT, 0);
                Some(p)
            }
            None => {
                self.pup = or_pu;
                None
            }
        }
    }

    fn bond_set(&mut self, l: &[u8], p: usize, n: i32) -> Option<usize> {
        let mut rs = [0u8; 16];
        if at(l, p) != b'{' {
            return None;
        }
        let p = parse_bonds(l, p + 1, &mut rs)?;
        if at(l, p) != b'}' {
            return None;
        }
        let slot = self.mem.sa(A_RULE_SETS) + 8 * n as i64;
        if self.r64(slot) == 0 {
            let a = self.malloc16();
            self.w64(slot, a);
        }
        let base = self.r64(slot);
        for (i, v) in rs.iter().enumerate() {
            self.wb(base + i as i64, *v);
        }
        Some(p + 1)
    }

    fn parse_rule_set(&mut self, l: &[u8], p: usize) -> Option<usize> {
        if self.seq_type != DNA {
            return None;
        }
        let (n, p) = rule_id(l, p)?;
        let p = ws(l, p);
        if at(l, p) != b'=' {
            return None;
        }
        self.bond_set(l, p + 1, n)
    }

    fn punit_parse(&mut self, l: &[u8], mut p: usize) -> Option<usize> {
        while let Some(p1) = self.parse_rule_set(l, p) {
            p = ws(l, p1);
        }
        let pu1 = self.pup;
        self.pup += SZ;
        let slot = ((pu1 - self.mem.sa(A_PU_S)) / SZ) as usize;
        if self.unit_offs.len() <= slot {
            self.unit_offs.resize(slot + 1, u32::MAX);
        }
        self.unit_offs[slot] = p as u32;
        if let Some((i, p1)) = name_assgn(l, p) {
            let slot = self.mem.sa(A_NAMES) + 8 * i as i64;
            if self.r64(slot) != 0 {
                return None; // the slot is not released, as in C
            }
            self.w64(slot, pu1);
            p = p1;
        }
        p = ws(l, p);
        let p1: usize;
        if at(l, p) == b'^' {
            self.w32(pu1 + O_TYPE, MATCH_START);
            p1 = p + 1;
        } else if at(l, p) == b'$' {
            self.w32(pu1 + O_TYPE, MATCH_END);
            p1 = p + 1;
        } else if let Some((i, j, q)) = range_pat(l, p) {
            self.w32(pu1 + O_TYPE, RANGE_PUNIT);
            self.w32(pu1 + O_U0, i);
            self.w32(pu1 + O_U4, j.wrapping_sub(i));
            p1 = q;
        } else if let Some((cm, q)) = self.any_pat(l, p) {
            self.w64(pu1 + O_U0, cm);
            self.w32(pu1 + O_TYPE, ANY_PUNIT);
            p1 = q;
        } else if let Some(q) = self.llim_pat(l, p, pu1) {
            self.w32(pu1 + O_TYPE, LLIM_PUNIT);
            p1 = q;
        } else if let Some((n, i, j, k, q)) = self.inv_rep_pat(l, p) {
            self.w32(pu1 + O_TYPE, INV_REP_PUNIT);
            self.w32(pu1 + O_U12, n);
            self.w32(pu1 + O_U8, i);
            self.w32(pu1 + O_U0, j);
            self.w32(pu1 + O_U4, k);
            p1 = q;
        } else if let Some((n, i, j, k, q)) = self.repeat_pat(l, p) {
            self.w32(pu1 + O_TYPE, REPEAT_PUNIT);
            self.w32(pu1 + O_U12, n);
            self.w32(pu1 + O_U8, i);
            self.w32(pu1 + O_U0, j);
            self.w32(pu1 + O_U4, k);
            p1 = q;
        } else if let Some((p2, (i, j, k, q))) = {
            let p2 = self.cvp;
            self.sim_pat(l, p).map(|r| (p2, r))
        } {
            if i == 0 && j == 0 && k == 0 {
                self.w32(pu1 + O_TYPE, EXACT_PUNIT);
                self.w64(pu1 + O_U8, p2);
                let len = self.strlen(p2);
                self.w32(pu1 + O_U0, len);
            } else {
                self.w32(pu1 + O_TYPE, SIM_PUNIT);
                self.w64(pu1 + O_U16, p2);
                let len = self.strlen(p2);
                self.w32(pu1 + O_U12, len);
                self.w32(pu1 + O_U8, i);
                self.w32(pu1 + O_U0, j);
                self.w32(pu1 + O_U4, k);
            }
            p1 = q;
        } else if let Some((rs, n, i, j, k, q)) = self.compl_pat(l, p) {
            self.w32(pu1 + O_TYPE, COMPL_PUNIT);
            self.w32(pu1 + O_U16, rs);
            self.w32(pu1 + O_U12, n);
            self.w32(pu1 + O_U8, i);
            self.w32(pu1 + O_U0, j);
            self.w32(pu1 + O_U4, k);
            p1 = q;
        } else if let Some((i1, (k, i, tupsz, j, q))) = {
            let i1 = self.ivp;
            self.wt_pat(l, p).map(|r| (i1, r))
        } {
            self.w32(pu1 + O_TYPE, WEIGHT_PUNIT);
            self.w64(pu1 + O_U8, i1);
            self.w32(pu1 + O_U0, j);
            self.w32(pu1 + O_U16, i);
            self.w32(pu1 + O_U24, k);
            self.w32(pu1 + O_U20, tupsz);
            p1 = q;
        } else {
            self.pup = pu1;
            return None;
        }
        Some(p1)
    }

    fn strlen(&self, mut a: i64) -> i32 {
        let mut n: i32 = 0;
        while self.rb(a) != 0 {
            n = n.wrapping_add(1);
            a += 1;
        }
        n
    }

    fn punit_list(&mut self, l: &[u8], p: usize) -> Option<usize> {
        let mut last: i64 = 0;
        let mut p = ws(l, p);
        let mut next = self.pup;
        while let Some(p1) = self.punit_parse(l, p) {
            self.w64(next + O_PREV, last);
            if last != 0 {
                self.w64(last + O_NXT, next);
            }
            last = next;
            next = self.pup;
            p = ws(l, p1);
        }
        if last == 0 {
            return None;
        }
        self.w64(last + O_NXT, 0);
        Some(p)
    }

    fn punit_cons(&mut self, l: &[u8], mut p: usize) -> Option<usize> {
        while let Some(p1) = self.parse_rule_set(l, p) {
            p = ws(l, p1);
        }
        if let Some(p1) = self.or_pat(l, p) {
            return Some(p1);
        }
        self.punit_list(l, p)
    }

    fn punit_cons_list(&mut self, l: &[u8], p: usize) -> Option<usize> {
        let mut last: i64 = 0;
        let mut p = ws(l, p);
        let mut next = self.pup;
        while let Some(p1) = self.punit_cons(l, p) {
            self.w64(next + O_PREV, last);
            if last != 0 {
                self.w64(last + O_NXT, next);
            }
            loop {
                let n = self.r64(next + O_NXT);
                if n == 0 {
                    break;
                }
                next = n;
            }
            last = next;
            next = self.pup;
            p = ws(l, p1);
        }
        if last == 0 {
            return None;
        }
        self.w64(last + O_NXT, 0);
        Some(p)
    }

    fn parser(&mut self, l: &[u8]) -> Option<i64> {
        let first = self.pup;
        self.w64(first + O_PREV, 0);
        let p = self.punit_cons_list(l, ws(l, 0))?;
        let p = ws(l, p);
        if at(l, p) == 0 { Some(first) } else { None }
    }

    fn set_anchors_on(&mut self, mut pu: i64, depth: u32) {
        if depth > 1_000_000 {
            segv();
        }
        loop {
            self.w32(pu + O_ANCH, 1);
            if self.r32(pu + O_TYPE) == OR_PUNIT {
                let (a, b) = (self.r64(pu + O_U0), self.r64(pu + O_U8));
                self.set_anchors_on(a, depth + 1);
                self.set_anchors_on(b, depth + 1);
            }
            pu = self.r64(pu + O_NXT);
            if pu == 0 {
                return;
            }
        }
    }

    fn set_anchors(&mut self, pu: i64, depth: u32) {
        if depth > 1_000_000 {
            segv();
        }
        self.w32(pu + O_ANCH, 0);
        if self.r32(pu + O_TYPE) == OR_PUNIT {
            let (a, b) = (self.r64(pu + O_U0), self.r64(pu + O_U8));
            self.set_anchors(a, depth + 1);
            self.set_anchors(b, depth + 1);
        }
        let n = self.r64(pu + O_NXT);
        if n != 0 {
            self.set_anchors_on(n, depth + 1);
        }
    }

    fn max_mat(&self, pu: i64, depth: u32) -> i32 {
        if depth > 1_000_000 {
            segv(); // unbounded recursion overflows the C stack
        }
        match self.r32(pu + O_TYPE) {
            EXACT_PUNIT => self.r32(pu + O_U0),
            LLIM_PUNIT => 0,
            RANGE_PUNIT => self.r32(pu + O_U0).wrapping_add(self.r32(pu + O_U4)),
            ANY_PUNIT => 1,
            COMPL_PUNIT | REPEAT_PUNIT | INV_REP_PUNIT => {
                let p1 = self.names(self.r32(pu + O_U12));
                self.max_mat(p1, depth + 1)
                    .wrapping_add(self.r32(pu + O_U0))
            }
            SIM_PUNIT => self.r32(pu + O_U12).wrapping_add(self.r32(pu + O_U0)),
            WEIGHT_PUNIT => self.r32(pu + O_U0),
            OR_PUNIT => {
                let a = self.max_mat(self.r64(pu + O_U0), depth + 1);
                let b = self.max_mat(self.r64(pu + O_U8), depth + 1);
                if a > b { a } else { b }
            }
            _ => 0,
        }
    }

    /// Every name and rule set that the pattern uses is defined, and no
    /// name refers to itself through other names (`p1=~p1`).  The C program
    /// crashes on these patterns.
    fn refs_ok(&self, mut pu: i64) -> bool {
        while pu != 0 {
            let ok = match self.r32(pu + O_TYPE) {
                COMPL_PUNIT => {
                    let rs = self.r32(pu + O_U16);
                    self.name_ok(self.r32(pu + O_U12))
                        && (rs == -1 || self.r64(self.mem.sa(A_RULE_SETS) + 8 * rs as i64) != 0)
                }
                REPEAT_PUNIT | INV_REP_PUNIT => self.name_ok(self.r32(pu + O_U12)),
                LLIM_PUNIT => {
                    let v = pu + O_U0;
                    (1..=self.r32(v)).all(|i| self.names(self.r32(v + 4 * i as i64)) != 0)
                }
                OR_PUNIT => self.refs_ok(self.r64(pu + O_U0)) && self.refs_ok(self.r64(pu + O_U8)),
                _ => true,
            };
            if !ok {
                return false;
            }
            pu = self.r64(pu + O_NXT);
        }
        true
    }

    /// Does the parsed pattern read a name that is not certainly matched
    /// before, in the same attempt?  (The name is defined later in the
    /// pattern, or only in one branch of an alternative.)  Such a unit uses
    /// the last match of the name, which can come from an earlier attempt
    /// or an earlier sequence, so the sequences are not independent.
    pub fn uses_earlier_state(&self) -> bool {
        let root = self.r64(self.mem.sa(A_AD_PU_S));
        if root == 0 {
            return false;
        }
        let mut stateful = false;
        self.walk_assigned(root, 0, &mut stateful);
        stateful
    }

    /// Names (bit n) that are certainly matched after the units from `pu`
    /// to the end of its list, given the names matched before.
    fn walk_assigned(&self, mut pu: i64, mut assigned: u64, stateful: &mut bool) -> u64 {
        let need = |n: i32, assigned: u64, stateful: &mut bool| {
            if !(0..=MAX_NAMES).contains(&n) || assigned & (1u64 << n) == 0 {
                *stateful = true;
            }
        };
        while pu != 0 {
            match self.r32(pu + O_TYPE) {
                COMPL_PUNIT | REPEAT_PUNIT | INV_REP_PUNIT => {
                    need(self.r32(pu + O_U12), assigned, stateful)
                }
                LLIM_PUNIT => {
                    let v = pu + O_U0;
                    for i in 1..=self.r32(v) {
                        need(self.r32(v + 4 * i as i64), assigned, stateful);
                    }
                }
                OR_PUNIT => {
                    let a = self.walk_assigned(self.r64(pu + O_U0), assigned, stateful);
                    let b = self.walk_assigned(self.r64(pu + O_U8), assigned, stateful);
                    assigned = a & b;
                }
                _ => {}
            }
            for n in 0..=MAX_NAMES {
                if self.names(n) == pu {
                    assigned |= 1u64 << n;
                }
            }
            pu = self.r64(pu + O_NXT);
        }
        assigned
    }

    /// Can one sequence be searched in pieces (`first_match_in`)?  The
    /// search must visit start positions in increasing order and the hits
    /// at one start must not depend on where the search began: not with
    /// matches of earlier sequences, a first unit that is an alternative
    /// (its second branch is not tried at later starts), or `^` (it
    /// matches where the search begins).
    pub fn can_split(&self) -> bool {
        let pu = self.r64(self.mem.sa(A_AD_PU_S));
        pu != 0
            && !self.uses_earlier_state()
            && self.r32(pu + O_TYPE) != OR_PUNIT
            && !self.has_match_start(pu, 0)
    }

    fn has_match_start(&self, mut pu: i64, depth: u32) -> bool {
        if depth > 1_000_000 {
            return true;
        }
        while pu != 0 {
            match self.r32(pu + O_TYPE) {
                MATCH_START => return true,
                OR_PUNIT
                    if self.has_match_start(self.r64(pu + O_U0), depth + 1)
                        || self.has_match_start(self.r64(pu + O_U8), depth + 1) =>
                {
                    return true;
                }
                _ => {}
            }
            pu = self.r64(pu + O_NXT);
        }
        false
    }

    /// Use `idx` as the k-mer index of the sequence loaded next (all
    /// engines searching pieces of one record strand share it); None: the
    /// engine builds its own.
    pub fn set_shared_index(&mut self, idx: Option<SharedIndex>) {
        self.shared_kidx = idx;
    }

    /// Threads that build a shared index (the other engines wait for it).
    pub fn set_index_threads(&mut self, n: usize) {
        self.index_threads = n.max(1);
    }

    /// Mark the range units whose following unit only checks a fixed string
    /// at one position (an exact word, an exact reverse complement or an
    /// exact repeat).  Trying a gap length where that unit fails changes
    /// nothing else, so such lengths can be skipped (`next_start`).  Not
    /// used when a unit may read a name matched in an earlier attempt.
    fn setup_skips(&mut self, mut pu: i64, on: bool, depth: u32) {
        if depth > 1_000_000 {
            return;
        }
        let m = self.mem;
        while pu != 0 {
            match self.r32(pu + O_TYPE) {
                RANGE_PUNIT => {
                    let n = Self::next_punit::<Checked>(m, pu);
                    let target = if on && n != 0 && self.skippable(n, pu, true) {
                        n
                    } else {
                        0
                    };
                    self.w64(pu + O_SKIP, target);
                    if target != 0 {
                        let chain = self.skip_chain(pu, target);
                        if !chain.is_empty() {
                            self.skip_chains.push((pu, chain));
                        }
                    }
                }
                OR_PUNIT => {
                    let (a, b) = (self.r64(pu + O_U0), self.r64(pu + O_U8));
                    self.setup_skips(a, on, depth + 1);
                    self.setup_skips(b, on, depth + 1);
                }
                _ => {}
            }
            pu = self.r64(pu + O_NXT);
        }
    }

    /// The units after `n` (the unit after the range `range`) that
    /// `next_start` checks together with `n`.  Each checks a fixed string
    /// right after the one before, and none refers to the range or to `n`
    /// or these units, so when one fails at a gap length the search goes
    /// on with the next length (there is no other choice to try).
    fn skip_chain(&self, range: i64, n: i64) -> Vec<i64> {
        let mut seen = vec![range, n];
        let mut chain = Vec::new();
        let mut u = n;
        while chain.len() < 16 {
            u = Self::next_punit::<Checked>(self.mem, u);
            if u == 0 || !self.skippable(u, range, false) {
                break;
            }
            let t = self.r32(u + O_TYPE);
            if (t == COMPL_PUNIT || t == REPEAT_PUNIT)
                && seen.contains(&self.names(self.r32(u + O_U12)))
            {
                break;
            }
            seen.push(u);
            chain.push(u);
        }
        chain
    }

    /// Can lengths of the range `range` be skipped when `n` follows it?
    /// `n` must check a fixed string: not one taken from the range itself
    /// (its text changes with the gap length).
    /// With `mismatches`, also a reverse complement with mismatches only
    /// (`~pN[m,0,0]`): it has a fixed length and no other choice to try.
    fn skippable(&self, n: i64, range: i64, mismatches: bool) -> bool {
        if self.r32(n + O_ANCH) == 0 {
            return false;
        }
        let t = self.r32(n + O_TYPE);
        if (t == COMPL_PUNIT || t == REPEAT_PUNIT) && self.names(self.r32(n + O_U12)) == range {
            return false;
        }
        let no_errors = |pu: i64| {
            self.r32(pu + O_U0) == 0 && self.r32(pu + O_U4) == 0 && self.r32(pu + O_U8) == 0
        };
        match self.r32(n + O_TYPE) {
            EXACT_PUNIT => self.r32(n + O_U0) >= 1,
            COMPL_PUNIT => {
                self.r32(n + O_U16) == -1
                    && (no_errors(n)
                        || (mismatches && self.r32(n + O_U0) == 0 && self.r32(n + O_U4) == 0))
            }
            REPEAT_PUNIT => no_errors(n),
            _ => false,
        }
    }

    /// The masks the unit `n` checks at its position: data code `d` (low
    /// nibble) matches mask `k` when `d` is one plain base and `d` is in
    /// the mask.  `Some(None)`: the unit matches anywhere (empty capture);
    /// `None`: it never matches.  Only the first `max` masks are made: a
    /// prefix is enough, the matcher checks the rest (and making all of
    /// them costs more than it saves when the unit is long).
    fn unit_masks(&mut self, n: i64, range: i64, max: usize) -> Option<Option<()>> {
        let mut masks = std::mem::take(&mut self.masks);
        masks.clear();
        let r = self.push_unit_masks(n, &mut masks, max);
        self.mask_mis = 0;
        self.mask_tol = 0;
        if self.r32(n + O_TYPE) == COMPL_PUNIT && self.r32(n + O_U8) != 0 {
            self.mask_mis = self.r32(n + O_U8);
            self.mask_tol = masks.len();
        }
        if r == Some(Some(())) {
            let chains = std::mem::take(&mut self.skip_chains);
            if let Some((_, chain)) = chains.iter().find(|c| c.0 == range) {
                for &u in chain {
                    // a prefix of the chain is enough
                    if masks.len() >= max || self.push_unit_masks(u, &mut masks, max).is_none() {
                        break;
                    }
                }
            }
            self.skip_chains = chains;
        }
        self.masks = masks;
        r
    }

    /// Append the masks of unit `n` (see `unit_masks`).
    fn push_unit_masks(&self, n: i64, masks: &mut Vec<u8>, max: usize) -> Option<Option<()>> {
        let m = self.mem;
        match self.r32(n + O_TYPE) {
            EXACT_PUNIT => {
                let code = self.r64(n + O_U8);
                for k in 0..self.r32(n + O_U0) as i64 {
                    if masks.len() >= max {
                        break;
                    }
                    masks.push(m.rb(code + k) & 15);
                }
                Some(Some(()))
            }
            COMPL_PUNIT => {
                let pu1 = m.names(self.r32(n + O_U12));
                let (p1, ln) = (self.r64(pu1 + O_HIT), self.r32(pu1 + O_MLEN) as i64);
                let mut ok = Some(Some(()));
                for k in 0..ln {
                    if masks.len() >= max {
                        break;
                    }
                    let c = m.rb(p1 + ln - 1 - k);
                    if KNOWN_CHAR[(c & 15) as usize] == 0 {
                        ok = None; // the C loop fails at this character
                        break;
                    }
                    masks.push((c >> 4) & 15);
                }
                if ln == 0 { Some(None) } else { ok }
            }
            _ => {
                // REPEAT
                let pu1 = m.names(self.r32(n + O_U12));
                let (p1, ln) = (self.r64(pu1 + O_HIT), self.r32(pu1 + O_MLEN) as i64);
                for k in 0..ln {
                    if masks.len() >= max {
                        break;
                    }
                    masks.push(m.rb(p1 + k) & 15);
                }
                if ln == 0 { Some(None) } else { Some(Some(())) }
            }
        }
    }

    /// Does the data at `x` match the masks (and end by `er`)?  The masks
    /// of the index key are not checked again: `mask_words` and
    /// `mask_bytes` (set by `set_checks`) cover the others that must
    /// match exactly, then the tolerant part is checked.
    #[inline(always)]
    fn masks_match_at(&self, x: i64, er: i64) -> bool {
        let n = self.masks.len() as i64;
        if x + n - 1 > er {
            return false;
        }
        // a plain base mask matches only the same code; x..=x+n-1 lies in
        // the coded sequence
        let words = self.mask_words.iter().all(|&(k, w)| {
            let d = unsafe { ((x as usize + k) as *const u64).read_unaligned() };
            d & 0x0f0f_0f0f_0f0f_0f0f == w
        });
        if !words {
            return false;
        }
        let m = self.mem;
        for &(a, b) in &self.mask_bytes {
            let ok = self.masks[a..b].iter().enumerate().all(|(k, &mk)| {
                let d = m.rb(x + (a + k) as i64) & 15;
                KNOWN_CHAR[d as usize] != 0 && (d & mk) == d
            });
            if !ok {
                return false;
            }
        }
        if self.mask_tol > 0 {
            // as loose_match without inserts and deletes: every base known,
            // at most mask_mis of them not in their mask
            let mut left = self.mask_mis;
            for (k, &mk) in self.masks[..self.mask_tol].iter().enumerate() {
                let d = m.rb(x + k as i64) & 15;
                if KNOWN_CHAR[d as usize] == 0 {
                    return false;
                }
                if (d & mk) != d {
                    left -= 1;
                    if left < 0 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Set how `masks_match_at` checks `masks[from..]`, except the index
    /// key `masks[key..key + KIDX_K]` (`key` is at least `from`): the part
    /// after the key in words when it is at least 8 plain bases.
    fn set_checks(&mut self, from: usize, key: Option<usize>) {
        self.mask_words.clear();
        let len = self.masks.len();
        let (a, b) = match key {
            Some(k) => ((from, k), (k + KIDX_K, len)),
            None => ((from, from), (from, len)),
        };
        self.mask_bytes = [a, (0, 0)];
        if b.1 < b.0 + 8 || self.masks[b.0..b.1].iter().any(|&mk| base2(mk).is_none()) {
            self.mask_bytes[1] = b;
            return;
        }
        let word = |k: usize| u64::from_le_bytes(self.masks[k..k + 8].try_into().unwrap());
        let mut k = b.0;
        while k + 8 < b.1 {
            self.mask_words.push((k, word(k)));
            k += 8;
        }
        self.mask_words.push((b.1 - 8, word(b.1 - 8)));
    }

    /// The first position in `xlo..=xhi` where the unit `n` (marked by
    /// setup_skips) can match: the same position the C search reaches by
    /// trying every gap length.
    #[inline(always)]
    fn next_start(&mut self, n: i64, range: i64, xlo: i64, xhi: i64, er: i64) -> Option<i64> {
        if xlo > xhi {
            return None;
        }
        if xhi - xlo < 8 {
            // a few lengths: making the masks costs more than letting the
            // matcher try each length (most fail at the first base)
            return Some(xlo);
        }
        self.next_start_far(n, range, xlo, xhi, er)
    }

    /// `next_start` for a window of at least 8 positions.
    #[inline(never)]
    fn next_start_far(&mut self, n: i64, range: i64, xlo: i64, xhi: i64, er: i64) -> Option<i64> {
        // masks made: the index key and 3 words to check, or 2 words when
        // every position of a short window is checked
        let max = if xhi - xlo >= SCAN_LIMIT {
            KIDX_K + 24
        } else {
            16
        };
        match self.unit_masks(n, range, max) {
            None => return None,
            Some(None) => return Some(xlo),
            Some(Some(())) => {}
        }
        let len = self.masks.len();
        let tol = self.mask_tol;
        // the index key: the first KIDX_K plain bases that must match
        let key_at = if xhi - xlo >= SCAN_LIMIT && len >= tol + KIDX_K {
            (tol..=len - KIDX_K).find(|&k| {
                self.masks[k..k + KIDX_K]
                    .iter()
                    .all(|&mk| base2(mk).is_some())
            })
        } else {
            None
        };
        let Some(ko) = key_at else {
            self.set_checks(tol, None);
            return (xlo..=xhi).find(|&x| self.masks_match_at(x, er));
        };
        let key = self.masks[ko..ko + KIDX_K]
            .iter()
            .fold(0u32, |key, &mk| key << 2 | base2(mk).unwrap());
        // the index holds the positions where the key matches
        self.set_checks(tol, Some(ko));
        let base = self.mem.cd as i64;
        if self.shared_kidx.is_none() {
            self.build_kidx(er);
        }
        let idx = match &self.shared_kidx {
            Some(sh) => sh.get_or_init(|| {
                let n = (er - base + 1).max(0) as usize;
                let seq = unsafe { std::slice::from_raw_parts(self.mem.cd as *const u8, n) };
                build_index_par(seq, self.index_threads)
            }),
            None => self.kidx.as_ref().unwrap(),
        };
        let list =
            &idx.pos[idx.starts[key as usize] as usize..idx.starts[key as usize + 1] as usize];
        // the key of a start x is at x + ko
        let lo = (xlo - base + ko as i64) as u32;
        let first = list.partition_point(|&p| p < lo);
        for &p in &list[first..] {
            let x = base + p as i64 - ko as i64;
            if x > xhi {
                break;
            }
            if self.masks_match_at(x, er) {
                return Some(x);
            }
        }
        None
    }

    /// Build the k-mer index of the current coded sequence (up to `er`).
    fn build_kidx(&mut self, er: i64) {
        if self
            .kidx
            .as_ref()
            .is_some_and(|k| k.generation == self.cd_gen)
        {
            return;
        }
        let m = self.mem;
        let base = m.cd as i64;
        let n = (er - base + 1).max(0) as usize;
        let seq = unsafe { std::slice::from_raw_parts(m.cd as *const u8, n) };
        let mut idx = build_index(seq, self.kidx.take());
        idx.generation = self.cd_gen;
        self.kidx = Some(idx);
    }

    /// Name `n` is defined, and following the names it refers to does not
    /// come back to a name already seen.
    fn name_ok(&self, n: i32) -> bool {
        let mut seen = [false; N_NAMES as usize];
        let mut n = n;
        loop {
            if seen[n as usize] {
                return false;
            }
            seen[n as usize] = true;
            let pu = self.names(n);
            if pu == 0 {
                return false;
            }
            match self.r32(pu + O_TYPE) {
                COMPL_PUNIT | REPEAT_PUNIT | INV_REP_PUNIT => n = self.r32(pu + O_U12),
                _ => return true,
            }
        }
    }

    fn max_mats(&self, mut pu: i64) -> i32 {
        let mut sum: i32 = 0;
        while pu != 0 {
            sum = sum.wrapping_add(self.max_mat(pu, 0));
            pu = self.r64(pu + O_NXT);
        }
        sum
    }

    /// parse_dna_cmd / parse_peptide_cmd.  Returns the C return value
    /// (0 means "failed to parse").
    pub fn parse_cmd(&mut self, line: &[u8], seq_type: i32) -> i32 {
        self.seq_type = seq_type;
        self.mem.pep = seq_type == PEPTIDE;
        self.w32(self.mem.sa(A_SEQ_TYPE), seq_type);
        if self.r32(self.mem.sa(A_INITIALIZED)) == 0 {
            self.build_conversion_tables();
        }
        for i in 0..N_NAMES {
            self.w64(self.mem.sa(A_NAMES) + 8 * i, 0);
        }
        self.ivp = self.mem.sa(A_IV);
        self.cvp = self.mem.sa(A_CV);
        self.pup = self.mem.sa(A_PU_S);
        self.w64(self.mem.sa(A_AD_PU_S), self.mem.sa(A_PU_S));
        if line.len() >= 1_000_000 {
            self.w64(self.mem.sa(A_AD_PU_S), 0);
            return 0;
        }
        match self.parser(line) {
            None => {
                self.w64(self.mem.sa(A_AD_PU_S), 0);
                0
            }
            Some(_) => {
                if !self.refs_ok(self.mem.sa(A_PU_S)) {
                    self.w64(self.mem.sa(A_AD_PU_S), 0);
                    return 0;
                }
                self.set_anchors(self.mem.sa(A_PU_S), 0);
                let skip = self.seq_type == DNA && !self.uses_earlier_state();
                self.skip_chains.clear();
                self.setup_skips(self.mem.sa(A_PU_S), skip, 0);
                self.w64(self.mem.sa(A_AD_PU_S), self.mem.sa(A_PU_S));
                self.max_mats(self.mem.sa(A_PU_S))
            }
        }
    }

    // ------------------------------------------------------------------
    // matcher

    /// The first `n` bytes of loose_match's first operand as a slice, when
    /// they all lie in one buffer.
    #[inline(always)]
    fn one_slice<'a>(&self, one: &One<'a>, n: usize) -> Option<&'a [u8]> {
        let m = self.mem;
        match one {
            One::Bytes(b) => {
                if b.len() >= n {
                    Some(&b[..n])
                } else {
                    None
                }
            }
            One::Mem(base) => {
                if m.cd_ok(*base, n as i64) {
                    Some(unsafe { std::slice::from_raw_parts(*base as *const u8, n) })
                } else {
                    let o = base.wrapping_sub(m.s as i64);
                    if o >= 0 && (o as usize) + n <= S_LEN {
                        Some(unsafe { std::slice::from_raw_parts(*base as *const u8, n) })
                    } else {
                        None
                    }
                }
            }
        }
    }

    /// loose_match with compl_flag for operands longer than MAX_CODES (the
    /// C program aborts there): the complemented operand is on the heap.
    #[cold]
    #[inline(never)]
    #[allow(clippy::too_many_arguments)]
    fn loose_match_long(
        &self,
        one_src: One,
        one_len: i32,
        two_start: i64,
        two_len: i32,
        max_ins: i32,
        max_del: i32,
        max_mis: i32,
        rule_set: i32,
    ) -> i32 {
        let m = self.mem;
        let n = one_len as usize;
        let mut v = vec![0u8; n];
        for k in 0..n {
            v[n - 1 - k] = compl_code(m.one_get(&one_src, k as i64) as u32, rule_set);
        }
        self.loose_match(
            One::Bytes(&v),
            one_len,
            two_start,
            two_len,
            max_ins,
            max_del,
            max_mis,
            rule_set,
            false,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn loose_match(
        &self,
        one_src: One,
        mut one_len: i32,
        two_start: i64,
        mut two_len: i32,
        mut max_ins: i32,
        mut max_del: i32,
        mut max_mis: i32,
        rule_set: i32,
        compl_flag: bool,
    ) -> i32 {
        let m = self.mem;
        // complemented operand (only the written part is read); longer
        // operands go through loose_match_long
        let mut result = [std::mem::MaybeUninit::<u8>::uninit(); MAX_CODES as usize];
        let one_src = if compl_flag {
            let mut i = 0;
            while i < one_len {
                if !m.known_char(m.one_get(&one_src, i as i64) & 15) {
                    return 0;
                }
                i += 1;
            }
            if one_len > MAX_CODES {
                return self.loose_match_long(
                    one_src, one_len, two_start, two_len, max_ins, max_del, max_mis, rule_set,
                );
            }
            let n = one_len.max(0) as usize;
            let mut len = one_len;
            let mut k: i64 = 0;
            while len > 0 {
                let i = m.one_get(&one_src, k) as u32;
                k += 1;
                len -= 1;
                result[len as usize] = std::mem::MaybeUninit::new(compl_code(i, rule_set));
            }
            One::Bytes(unsafe { std::slice::from_raw_parts(result.as_ptr() as *const u8, n) })
        } else {
            one_src
        };
        let mut one: i64 = 0;
        let mut two: i64 = two_start;

        if max_ins == 0 && max_del == 0 {
            if one_len > two_len {
                return 0;
            }
            // same loop on plain slices when every byte it may read is in
            // a known buffer
            if one_len >= 1 && m.cd_ok(two_start, one_len as i64) {
                let n = one_len as usize;
                if let Some(ob) = self.one_slice(&one_src, n) {
                    let tb = unsafe { std::slice::from_raw_parts(two_start as *const u8, n) };
                    for k in 0..n {
                        let t = tb[k];
                        if !m.known_char(t & 15)
                            || (!m.ex_matches(rule_set, t, ob[k]) && {
                                max_mis = max_mis.wrapping_sub(1);
                                max_mis < 0
                            })
                        {
                            return 0;
                        }
                    }
                    return one_len + 1;
                }
            }
            let mut i = one_len;
            while i >= 1 {
                let t = m.rb(two);
                if !m.known_char(t & 15)
                    || (!m.ex_matches(rule_set, t, m.one_get(&one_src, one)) && {
                        max_mis = max_mis.wrapping_sub(1);
                        max_mis < 0
                    })
                {
                    return 0;
                }
                two += 1;
                one += 1;
                i -= 1;
            }
            return (two - two_start + 1) as i32;
        }

        // choice points: the first STACK_INLINE on the stack, the rest on
        // the heap (the C program has a fixed stack[100])
        const STACK_INLINE: usize = 64;
        let mut stack = [std::mem::MaybeUninit::<StackEnt>::uninit(); STACK_INLINE];
        let mut stack_heap: Vec<StackEnt> = Vec::new();
        let mut nxtent: i32 = 0;
        // every entry below nxtent was written by push
        macro_rules! ent {
            ($k:expr) => {{
                let k = $k as usize;
                if k < STACK_INLINE {
                    unsafe { &mut *stack[k].as_mut_ptr() }
                } else {
                    &mut stack_heap[k - STACK_INLINE]
                }
            }};
        }
        macro_rules! ret {
            ($v:expr) => {{
                return $v;
            }};
        }
        macro_rules! push {
            ($n:expr) => {{
                let k = nxtent as usize;
                let e = StackEnt {
                    p1: one,
                    p2: two,
                    n1: one_len,
                    n2: two_len,
                    mis: max_mis,
                    ins: max_ins,
                    del: max_del,
                    next_choice: $n,
                };
                if k < STACK_INLINE {
                    stack[k] = std::mem::MaybeUninit::new(e);
                } else if k - STACK_INLINE < stack_heap.len() {
                    stack_heap[k - STACK_INLINE] = e;
                } else {
                    stack_heap.push(e);
                }
                nxtent += 1;
            }};
        }
        let mut lbl = Lm::Top;
        loop {
            match lbl {
                Lm::Top => {
                    if !(two_len != 0 || nxtent != 0) {
                        ret!(0);
                    }
                    if two_len != 0 && one_len != 0 && {
                        let t = m.rb(two);
                        m.known_char(t & 15) && m.ex_matches(rule_set, t, m.one_get(&one_src, one))
                    } {
                        two += 1;
                        one += 1;
                        two_len = two_len.wrapping_sub(1);
                        one_len = one_len.wrapping_sub(1);
                        if one_len == 0 {
                            ret!((two - two_start + 1) as i32);
                        }
                    } else if max_mis != 0 && one_len >= 1 && two_len >= 1 {
                        if max_ins != 0 {
                            push!(1);
                        } else if max_del != 0 {
                            push!(2);
                        }
                        max_mis = max_mis.wrapping_sub(1);
                        one += 1;
                        two += 1;
                        one_len -= 1;
                        two_len -= 1;
                        if one_len == 0 {
                            ret!((two - two_start + 1) as i32);
                        }
                    } else if max_ins != 0 && one_len >= 1 {
                        if max_del != 0 && two_len >= 1 {
                            push!(2);
                        }
                        lbl = Lm::Ins;
                    } else if max_del != 0 && two_len >= 1 {
                        lbl = Lm::Del;
                    } else if nxtent != 0 {
                        let e = *ent!(nxtent - 1);
                        one = e.p1;
                        two = e.p2;
                        one_len = e.n1;
                        two_len = e.n2;
                        max_mis = e.mis;
                        max_ins = e.ins;
                        max_del = e.del;
                        if e.next_choice == 1 {
                            if max_del != 0 {
                                // the "delete" choice of this point is still
                                // to come: the entry stays on the stack (the
                                // C code released it and lost that choice)
                                ent!(nxtent - 1).next_choice = 2;
                            } else {
                                nxtent -= 1;
                            }
                            lbl = Lm::Ins;
                        } else {
                            nxtent -= 1;
                            lbl = Lm::Del;
                        }
                    } else {
                        ret!(0);
                    }
                }
                Lm::Ins => {
                    max_ins = max_ins.wrapping_sub(1);
                    one += 1;
                    one_len = one_len.wrapping_sub(1);
                    if one_len == 0 {
                        ret!((two - two_start + 1) as i32);
                    }
                    lbl = Lm::Top;
                }
                Lm::Del => {
                    max_del = max_del.wrapping_sub(1);
                    two += 1;
                    two_len = two_len.wrapping_sub(1);
                    if one_len == 0 {
                        ret!((two - two_start + 1) as i32);
                    }
                    lbl = Lm::Top;
                }
            }
        }
    }

    /// next_punit(); in fast mode every pointer is checked before it is
    /// followed, and the checked version takes over if one is outside.
    #[inline(always)]
    fn next_punit<P: Pol>(m: M, pu: i64) -> i64 {
        if !P::CHECKED
            && let Some(n) = Self::next_punit_fast(m, pu)
        {
            return n;
        }
        let n = m.r64(pu + O_NXT);
        if n != 0 {
            return n;
        }
        let mut pu1 = pu;
        let mut pu2 = m.r64(pu1 + O_PREV);
        while pu2 != 0 {
            let n = m.r64(pu2 + O_NXT);
            if n == pu1 || n == 0 {
                pu1 = pu2;
                pu2 = m.r64(pu1 + O_PREV);
            } else {
                break;
            }
        }
        if pu2 != 0 && m.r64(pu2 + O_NXT) != pu1 {
            return m.r64(pu2 + O_NXT);
        }
        0
    }

    #[inline(always)]
    fn next_punit_fast(m: M, pu: i64) -> Option<i64> {
        let n = hr64(m.up(pu), O_NXT);
        if n != 0 {
            return Some(n);
        }
        let mut pu1 = pu;
        let mut pu2 = hr64(m.up(pu1), O_PREV);
        while pu2 != 0 {
            if !m.pu_ok(pu2) {
                return None;
            }
            let n = hr64(m.up(pu2), O_NXT);
            if n == pu1 || n == 0 {
                pu1 = pu2;
                pu2 = hr64(m.up(pu1), O_PREV);
            } else {
                break;
            }
        }
        if pu2 != 0 {
            let n = hr64(m.up(pu2), O_NXT);
            if n != pu1 {
                return Some(n);
            }
        }
        Some(0)
    }

    fn collect_hits(&self, pu: i64, revhits: &mut Vec<i64>) {
        let m = self.mem;
        if !Self::collect_hits_fast(m, pu, revhits) {
            Self::collect_hits_checked(m, pu, revhits);
        }
    }

    /// collect_hits() while every unit pointer lies in the static block;
    /// false if one does not (the caller then redoes it checked).
    fn collect_hits_fast(m: M, pu: i64, revhits: &mut Vec<i64>) -> bool {
        revhits.clear();
        if revhits.capacity() < 256 {
            revhits.reserve(256);
        }
        let mut cap = revhits.capacity();
        let mut out = revhits.as_mut_ptr();
        let mut n = 0usize;
        let mut last: i64 = 0;
        let mut pu = pu;
        while pu != 0 {
            if !m.pu_ok(pu) {
                return false;
            }
            let up = m.up(pu);
            if hr32(up, O_TYPE) != OR_PUNIT {
                if n == cap {
                    if n > 1 << 20 {
                        return false; // let the checked walk run off the stack
                    }
                    unsafe { revhits.set_len(n) };
                    revhits.reserve(n);
                    cap = revhits.capacity();
                    out = revhits.as_mut_ptr();
                }
                unsafe { *out.add(n) = hr64(up, O_HIT) };
                n += 1;
                last = pu;
                pu = hr64(up, O_PREV);
            } else if hr64(up, O_NXT) == last {
                pu = if hr32(up, O_U24) == 1 {
                    hr64(up, O_U0)
                } else {
                    hr64(up, O_U8)
                };
                loop {
                    if !m.pu_ok(pu) {
                        return false;
                    }
                    let nx = hr64(m.up(pu), O_NXT);
                    if nx == 0 {
                        break;
                    }
                    pu = nx;
                }
                last = 0;
            } else {
                last = pu;
                pu = hr64(up, O_PREV);
            }
        }
        unsafe { revhits.set_len(n) };
        true
    }

    /// The unit slot of each hit entry, in the order of the entries: the
    /// walk of collect_hits_checked, recording the unit instead of its hit.
    fn collect_slots(&mut self, pu: i64) {
        let m = self.mem;
        let base = m.sa(A_PU_S);
        self.hit_slots.clear();
        let mut last: i64 = 0;
        let mut pu = pu;
        while pu != 0 {
            if m.r32(pu + O_TYPE) != OR_PUNIT {
                let k = pu.wrapping_sub(base);
                self.hit_slots
                    .push(if k >= 0 && k % SZ == 0 && k / SZ < N_SLOTS {
                        (k / SZ) as u32
                    } else {
                        u32::MAX
                    });
                last = pu;
                pu = m.r64(pu + O_PREV);
            } else if m.r64(pu + O_NXT) == last {
                pu = if m.r32(pu + O_U24) == 1 {
                    m.r64(pu + O_U0)
                } else {
                    m.r64(pu + O_U8)
                };
                loop {
                    let n = m.r64(pu + O_NXT);
                    if n == 0 {
                        break;
                    }
                    pu = n;
                }
                last = 0;
            } else {
                last = pu;
                pu = m.r64(pu + O_PREV);
            }
        }
        self.hit_slots.reverse();
    }

    /// Record which unit gave each entry of a hit (see `hit_slots`).
    pub fn set_track_units(&mut self, on: bool) {
        self.track_units = on;
    }

    /// The unit slot of each entry of the last hit (`set_track_units`);
    /// u32::MAX for an entry that is not a unit of the pattern.
    pub fn hit_slots(&self) -> &[u32] {
        &self.hit_slots
    }

    /// (slot, byte offset in the pattern line) of every unit of the
    /// parsed pattern except the alternatives `( | )` themselves.
    pub fn unit_offsets(&self) -> Vec<(u32, usize)> {
        let base = self.mem.sa(A_PU_S);
        let n = ((self.pup - base) / SZ) as usize;
        (0..n)
            .filter(|&k| self.r32(base + k as i64 * SZ + O_TYPE) != OR_PUNIT)
            .filter_map(|k| {
                let off = *self.unit_offs.get(k)?;
                (off != u32::MAX).then_some((k as u32, off as usize))
            })
            .collect()
    }

    fn collect_hits_checked(m: M, pu: i64, revhits: &mut Vec<i64>) {
        revhits.clear();
        let mut last: i64 = 0;
        let mut pu = pu;
        while pu != 0 {
            if m.r32(pu + O_TYPE) != OR_PUNIT {
                revhits.push(m.r64(pu + O_HIT));
                if revhits.len() > 1 << 20 {
                    segv(); // runs off the top of the C stack
                }
                last = pu;
                pu = m.r64(pu + O_PREV);
            } else if m.r64(pu + O_NXT) == last {
                pu = if m.r32(pu + O_U24) == 1 {
                    m.r64(pu + O_U0)
                } else {
                    m.r64(pu + O_U8)
                };
                loop {
                    let n = m.r64(pu + O_NXT);
                    if n == 0 {
                        break;
                    }
                    pu = n;
                }
                last = 0;
            } else {
                last = pu;
                pu = m.r64(pu + O_PREV);
            }
        }
    }

    fn pattern_match(
        &mut self,
        pu: i64,
        start: i64,
        end: i64,
        hits: &mut Vec<i64>,
        first: bool,
    ) -> i32 {
        let br = self.br1;
        match self.pm::<Fast>(start, end, hits, !first, pu, start, br) {
            Out::Done(n) => n,
            Out::Bail { back, cr, sr, br } => {
                match self.pm::<Checked>(start, end, hits, back, cr, sr, br) {
                    Out::Done(n) => n,
                    Out::Bail { .. } => unreachable!(),
                }
            }
        }
    }

    /// The body of pattern_match() from label TRY (or BACKTRACK when
    /// `back`) with the given registers.
    #[allow(clippy::too_many_arguments)]
    #[inline(never)]
    fn pm<P: Pol>(
        &mut self,
        start: i64,
        end: i64,
        hits: &mut Vec<i64>,
        back: bool,
        cr0: i64,
        sr0: i64,
        br0: i64,
    ) -> Out {
        let m = self.mem;
        let mut br = br0;
        let mut sr = sr0;
        let er = end;
        let mut cr = cr0;
        let mut skip_try = back;
        // host pointer to the current unit (fast mode)
        let mut cp: *mut u8 = m.s;

        // fields of the current unit
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r32(cr + $o)
                } else {
                    hr32(cp, $o)
                }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r64(cr + $o)
                } else {
                    hr64(cp, $o)
                }
            };
        }
        macro_rules! s32 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w32(cr + $o, $v)
                } else {
                    hw32(cp, $o, $v)
                }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w64(cr + $o, $v)
                } else {
                    hw64(cp, $o, $v)
                }
            };
        }

        // SUCCESS: record the match length, go on with the next unit (TRY)
        // or, after the last unit, LEAVE
        macro_rules! success {
            ($try:lifetime, $leave:lifetime) => {{
                let h = g64!(O_HIT);
                s32!(O_MLEN, (sr - h) as i32);
                let n = Self::next_punit::<P>(m, cr);
                if n != 0 {
                    cr = n;
                    continue $try;
                } else {
                    break $leave;
                }
            }};
        }
        macro_rules! push_br {
            () => {{
                s64!(O_BR, br);
                br = cr;
            }};
        }

        'main: loop {
            if !skip_try {
                // TRY
                'tryl: loop {
                    if !P::CHECKED {
                        if !m.pu_ok(cr) {
                            return Out::Bail {
                                back: false,
                                cr,
                                sr,
                                br,
                            };
                        }
                        cp = m.up(cr);
                    }
                    match g32!(O_TYPE) {
                        MATCH_START => {
                            if sr == start {
                                s64!(O_HIT, sr);
                                success!('tryl, 'main);
                            } else {
                                break 'tryl;
                            }
                        }
                        MATCH_END => {
                            if sr == end + 1 {
                                s64!(O_HIT, sr);
                                success!('tryl, 'main);
                            } else {
                                break 'tryl;
                            }
                        }
                        ANY_PUNIT => {
                            let (n, c) = self.t_any::<P>(cr, sr, er);
                            if c == 0 {
                                break 'tryl;
                            }
                            if c == 2 {
                                push_br!();
                            }
                            sr = n;
                            success!('tryl, 'main);
                        }
                        LLIM_PUNIT => {
                            let (n, c) = self.t_llim::<P>(cr, sr, er);
                            if c == 0 {
                                break 'tryl;
                            }
                            if c == 2 {
                                push_br!();
                            }
                            sr = n;
                            success!('tryl, 'main);
                        }
                        RANGE_PUNIT => {
                            let min = g32!(O_U0);
                            let i = min.wrapping_sub(1);
                            if sr + i as i64 <= er {
                                s64!(O_HIT, sr);
                                sr += i as i64 + 1;
                                if (sr <= er && g32!(O_U4) != 0) || (g32!(O_ANCH) == 0 && sr <= er)
                                {
                                    push_br!();
                                    let mn = g32!(O_U0);
                                    s32!(O_U8, mn.wrapping_add(1));
                                    let n = if P::CHECKED { 0 } else { g64!(O_SKIP) };
                                    if n != 0 {
                                        // go straight to the first gap length
                                        // where the next unit can match
                                        let hit = g64!(O_HIT);
                                        // the first length (min) is always
                                        // tried, then min+1 ..= min+width
                                        let lmax = (mn as i64 + g32!(O_U4) as i64)
                                            .min(er - hit + 1)
                                            .max(mn as i64);
                                        match self.next_start(n, cr, sr, hit + lmax, er) {
                                            Some(x) => {
                                                s32!(O_U8, (x - hit + 1) as i32);
                                                sr = x;
                                            }
                                            None => {
                                                // every length fails: as after trying them all
                                                s32!(
                                                    O_U8,
                                                    (mn as i64 + g32!(O_U4) as i64 + 1) as i32
                                                );
                                                break 'tryl;
                                            }
                                        }
                                    }
                                }
                                success!('tryl, 'main);
                            } else {
                                break 'tryl;
                            }
                        }
                        EXACT_PUNIT => {
                            let len = g32!(O_U0);
                            let mut last = er + 1 - len as i64;
                            if last > sr && g32!(O_ANCH) != 0 {
                                last = sr;
                            }
                            if last > self.start_lim && g32!(O_ANCH) == 0 {
                                last = self.start_lim;
                            }
                            let p1 = g64!(O_U8);
                            let ln = len.wrapping_sub(1);
                            let c0 = m.rb(p1);
                            while sr <= last {
                                if m.matches(m.rb(sr), c0) {
                                    let mut p2 = sr + 1;
                                    let mut p3 = p1 + 1;
                                    let mut i = ln;
                                    while i != 0 && m.matches(m.rb(p2), m.rb(p3)) {
                                        i = i.wrapping_sub(1);
                                        p3 += 1;
                                        p2 += 1;
                                    }
                                    if i == 0 {
                                        break;
                                    }
                                }
                                sr += 1;
                            }
                            if sr > last {
                                break 'tryl;
                            }
                            s64!(O_HIT, sr);
                            if sr < last {
                                push_br!();
                            }
                            sr += g32!(O_U0) as i64;
                            success!('tryl, 'main);
                        }
                        COMPL_PUNIT => {
                            let pu1 = m.names_p::<P>(g32!(O_U12));
                            let pp = if !P::CHECKED && m.pu_ok(pu1) {
                                m.up(pu1)
                            } else {
                                std::ptr::null_mut()
                            };
                            let pu1_hit = |m: M| {
                                if pp.is_null() {
                                    m.r64(pu1 + O_HIT)
                                } else {
                                    hr64(pp, O_HIT)
                                }
                            };
                            let p1 = pu1_hit(m);
                            let mut ln = if pp.is_null() {
                                m.r32(pu1 + O_MLEN)
                            } else {
                                hr32(pp, O_MLEN)
                            };
                            let rs = g32!(O_U16);
                            let (ins, del, mis) = (g32!(O_U0), g32!(O_U4), g32!(O_U8));
                            if rs != -1 || ins != 0 || del != 0 || mis != 0 {
                                let i = self.loose_match(
                                    One::Mem(p1),
                                    ln,
                                    sr,
                                    (er + 1 - sr) as i32,
                                    ins,
                                    del,
                                    mis,
                                    rs,
                                    true,
                                );
                                if i != 0 {
                                    let i = i - 1;
                                    s64!(O_HIT, sr);
                                    sr += i as i64;
                                    success!('tryl, 'main);
                                } else {
                                    break 'tryl;
                                }
                            } else if er - sr >= ln as i64 - 1 {
                                s64!(O_HIT, sr);
                                // pu1->hit is read again, as in C (pu1 may be this unit)
                                let mut q = pu1_hit(m) + (ln as i64 - 1);
                                let mut ok = true;
                                if !P::CHECKED
                                    && !m.pep
                                    && ln > 0
                                    && m.cd_ok(q - (ln as i64 - 1), ln as i64)
                                    && m.cd_ok(sr, ln as i64)
                                {
                                    // all reads are inside the sequence buffer
                                    let qp = q as *const u8;
                                    let sp = sr as *const u8;
                                    let n = ln as usize;
                                    let mut k = 0usize;
                                    while k < n {
                                        let c = unsafe { *qp.wrapping_sub(k) };
                                        if KNOWN_CHAR[(c & 15) as usize] == 0
                                            || ((c >> 4) & 15) != (unsafe { *sp.add(k) } & 15)
                                        {
                                            ok = false;
                                            break;
                                        }
                                        k += 1;
                                    }
                                    ln = 0;
                                    q -= k as i64;
                                    sr += k as i64;
                                    let _ = q;
                                }
                                while ln != 0 {
                                    ln = ln.wrapping_sub(1);
                                    let c = m.rb(q);
                                    if !m.known_char(c & 15) || ((c >> 4) & 15) != (m.rb(sr) & 15) {
                                        ok = false;
                                        break;
                                    }
                                    q -= 1;
                                    sr += 1;
                                }
                                if ok {
                                    success!('tryl, 'main);
                                } else {
                                    break 'tryl;
                                }
                            } else {
                                break 'tryl;
                            }
                        }
                        REPEAT_PUNIT => {
                            let (n, c) = self.t_repeat::<P>(cr, sr, er);
                            if c == 0 {
                                break 'tryl;
                            }
                            if c == 2 {
                                push_br!();
                            }
                            sr = n;
                            success!('tryl, 'main);
                        }
                        INV_REP_PUNIT => {
                            let (n, c) = self.t_inv_rep::<P>(cr, sr, er);
                            if c == 0 {
                                break 'tryl;
                            }
                            if c == 2 {
                                push_br!();
                            }
                            sr = n;
                            success!('tryl, 'main);
                        }
                        SIM_PUNIT => {
                            let (n, c) = self.t_sim::<P>(cr, sr, er);
                            if c == 0 {
                                break 'tryl;
                            }
                            if c == 2 {
                                push_br!();
                            }
                            sr = n;
                            success!('tryl, 'main);
                        }
                        WEIGHT_PUNIT => {
                            let (n, c) = self.t_weight::<P>(cr, sr, er);
                            if c == 0 {
                                break 'tryl;
                            }
                            if c == 2 {
                                push_br!();
                            }
                            sr = n;
                            success!('tryl, 'main);
                        }
                        OR_PUNIT => {
                            push_br!();
                            s64!(O_U16, sr);
                            s64!(O_HIT, sr);
                            s32!(O_U24, 1);
                            cr = g64!(O_U0);
                            continue 'tryl;
                        }
                        _ => break 'tryl,
                    }
                }
            }
            skip_try = false;
            // BACKTRACK
            'backl: loop {
                if br == 0 {
                    return Out::Done(0);
                }
                if !P::CHECKED {
                    if !m.pu_ok(br) {
                        return Out::Bail {
                            back: true,
                            cr,
                            sr,
                            br,
                        };
                    }
                    cp = m.up(br);
                }
                cr = br;
                br = g64!(O_BR);
                sr = g64!(O_HIT);
                match g32!(O_TYPE) {
                    RANGE_PUNIT => {
                        let nx = g32!(O_U8);
                        let min = g32!(O_U0);
                        let width = g32!(O_U4);
                        let n = if P::CHECKED { 0 } else { g64!(O_SKIP) };
                        if n != 0 {
                            // same as below, but only gap lengths where the
                            // next unit can match are tried
                            let mut hit = sr;
                            let mut from = nx as i64;
                            loop {
                                let mut lmax = (min as i64 + width as i64).min(er - hit + 1);
                                if from == min as i64 {
                                    // a new start: its first length is
                                    // always tried
                                    lmax = lmax.max(min as i64);
                                }
                                if let Some(x) = self.next_start(n, cr, hit + from, hit + lmax, er)
                                {
                                    s64!(O_HIT, hit);
                                    s32!(O_U8, (x - hit + 1) as i32);
                                    sr = x;
                                    br = cr;
                                    success!('main, 'main);
                                }
                                // all lengths fail: next start, if unanchored
                                hit += 1;
                                s64!(O_HIT, hit);
                                if hit + min as i64 - 1 <= er
                                    && g32!(O_ANCH) == 0
                                    && hit <= self.start_lim
                                {
                                    from = min as i64;
                                } else {
                                    continue 'backl;
                                }
                            }
                        }
                        if nx <= min.wrapping_add(width) && sr + nx as i64 - 1 <= er {
                            s32!(O_U8, nx.wrapping_add(1));
                            sr += nx as i64;
                            br = cr;
                            success!('main, 'main);
                        } else {
                            let h = g64!(O_HIT) + 1;
                            s64!(O_HIT, h);
                            if h + g32!(O_U0) as i64 - 1 <= er
                                && g32!(O_ANCH) == 0
                                && h <= self.start_lim
                            {
                                let mn = g32!(O_U0);
                                s32!(O_U8, mn.wrapping_add(1));
                                sr = g64!(O_HIT) + g32!(O_U0) as i64;
                                br = cr;
                                success!('main, 'main);
                            } else {
                                continue 'backl;
                            }
                        }
                    }
                    LLIM_PUNIT | ANY_PUNIT | EXACT_PUNIT | SIM_PUNIT | WEIGHT_PUNIT => {
                        sr += 1;
                        continue 'main;
                    }
                    OR_PUNIT => {
                        if g32!(O_U24) == 1 {
                            s32!(O_U24, 2);
                            cr = g64!(O_U8);
                            continue 'main;
                        } else if g32!(O_ANCH) == 0 && g32!(O_U24) == 2 {
                            sr += 1;
                            continue 'main;
                        } else {
                            continue 'backl;
                        }
                    }
                    _ => break 'main,
                }
            }
        }
        // LEAVE
        {
            let mut revhits = std::mem::take(&mut self.revhits);
            self.collect_hits(cr, &mut revhits);
            if self.track_units {
                self.collect_slots(cr);
            }
            let n = revhits.len();
            if hits.len() < n + 1 {
                hits.resize(n + 1, 0);
            }
            hits[..n].copy_from_slice(&revhits[..n]);
            hits[..n].reverse();
            hits[n] = sr;
            self.revhits = revhits;
            self.br1 = br;
            Out::Done(n as i32)
        }
    }

    #[inline(never)]
    #[allow(unused_macros, unused_mut, unused_variables)]
    fn t_any<P: Pol>(&self, cr: i64, mut sr: i64, er: i64) -> (i64, u8) {
        let m = self.mem;
        let mut pushed = false;
        let cp = m.up(cr);
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r32(cr + $o)
                } else {
                    hr32(cp, $o)
                }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r64(cr + $o)
                } else {
                    hr64(cp, $o)
                }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w64(cr + $o, $v)
                } else {
                    hw64(cp, $o, $v)
                }
            };
        }
        let mut last = er;
        if last > sr && g32!(O_ANCH) != 0 {
            last = sr;
        }
        if last > self.start_lim && g32!(O_ANCH) == 0 {
            last = self.start_lim;
        }
        let cm = g64!(O_U0);
        while sr <= last {
            let i = m.rb(sr) as i8 as i32;
            if i >= b'A' as i32
                && i <= b'Z' as i32
                && ((1i32 << (i - b'A' as i32)) as i64 & cm) != 0
            {
                break;
            }
            sr += 1;
        }
        if sr > last {
            (0, 0)
        } else {
            s64!(O_HIT, sr);
            if sr < last {
                {
                    pushed = true;
                }
            }
            sr += 1;
            (sr, if pushed { 2 } else { 1 })
        }
    }

    #[inline(never)]
    #[allow(unused_macros, unused_mut, unused_variables)]
    fn t_llim<P: Pol>(&self, cr: i64, mut sr: i64, er: i64) -> (i64, u8) {
        let m = self.mem;
        let mut pushed = false;
        let cp = m.up(cr);
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r32(cr + $o)
                } else {
                    hr32(cp, $o)
                }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r64(cr + $o)
                } else {
                    hr64(cp, $o)
                }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w64(cr + $o, $v)
                } else {
                    hw64(cp, $o, $v)
                }
            };
        }
        let v = cr + O_U0;
        let mut ln: i32 = 0;
        let mut i = m.r32(v);
        while i != 0 {
            let nm = m.r32(v + 4 * i as i64);
            i = i.wrapping_sub(1);
            let pu1 = m.names_p::<P>(nm);
            ln = ln.wrapping_add(m.pu_r32::<P>(pu1, O_MLEN));
        }
        if ln < g32!(O_LLIM_BOUND) {
            s64!(O_HIT, sr);
            (sr, if pushed { 2 } else { 1 })
        } else {
            (0, 0)
        }
    }

    #[inline(never)]
    #[allow(unused_macros, unused_mut, unused_variables)]
    fn t_repeat<P: Pol>(&self, cr: i64, mut sr: i64, er: i64) -> (i64, u8) {
        let m = self.mem;
        let mut pushed = false;
        let cp = m.up(cr);
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r32(cr + $o)
                } else {
                    hr32(cp, $o)
                }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r64(cr + $o)
                } else {
                    hr64(cp, $o)
                }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w64(cr + $o, $v)
                } else {
                    hw64(cp, $o, $v)
                }
            };
        }
        let pu1 = m.names_p::<P>(g32!(O_U12));
        let p1 = m.pu_r64::<P>(pu1, O_HIT);
        let ln = m.pu_r32::<P>(pu1, O_MLEN);
        if ln == 0 {
            s64!(O_HIT, sr);
            (sr, if pushed { 2 } else { 1 })
        } else {
            let (ins, del, mis) = (g32!(O_U0), g32!(O_U4), g32!(O_U8));
            let i = self.loose_match(
                One::Mem(p1),
                ln,
                sr,
                (er + 1 - sr) as i32,
                ins,
                del,
                mis,
                -1,
                false,
            );
            if i != 0 {
                let i = i - 1;
                s64!(O_HIT, sr);
                sr += i as i64;
                (sr, if pushed { 2 } else { 1 })
            } else {
                (0, 0)
            }
        }
    }

    #[inline(never)]
    #[allow(unused_macros, unused_mut, unused_variables)]
    fn t_inv_rep<P: Pol>(&self, cr: i64, mut sr: i64, er: i64) -> (i64, u8) {
        let m = self.mem;
        let mut pushed = false;
        let cp = m.up(cr);
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r32(cr + $o)
                } else {
                    hr32(cp, $o)
                }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r64(cr + $o)
                } else {
                    hr64(cp, $o)
                }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w64(cr + $o, $v)
                } else {
                    hw64(cp, $o, $v)
                }
            };
        }
        let pu1 = m.names_p::<P>(g32!(O_U12));
        let p1 = m.pu_r64::<P>(pu1, O_HIT);
        let ln = m.pu_r32::<P>(pu1, O_MLEN);
        if ln == 0 {
            s64!(O_HIT, sr);
            (sr, if pushed { 2 } else { 1 })
        } else {
            let n = ln.max(0) as usize;
            // C uses char scratch[4000] and mallocs only above that
            let mut scratch = [std::mem::MaybeUninit::<u8>::uninit(); 4000];
            let mut heap_copy: Vec<u8> = Vec::new();
            let p3: &[u8] = if n <= 4000 {
                for (i, x) in scratch[..n].iter_mut().enumerate() {
                    *x = std::mem::MaybeUninit::new(m.rb(p1 + (n - 1 - i) as i64));
                }
                unsafe { std::slice::from_raw_parts(scratch.as_ptr() as *const u8, n) }
            } else {
                heap_copy.extend((0..n).map(|i| m.rb(p1 + (n - 1 - i) as i64)));
                &heap_copy[..]
            };
            let (ins, del, mis) = (g32!(O_U0), g32!(O_U4), g32!(O_U8));
            let i = self.loose_match(
                One::Bytes(p3),
                ln,
                sr,
                (er + 1 - sr) as i32,
                ins,
                del,
                mis,
                -1,
                false,
            );
            if i != 0 {
                let i = i - 1;
                s64!(O_HIT, sr);
                sr += i as i64;
                (sr, if pushed { 2 } else { 1 })
            } else {
                (0, 0)
            }
        }
    }

    #[inline(never)]
    #[allow(unused_macros, unused_mut, unused_variables)]
    fn t_sim<P: Pol>(&self, cr: i64, mut sr: i64, er: i64) -> (i64, u8) {
        let m = self.mem;
        let mut pushed = false;
        let cp = m.up(cr);
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r32(cr + $o)
                } else {
                    hr32(cp, $o)
                }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r64(cr + $o)
                } else {
                    hr64(cp, $o)
                }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w64(cr + $o, $v)
                } else {
                    hw64(cp, $o, $v)
                }
            };
        }
        let len = g32!(O_U12);
        let ins = g32!(O_U0);
        let mut last = er + 1 + ins as i64 - len as i64;
        // more inserts than letters: never start past the end of the
        // sequence (C reads outside the buffer there)
        if last > er + 1 {
            last = er + 1;
        }
        if last > sr && g32!(O_ANCH) != 0 {
            last = sr;
        }
        if last > self.start_lim && g32!(O_ANCH) == 0 {
            last = self.start_lim;
        }
        let mut found = false;
        while sr <= last {
            let code = g64!(O_U16);
            let (len, ins, del, mis) = (g32!(O_U12), g32!(O_U0), g32!(O_U4), g32!(O_U8));
            let i = self.loose_match(
                One::Mem(code),
                len,
                sr,
                (er + 1 - sr) as i32,
                ins,
                del,
                mis,
                -1,
                false,
            );
            if i != 0 {
                let i = i - 1;
                s64!(O_HIT, sr);
                if sr < last {
                    {
                        pushed = true;
                    }
                }
                sr += i as i64;
                found = true;
                break;
            }
            sr += 1;
        }
        if found {
            (sr, if pushed { 2 } else { 1 })
        } else {
            (0, 0)
        }
    }

    #[inline(never)]
    #[allow(unused_macros, unused_mut, unused_variables)]
    fn t_weight<P: Pol>(&self, cr: i64, mut sr: i64, er: i64) -> (i64, u8) {
        let m = self.mem;
        let mut pushed = false;
        let cp = m.up(cr);
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r32(cr + $o)
                } else {
                    hr32(cp, $o)
                }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED {
                    m.r64(cr + $o)
                } else {
                    hr64(cp, $o)
                }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED {
                    m.w64(cr + $o, $v)
                } else {
                    hw64(cp, $o, $v)
                }
            };
        }
        let wlen = g32!(O_U0);
        let mut last = er + 1 - wlen as i64;
        if last > sr && g32!(O_ANCH) != 0 {
            last = sr;
        }
        if last > self.start_lim && g32!(O_ANCH) == 0 {
            last = self.start_lim;
        }
        let pv1 = g64!(O_U8);
        let tupsz = g32!(O_U20);
        let cutoff = g32!(O_U16);
        let maxwt = g32!(O_U24);
        while sr <= last {
            let mut wval: i32 = 0;
            let mut p1 = sr;
            let mut i = wlen;
            if tupsz == 4 {
                let mut pv = pv1;
                while i != 0 {
                    let c = m.rb(p1) & 15;
                    p1 += 1;
                    let aw = || m.r32(pv);
                    let cw = || m.r32(pv + 4);
                    let gw = || m.r32(pv + 8);
                    let tw = || m.r32(pv + 12);
                    let add = match c {
                        0x1 => aw(),
                        0x2 => cw(),
                        0x4 => gw(),
                        0x8 => tw(),
                        0x3 => (aw() >> 1).wrapping_add(cw() >> 1),
                        0x5 => (aw() >> 1).wrapping_add(gw() >> 1),
                        0x9 => (aw() >> 1).wrapping_add(tw() >> 1),
                        0x6 => (cw() >> 1).wrapping_add(gw() >> 1),
                        0xA => (cw() >> 1).wrapping_add(tw() >> 1),
                        0xC => (gw() >> 1).wrapping_add(tw() >> 1),
                        0xE => (cw() / 3).wrapping_add(gw() / 3).wrapping_add(tw() / 3),
                        0xD => (aw() / 3).wrapping_add(gw() / 3).wrapping_add(tw() / 3),
                        0xB => (aw() / 3).wrapping_add(cw() / 3).wrapping_add(tw() / 3),
                        0x7 => (aw() / 3).wrapping_add(cw() / 3).wrapping_add(gw() / 3),
                        0xF => (aw() >> 2)
                            .wrapping_add(cw() >> 2)
                            .wrapping_add(gw() >> 2)
                            .wrapping_add(tw() >> 2),
                        _ => 0,
                    };
                    wval = wval.wrapping_add(add);
                    i = i.wrapping_sub(1);
                    pv += 16;
                }
            } else {
                let mut pv = pv1;
                while i != 0 {
                    let k: i64 = match m.rb(p1) {
                        b'A' => 0,
                        b'C' => 1,
                        b'D' => 2,
                        b'E' => 3,
                        b'F' => 4,
                        b'G' => 5,
                        b'H' => 6,
                        b'I' => 7,
                        b'K' => 8,
                        b'L' => 9,
                        b'M' => 10,
                        b'N' => 11,
                        b'P' => 12,
                        b'Q' => 13,
                        b'R' => 14,
                        b'S' => 15,
                        b'T' => 16,
                        b'V' => 17,
                        b'W' => 18,
                        b'Y' => 19,
                        _ => {
                            if tupsz > 20 {
                                20
                            } else {
                                -1
                            }
                        }
                    };
                    p1 += 1;
                    if k >= 0 {
                        wval = wval.wrapping_add(m.r32(pv + 4 * k));
                    }
                    i = i.wrapping_sub(1);
                    pv += 4 * tupsz as i64;
                }
            }
            if wval > cutoff && wval < maxwt {
                break;
            }
            sr += 1;
        }
        if sr > last {
            (0, 0)
        } else {
            s64!(O_HIT, sr);
            if sr < last {
                {
                    pushed = true;
                }
            }
            sr += g32!(O_U0) as i64;
            (sr, if pushed { 2 } else { 1 })
        }
    }

    pub fn first_match(&mut self, len: i32, hits: &mut Vec<i64>) -> i32 {
        self.start_lim = i64::MAX;
        self.first_match_at(len, 0, hits)
    }

    /// Search one piece of the sequence: only hits that start at offsets
    /// `from..=last` (they may end past `last`).  next_match and
    /// cont_match then go on inside the same piece.  Only for patterns
    /// where `can_split()` holds.
    pub fn first_match_in(&mut self, len: i32, from: i64, last: i64, hits: &mut Vec<i64>) -> i32 {
        self.start_lim = self.mem.cd as i64 + last;
        self.first_match_at(len, from, hits)
    }

    fn first_match_at(&mut self, len: i32, from: i64, hits: &mut Vec<i64>) -> i32 {
        let m = self.mem;
        let start = m.cd as i64 + from;
        m.w64(self.mem.sa(A_START_SRCH), start);
        self.end_srch = m.cd as i64 + (len as i64 - 1);
        self.br1 = 0;
        let pu = m.r64(self.mem.sa(A_AD_PU_S));
        let i = self.pattern_match(pu, start, self.end_srch, hits, true);
        let v = hits[i as usize];
        m.w64(self.mem.sa(A_PAST_LAST), v);
        i
    }

    pub fn next_match(&mut self, hits: &mut Vec<i64>) -> i32 {
        let m = self.mem;
        let pu = m.r64(self.mem.sa(A_AD_PU_S));
        let (s, e) = (self.start_srch(), self.end_srch);
        let i = self.pattern_match(pu, s, e, hits, false);
        let v = hits[i as usize];
        m.w64(self.mem.sa(A_PAST_LAST), v);
        i
    }

    pub fn cont_match(&mut self, hits: &mut Vec<i64>) -> i32 {
        let m = self.mem;
        let past_last1 = self.past_last();
        let mut i;
        loop {
            i = self.next_match(hits);
            if !(i > 0 && hits[0] < past_last1) {
                break;
            }
        }
        let v = hits[i as usize];
        m.w64(self.mem.sa(A_PAST_LAST), v);
        i
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_with_threads_is_the_same() {
        // coded bases (1 2 4 8) with runs of other codes (N etc.)
        let mut x: u64 = 12345;
        let mut seq = Vec::new();
        for _ in 0..(3 << 20) + 777 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            seq.push(match x % 50 {
                0 => 15,
                1 => 0,
                r => [1, 2, 4, 8][(r % 4) as usize],
            });
        }
        let one = build_index(&seq, None);
        for t in [1, 2, 3, 7, 64] {
            let par = build_index_par(&seq, t);
            assert!(
                one.starts == par.starts && one.pos == par.pos,
                "threads {t}"
            );
        }
        // a part boundary inside a run of plain bases and next to an N
        let mut s2 = vec![1u8; 2 << 20];
        s2[1 << 20] = 15;
        s2[(1 << 20) + 2] = 15;
        let one = build_index(&s2, None);
        let par = build_index_par(&s2, 2);
        assert!(one.starts == par.starts && one.pos == par.pos);
    }

    /// The search of loose_match written as plain recursion: at each step
    /// a matching character is always taken; otherwise the choices are
    /// tried in the order mismatch, insert, delete (insert/delete only as
    /// the C code offers them).  `pending` tells if an earlier choice is
    /// still open (C: `nxtent != 0`), because the C loop stops at the end
    /// of the sequence when nothing is open.  Returns the number of
    /// sequence characters used, or None.
    #[allow(clippy::too_many_arguments)]
    fn reference(
        one: &[u8],
        two: &[u8],
        mut i: usize,
        mut j: usize,
        mis: i32,
        ins: i32,
        del: i32,
        pending: bool,
    ) -> Option<usize> {
        let known = |c: u8| KNOWN_CHAR[(c & 15) as usize] != 0;
        let matches = |c1: u8, c2: u8| {
            let a = c1 & 15;
            known(a) && (a & (c2 & 15)) == a
        };
        loop {
            let one_len = one.len() - i;
            let two_len = two.len() - j;
            if two_len == 0 && !pending {
                return None;
            }
            if two_len > 0 && one_len > 0 && matches(two[j], one[i]) {
                i += 1;
                j += 1;
                if one.len() == i {
                    return Some(j);
                }
                continue;
            }
            // (choice, is it the last choice of this point)
            let mut choices: Vec<char> = Vec::new();
            if mis > 0 && one_len >= 1 && two_len >= 1 {
                choices.push('m');
                if ins > 0 {
                    choices.push('i');
                    if del > 0 {
                        choices.push('d');
                    }
                } else if del > 0 {
                    choices.push('d');
                }
            } else if ins > 0 && one_len >= 1 {
                choices.push('i');
                if del > 0 && two_len >= 1 {
                    choices.push('d');
                }
            } else if del > 0 && two_len >= 1 {
                // no other choice: continue in this loop
                j += 1;
                let del2 = del - 1;
                if one.len() == i {
                    return Some(j);
                }
                return reference(one, two, i, j, mis, ins, del2, pending);
            } else {
                return None;
            }
            let n = choices.len();
            for (k, c) in choices.into_iter().enumerate() {
                let p = pending || k + 1 < n;
                let r = match c {
                    'm' => {
                        if one.len() == i + 1 {
                            return Some(j + 1);
                        }
                        reference(one, two, i + 1, j + 1, mis - 1, ins, del, p)
                    }
                    'i' => {
                        if one.len() == i + 1 {
                            return Some(j);
                        }
                        reference(one, two, i + 1, j, mis, ins - 1, del, p)
                    }
                    _ => {
                        if one.len() == i {
                            return Some(j + 1);
                        }
                        reference(one, two, i, j + 1, mis, ins, del - 1, p)
                    }
                };
                if r.is_some() {
                    return r;
                }
            }
            return None;
        }
    }

    fn parse(e: &mut Engine, pat: &str) -> bool {
        e.parse_cmd(pat.as_bytes(), DNA) != 0
    }

    #[test]
    fn earlier_state_detection() {
        let cases = [
            ("p1=3...3 2...4 ~p1", false),
            ("p1=2...2 (~p1 | p1) p2=1...1 <p2", false),
            ("p1=3...3 p2=1...2 length(p1+p2) < 5", false),
            ("~p1 p1=3...3", true),
            ("(p1=3...3 | AC) ~p1", true),
            ("(p1=3...3 AA | p2=2...2) length(p1) < 5", true),
            ("length(p1+p2) < 5 p1=2...2 p2=2...2", true),
            ("(p1=2...2 | p1x) 1...2", false),
        ];
        for (pat, want) in cases {
            let mut e = Engine::new();
            if !parse(&mut e, pat) {
                continue;
            }
            assert_eq!(e.uses_earlier_state(), want, "{pat}");
        }
    }

    /// Engines in several threads at once give the same hits as one engine.
    #[test]
    fn engines_in_parallel_threads() {
        let pat = "p1=4...6 2...6 ~p1[1,0,0] 0...3 p2=3...3 p2";
        let seqs: Vec<Vec<u8>> = (0..8u64)
            .map(|k| {
                let mut x = 1 + k;
                (0..20_000)
                    .map(|_| {
                        x ^= x << 13;
                        x ^= x >> 7;
                        x ^= x << 17;
                        b"acgt"[(x % 4) as usize]
                    })
                    .collect()
            })
            .collect();
        let scan = |seq: &[u8]| -> Vec<Vec<i64>> {
            let mut e = Engine::new();
            assert!(parse(&mut e, pat));
            let mut data = Buf::new();
            data.store(seq);
            e.comp_data(&data);
            let base = e.cdata_base();
            let mut hits = vec![0i64; 2000];
            let mut out = Vec::new();
            let mut i = e.first_match(seq.len() as i32, &mut hits);
            while i > 0 {
                out.push(hits[..=i as usize].iter().map(|h| h - base).collect());
                i = e.cont_match(&mut hits);
            }
            out
        };
        let one: Vec<_> = seqs.iter().map(|s| scan(s)).collect();
        assert!(one.iter().all(|h| !h.is_empty()));
        let many: Vec<_> = std::thread::scope(|sc| {
            let hs: Vec<_> = seqs.iter().map(|s| sc.spawn(move || scan(s))).collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        assert_eq!(one, many);
    }

    #[test]
    fn loose_match_explores_every_choice() {
        let mut e = Engine::new();
        assert!(e.parse_cmd(b"A", DNA) != 0);
        let mut seed: u64 = 12345;
        let mut rnd = |n: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % n
        };
        let mut data = Buf::new();
        for case in 0..200_000 {
            let one_len = 1 + rnd(6) as usize;
            let two_len = rnd(9) as usize;
            let bases = b"acgtn";
            let one: Vec<u8> = (0..one_len)
                .map(|_| e.p2c(bases[rnd(4) as usize]))
                .collect();
            let two_chars: Vec<u8> = (0..two_len)
                .map(|_| {
                    let k = if rnd(8) == 0 { 5 } else { 4 };
                    bases[rnd(k) as usize]
                })
                .collect();
            data.store(&two_chars);
            e.comp_data(&data);
            let two: Vec<u8> = two_chars.iter().map(|&c| e.p2c(c)).collect();
            let (mis, ins, del) = (rnd(3) as i32, rnd(3) as i32, rnd(3) as i32);
            if ins == 0 && del == 0 {
                continue;
            }
            let got = e.loose_match(
                One::Bytes(&one),
                one_len as i32,
                e.cdata_base(),
                two_len as i32,
                ins,
                del,
                mis,
                -1,
                false,
            );
            let want =
                reference(&one, &two, 0, 0, mis, ins, del, false).map_or(0, |j| j as i32 + 1);
            assert_eq!(
                got, want,
                "case {case}: one {one:?} two {two:?} mis {mis} ins {ins} del {del}"
            );
        }
    }
}
