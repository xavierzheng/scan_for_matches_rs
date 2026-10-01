//! Port of ggpunit.c: the pattern parser and the backtracking matcher.
//!
//! The code follows the C original statement by statement.
//!
//! The C program keeps its pattern in fixed-size static arrays (100 pattern
//! units, 600 code bytes, 10500 weights) and its state in globals.  When a
//! pattern is too big for these arrays the C code silently writes past them
//! into the neighbouring arrays, which changes its results.  To give the same
//! results in those cases too, all of that state lives here in one byte
//! array laid out exactly like the static data of the reference build
//! (`gcc -std=gnu89 -O2` with Apple clang on macOS).  Pointers are stored in
//! it as 8-byte addresses, so overlapping writes behave as they do in C.
//!
//! Address map (same numbers as the reference binary):
//!   0x1_0000_c000 .. 0x1_0002_8000   __DATA segment (globals, pu_s, cv, iv)
//!   0x1_0002_8000 .. 0x1_0002_c000   __LINKEDIT (read only)
//!   CDATA_BASE ..                    the coded sequence buffer (malloc)
//!   HEAP_BASE ..                     malloc(16) blocks for rule sets

use crate::sys::{abort, sigbus, segv};

pub const PEPTIDE: i32 = 1;
pub const DNA: i32 = 2;

#[allow(dead_code)]
pub const MAX_SEQ_LEN: usize = 250_000_000;
/// malloc(MAX_SEQ_LEN+1) is rounded up to whole 16 KiB pages; bytes up to
/// this size are addressable, the next one faults.
pub const ALLOC_LEN: i64 = 250_003_456;
/// Bytes just before a buffer are mapped (they read as zero in the
/// reference build); the second buffer starts this far after the first.
const LOW_SLACK: i64 = 1_671_168;
const MAX_NAMES: i32 = 50;

// ---- address map ---------------------------------------------------------
const S_BASE: i64 = 0x1_0000_c000;
const S_END: i64 = 0x1_0002_8000;
const LINKEDIT_END: i64 = 0x1_0002_c000;
const S_LEN: usize = (S_END - S_BASE) as usize;

const A_KNOWN_CHAR: i64 = 0x1_0000_c000;
const A_KNOWN_CHAR_INDEX: i64 = 0x1_0000_c010;
const A_INITIALIZED: i64 = 0x1_0000_c100;
const A_AD_PU_S: i64 = 0x1_0000_c108;
const A_CODE_TO_PUNIT: i64 = 0x1_0000_c200;
const A_NAMES: i64 = 0x1_0000_c310;
const A_PAST_LAST: i64 = 0x1_0000_c4a0;
const A_SEQ_TYPE: i64 = 0x1_0000_c4a8;
const A_P2C: i64 = 0x1_0000_c500;
const A_RULE_SETS: i64 = 0x1_0000_c600;
const A_START_SRCH: i64 = 0x1_0000_c790;
const A_PU_S: i64 = 0x1_0000_c798;
const A_CV: i64 = 0x1_0001_2eb8;
const A_IV: i64 = 0x1_0001_3110;

pub const CDATA_BASE: i64 = 0x3_1040_8000; // low 24 bits as in the reference run
const HEAP_BASE: i64 = 0x6000_0000_0000;
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

/// Fixed stack arrays of the C code; overrunning them trips the stack
/// protector (abort) in the reference build.
const MAX_PUNITS: usize = 100; // revhits[] in pattern_match
const MAX_CODES: i32 = 100; // result[] in loose_match
const LOOSE_STACK_SAFE: i32 = 102; // stack[102] reaches the canary

/// A byte buffer that the C code mallocs once (MAX_SEQ_LEN+1 bytes, which
/// the allocator rounds up to ALLOC_LEN) and reuses for every sequence:
/// bytes past the current sequence keep what an earlier, longer sequence
/// left there.  The zeroed allocation only uses memory for touched pages.
pub struct Buf {
    pub v: Vec<u8>,
}

impl Buf {
    pub fn new() -> Buf {
        Buf { v: vec![0u8; ALLOC_LEN as usize] }
    }

    /// None when the C program would fault on this read.
    #[inline(always)]
    pub fn try_get(&self, off: i64) -> Option<u8> {
        if (off as u64) < ALLOC_LEN as u64 {
            Some(self.v[off as usize])
        } else if (-LOW_SLACK..0).contains(&off) {
            Some(0)
        } else {
            None
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
fn pu_ok(p: i64) -> bool {
    (p.wrapping_sub(S_BASE) as u64) <= (S_LEN as i64 - SZ) as u64
}

/// How the fast matcher stopped.
enum Out {
    Done(i32),
    /// a unit pointer outside the static block: continue in checked mode
    Bail { back: bool, cr: i64, sr: i64, br: i64 },
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
    let o = a.wrapping_sub(CDATA_BASE) as u64;
    if o < ALLOC_LEN as u64 {
        return unsafe { *sm.cd.add(o as usize) };
    }
    let o = a.wrapping_sub(S_BASE) as u64;
    if o < S_LEN as u64 {
        return unsafe { *sm.s.add(o as usize) };
    }
    let o = a.wrapping_sub(CDATA_BASE);
    if (-LOW_SLACK..0).contains(&o) {
        return 0;
    }
    if (S_END..LINKEDIT_END).contains(&a) {
        return 0;
    }
    let o = a.wrapping_sub(HEAP_BASE);
    if o >= 0 && (o as usize) < sm.heap_len {
        return unsafe { *sm.heap.add(o as usize) };
    }
    segv()
}

#[inline(always)]
fn wb_full(sm: &SlowMem, a: i64, v: u8) {
    let o = a.wrapping_sub(S_BASE) as u64;
    if o < S_LEN as u64 {
        unsafe { *sm.s.add(o as usize) = v };
        return;
    }
    if (S_END..LINKEDIT_END).contains(&a) {
        sigbus(); // read-only segment
    }
    let o = a.wrapping_sub(HEAP_BASE);
    if o >= 0 && (o as usize) < sm.heap_len {
        unsafe { *sm.heap.add(o as usize) = v };
        return;
    }
    let o = a.wrapping_sub(CDATA_BASE) as u64;
    if o < ALLOC_LEN as u64 {
        unsafe { *sm.cd.add(o as usize) = v };
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
    i32::from_le_bytes([rb_full(sm, a), rb_full(sm, a + 1), rb_full(sm, a + 2), rb_full(sm, a + 3)])
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
    #[inline(always)]
    fn rb(self, a: i64) -> u8 {
        let o = a.wrapping_sub(CDATA_BASE) as u64;
        if o < ALLOC_LEN as u64 {
            return unsafe { *self.cd.add(o as usize) };
        }
        let o = a.wrapping_sub(S_BASE) as u64;
        if o < S_LEN as u64 {
            return unsafe { *self.s.add(o as usize) };
        }
        rb_slow(self.sm, a)
    }

    #[inline(always)]
    fn wb(self, a: i64, v: u8) {
        let o = a.wrapping_sub(S_BASE) as u64;
        if o < S_LEN as u64 {
            unsafe { *self.s.add(o as usize) = v };
            return;
        }
        wb_slow(self.sm, a, v)
    }

    #[inline(always)]
    fn r32(self, a: i64) -> i32 {
        let o = a.wrapping_sub(S_BASE) as u64;
        if o + 4 <= S_LEN as u64 {
            return unsafe { (self.s.add(o as usize) as *const i32).read_unaligned() };
        }
        r32_slow(self.sm, a)
    }

    #[inline(always)]
    fn r64(self, a: i64) -> i64 {
        let o = a.wrapping_sub(S_BASE) as u64;
        if o + 8 <= S_LEN as u64 {
            return unsafe { (self.s.add(o as usize) as *const i64).read_unaligned() };
        }
        r64_slow(self.sm, a)
    }

    #[inline(always)]
    fn w32(self, a: i64, v: i32) {
        let o = a.wrapping_sub(S_BASE) as u64;
        if o + 4 <= S_LEN as u64 {
            unsafe { (self.s.add(o as usize) as *mut i32).write_unaligned(v) };
            return;
        }
        wn_slow(self.sm, a, v as u32 as u64, 4);
    }

    #[inline(always)]
    fn w64(self, a: i64, v: i64) {
        let o = a.wrapping_sub(S_BASE) as u64;
        if o + 8 <= S_LEN as u64 {
            unsafe { (self.s.add(o as usize) as *mut i64).write_unaligned(v) };
            return;
        }
        wn_slow(self.sm, a, v as u64, 8);
    }

    #[inline(always)]
    fn pr32<P: Pol>(self, a: i64) -> i32 {
        if P::CHECKED {
            self.r32(a)
        } else {
            unsafe { (self.s.add(a.wrapping_sub(S_BASE) as usize) as *const i32).read_unaligned() }
        }
    }

    #[inline(always)]
    fn pr64<P: Pol>(self, a: i64) -> i64 {
        if P::CHECKED {
            self.r64(a)
        } else {
            unsafe { (self.s.add(a.wrapping_sub(S_BASE) as usize) as *const i64).read_unaligned() }
        }
    }

    /// Host pointer to the unit at `pu` (only valid when pu_ok(pu)).
    #[inline(always)]
    fn up(self, pu: i64) -> *mut u8 {
        unsafe { self.s.add(pu.wrapping_sub(S_BASE) as usize) }
    }

    /// names[i] (static slot unless `i` is out of range)
    #[inline(always)]
    fn names_p<P: Pol>(self, i: i32) -> i64 {
        if !P::CHECKED && (0..=MAX_NAMES).contains(&i) {
            self.pr64::<Fast>(A_NAMES + 8 * i as i64)
        } else {
            self.names(i)
        }
    }

    /// a field of a unit reached through a pointer that has not been checked
    #[inline(always)]
    fn pu_r64<P: Pol>(self, pu: i64, off: i64) -> i64 {
        if !P::CHECKED && pu_ok(pu) { self.pr64::<Fast>(pu + off) } else { self.r64(pu + off) }
    }

    #[inline(always)]
    fn pu_r32<P: Pol>(self, pu: i64, off: i64) -> i32 {
        if !P::CHECKED && pu_ok(pu) { self.pr32::<Fast>(pu + off) } else { self.r32(pu + off) }
    }

    /// Are the `n` bytes from `a` inside the coded sequence buffer?
    #[inline(always)]
    fn cd_ok(a: i64, n: i64) -> bool {
        let o = a.wrapping_sub(CDATA_BASE);
        o >= 0 && n >= 0 && o.wrapping_add(n) <= ALLOC_LEN
    }

    #[inline(always)]
    fn names(self, i: i32) -> i64 {
        self.r64(A_NAMES + 8 * i as i64)
    }

    /// KnownChar(C)
    #[inline(always)]
    fn known_char(self, c: u8) -> bool {
        if self.pep { true } else { KNOWN_CHAR[c as usize] != 0 }
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
            let base = self.r64(A_RULE_SETS + 8 * rule_set as i64);
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
        Some((n, p)) if n >= 0 && n <= MAX_NAMES => Some((n, p)),
        _ => None,
    }
}

fn name_assgn(l: &[u8], p: usize) -> Option<(i32, usize)> {
    let (n, p) = name_id(l, p)?;
    let p = ws(l, p);
    if at(l, p) == b'=' { Some((n, p + 1)) } else { None }
}

fn rule_id(l: &[u8], p: usize) -> Option<(i32, usize)> {
    if at(l, p) != b'r' {
        return None;
    }
    match num(l, p + 1) {
        Some((n, p)) if n >= 0 && n <= MAX_NAMES => Some((n, p)),
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
    Some((mis, ins, del, p + 1))
}

fn elipses(l: &[u8], p: usize) -> Option<usize> {
    if at(l, p) == b'.' && at(l, p + 1) == b'.' && at(l, p + 2) == b'.' { Some(p + 3) } else { None }
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
        match bond(l, p + 1, rs) {
            Some(q) => p = q,
            None => return None,
        }
    }
}

/// C `compl()` from scan_for_matches.c (note: 's' maps to 'S').
pub fn compl(c: u8) -> u8 {
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
        b's' => b'S',
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
    /// static data of the C program (accessed through `mem`)
    _s: Vec<u8>,
    /// malloc'd rule-set blocks (accessed through `mem`)
    _heap: Vec<u8>,
    slow: Box<SlowMem>,
    /// the coded sequence, cdata (accessed through `mem`)
    _cdata: Buf,
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
    /// punit_to_code for indexes 0..127 (never overwritten)
    p2c_lo: [u8; 128],
}

impl Engine {
    pub fn new() -> Engine {
        let mut s = vec![0u8; S_LEN];
        let mut heap = vec![0u8; HEAP_CAP];
        let mut cdata = Buf::new();
        let slow = Box::new(SlowMem { s: s.as_mut_ptr(), cd: cdata.v.as_mut_ptr(), heap: heap.as_mut_ptr(), heap_len: 0 });
        let mem = M { s: slow.s, cd: slow.cd, sm: &*slow as *const SlowMem, pep: false };
        let mut e = Engine {
            _s: s,
            _heap: heap,
            slow,
            _cdata: cdata,
            mem,
            seq_type: 0,
            pup: A_PU_S,
            cvp: A_CV,
            ivp: A_IV,
            br1: 0,
            end_srch: 0,
            revhits: Vec::with_capacity(256),
            p2c_lo: [0; 128],
        };
        for (i, v) in KNOWN_CHAR.iter().enumerate() {
            e.wb(A_KNOWN_CHAR + i as i64, *v);
        }
        for (i, v) in KNOWN_CHAR_INDEX.iter().enumerate() {
            e.wb(A_KNOWN_CHAR_INDEX + i as i64, *v as u8);
        }
        e
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
        let a = HEAP_BASE + self.slow.heap_len as i64;
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
        self.r64(A_PAST_LAST)
    }

    #[inline(always)]
    fn start_srch(&self) -> i64 {
        self.r64(A_START_SRCH)
    }

    /// punit_to_code[c] where C subscripts with a signed char.
    #[inline(always)]
    fn p2c(&self, c: u8) -> u8 {
        if c < 0x80 {
            self.p2c_lo[c as usize]
        } else {
            self.rb(A_P2C + (c as i8) as i64)
        }
    }

    fn build_conversion_tables(&mut self) {
        for the_char in 0..256i64 {
            let lc = if the_char < 128 { (the_char as u8).to_ascii_lowercase() } else { the_char as u8 };
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
            self.wb(A_P2C + the_char, v);
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
                self.wb(A_CODE_TO_PUNIT + the_char, c);
            }
        }
        self.w32(A_INITIALIZED, 1);
    }

    /// comp_data(data, cdata): translate characters to codes, stop at NUL.
    pub fn comp_data(&mut self, data: &Buf) {
        let mut k = 0usize;
        loop {
            let c = data.v[k];
            if c == 0 {
                break;
            }
            let code = self.p2c(c);
            if k as i64 >= ALLOC_LEN - 1 {
                segv();
            }
            unsafe { *self.mem.cd.add(k) = code };
            k += 1;
        }
        unsafe { *self.mem.cd.add(k) = 0 };
    }

    /// strcpy(cdata, data)
    pub fn copy_data(&mut self, data: &Buf) {
        let mut k = 0usize;
        loop {
            let c = data.v[k];
            unsafe { *self.mem.cd.add(k) = c };
            if c == 0 {
                break;
            }
            k += 1;
        }
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
        if self.seq_type == DNA {
            if let Some((mis, ins, del, p1)) = misinsdel(l, p) {
                return Some((n, mis, ins, del, p1));
            }
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
        if n >= 0 && n <= MAX_NAMES { Some(q) } else { None }
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
        let slot = A_RULE_SETS + 8 * n as i64;
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
        if let Some((i, p1)) = name_assgn(l, p) {
            let slot = A_NAMES + 8 * i as i64;
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
                self.max_mat(p1, depth + 1).wrapping_add(self.r32(pu + O_U0))
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
        self.w32(A_SEQ_TYPE, seq_type);
        if self.r32(A_INITIALIZED) == 0 {
            self.build_conversion_tables();
        }
        for i in 0..MAX_NAMES as i64 {
            self.w64(A_NAMES + 8 * i, 0);
        }
        self.ivp = A_IV;
        self.cvp = A_CV;
        self.pup = A_PU_S;
        self.w64(A_AD_PU_S, A_PU_S);
        if line.len() >= 1_000_000 {
            self.w64(A_AD_PU_S, 0);
            return 0;
        }
        match self.parser(line) {
            None => {
                self.w64(A_AD_PU_S, 0);
                0
            }
            Some(_) => {
                self.set_anchors(A_PU_S, 0);
                self.w64(A_AD_PU_S, A_PU_S);
                self.max_mats(A_PU_S)
            }
        }
    }

    // ------------------------------------------------------------------
    // matcher



    /// ExMatches(RuleSet, C1, C2)


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
                if M::cd_ok(*base, n as i64) {
                    Some(unsafe { std::slice::from_raw_parts(m.cd.add(base.wrapping_sub(CDATA_BASE) as usize), n) })
                } else {
                    let o = base.wrapping_sub(S_BASE);
                    if o >= 0 && (o as usize) + n <= S_LEN {
                        Some(unsafe { std::slice::from_raw_parts(m.s.add(o as usize), n) })
                    } else {
                        None
                    }
                }
            }
        }
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
        // like the C locals: not initialised, only written parts are read
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
                // result[MAX_CODES] overflows into the stack guard
                abort();
            }
            let n = one_len.max(0) as usize;
            let mut len = one_len;
            let mut k: i64 = 0;
            while len > 0 {
                let i = m.one_get(&one_src, k) as u32;
                k += 1;
                len -= 1;
                result[len as usize] = std::mem::MaybeUninit::new(if rule_set == -1 {
                    ((i >> 4) & 15) as u8
                } else {
                    ((KNOWN_CHAR_INDEX[(i & 15) as usize] as i32) << 2) as u8
                });
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
            if one_len >= 1 && M::cd_ok(two_start, one_len as i64) {
                let n = one_len as usize;
                if let Some(ob) = self.one_slice(&one_src, n) {
                    let tb = unsafe { std::slice::from_raw_parts(m.cd.add(two_start.wrapping_sub(CDATA_BASE) as usize), n) };
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

        let mut stack = [std::mem::MaybeUninit::<StackEnt>::uninit(); LOOSE_STACK_SAFE as usize];
        let mut nxtent: i32 = 0;
        macro_rules! ret {
            ($v:expr) => {{
                return $v;
            }};
        }
        macro_rules! push {
            ($n:expr) => {{
                let k = nxtent as usize;
                if nxtent >= LOOSE_STACK_SAFE {
                    // stack[] has run into the stack guard: the C function
                    // aborts when it returns (it has no other side effects)
                    abort();
                }
                stack[k] = std::mem::MaybeUninit::new(StackEnt {
                    p1: one,
                    p2: two,
                    n1: one_len,
                    n2: two_len,
                    mis: max_mis,
                    ins: max_ins,
                    del: max_del,
                    next_choice: $n,
                });
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
                        nxtent -= 1;
                        // entries below nxtent were all written by push
                        let e = unsafe { stack[nxtent as usize].assume_init() };
                        one = e.p1;
                        two = e.p2;
                        one_len = e.n1;
                        two_len = e.n2;
                        max_mis = e.mis;
                        max_ins = e.ins;
                        max_del = e.del;
                        if e.next_choice == 1 {
                            if max_del != 0 {
                                unsafe { (*stack[nxtent as usize].as_mut_ptr()).next_choice = 2 };
                            }
                            lbl = Lm::Ins;
                        } else {
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
        if !P::CHECKED {
            if let Some(n) = Self::next_punit_fast(m, pu) {
                return n;
            }
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
            if !pu_ok(pu2) {
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
            if !pu_ok(pu) {
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
                pu = if hr32(up, O_U24) == 1 { hr64(up, O_U0) } else { hr64(up, O_U8) };
                loop {
                    if !pu_ok(pu) {
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
                pu = if m.r32(pu + O_U24) == 1 { m.r64(pu + O_U0) } else { m.r64(pu + O_U8) };
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

    fn pattern_match(&mut self, pu: i64, start: i64, end: i64, hits: &mut Vec<i64>, first: bool) -> i32 {
        let br = self.br1;
        match self.pm::<Fast>(start, end, hits, !first, pu, start, br) {
            Out::Done(n) => n,
            Out::Bail { back, cr, sr, br } => match self.pm::<Checked>(start, end, hits, back, cr, sr, br) {
                Out::Done(n) => n,
                Out::Bail { .. } => unreachable!(),
            },
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
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s32 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w32(cr + $o, $v) } else { hw32(cp, $o, $v) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
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
                        if !pu_ok(cr) {
                            return Out::Bail { back: false, cr, sr, br };
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
                            if (sr <= er && g32!(O_U4) != 0) || (g32!(O_ANCH) == 0 && sr <= er) {
                                push_br!();
                                let mn = g32!(O_U0);
                                s32!(O_U8, mn.wrapping_add(1));
                            }
                            success!('tryl, 'main);
                        } else {
                            break 'tryl;
                        }
                    }
                    EXACT_PUNIT => {
                        let (n, c) = self.t_exact::<P>(cr, sr, er);
                        if c == 0 {
                            break 'tryl;
                        }
                        if c == 2 {
                            push_br!();
                        }
                        sr = n;
                        success!('tryl, 'main);
                    }
                    COMPL_PUNIT => {
                        let pu1 = m.names_p::<P>(g32!(O_U12));
                        let pp = if !P::CHECKED && pu_ok(pu1) { m.up(pu1) } else { std::ptr::null_mut() };
                        let pu1_hit = |m: M| if pp.is_null() { m.r64(pu1 + O_HIT) } else { hr64(pp, O_HIT) };
                        let p1 = pu1_hit(m);
                        let mut ln = if pp.is_null() { m.r32(pu1 + O_MLEN) } else { hr32(pp, O_MLEN) };
                        let rs = g32!(O_U16);
                        let (ins, del, mis) = (g32!(O_U0), g32!(O_U4), g32!(O_U8));
                        if rs != -1 || ins != 0 || del != 0 || mis != 0 {
                            let i = self.loose_match(One::Mem(p1), ln, sr, (er + 1 - sr) as i32, ins, del, mis, rs, true);
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
                            if !P::CHECKED && !m.pep && ln > 0 && M::cd_ok(q - (ln as i64 - 1), ln as i64) && M::cd_ok(sr, ln as i64) {
                                // all reads are inside the sequence buffer
                                let qp = m.cd.wrapping_add(q.wrapping_sub(CDATA_BASE) as usize);
                                let sp = m.cd.wrapping_add(sr.wrapping_sub(CDATA_BASE) as usize);
                                let n = ln as usize;
                                let mut k = 0usize;
                                while k < n {
                                    let c = unsafe { *qp.wrapping_sub(k) };
                                    if KNOWN_CHAR[(c & 15) as usize] == 0 || ((c >> 4) & 15) != (unsafe { *sp.add(k) } & 15) {
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
                        if !pu_ok(br) {
                            return Out::Bail { back: true, cr, sr, br };
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
                            if nx <= min.wrapping_add(width) && sr + nx as i64 - 1 <= er {
                                s32!(O_U8, nx.wrapping_add(1));
                                sr += nx as i64;
                                br = cr;
                                success!('main, 'main);
                            } else {
                                let h = g64!(O_HIT) + 1;
                                s64!(O_HIT, h);
                                if h + g32!(O_U0) as i64 - 1 <= er && g32!(O_ANCH) == 0 {
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
                    let n = revhits.len();
                    if hits.len() < n + 1 {
                        hits.resize(n + 1, 0);
                    }
                    hits[..n].copy_from_slice(&revhits[..n]);
                    hits[..n].reverse();
                    hits[n] = sr;
                    self.revhits = revhits;
                    self.br1 = br;
                    if n > MAX_PUNITS {
                        // revhits[MAX_PUNITS] overflows into the stack guard
                        abort();
                    }
                    return Out::Done(n as i32);
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
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
            };
        }
        let mut last = er;
        if last > sr && g32!(O_ANCH) != 0 {
            last = sr;
        }
        let cm = g64!(O_U0);
        while sr <= last {
            let i = m.rb(sr) as i8 as i32;
            if i >= b'A' as i32 && i <= b'Z' as i32 && ((1i32 << (i - b'A' as i32)) as i64 & cm) != 0 {
                break;
            }
            sr += 1;
        }
        if sr > last {
            return (0, 0);
        } else {
            s64!(O_HIT, sr);
            if sr < last {
                {
            pushed = true;
        }
            }
            sr += 1;
            return (sr, if pushed { 2 } else { 1 });
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
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
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
            return (sr, if pushed { 2 } else { 1 });
        } else {
            return (0, 0);
        }
    }

    #[inline(never)]
    #[allow(unused_macros, unused_mut, unused_variables)]
    fn t_exact<P: Pol>(&self, cr: i64, mut sr: i64, er: i64) -> (i64, u8) {
        let m = self.mem;
        let mut pushed = false;
        let cp = m.up(cr);
        macro_rules! g32 {
            ($o:expr) => {
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
            };
        }
        let len = g32!(O_U0);
        let mut last = er + 1 - len as i64;
        if last > sr && g32!(O_ANCH) != 0 {
            last = sr;
        }
        let p1 = g64!(O_U8);
        let ln = len.wrapping_sub(1);
        while sr <= last {
            if m.matches(m.rb(sr), m.rb(p1)) {
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
            return (0, 0);
        } else {
            s64!(O_HIT, sr);
            if sr < last {
                {
            pushed = true;
        }
            }
            sr += g32!(O_U0) as i64;
            return (sr, if pushed { 2 } else { 1 });
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
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
            };
        }
        let pu1 = m.names_p::<P>(g32!(O_U12));
        let p1 = m.pu_r64::<P>(pu1, O_HIT);
        let ln = m.pu_r32::<P>(pu1, O_MLEN);
        if ln == 0 {
            s64!(O_HIT, sr);
            return (sr, if pushed { 2 } else { 1 });
        } else {
            let (ins, del, mis) = (g32!(O_U0), g32!(O_U4), g32!(O_U8));
            let i = self.loose_match(One::Mem(p1), ln, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
            if i != 0 {
                let i = i - 1;
                s64!(O_HIT, sr);
                sr += i as i64;
                return (sr, if pushed { 2 } else { 1 });
            } else {
                return (0, 0);
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
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
            };
        }
        let pu1 = m.names_p::<P>(g32!(O_U12));
        let p1 = m.pu_r64::<P>(pu1, O_HIT);
        let ln = m.pu_r32::<P>(pu1, O_MLEN);
        if ln == 0 {
            s64!(O_HIT, sr);
            return (sr, if pushed { 2 } else { 1 });
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
            let i = self.loose_match(One::Bytes(p3), ln, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
            if i != 0 {
                let i = i - 1;
                s64!(O_HIT, sr);
                sr += i as i64;
                return (sr, if pushed { 2 } else { 1 });
            } else {
                return (0, 0);
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
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
            };
        }
        let len = g32!(O_U12);
        let ins = g32!(O_U0);
        let mut last = er + 1 + ins as i64 - len as i64;
        if last > sr && g32!(O_ANCH) != 0 {
            last = sr;
        }
        let mut found = false;
        while sr <= last {
            let code = g64!(O_U16);
            let (len, ins, del, mis) =
                (g32!(O_U12), g32!(O_U0), g32!(O_U4), g32!(O_U8));
            let i = self.loose_match(One::Mem(code), len, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
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
            return (sr, if pushed { 2 } else { 1 });
        } else {
            return (0, 0);
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
                if P::CHECKED { m.r32(cr + $o) } else { hr32(cp, $o) }
            };
        }
        macro_rules! g64 {
            ($o:expr) => {
                if P::CHECKED { m.r64(cr + $o) } else { hr64(cp, $o) }
            };
        }
        macro_rules! s64 {
            ($o:expr, $v:expr) => {
                if P::CHECKED { m.w64(cr + $o, $v) } else { hw64(cp, $o, $v) }
            };
        }
        let wlen = g32!(O_U0);
        let mut last = er + 1 - wlen as i64;
        if last > sr && g32!(O_ANCH) != 0 {
            last = sr;
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
            return (0, 0);
        } else {
            s64!(O_HIT, sr);
            if sr < last {
                {
            pushed = true;
        }
            }
            sr += g32!(O_U0) as i64;
            return (sr, if pushed { 2 } else { 1 });
        }
    }

    pub fn first_match(&mut self, len: i32, hits: &mut Vec<i64>) -> i32 {
        let m = self.mem;
        let start = CDATA_BASE;
        m.w64(A_START_SRCH, start);
        self.end_srch = start + (len as i64 - 1);
        self.br1 = 0;
        let pu = m.r64(A_AD_PU_S);
        let i = self.pattern_match(pu, start, self.end_srch, hits, true);
        let v = hits[i as usize];
        m.w64(A_PAST_LAST, v);
        i
    }

    pub fn next_match(&mut self, hits: &mut Vec<i64>) -> i32 {
        let m = self.mem;
        let pu = m.r64(A_AD_PU_S);
        let (s, e) = (self.start_srch(), self.end_srch);
        let i = self.pattern_match(pu, s, e, hits, false);
        let v = hits[i as usize];
        m.w64(A_PAST_LAST, v);
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
        m.w64(A_PAST_LAST, v);
        i
    }
}
