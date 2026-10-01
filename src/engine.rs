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

/// A byte buffer that the C code mallocs once (MAX_SEQ_LEN+1 bytes) and
/// reuses for every sequence: bytes past the current sequence keep what an
/// earlier, longer sequence left there; never-written bytes read as zero.
pub struct Buf {
    pub v: Vec<u8>,
}

impl Buf {
    pub fn new() -> Buf {
        Buf { v: Vec::new() }
    }

    /// Read at an offset from the start of the buffer.
    #[inline(always)]
    pub fn get(&self, off: i64) -> u8 {
        match self.try_get(off) {
            Some(b) => b,
            None => segv(),
        }
    }

    /// None when the C program would fault on this read.
    #[inline(always)]
    pub fn try_get(&self, off: i64) -> Option<u8> {
        if off >= 0 && (off as usize) < self.v.len() {
            Some(unsafe { *self.v.get_unchecked(off as usize) })
        } else if (-LOW_SLACK..ALLOC_LEN).contains(&off) {
            Some(0)
        } else {
            None
        }
    }

    #[inline(always)]
    pub fn set(&mut self, off: usize, b: u8) {
        if off >= self.v.len() {
            self.v.resize(off + 1, 0);
        }
        self.v[off] = b;
    }

    pub fn reserve_len(&mut self, n: usize) {
        if self.v.len() < n {
            self.v.resize(n, 0);
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

enum L {
    Try,
    Back,
    Leave,
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
    /// static data of the C program
    s: Vec<u8>,
    /// malloc'd rule-set blocks
    heap: Vec<u8>,
    /// the coded sequence (cdata)
    pub cdata: Buf,
    /// punit_sequence_type (also stored in `s`)
    seq_type: i32,
    /// parse cursors (C: *pu_s, *cv, *iv)
    pup: i64,
    cvp: i64,
    ivp: i64,
    br1: i64,
    end_srch: i64,
    /// punit_to_code for indexes 0..127 (never overwritten)
    p2c_lo: [u8; 128],
}

impl Engine {
    pub fn new() -> Engine {
        let mut e = Engine {
            s: vec![0u8; S_LEN],
            heap: Vec::new(),
            cdata: Buf::new(),
            seq_type: 0,
            pup: A_PU_S,
            cvp: A_CV,
            ivp: A_IV,
            br1: 0,
            end_srch: 0,
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
        let o = a.wrapping_sub(CDATA_BASE);
        if (-LOW_SLACK..ALLOC_LEN).contains(&o) {
            return self.cdata.get(o);
        }
        let o = a.wrapping_sub(S_BASE);
        if o >= 0 && (o as usize) < S_LEN {
            return unsafe { *self.s.get_unchecked(o as usize) };
        }
        self.rb_slow(a)
    }

    #[cold]
    fn rb_slow(&self, a: i64) -> u8 {
        if (S_END..LINKEDIT_END).contains(&a) {
            return 0;
        }
        let o = a.wrapping_sub(HEAP_BASE);
        if o >= 0 && (o as usize) < self.heap.len() {
            return self.heap[o as usize];
        }
        segv()
    }

    #[inline(always)]
    fn wb(&mut self, a: i64, v: u8) {
        let o = a.wrapping_sub(S_BASE);
        if o >= 0 && (o as usize) < S_LEN {
            self.s[o as usize] = v;
            return;
        }
        self.wb_slow(a, v)
    }

    #[cold]
    fn wb_slow(&mut self, a: i64, v: u8) {
        if (S_END..LINKEDIT_END).contains(&a) {
            sigbus(); // read-only segment
        }
        let o = a.wrapping_sub(HEAP_BASE);
        if o >= 0 && (o as usize) < self.heap.len() {
            self.heap[o as usize] = v;
            return;
        }
        let o = a.wrapping_sub(CDATA_BASE);
        if (0..ALLOC_LEN).contains(&o) {
            self.cdata.set(o as usize, v);
            return;
        }
        segv()
    }

    #[inline(always)]
    fn r32(&self, a: i64) -> i32 {
        let o = a.wrapping_sub(S_BASE);
        if o >= 0 && (o as usize) + 4 <= S_LEN {
            let o = o as usize;
            return i32::from_le_bytes([self.s[o], self.s[o + 1], self.s[o + 2], self.s[o + 3]]);
        }
        i32::from_le_bytes([self.rb(a), self.rb(a + 1), self.rb(a + 2), self.rb(a + 3)])
    }

    #[inline(always)]
    fn r64(&self, a: i64) -> i64 {
        let o = a.wrapping_sub(S_BASE);
        if o >= 0 && (o as usize) + 8 <= S_LEN {
            let o = o as usize;
            let mut b = [0u8; 8];
            b.copy_from_slice(&self.s[o..o + 8]);
            return i64::from_le_bytes(b);
        }
        let mut b = [0u8; 8];
        for (k, x) in b.iter_mut().enumerate() {
            *x = self.rb(a + k as i64);
        }
        i64::from_le_bytes(b)
    }

    #[inline(always)]
    fn w32(&mut self, a: i64, v: i32) {
        for (k, x) in v.to_le_bytes().iter().enumerate() {
            self.wb(a + k as i64, *x);
        }
    }

    #[inline(always)]
    fn w64(&mut self, a: i64, v: i64) {
        for (k, x) in v.to_le_bytes().iter().enumerate() {
            self.wb(a + k as i64, *x);
        }
    }

    fn malloc16(&mut self) -> i64 {
        let a = HEAP_BASE + self.heap.len() as i64;
        self.heap.resize(self.heap.len() + 16, 0);
        a
    }

    // globals that other C objects alias
    #[inline(always)]
    fn names(&self, i: i32) -> i64 {
        self.r64(A_NAMES + 8 * i as i64)
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
            let c = data.get(k as i64);
            if c == 0 {
                break;
            }
            let code = self.p2c(c);
            self.cdata.set(k, code);
            k += 1;
        }
        self.cdata.set(k, 0);
    }

    /// strcpy(cdata, data)
    pub fn copy_data(&mut self, data: &Buf) {
        let mut k = 0usize;
        loop {
            let c = data.get(k as i64);
            self.cdata.set(k, c);
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

    #[inline(always)]
    fn known_char(&self, c: u8) -> bool {
        if self.seq_type == PEPTIDE { true } else { KNOWN_CHAR[c as usize] != 0 }
    }

    #[inline(always)]
    fn matches(&self, c1: u8, c2: u8) -> bool {
        if self.seq_type == PEPTIDE {
            c1 == c2 || c2 == b'X'
        } else {
            let a = c1 & 15;
            KNOWN_CHAR[a as usize] != 0 && (a & (c2 & 15)) == a
        }
    }

    /// ExMatches(RuleSet, C1, C2)
    #[inline(always)]
    fn ex_matches(&self, rule_set: i32, c1: u8, c2: u8) -> bool {
        if rule_set == -1 {
            self.matches(c1, c2)
        } else {
            let base = self.r64(A_RULE_SETS + 8 * rule_set as i64);
            let idx = c2 as i64 + KNOWN_CHAR_INDEX[(c1 & 15) as usize] as i64;
            self.rb(base + idx) != 0
        }
    }

    #[inline(always)]
    fn one_get(&self, one: &One, i: i64) -> u8 {
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
        let result: Vec<u8>;
        let one_src = if compl_flag {
            let mut i = 0;
            while i < one_len {
                if !self.known_char(self.one_get(&one_src, i as i64) & 15) {
                    return 0;
                }
                i += 1;
            }
            if one_len > MAX_CODES {
                // result[MAX_CODES] overflows into the stack guard
                abort();
            }
            let n = one_len.max(0) as usize;
            let mut r = vec![0u8; n];
            let mut len = one_len;
            let mut k: i64 = 0;
            while len > 0 {
                let i = self.one_get(&one_src, k) as u32;
                k += 1;
                len -= 1;
                r[len as usize] = if rule_set == -1 {
                    ((i >> 4) & 15) as u8
                } else {
                    ((KNOWN_CHAR_INDEX[(i & 15) as usize] as i32) << 2) as u8
                };
            }
            result = r;
            One::Bytes(&result)
        } else {
            one_src
        };
        let mut one: i64 = 0;
        let mut two: i64 = two_start;

        if max_ins == 0 && max_del == 0 {
            if one_len > two_len {
                return 0;
            }
            let mut i = one_len;
            while i >= 1 {
                let t = self.rb(two);
                if !self.known_char(t & 15)
                    || (!self.ex_matches(rule_set, t, self.one_get(&one_src, one)) && {
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

        let mut stack: Vec<StackEnt> = Vec::new();
        let mut nxtent: i32 = 0;
        let mut smashed = false;
        macro_rules! ret {
            ($v:expr) => {{
                if smashed {
                    abort();
                }
                return $v;
            }};
        }
        macro_rules! push {
            ($n:expr) => {{
                let k = nxtent as usize;
                if nxtent >= LOOSE_STACK_SAFE {
                    smashed = true;
                }
                if stack.len() <= k {
                    stack.resize(k + 1, StackEnt::default());
                }
                stack[k] = StackEnt {
                    p1: one,
                    p2: two,
                    n1: one_len,
                    n2: two_len,
                    mis: max_mis,
                    ins: max_ins,
                    del: max_del,
                    next_choice: $n,
                };
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
                        let t = self.rb(two);
                        self.known_char(t & 15) && self.ex_matches(rule_set, t, self.one_get(&one_src, one))
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
                        let e = stack[nxtent as usize];
                        one = e.p1;
                        two = e.p2;
                        one_len = e.n1;
                        two_len = e.n2;
                        max_mis = e.mis;
                        max_ins = e.ins;
                        max_del = e.del;
                        if e.next_choice == 1 {
                            if max_del != 0 {
                                stack[nxtent as usize].next_choice = 2;
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

    fn next_punit(&self, pu: i64) -> i64 {
        let n = self.r64(pu + O_NXT);
        if n != 0 {
            return n;
        }
        let mut pu1 = pu;
        let mut pu2 = self.r64(pu1 + O_PREV);
        while pu2 != 0 {
            let n = self.r64(pu2 + O_NXT);
            if n == pu1 || n == 0 {
                pu1 = pu2;
                pu2 = self.r64(pu1 + O_PREV);
            } else {
                break;
            }
        }
        if pu2 != 0 && self.r64(pu2 + O_NXT) != pu1 {
            return self.r64(pu2 + O_NXT);
        }
        0
    }

    fn collect_hits(&self, pu: i64, revhits: &mut Vec<i64>) {
        revhits.clear();
        let mut last: i64 = 0;
        let mut pu = pu;
        while pu != 0 {
            if self.r32(pu + O_TYPE) != OR_PUNIT {
                revhits.push(self.r64(pu + O_HIT));
                if revhits.len() > 1 << 20 {
                    segv(); // runs off the top of the C stack
                }
                last = pu;
                pu = self.r64(pu + O_PREV);
            } else if self.r64(pu + O_NXT) == last {
                pu = if self.r32(pu + O_U24) == 1 { self.r64(pu + O_U0) } else { self.r64(pu + O_U8) };
                loop {
                    let n = self.r64(pu + O_NXT);
                    if n == 0 {
                        break;
                    }
                    pu = n;
                }
                last = 0;
            } else {
                last = pu;
                pu = self.r64(pu + O_PREV);
            }
        }
    }

    fn pattern_match(&mut self, pu: i64, start: i64, end: i64, hits: &mut Vec<i64>, first: bool) -> i32 {
        let mut br = self.br1;
        let mut sr = start;
        let er = end;
        let mut cr = pu;
        let mut lbl = if first { L::Try } else { L::Back };

        macro_rules! success {
            () => {{
                let h = self.r64(cr + O_HIT);
                self.w32(cr + O_MLEN, (sr - h) as i32);
                let n = self.next_punit(cr);
                if n != 0 {
                    cr = n;
                    lbl = L::Try;
                } else {
                    lbl = L::Leave;
                }
            }};
        }
        macro_rules! push_br {
            () => {{
                self.w64(cr + O_BR, br);
                br = cr;
            }};
        }

        loop {
            match lbl {
                L::Try => match self.r32(cr + O_TYPE) {
                    MATCH_START => {
                        if sr == start {
                            self.w64(cr + O_HIT, sr);
                            success!();
                        } else {
                            lbl = L::Back;
                        }
                    }
                    MATCH_END => {
                        if sr == end + 1 {
                            self.w64(cr + O_HIT, sr);
                            success!();
                        } else {
                            lbl = L::Back;
                        }
                    }
                    ANY_PUNIT => {
                        let mut last = er;
                        if last > sr && self.r32(cr + O_ANCH) != 0 {
                            last = sr;
                        }
                        let cm = self.r64(cr + O_U0);
                        while sr <= last {
                            let i = self.rb(sr) as i8 as i32;
                            if i >= b'A' as i32 && i <= b'Z' as i32 && ((1i32 << (i - b'A' as i32)) as i64 & cm) != 0 {
                                break;
                            }
                            sr += 1;
                        }
                        if sr > last {
                            lbl = L::Back;
                        } else {
                            self.w64(cr + O_HIT, sr);
                            if sr < last {
                                push_br!();
                            }
                            sr += 1;
                            success!();
                        }
                    }
                    LLIM_PUNIT => {
                        let v = cr + O_U0;
                        let mut ln: i32 = 0;
                        let mut i = self.r32(v);
                        while i != 0 {
                            let nm = self.r32(v + 4 * i as i64);
                            i = i.wrapping_sub(1);
                            let pu1 = self.names(nm);
                            ln = ln.wrapping_add(self.r32(pu1 + O_MLEN));
                        }
                        if ln < self.r32(cr + O_LLIM_BOUND) {
                            self.w64(cr + O_HIT, sr);
                            success!();
                        } else {
                            lbl = L::Back;
                        }
                    }
                    RANGE_PUNIT => {
                        let min = self.r32(cr + O_U0);
                        let i = min.wrapping_sub(1);
                        if sr + i as i64 <= er {
                            self.w64(cr + O_HIT, sr);
                            sr += i as i64 + 1;
                            if (sr <= er && self.r32(cr + O_U4) != 0) || (self.r32(cr + O_ANCH) == 0 && sr <= er) {
                                push_br!();
                                let m = self.r32(cr + O_U0);
                                self.w32(cr + O_U8, m.wrapping_add(1));
                            }
                            success!();
                        } else {
                            lbl = L::Back;
                        }
                    }
                    EXACT_PUNIT => {
                        let len = self.r32(cr + O_U0);
                        let mut last = er + 1 - len as i64;
                        if last > sr && self.r32(cr + O_ANCH) != 0 {
                            last = sr;
                        }
                        let p1 = self.r64(cr + O_U8);
                        let ln = len.wrapping_sub(1);
                        while sr <= last {
                            if self.matches(self.rb(sr), self.rb(p1)) {
                                let mut p2 = sr + 1;
                                let mut p3 = p1 + 1;
                                let mut i = ln;
                                while i != 0 && self.matches(self.rb(p2), self.rb(p3)) {
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
                            lbl = L::Back;
                        } else {
                            self.w64(cr + O_HIT, sr);
                            if sr < last {
                                push_br!();
                            }
                            sr += self.r32(cr + O_U0) as i64;
                            success!();
                        }
                    }
                    COMPL_PUNIT => {
                        let pu1 = self.names(self.r32(cr + O_U12));
                        let p1 = self.r64(pu1 + O_HIT);
                        let mut ln = self.r32(pu1 + O_MLEN);
                        let rs = self.r32(cr + O_U16);
                        let (ins, del, mis) = (self.r32(cr + O_U0), self.r32(cr + O_U4), self.r32(cr + O_U8));
                        if rs != -1 || ins != 0 || del != 0 || mis != 0 {
                            let i = self.loose_match(One::Mem(p1), ln, sr, (er + 1 - sr) as i32, ins, del, mis, rs, true);
                            if i != 0 {
                                let i = i - 1;
                                self.w64(cr + O_HIT, sr);
                                sr += i as i64;
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        } else if er - sr >= ln as i64 - 1 {
                            self.w64(cr + O_HIT, sr);
                            let mut q = self.r64(pu1 + O_HIT) + (ln as i64 - 1);
                            let mut ok = true;
                            while ln != 0 {
                                ln = ln.wrapping_sub(1);
                                let c = self.rb(q);
                                if !self.known_char(c & 15) || ((c >> 4) & 15) != (self.rb(sr) & 15) {
                                    ok = false;
                                    break;
                                }
                                q -= 1;
                                sr += 1;
                            }
                            if ok {
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        } else {
                            lbl = L::Back;
                        }
                    }
                    REPEAT_PUNIT => {
                        let pu1 = self.names(self.r32(cr + O_U12));
                        let p1 = self.r64(pu1 + O_HIT);
                        let ln = self.r32(pu1 + O_MLEN);
                        if ln == 0 {
                            self.w64(cr + O_HIT, sr);
                            success!();
                        } else {
                            let (ins, del, mis) = (self.r32(cr + O_U0), self.r32(cr + O_U4), self.r32(cr + O_U8));
                            let i = self.loose_match(One::Mem(p1), ln, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
                            if i != 0 {
                                let i = i - 1;
                                self.w64(cr + O_HIT, sr);
                                sr += i as i64;
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        }
                    }
                    INV_REP_PUNIT => {
                        let pu1 = self.names(self.r32(cr + O_U12));
                        let p1 = self.r64(pu1 + O_HIT);
                        let ln = self.r32(pu1 + O_MLEN);
                        if ln == 0 {
                            self.w64(cr + O_HIT, sr);
                            success!();
                        } else {
                            let n = ln.max(0) as usize;
                            let mut p3 = vec![0u8; n];
                            for (i, x) in p3.iter_mut().enumerate() {
                                *x = self.rb(p1 + (n - 1 - i) as i64);
                            }
                            let (ins, del, mis) = (self.r32(cr + O_U0), self.r32(cr + O_U4), self.r32(cr + O_U8));
                            let i = self.loose_match(One::Bytes(&p3), ln, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
                            if i != 0 {
                                let i = i - 1;
                                self.w64(cr + O_HIT, sr);
                                sr += i as i64;
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        }
                    }
                    SIM_PUNIT => {
                        let len = self.r32(cr + O_U12);
                        let ins = self.r32(cr + O_U0);
                        let mut last = er + 1 + ins as i64 - len as i64;
                        if last > sr && self.r32(cr + O_ANCH) != 0 {
                            last = sr;
                        }
                        let mut found = false;
                        while sr <= last {
                            let code = self.r64(cr + O_U16);
                            let (len, ins, del, mis) =
                                (self.r32(cr + O_U12), self.r32(cr + O_U0), self.r32(cr + O_U4), self.r32(cr + O_U8));
                            let i = self.loose_match(One::Mem(code), len, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
                            if i != 0 {
                                let i = i - 1;
                                self.w64(cr + O_HIT, sr);
                                if sr < last {
                                    push_br!();
                                }
                                sr += i as i64;
                                found = true;
                                break;
                            }
                            sr += 1;
                        }
                        if found {
                            success!();
                        } else {
                            lbl = L::Back;
                        }
                    }
                    WEIGHT_PUNIT => {
                        let wlen = self.r32(cr + O_U0);
                        let mut last = er + 1 - wlen as i64;
                        if last > sr && self.r32(cr + O_ANCH) != 0 {
                            last = sr;
                        }
                        let pv1 = self.r64(cr + O_U8);
                        let tupsz = self.r32(cr + O_U20);
                        let cutoff = self.r32(cr + O_U16);
                        let maxwt = self.r32(cr + O_U24);
                        while sr <= last {
                            let mut wval: i32 = 0;
                            let mut p1 = sr;
                            let mut i = wlen;
                            if tupsz == 4 {
                                let mut pv = pv1;
                                while i != 0 {
                                    let c = self.rb(p1) & 15;
                                    p1 += 1;
                                    let aw = || self.r32(pv);
                                    let cw = || self.r32(pv + 4);
                                    let gw = || self.r32(pv + 8);
                                    let tw = || self.r32(pv + 12);
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
                                    let k: i64 = match self.rb(p1) {
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
                                        wval = wval.wrapping_add(self.r32(pv + 4 * k));
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
                            lbl = L::Back;
                        } else {
                            self.w64(cr + O_HIT, sr);
                            if sr < last {
                                push_br!();
                            }
                            sr += self.r32(cr + O_U0) as i64;
                            success!();
                        }
                    }
                    OR_PUNIT => {
                        push_br!();
                        self.w64(cr + O_U16, sr);
                        self.w64(cr + O_HIT, sr);
                        self.w32(cr + O_U24, 1);
                        cr = self.r64(cr + O_U0);
                        lbl = L::Try;
                    }
                    _ => lbl = L::Back,
                },
                L::Back => {
                    if br == 0 {
                        return 0;
                    }
                    cr = br;
                    br = self.r64(cr + O_BR);
                    sr = self.r64(cr + O_HIT);
                    match self.r32(cr + O_TYPE) {
                        RANGE_PUNIT => {
                            let nx = self.r32(cr + O_U8);
                            let min = self.r32(cr + O_U0);
                            let width = self.r32(cr + O_U4);
                            if nx <= min.wrapping_add(width) && sr + nx as i64 - 1 <= er {
                                self.w32(cr + O_U8, nx.wrapping_add(1));
                                sr += nx as i64;
                                br = cr;
                                success!();
                            } else {
                                let h = self.r64(cr + O_HIT) + 1;
                                self.w64(cr + O_HIT, h);
                                if h + self.r32(cr + O_U0) as i64 - 1 <= er && self.r32(cr + O_ANCH) == 0 {
                                    let m = self.r32(cr + O_U0);
                                    self.w32(cr + O_U8, m.wrapping_add(1));
                                    sr = self.r64(cr + O_HIT) + self.r32(cr + O_U0) as i64;
                                    br = cr;
                                    success!();
                                } else {
                                    lbl = L::Back;
                                }
                            }
                        }
                        LLIM_PUNIT | ANY_PUNIT | EXACT_PUNIT | SIM_PUNIT | WEIGHT_PUNIT => {
                            sr += 1;
                            lbl = L::Try;
                        }
                        OR_PUNIT => {
                            if self.r32(cr + O_U24) == 1 {
                                self.w32(cr + O_U24, 2);
                                cr = self.r64(cr + O_U8);
                                lbl = L::Try;
                            } else if self.r32(cr + O_ANCH) == 0 && self.r32(cr + O_U24) == 2 {
                                sr += 1;
                                lbl = L::Try;
                            } else {
                                lbl = L::Back;
                            }
                        }
                        _ => lbl = L::Leave,
                    }
                }
                L::Leave => {
                    let mut revhits = Vec::new();
                    self.collect_hits(cr, &mut revhits);
                    let n = revhits.len();
                    if hits.len() < n + 1 {
                        hits.resize(n + 1, 0);
                    }
                    for (j, k) in (0..n).rev().enumerate() {
                        hits[j] = revhits[k];
                    }
                    hits[n] = sr;
                    self.br1 = br;
                    if n > MAX_PUNITS {
                        // revhits[MAX_PUNITS] overflows into the stack guard
                        abort();
                    }
                    return n as i32;
                }
            }
        }
    }

    pub fn first_match(&mut self, len: i32, hits: &mut Vec<i64>) -> i32 {
        let start = CDATA_BASE;
        self.w64(A_START_SRCH, start);
        self.end_srch = start + (len as i64 - 1);
        self.br1 = 0;
        let pu = self.r64(A_AD_PU_S);
        let i = self.pattern_match(pu, start, self.end_srch, hits, true);
        let v = hits[i as usize];
        self.w64(A_PAST_LAST, v);
        i
    }

    pub fn next_match(&mut self, hits: &mut Vec<i64>) -> i32 {
        let pu = self.r64(A_AD_PU_S);
        let (s, e) = (self.start_srch(), self.end_srch);
        let i = self.pattern_match(pu, s, e, hits, false);
        let v = hits[i as usize];
        self.w64(A_PAST_LAST, v);
        i
    }

    pub fn cont_match(&mut self, hits: &mut Vec<i64>) -> i32 {
        let past_last1 = self.past_last();
        let mut i;
        loop {
            i = self.next_match(hits);
            if !(i > 0 && hits[0] < past_last1) {
                break;
            }
        }
        let v = hits[i as usize];
        self.w64(A_PAST_LAST, v);
        i
    }
}
