//! Port of ggpunit.c: the pattern parser and the backtracking matcher.
//!
//! The code follows the C original statement by statement.  Pointers into
//! the sequence buffer are kept as `isize` offsets from the start of the
//! buffer, and pattern units live in an arena indexed by `usize`.  All
//! state that the C code keeps in globals (names, rule sets, the
//! backtrack register, last hit position, per-unit hit/mlen) persists
//! across calls just as it does in C.

use crate::sys::segv;

pub const PEPTIDE: i32 = 1;
pub const DNA: i32 = 2;

pub const MAX_SEQ_LEN: usize = 250_000_000;
const MAX_NAMES: i32 = 50;

/// Offset used for a NULL pointer (never a valid buffer offset).
pub const NULLP: isize = isize::MIN / 4;

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

/// A sequence buffer.  The C program mallocs MAX_SEQ_LEN+1 bytes once and
/// reuses them for every sequence, so bytes past the current sequence keep
/// whatever an earlier (longer) sequence left there.  Never-written bytes
/// read as zero.
pub struct Buf {
    pub v: Vec<u8>,
}

impl Buf {
    pub fn new() -> Buf {
        Buf { v: Vec::new() }
    }

    #[inline(always)]
    pub fn get(&self, off: isize) -> u8 {
        if off >= 0 && (off as usize) < self.v.len() {
            unsafe { *self.v.get_unchecked(off as usize) }
        } else if off >= 0 && (off as usize) <= MAX_SEQ_LEN {
            0
        } else {
            segv()
        }
    }

    #[inline(always)]
    pub fn set(&mut self, off: usize, b: u8) {
        if off >= self.v.len() {
            self.v.resize(off + 1, 0);
        }
        self.v[off] = b;
    }

    /// Make sure `n` bytes are addressable.
    pub fn reserve_len(&mut self, n: usize) {
        if self.v.len() < n {
            self.v.resize(n, 0);
        }
    }
}

#[derive(Clone)]
struct Punit {
    typ: i32,
    nxt: Option<usize>,
    prev: Option<usize>,
    br: Option<usize>,
    anchored: i32,
    hit: isize,
    mlen: i32,
    // or
    or1: usize,
    or2: usize,
    alt: i32,
    // exact / sim
    len: i32,
    code: usize,
    // any
    code_matrix: i64,
    // range
    min: i32,
    width: i32,
    rnxt: i32,
    // compl / repeat / sim
    ins: i32,
    del: i32,
    mis: i32,
    of: i32,
    rule_set: i32,
    // llim
    llim: Vec<i32>,
    bound: i32,
    // weight
    wlen: i32,
    vec: usize,
    cutoff: i32,
    tupsz: i32,
    maxwt: i32,
}

impl Default for Punit {
    fn default() -> Punit {
        Punit {
            typ: 0,
            nxt: None,
            prev: None,
            br: None,
            anchored: 0,
            hit: NULLP,
            mlen: 0,
            or1: 0,
            or2: 0,
            alt: 0,
            len: 0,
            code: 0,
            code_matrix: 0,
            min: 0,
            width: 0,
            rnxt: 0,
            ins: 0,
            del: 0,
            mis: 0,
            of: 0,
            rule_set: 0,
            llim: Vec::new(),
            bound: 0,
            wlen: 0,
            vec: 0,
            cutoff: 0,
            tupsz: 0,
            maxwt: 0,
        }
    }
}

/// Where the "one" side of loose_match reads from.
enum One<'a> {
    /// the sequence code buffer at an offset
    Cd(isize),
    /// a private byte array (pattern codes, reversed or complemented copy)
    Bytes(&'a [u8], usize),
}

pub struct Engine {
    pu: Vec<Punit>,
    pup: usize,
    names: [Option<usize>; 51],
    rule_sets: [Option<[u8; 16]>; 51],
    cv: Vec<u8>,
    cvp: usize,
    iv: Vec<i32>,
    ivp: usize,
    pub seq_type: i32,
    root: usize,
    br1: Option<usize>,
    start_srch: isize,
    end_srch: isize,
    past_last: isize,
    pub punit_to_code: [u8; 256],
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

fn build_punit_to_code() -> [u8; 256] {
    let mut t = [0u8; 256];
    for c in 0..256usize {
        let lc = if c < 128 { (c as u8).to_ascii_lowercase() } else { c as u8 };
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
        t[c] = v;
    }
    t
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
    p1: usize,
    p2: isize,
    n1: i32,
    n2: i32,
    mis: i32,
    ins: i32,
    del: i32,
    next_choice: i32,
}

impl Engine {
    pub fn new() -> Engine {
        Engine {
            pu: Vec::new(),
            pup: 0,
            names: [None; 51],
            rule_sets: [None; 51],
            cv: Vec::new(),
            cvp: 0,
            iv: Vec::new(),
            ivp: 0,
            seq_type: DNA,
            root: 0,
            br1: None,
            start_srch: 0,
            end_srch: 0,
            past_last: 0,
            punit_to_code: build_punit_to_code(),
        }
    }

    /// `punit_to_code[c]` where C indexes with a (signed) `char`.  Bytes
    /// >= 0x80 give a negative subscript and read the globals stored just
    /// before the table.  With the reference build (clang, macOS) 0xA8 lands
    /// on the low byte of `punit_sequence_type`; the other bytes there are
    /// zero padding, unset name slots or pointer bytes, and read as 0 here.
    #[inline(always)]
    fn p2c(&self, c: u8) -> u8 {
        if c < 0x80 {
            self.punit_to_code[c as usize]
        } else if c == 0xA8 {
            self.seq_type as u8
        } else {
            0
        }
    }

    /// comp_data(): translate characters to nucleotide codes (stops at NUL).
    pub fn comp_data(&self, data: &Buf, cdata: &mut Buf) {
        let mut k = 0usize;
        loop {
            let c = data.get(k as isize);
            if c == 0 {
                break;
            }
            cdata.set(k, self.p2c(c));
            k += 1;
        }
        cdata.set(k, 0);
    }

    // ------------------------------------------------------------------
    // arena helpers

    fn slot(&mut self, i: usize) -> &mut Punit {
        if self.pu.len() <= i {
            self.pu.resize(i + 1, Punit::default());
        }
        &mut self.pu[i]
    }

    fn alloc(&mut self) -> usize {
        let i = self.pup;
        self.pup += 1;
        self.slot(i);
        i
    }

    fn cv_set(&mut self, i: usize, b: u8) {
        if self.cv.len() <= i {
            self.cv.resize(i + 1, 0);
        }
        self.cv[i] = b;
    }

    #[inline(always)]
    fn cv_get(&self, i: usize) -> u8 {
        if i < self.cv.len() { self.cv[i] } else { 0 }
    }

    fn iv_set(&mut self, i: usize, v: i32) {
        if self.iv.len() <= i {
            self.iv.resize(i + 1, 0);
        }
        self.iv[i] = v;
    }

    #[inline(always)]
    fn iv_get(&self, i: usize) -> i32 {
        if i < self.iv.len() { self.iv[i] } else { 0 }
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
            self.cv_set(p1, code);
            if code != 0 {
                p1 += 1;
            } else {
                return None;
            }
            p += 1;
        }
        if p1 > start {
            self.cv_set(p1, 0);
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
                let c = at(l, p).to_ascii_uppercase();
                self.cv_set(p1, c);
                p1 += 1;
                p += 1;
            }
            self.cv_set(p1, 0);
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

    /// n_tuple(): writes the numbers at iv[t..]; returns (count, p)
    fn n_tuple(&mut self, l: &[u8], p: usize, t: usize) -> Option<(i32, usize)> {
        if at(l, p) != b'(' {
            return None;
        }
        let p = ws(l, p + 1);
        let (v, mut p) = num(l, p)?;
        self.iv_set(t, v);
        let mut n: i32 = 1;
        loop {
            if at(l, p) == b',' {
                p += 1;
                p = ws(l, p);
                match num(l, p) {
                    Some((v, q)) => {
                        self.iv_set(t + n as usize, v);
                        n += 1;
                        p = ws(l, q);
                    }
                    None => return None,
                }
            } else {
                let c = at(l, p);
                p += 1;
                if c == b')' {
                    return Some((n, p));
                } else {
                    return None;
                }
            }
        }
    }

    /// wt_template(): returns (tupsz, n, p).  Advances the iv cursor even on
    /// failure, as the C code does.
    fn wt_template(&mut self, l: &[u8], p: usize) -> Option<(i32, i32, usize)> {
        let mut n: i32 = 0;
        if at(l, p) != b'{' {
            return None;
        }
        let p = ws(l, p + 1);
        let (tupsz, mut p) = self.n_tuple(l, p, self.ivp)?;
        self.ivp = (self.ivp as isize + tupsz as isize) as usize;
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
                    self.ivp = (self.ivp as isize + sz as isize) as usize;
                    p = q;
                }
                _ => return None,
            }
        }
        let c = at(l, p);
        if c == b'}' && (tupsz == 4 || tupsz == 20 || tupsz == 21) {
            Some((tupsz, n, p + 1))
        } else {
            None
        }
    }

    /// returns (maxwt, cutoff, tupsz, n, p)
    fn wt_pat(&mut self, l: &[u8], p: usize) -> Option<(i32, i32, i32, i32, usize)> {
        // form: maxwt > {template} > cutoff
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

    /// llim_pat(): writes the name vector into `v` (v[0] = count)
    fn llim_pat(&self, l: &[u8], p: usize, v: &mut Vec<i32>) -> Option<(i32, usize)> {
        v.clear();
        v.resize(52, 0);
        v[0] = 1;
        let p = word(b"length(", l, p)?;
        let (n, mut p) = name_id(l, p)?;
        v[1] = n;
        loop {
            let p1 = ws(l, p);
            if at(l, p1) != b'+' {
                break;
            }
            let p1 = ws(l, p1 + 1);
            let idx = 1 + v[0] as usize;
            v[0] += 1;
            match name_id(l, p1) {
                Some((n, q)) => {
                    if idx >= v.len() {
                        v.resize(idx + 1, 0);
                    }
                    v[idx] = n;
                    p = q;
                }
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
        Some((max, p))
    }

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
        let or_pu = self.alloc();
        let ok = (|| {
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
        match ok {
            Some((p1, p2, p)) => {
                {
                    let o = self.slot(or_pu);
                    o.typ = OR_PUNIT;
                    o.or1 = p1;
                    o.or2 = p2;
                }
                self.slot(p1).prev = Some(or_pu);
                self.slot(p2).prev = Some(or_pu);
                let o = self.slot(or_pu);
                o.nxt = None;
                o.prev = None;
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
        self.rule_sets[n as usize] = Some(rs);
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
        let pu1 = self.alloc();
        if let Some((i, p1)) = name_assgn(l, p) {
            if self.names[i as usize].is_some() {
                return None; // slot is not released, as in C
            }
            self.names[i as usize] = Some(pu1);
            p = p1;
        }
        p = ws(l, p);
        let p1: usize;
        if at(l, p) == b'^' {
            self.pu[pu1].typ = MATCH_START;
            p1 = p + 1;
        } else if at(l, p) == b'$' {
            self.pu[pu1].typ = MATCH_END;
            p1 = p + 1;
        } else if let Some((i, j, q)) = range_pat(l, p) {
            let u = &mut self.pu[pu1];
            u.typ = RANGE_PUNIT;
            u.min = i;
            u.width = j.wrapping_sub(i);
            p1 = q;
        } else if let Some((cm, q)) = self.any_pat(l, p) {
            let u = &mut self.pu[pu1];
            u.code_matrix = cm;
            u.typ = ANY_PUNIT;
            p1 = q;
        } else if let Some((bound, q)) = {
            let mut v = std::mem::take(&mut self.pu[pu1].llim);
            let r = self.llim_pat(l, p, &mut v);
            self.pu[pu1].llim = v;
            r
        } {
            let u = &mut self.pu[pu1];
            u.bound = bound;
            u.typ = LLIM_PUNIT;
            p1 = q;
        } else if let Some((n, i, j, k, q)) = self.inv_rep_pat(l, p) {
            let u = &mut self.pu[pu1];
            u.typ = INV_REP_PUNIT;
            u.of = n;
            u.mis = i;
            u.ins = j;
            u.del = k;
            p1 = q;
        } else if let Some((n, i, j, k, q)) = self.repeat_pat(l, p) {
            let u = &mut self.pu[pu1];
            u.typ = REPEAT_PUNIT;
            u.of = n;
            u.mis = i;
            u.ins = j;
            u.del = k;
            p1 = q;
        } else if let Some((p2, (i, j, k, q))) = {
            let p2 = self.cvp;
            self.sim_pat(l, p).map(|r| (p2, r))
        } {
            let mut len = 0usize;
            while self.cv_get(p2 + len) != 0 {
                len += 1;
            }
            let u = &mut self.pu[pu1];
            u.code = p2;
            u.len = len as i32;
            if i == 0 && j == 0 && k == 0 {
                u.typ = EXACT_PUNIT;
            } else {
                u.typ = SIM_PUNIT;
                u.mis = i;
                u.ins = j;
                u.del = k;
            }
            p1 = q;
        } else if let Some((rs, n, i, j, k, q)) = self.compl_pat(l, p) {
            let u = &mut self.pu[pu1];
            u.typ = COMPL_PUNIT;
            u.rule_set = rs;
            u.of = n;
            u.mis = i;
            u.ins = j;
            u.del = k;
            p1 = q;
        } else if let Some((i1, (k, i, tupsz, j, q))) = {
            let i1 = self.ivp;
            self.wt_pat(l, p).map(|r| (i1, r))
        } {
            let u = &mut self.pu[pu1];
            u.typ = WEIGHT_PUNIT;
            u.vec = i1;
            u.wlen = j;
            u.cutoff = i;
            u.maxwt = k;
            u.tupsz = tupsz;
            p1 = q;
        } else {
            self.pup = pu1;
            return None;
        }
        Some(p1)
    }

    fn punit_list(&mut self, l: &[u8], p: usize) -> Option<usize> {
        let mut last: Option<usize> = None;
        let mut p = ws(l, p);
        let mut next = self.pup;
        while let Some(p1) = self.punit_parse(l, p) {
            self.slot(next).prev = last;
            if let Some(la) = last {
                self.slot(la).nxt = Some(next);
            }
            last = Some(next);
            next = self.pup;
            p = ws(l, p1);
        }
        let la = last?;
        self.slot(la).nxt = None;
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
        let mut last: Option<usize> = None;
        let mut p = ws(l, p);
        let mut next = self.pup;
        while let Some(p1) = self.punit_cons(l, p) {
            self.slot(next).prev = last;
            if let Some(la) = last {
                self.slot(la).nxt = Some(next);
            }
            while let Some(n) = self.slot(next).nxt {
                next = n;
            }
            last = Some(next);
            next = self.pup;
            p = ws(l, p1);
        }
        let la = last?;
        self.slot(la).nxt = None;
        Some(p)
    }

    fn parser(&mut self, l: &[u8]) -> Option<usize> {
        let first = self.pup;
        self.slot(first).prev = None;
        let p = self.punit_cons_list(l, ws(l, 0))?;
        let p = ws(l, p);
        if at(l, p) == 0 { Some(first) } else { None }
    }

    fn set_anchors_on(&mut self, pu: usize) {
        let mut cur = Some(pu);
        while let Some(c) = cur {
            self.pu[c].anchored = 1;
            if self.pu[c].typ == OR_PUNIT {
                let (a, b) = (self.pu[c].or1, self.pu[c].or2);
                self.set_anchors_on(a);
                self.set_anchors_on(b);
            }
            cur = self.pu[c].nxt;
        }
    }

    fn set_anchors(&mut self, pu: usize) {
        self.pu[pu].anchored = 0;
        if self.pu[pu].typ == OR_PUNIT {
            let (a, b) = (self.pu[pu].or1, self.pu[pu].or2);
            self.set_anchors(a);
            self.set_anchors(b);
        }
        if let Some(n) = self.pu[pu].nxt {
            self.set_anchors_on(n);
        }
    }

    fn name(&self, i: i32) -> usize {
        match self.names.get(i as usize).copied().flatten() {
            Some(p) => p,
            None => segv(),
        }
    }

    fn max_mat(&self, pu: usize, depth: u32) -> i32 {
        if depth > 1_000_000 {
            segv(); // unbounded recursion overflows the C stack
        }
        let u = &self.pu[pu];
        match u.typ {
            EXACT_PUNIT => u.len,
            LLIM_PUNIT => 0,
            RANGE_PUNIT => u.min.wrapping_add(u.width),
            ANY_PUNIT => 1,
            COMPL_PUNIT | REPEAT_PUNIT | INV_REP_PUNIT => {
                let p1 = self.name(u.of);
                self.max_mat(p1, depth + 1).wrapping_add(u.ins)
            }
            SIM_PUNIT => u.len.wrapping_add(u.ins),
            WEIGHT_PUNIT => u.wlen,
            OR_PUNIT => {
                let a = self.max_mat(u.or1, depth + 1);
                let b = self.max_mat(u.or2, depth + 1);
                if a > b { a } else { b }
            }
            _ => 0,
        }
    }

    fn max_mats(&self, pu: usize) -> i32 {
        let mut sum: i32 = 0;
        let mut cur = Some(pu);
        while let Some(c) = cur {
            sum = sum.wrapping_add(self.max_mat(c, 0));
            cur = self.pu[c].nxt;
        }
        sum
    }

    /// parse_dna_cmd / parse_peptide_cmd.  Returns the C return value
    /// (0 means "failed to parse").
    pub fn parse_cmd(&mut self, line: &[u8], seq_type: i32) -> i32 {
        self.seq_type = seq_type;
        for i in 0..MAX_NAMES as usize {
            self.names[i] = None;
        }
        self.pu.clear();
        self.pu.resize(100, Punit::default());
        self.pup = 0;
        self.cvp = 0;
        self.ivp = 0;
        if line.len() >= 1_000_000 {
            return 0;
        }
        match self.parser(line) {
            None => 0,
            Some(first) => {
                self.root = first;
                self.set_anchors(first);
                self.max_mats(first)
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

    #[inline(always)]
    fn ex_matches(&self, rule_set: i32, c1: u8, c2: u8, cd: &Buf) -> bool {
        if rule_set == -1 {
            self.matches(c1, c2)
        } else {
            let idx = c2 as isize + KNOWN_CHAR_INDEX[(c1 & 15) as usize] as isize;
            if rule_set == 50 {
                // rule_sets[50] lies past the end of the C array and aliases
                // start_srch, i.e. the start of the sequence code buffer.
                return cd.get(idx) != 0;
            }
            match self.rule_sets.get(rule_set as usize).copied().flatten() {
                Some(rs) => {
                    if (0..16).contains(&idx) {
                        rs[idx as usize] != 0
                    } else {
                        false
                    }
                }
                None => segv(),
            }
        }
    }

    #[inline(always)]
    fn one_get(&self, one: &One, i: isize, cd: &Buf) -> u8 {
        match one {
            One::Cd(base) => cd.get(base + i),
            One::Bytes(b, base) => {
                let k = *base as isize + i;
                if k >= 0 && (k as usize) < b.len() { b[k as usize] } else { 0 }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn loose_match(
        &self,
        cd: &Buf,
        one_src: One,
        mut one_len: i32,
        two_start: isize,
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
                if !self.known_char(self.one_get(&one_src, i as isize, cd) & 15) {
                    return 0;
                }
                i += 1;
            }
            let n = if one_len > 0 { one_len as usize } else { 0 };
            let mut r = vec![0u8; n];
            let mut len = one_len;
            let mut k: isize = 0;
            while len > 0 {
                let i = self.one_get(&one_src, k, cd) as u32;
                k += 1;
                len -= 1;
                r[len as usize] = if rule_set == -1 {
                    ((i >> 4) & 15) as u8
                } else {
                    ((KNOWN_CHAR_INDEX[(i & 15) as usize] as i32) << 2) as u8
                };
            }
            result = r;
            One::Bytes(&result, 0)
        } else {
            one_src
        };
        // position of "one" relative to its source
        let mut one: isize = 0;
        let mut two: isize = two_start;

        if max_ins == 0 && max_del == 0 {
            if one_len > two_len {
                return 0;
            }
            let mut i = one_len;
            while i >= 1 {
                let t = cd.get(two);
                if !self.known_char(t & 15)
                    || (!self.ex_matches(rule_set, t, self.one_get(&one_src, one, cd), cd) && {
                        max_mis -= 1;
                        max_mis < 0
                    })
                {
                    return 0;
                } else {
                    two += 1;
                    one += 1;
                }
                i -= 1;
            }
            return (two - two_start + 1) as i32;
        }

        let mut stack: Vec<StackEnt> = Vec::new();
        let mut nxtent: i32 = 0;
        macro_rules! push {
            ($n:expr) => {{
                let e = StackEnt {
                    p1: one as usize,
                    p2: two,
                    n1: one_len,
                    n2: two_len,
                    mis: max_mis,
                    ins: max_ins,
                    del: max_del,
                    next_choice: $n,
                };
                let k = nxtent as usize;
                if stack.len() <= k {
                    stack.resize(k + 1, StackEnt::default());
                }
                stack[k] = e;
                nxtent += 1;
            }};
        }
        let mut lbl = Lm::Top;
        loop {
            match lbl {
                Lm::Top => {
                    if !(two_len != 0 || nxtent != 0) {
                        return 0;
                    }
                    if two_len != 0 && one_len != 0 && {
                        let t = cd.get(two);
                        self.known_char(t & 15) && self.ex_matches(rule_set, t, self.one_get(&one_src, one, cd), cd)
                    } {
                        two += 1;
                        one += 1;
                        two_len -= 1;
                        one_len -= 1;
                        if one_len == 0 {
                            return (two - two_start + 1) as i32;
                        }
                    } else if max_mis != 0 && one_len >= 1 && two_len >= 1 {
                        if max_ins != 0 {
                            push!(1);
                        } else if max_del != 0 {
                            push!(2);
                        }
                        max_mis -= 1;
                        one += 1;
                        two += 1;
                        one_len -= 1;
                        two_len -= 1;
                        if one_len == 0 {
                            return (two - two_start + 1) as i32;
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
                        one = e.p1 as isize;
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
                        return 0;
                    }
                }
                Lm::Ins => {
                    max_ins -= 1;
                    one += 1;
                    one_len -= 1;
                    if one_len == 0 {
                        return (two - two_start + 1) as i32;
                    }
                    lbl = Lm::Top;
                }
                Lm::Del => {
                    max_del -= 1;
                    two += 1;
                    two_len -= 1;
                    if one_len == 0 {
                        return (two - two_start + 1) as i32;
                    }
                    lbl = Lm::Top;
                }
            }
        }
    }

    fn next_punit(&self, pu: usize) -> Option<usize> {
        if let Some(n) = self.pu[pu].nxt {
            return Some(n);
        }
        let mut pu1 = pu;
        let mut pu2 = self.pu[pu1].prev;
        while let Some(p2) = pu2 {
            let n = self.pu[p2].nxt;
            if n == Some(pu1) || n.is_none() {
                pu1 = p2;
                pu2 = self.pu[pu1].prev;
            } else {
                break;
            }
        }
        if let Some(p2) = pu2 {
            if self.pu[p2].nxt != Some(pu1) {
                return self.pu[p2].nxt;
            }
        }
        None
    }

    fn collect_hits(&self, pu: usize, revhits: &mut Vec<isize>) {
        revhits.clear();
        let mut last: Option<usize> = None;
        let mut cur = Some(pu);
        while let Some(p) = cur {
            let u = &self.pu[p];
            if u.typ != OR_PUNIT {
                revhits.push(u.hit);
                last = Some(p);
                cur = u.prev;
            } else if u.nxt == last {
                let mut q = if u.alt == 1 { u.or1 } else { u.or2 };
                while let Some(n) = self.pu[q].nxt {
                    q = n;
                }
                cur = Some(q);
                last = None;
            } else {
                last = Some(p);
                cur = u.prev;
            }
        }
    }

    fn pattern_match(&mut self, cd: &Buf, start: isize, end: isize, hits: &mut Vec<isize>, first: bool) -> i32 {
        let mut br = self.br1;
        let mut sr = start;
        let er = end;
        let mut cr = self.root;
        let mut lbl = if first { L::Try } else { L::Back };

        macro_rules! success {
            () => {{
                self.pu[cr].mlen = (sr - self.pu[cr].hit) as i32;
                match self.next_punit(cr) {
                    Some(n) => {
                        cr = n;
                        lbl = L::Try;
                    }
                    None => lbl = L::Leave,
                }
            }};
        }
        macro_rules! push_br {
            () => {{
                self.pu[cr].br = br;
                br = Some(cr);
            }};
        }

        loop {
            match lbl {
                L::Try => {
                    let typ = self.pu[cr].typ;
                    match typ {
                        MATCH_START => {
                            if sr == start {
                                self.pu[cr].hit = sr;
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        }
                        MATCH_END => {
                            if sr == end + 1 {
                                self.pu[cr].hit = sr;
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        }
                        ANY_PUNIT => {
                            let mut last = er;
                            if last > sr && self.pu[cr].anchored != 0 {
                                last = sr;
                            }
                            let cm = self.pu[cr].code_matrix;
                            while sr <= last {
                                let i = cd.get(sr) as i8 as i32;
                                if i >= b'A' as i32 && i <= b'Z' as i32 && ((1i32 << (i - b'A' as i32)) as i64 & cm) != 0 {
                                    break;
                                }
                                sr += 1;
                            }
                            if sr > last {
                                lbl = L::Back;
                            } else {
                                self.pu[cr].hit = sr;
                                if sr < last {
                                    push_br!();
                                }
                                sr += 1;
                                success!();
                            }
                        }
                        LLIM_PUNIT => {
                            let mut ln: i32 = 0;
                            let mut i = self.pu[cr].llim.first().copied().unwrap_or(0);
                            while i != 0 {
                                let nm = self.pu[cr].llim.get(i as usize).copied().unwrap_or(0);
                                i -= 1;
                                let pu1 = self.name(nm);
                                ln = ln.wrapping_add(self.pu[pu1].mlen);
                            }
                            if ln < self.pu[cr].bound {
                                self.pu[cr].hit = sr;
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        }
                        RANGE_PUNIT => {
                            let min = self.pu[cr].min;
                            let i = min.wrapping_sub(1);
                            if sr + i as isize <= er {
                                self.pu[cr].hit = sr;
                                sr += i as isize + 1;
                                if (sr <= er && self.pu[cr].width != 0) || (self.pu[cr].anchored == 0 && sr <= er) {
                                    push_br!();
                                    self.pu[cr].rnxt = min.wrapping_add(1);
                                }
                                success!();
                            } else {
                                lbl = L::Back;
                            }
                        }
                        EXACT_PUNIT => {
                            let len = self.pu[cr].len;
                            let mut last = er + 1 - len as isize;
                            if last > sr && self.pu[cr].anchored != 0 {
                                last = sr;
                            }
                            let p1 = self.pu[cr].code;
                            let ln = len - 1;
                            let c0 = self.cv_get(p1);
                            while sr <= last {
                                if self.matches(cd.get(sr), c0) {
                                    let mut p2 = sr + 1;
                                    let mut p3 = p1 + 1;
                                    let mut i = ln;
                                    while i != 0 && self.matches(cd.get(p2), self.cv_get(p3)) {
                                        i -= 1;
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
                                self.pu[cr].hit = sr;
                                if sr < last {
                                    push_br!();
                                }
                                sr += len as isize;
                                success!();
                            }
                        }
                        COMPL_PUNIT => {
                            let u = &self.pu[cr];
                            let (rs, ins, del, mis, of) = (u.rule_set, u.ins, u.del, u.mis, u.of);
                            let pu1 = self.name(of);
                            let p1 = self.pu[pu1].hit;
                            let mut ln = self.pu[pu1].mlen;
                            if rs != -1 || ins != 0 || del != 0 || mis != 0 {
                                let i = self.loose_match(cd, One::Cd(p1), ln, sr, (er + 1 - sr) as i32, ins, del, mis, rs, true);
                                if i != 0 {
                                    let i = i - 1;
                                    self.pu[cr].hit = sr;
                                    sr += i as isize;
                                    success!();
                                } else {
                                    lbl = L::Back;
                                }
                            } else if er - sr >= (ln as isize) - 1 {
                                self.pu[cr].hit = sr;
                                let mut q = p1 + (ln as isize - 1);
                                let mut ok = true;
                                while ln != 0 {
                                    ln = ln.wrapping_sub(1);
                                    let c = cd.get(q);
                                    if !self.known_char(c & 15) || ((c >> 4) & 15) != (cd.get(sr) & 15) {
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
                            let u = &self.pu[cr];
                            let (ins, del, mis, of) = (u.ins, u.del, u.mis, u.of);
                            let pu1 = self.name(of);
                            let p1 = self.pu[pu1].hit;
                            let ln = self.pu[pu1].mlen;
                            if ln == 0 {
                                self.pu[cr].hit = sr;
                                success!();
                            } else {
                                let i = self.loose_match(cd, One::Cd(p1), ln, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
                                if i != 0 {
                                    let i = i - 1;
                                    self.pu[cr].hit = sr;
                                    sr += i as isize;
                                    success!();
                                } else {
                                    lbl = L::Back;
                                }
                            }
                        }
                        INV_REP_PUNIT => {
                            let u = &self.pu[cr];
                            let (ins, del, mis, of) = (u.ins, u.del, u.mis, u.of);
                            let pu1 = self.name(of);
                            let p1 = self.pu[pu1].hit;
                            let ln = self.pu[pu1].mlen;
                            if ln == 0 {
                                self.pu[cr].hit = sr;
                                success!();
                            } else {
                                let n = if ln > 0 { ln as usize } else { 0 };
                                let mut p3 = vec![0u8; n];
                                for i in 0..n {
                                    p3[i] = cd.get(p1 + (n - 1 - i) as isize);
                                }
                                let i = self.loose_match(cd, One::Bytes(&p3, 0), ln, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
                                if i != 0 {
                                    let i = i - 1;
                                    self.pu[cr].hit = sr;
                                    sr += i as isize;
                                    success!();
                                } else {
                                    lbl = L::Back;
                                }
                            }
                        }
                        SIM_PUNIT => {
                            let u = &self.pu[cr];
                            let (ins, del, mis, len, code) = (u.ins, u.del, u.mis, u.len, u.code);
                            let mut last = er + 1 + ins as isize - len as isize;
                            if last > sr && u.anchored != 0 {
                                last = sr;
                            }
                            let mut found = false;
                            while sr <= last {
                                let i = self.loose_match(cd, One::Bytes(&self.cv, code), len, sr, (er + 1 - sr) as i32, ins, del, mis, -1, false);
                                if i != 0 {
                                    let i = i - 1;
                                    self.pu[cr].hit = sr;
                                    if sr < last {
                                        push_br!();
                                    }
                                    sr += i as isize;
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
                            let u = &self.pu[cr];
                            let (wlen, vec, tupsz, cutoff, maxwt) = (u.wlen, u.vec, u.tupsz, u.cutoff, u.maxwt);
                            let mut last = er + 1 - wlen as isize;
                            if last > sr && u.anchored != 0 {
                                last = sr;
                            }
                            while sr <= last {
                                let mut wval: i32 = 0;
                                let mut p1 = sr;
                                if tupsz == 4 {
                                    let mut pv = vec;
                                    let mut i = wlen;
                                    while i != 0 {
                                        let aw = self.iv_get(pv);
                                        let cw = self.iv_get(pv + 1);
                                        let gw = self.iv_get(pv + 2);
                                        let tw = self.iv_get(pv + 3);
                                        let add = match cd.get(p1) & 15 {
                                            0x1 => aw,
                                            0x2 => cw,
                                            0x4 => gw,
                                            0x8 => tw,
                                            0x3 => (aw >> 1).wrapping_add(cw >> 1),
                                            0x5 => (aw >> 1).wrapping_add(gw >> 1),
                                            0x9 => (aw >> 1).wrapping_add(tw >> 1),
                                            0x6 => (cw >> 1).wrapping_add(gw >> 1),
                                            0xA => (cw >> 1).wrapping_add(tw >> 1),
                                            0xC => (gw >> 1).wrapping_add(tw >> 1),
                                            0xE => (cw / 3).wrapping_add(gw / 3).wrapping_add(tw / 3),
                                            0xD => (aw / 3).wrapping_add(gw / 3).wrapping_add(tw / 3),
                                            0xB => (aw / 3).wrapping_add(cw / 3).wrapping_add(tw / 3),
                                            0x7 => (aw / 3).wrapping_add(cw / 3).wrapping_add(gw / 3),
                                            0xF => (aw >> 2).wrapping_add(cw >> 2).wrapping_add(gw >> 2).wrapping_add(tw >> 2),
                                            _ => 0,
                                        };
                                        wval = wval.wrapping_add(add);
                                        p1 += 1;
                                        i -= 1;
                                        pv += 4;
                                    }
                                } else {
                                    let mut pv = vec;
                                    let mut i = wlen;
                                    while i != 0 {
                                        let k: Option<usize> = match cd.get(p1) {
                                            b'A' => Some(0),
                                            b'C' => Some(1),
                                            b'D' => Some(2),
                                            b'E' => Some(3),
                                            b'F' => Some(4),
                                            b'G' => Some(5),
                                            b'H' => Some(6),
                                            b'I' => Some(7),
                                            b'K' => Some(8),
                                            b'L' => Some(9),
                                            b'M' => Some(10),
                                            b'N' => Some(11),
                                            b'P' => Some(12),
                                            b'Q' => Some(13),
                                            b'R' => Some(14),
                                            b'S' => Some(15),
                                            b'T' => Some(16),
                                            b'V' => Some(17),
                                            b'W' => Some(18),
                                            b'Y' => Some(19),
                                            _ => {
                                                if tupsz > 20 {
                                                    Some(20)
                                                } else {
                                                    None
                                                }
                                            }
                                        };
                                        if let Some(k) = k {
                                            wval = wval.wrapping_add(self.iv_get(pv + k));
                                        }
                                        p1 += 1;
                                        i -= 1;
                                        pv = (pv as isize + tupsz as isize) as usize;
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
                                self.pu[cr].hit = sr;
                                if sr < last {
                                    push_br!();
                                }
                                sr += wlen as isize;
                                success!();
                            }
                        }
                        OR_PUNIT => {
                            push_br!();
                            self.pu[cr].hit = sr;
                            self.pu[cr].alt = 1;
                            cr = self.pu[cr].or1;
                            lbl = L::Try;
                        }
                        _ => {
                            lbl = L::Back;
                        }
                    }
                }
                L::Back => {
                    let c = match br {
                        None => return 0,
                        Some(c) => c,
                    };
                    cr = c;
                    br = self.pu[cr].br;
                    sr = self.pu[cr].hit;
                    match self.pu[cr].typ {
                        RANGE_PUNIT => {
                            let u = &self.pu[cr];
                            let (nx, min, width, anchored) = (u.rnxt, u.min, u.width, u.anchored);
                            if nx <= min.wrapping_add(width) && sr + nx as isize - 1 <= er {
                                sr += nx as isize;
                                self.pu[cr].rnxt = nx.wrapping_add(1);
                                br = Some(cr);
                                success!();
                            } else {
                                self.pu[cr].hit += 1;
                                let h = self.pu[cr].hit;
                                if h + min as isize - 1 <= er && anchored == 0 {
                                    self.pu[cr].rnxt = min.wrapping_add(1);
                                    sr = h + min as isize;
                                    br = Some(cr);
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
                            if self.pu[cr].alt == 1 {
                                self.pu[cr].alt = 2;
                                cr = self.pu[cr].or2;
                                lbl = L::Try;
                            } else if self.pu[cr].anchored == 0 && self.pu[cr].alt == 2 {
                                sr += 1;
                                lbl = L::Try;
                            } else {
                                lbl = L::Back;
                            }
                        }
                        _ => {
                            lbl = L::Leave;
                        }
                    }
                }
                L::Leave => {
                    let mut revhits = Vec::new();
                    self.collect_hits(cr, &mut revhits);
                    let n = revhits.len();
                    if hits.len() < n + 1 {
                        hits.resize(n + 1, NULLP);
                    }
                    let mut j = 0usize;
                    for k in (0..n).rev() {
                        hits[j] = revhits[k];
                        j += 1;
                    }
                    hits[j] = sr;
                    self.br1 = br;
                    return j as i32;
                }
            }
        }
    }

    pub fn first_match(&mut self, cd: &Buf, len: i32, hits: &mut Vec<isize>) -> i32 {
        self.start_srch = 0;
        self.end_srch = len as isize - 1;
        self.br1 = None;
        let i = self.pattern_match(cd, self.start_srch, self.end_srch, hits, true);
        self.past_last = hits[i as usize];
        i
    }

    pub fn next_match(&mut self, cd: &Buf, hits: &mut Vec<isize>) -> i32 {
        let i = self.pattern_match(cd, self.start_srch, self.end_srch, hits, false);
        self.past_last = hits[i as usize];
        i
    }

    pub fn cont_match(&mut self, cd: &Buf, hits: &mut Vec<isize>) -> i32 {
        let past_last1 = self.past_last;
        let mut i;
        loop {
            i = self.next_match(cd, hits);
            if !(i > 0 && hits[0] < past_last1) {
                break;
            }
        }
        self.past_last = hits[i as usize];
        i
    }
}
