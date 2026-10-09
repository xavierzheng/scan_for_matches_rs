//! `--lint` and `--explain`: check a pattern for traps, and describe each
//! unit in plain words.  The engine is not changed; this module reads the
//! pattern text again (same grammar as `engine.rs`).

use crate::engine::compl;
use crate::fmt;
use std::cmp::Ordering;
use std::fmt::Write as _;

/// The pattern as read from the file.
pub struct Source<'a> {
    /// file name (for messages)
    pub file: &'a [u8],
    /// the joined pattern line, as given to the engine
    pub line: &'a [u8],
    /// source line (1-based) of each byte of `line`
    pub src_line: &'a [u32],
    /// (source line, text after '%') of each comment
    pub comments: &'a [(u32, Vec<u8>)],
    /// the file was cut at 31 999 bytes and text was lost
    pub truncated: bool,
}

/// What the engine found when it parsed the pattern.
pub struct Facts {
    /// the engine parsed the pattern (false: "failed to parse pattern")
    pub parsed: bool,
    /// `Engine::uses_earlier_state`
    pub uses_earlier_state: bool,
    /// `Engine::can_split`
    pub can_split: bool,
    /// `Engine::unit_offsets` (empty when not parsed)
    pub units: Vec<(u32, usize)>,
}

/// The options of the run.
pub struct RunOpts {
    pub protein: bool,
    pub complements: bool,
    pub dedup: bool,
    /// `--format` was given
    pub format: bool,
    /// the value given to `-o` (None: no `-o`)
    pub o_value: Option<Vec<u8>>,
    /// `--strict-n`
    pub strict_n: bool,
}

const MAX_NAMES: i32 = 50;
const N_NAMES: usize = MAX_NAMES as usize + 1;

// ---------------------------------------------------------------------
// the pattern units

/// Mismatches, deletions and insertions of an inexact unit (`[m,d,i]`,
/// as in README.original): a deletion is a letter of the pattern missing
/// in the data, an insertion an extra letter in the data.
#[derive(Clone, Copy, Default, PartialEq)]
struct Err3 {
    mis: i32,
    del: i32,
    ins: i32,
}

impl Err3 {
    fn exact(self) -> bool {
        self == Err3::default()
    }
}

enum Kind {
    Start,
    End,
    /// min...max
    Range(i32, i32),
    /// a word (upper case) with its errors; exact when all are 0
    Word(Vec<u8>, Err3),
    /// `pN`
    Repeat(i32, Err3),
    /// `~pN` or `rK~pN`
    Compl(Option<i32>, i32, Err3),
    /// `<pN`
    Inv(i32, Err3),
    /// `max > {(..),..} > cut` (max None: not given)
    Weight(Option<i32>, Vec<Vec<i32>>, i32),
    /// `length(pA+pB) < bound`; a name is None when the text after `+`
    /// is not a name (the engine then reads unset memory)
    Length(Vec<Option<i32>>, i32),
    /// `any(..)` (false) or `notany(..)` (true)
    Any(Vec<u8>, bool),
    /// `( A | B )`: the branches and the offsets of `|` and `)`
    Or(Vec<Unit>, Vec<Unit>, usize, usize),
}

struct Unit {
    kind: Kind,
    /// `pN=` before the unit
    def: Option<i32>,
    /// slot number in the engine (pre-order, alternatives included)
    slot: u32,
    /// byte offset where the unit starts (as `Engine::unit_offsets`)
    off: usize,
    /// byte offset after the unit
    end: usize,
}

/// A rule set `rK={au,ua,...}`.
struct RuleSet {
    id: i32,
    pairs: Vec<[u8; 2]>,
    off: usize,
}

// ---------------------------------------------------------------------
// the reader (each function follows the one of the same name in engine.rs)

fn isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn at(l: &[u8], p: usize) -> u8 {
    l.get(p).copied().unwrap_or(0)
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
    if !at(l, p).is_ascii_digit() {
        return None;
    }
    let mut i: i32 = 0;
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
    num(l, p + 1).filter(|(n, _)| (0..=MAX_NAMES).contains(n))
}

fn name_assgn(l: &[u8], p: usize) -> Option<(i32, usize)> {
    let (n, p) = name_id(l, p)?;
    let p = ws(l, p);
    (at(l, p) == b'=').then_some((n, p + 1))
}

fn rule_id(l: &[u8], p: usize) -> Option<(i32, usize)> {
    if at(l, p) != b'r' {
        return None;
    }
    num(l, p + 1).filter(|(n, _)| (0..=MAX_NAMES).contains(n))
}

fn misinsdel(l: &[u8], p: usize) -> Option<(Err3, usize)> {
    if at(l, p) != b'[' {
        return None;
    }
    let (mis, p) = num(l, p + 1)?;
    if at(l, p) != b',' {
        return None;
    }
    let (del, p) = num(l, p + 1)?;
    if at(l, p) != b',' {
        return None;
    }
    let (ins, p) = num(l, p + 1)?;
    if at(l, p) != b']' || mis < 0 || ins < 0 || del < 0 {
        return None;
    }
    Some((Err3 { mis, del, ins }, p + 1))
}

fn opt_errs(l: &[u8], p: usize) -> (Err3, usize) {
    misinsdel(l, p).unwrap_or((Err3::default(), p))
}

fn starts(l: &[u8], p: usize, s: &[u8]) -> bool {
    l.len() >= p + s.len() && &l[p..p + s.len()] == s
}

fn range_pat(l: &[u8], p: usize) -> Option<(i32, i32, usize)> {
    let (min, p) = num(l, p)?;
    let p = ws(l, p);
    if !starts(l, p, b"...") {
        return None;
    }
    let p = ws(l, p + 3);
    let (max, p) = num(l, p)?;
    (min >= 0 && max >= 0).then_some((min, max, p))
}

/// A pattern letter in DNA mode (the IUPAC codes, U = T).
fn dna_letter(c: u8) -> bool {
    c < 0x80 && b"acgtumrwsykbdhvn".contains(&c.to_ascii_lowercase())
}

fn dna_char(c: u8) -> Option<u8> {
    match c.to_ascii_lowercase() {
        b'u' => Some(b't'),
        c @ (b'a' | b'c' | b'g' | b't') => Some(c),
        _ => None,
    }
}

fn n_tuple(l: &[u8], p: usize) -> Option<(Vec<i32>, usize)> {
    if at(l, p) != b'(' {
        return None;
    }
    let (v, mut p) = num(l, ws(l, p + 1))?;
    let mut t = vec![v];
    loop {
        if at(l, p) == b',' {
            let (v, q) = num(l, ws(l, p + 1))?;
            t.push(v);
            p = ws(l, q);
        } else {
            return (at(l, p) == b')').then_some((t, p + 1));
        }
    }
}

fn wt_template(l: &[u8], p: usize) -> Option<(Vec<Vec<i32>>, usize)> {
    if at(l, p) != b'{' {
        return None;
    }
    let (first, mut p) = n_tuple(l, ws(l, p + 1))?;
    let mut rows = vec![first];
    loop {
        p = ws(l, p);
        if at(l, p) != b',' {
            break;
        }
        match n_tuple(l, ws(l, p + 1)) {
            Some((t, q)) if t.len() == rows[0].len() => {
                rows.push(t);
                p = q;
            }
            _ => return None,
        }
    }
    (at(l, p) == b'}' && matches!(rows[0].len(), 4 | 20 | 21)).then_some((rows, p + 1))
}

fn wt_pat(l: &[u8], p: usize) -> Option<(Kind, usize)> {
    if let Some((max, p1)) = num(l, p) {
        let p1 = ws(l, p1);
        if at(l, p1) == b'>'
            && let Some((rows, p1)) = wt_template(l, ws(l, p1 + 1))
        {
            let p1 = ws(l, p1);
            if at(l, p1) == b'>'
                && let Some((cut, p1)) = num(l, ws(l, p1 + 1))
            {
                return Some((Kind::Weight(Some(max), rows, cut), p1));
            }
        }
    }
    let (rows, p) = wt_template(l, p)?;
    let p = ws(l, p);
    if at(l, p) != b'>' {
        return None;
    }
    let (cut, p) = num(l, ws(l, p + 1))?;
    Some((Kind::Weight(None, rows, cut), p))
}

fn llim_pat(l: &[u8], p: usize) -> Option<(Kind, usize)> {
    if !starts(l, p, b"length(") {
        return None;
    }
    let (n, mut p) = name_id(l, p + 7)?;
    let mut names = vec![Some(n)];
    loop {
        let p1 = ws(l, p);
        if at(l, p1) != b'+' {
            break;
        }
        let p1 = ws(l, p1 + 1);
        let stored = if at(l, p1) == b'p' {
            num(l, p1 + 1)
        } else {
            None
        };
        names.push(stored.map(|(n, _)| n));
        match stored {
            Some((n, q)) if (0..=MAX_NAMES).contains(&n) => p = q,
            _ => break,
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
    let (bound, p) = num(l, ws(l, p + 1))?;
    let names = names
        .into_iter()
        .map(|n| n.filter(|n| (0..=MAX_NAMES).contains(n)))
        .collect();
    Some((Kind::Length(names, bound), p))
}

/// A `length()` whose names the engine cannot read safely (a `+` not
/// followed by a name, or more than 50 names).
fn length_bad(names: &[Option<i32>]) -> bool {
    names.len() > MAX_NAMES as usize || names.iter().any(|n| n.is_none())
}

fn any_pat(l: &[u8], p: usize, pep: bool) -> Option<(Kind, usize)> {
    if !pep {
        return None;
    }
    let (mut p, not) = if starts(l, p, b"any(") {
        (p + 4, false)
    } else if starts(l, p, b"notany(") {
        (p + 7, true)
    } else {
        return None;
    };
    let s = p;
    while at(l, p).is_ascii_alphabetic() {
        p += 1;
    }
    let set = l[s..p].to_ascii_uppercase();
    (at(l, p) == b')').then_some((Kind::Any(set, not), p + 1))
}

/// The pattern as read, also when reading stopped early.
struct Parsed {
    pep: bool,
    /// the top units (only those read before the stop when `!ok`)
    units: Vec<Unit>,
    /// the whole line was read
    ok: bool,
    /// the farthest place where reading failed, and why
    stop: Option<(usize, &'static str)>,
    /// slot of the unit that defines each name
    names: [Option<u32>; N_NAMES],
    /// rule sets as defined (also inside parts read again later)
    rules: Vec<RuleSet>,
    rule_def: [bool; N_NAMES],
    /// (name, offset) of each name defined a second time
    dups: Vec<(i32, usize)>,
    /// next free slot
    pup: u32,
}

impl Parsed {
    fn new(l: &[u8], pep: bool) -> Parsed {
        let mut ps = Parsed {
            pep,
            units: Vec::new(),
            ok: false,
            stop: None,
            names: [None; N_NAMES],
            rules: Vec::new(),
            rule_def: [false; N_NAMES],
            dups: Vec::new(),
            pup: 0,
        };
        let p = ws(l, 0);
        match ps.cons_list(l, p) {
            Some((units, p)) => {
                ps.units = units;
                let p = ws(l, p);
                ps.ok = at(l, p) == 0;
                if !ps.ok {
                    ps.fail(p, "this text cannot be read");
                }
            }
            None => ps.fail(p, "no pattern unit can be read here"),
        }
        ps
    }

    fn fail(&mut self, p: usize, why: &'static str) {
        if self.stop.is_none_or(|(q, _)| p > q) {
            self.stop = Some((p, why));
        }
    }

    fn rule_set(&mut self, l: &[u8], p: usize) -> Option<usize> {
        if self.pep {
            return None;
        }
        let off = p;
        let (id, p) = rule_id(l, p)?;
        let p = ws(l, p);
        if at(l, p) != b'=' || at(l, p + 1) != b'{' {
            return None;
        }
        let mut pairs = Vec::new();
        let mut q = p + 2;
        loop {
            let a = dna_char(at(l, q))?;
            let b = dna_char(at(l, q + 1))?;
            pairs.push([a, b]);
            q += 2;
            if at(l, q) != b',' {
                break;
            }
            q += 1;
        }
        if at(l, q) != b'}' {
            return None;
        }
        self.rule_def[id as usize] = true;
        self.rules.push(RuleSet { id, pairs, off });
        Some(q + 1)
    }

    fn word(&self, l: &[u8], p: usize) -> Option<(Vec<u8>, usize)> {
        let mut q = p;
        if self.pep {
            while at(l, q).is_ascii_alphabetic() {
                q += 1;
            }
        } else {
            loop {
                let c = at(l, q);
                if matches!(c, 0 | b' ' | b'\t' | b'[' | b'\n' | b')') {
                    break;
                }
                if !dna_letter(c) {
                    return None;
                }
                q += 1;
            }
        }
        (q > p).then(|| (l[p..q].to_ascii_uppercase(), q))
    }

    fn compl_pat(&self, l: &[u8], p: usize) -> Option<(Kind, usize)> {
        if self.pep {
            return None;
        }
        let (rs, p) = if at(l, p) == b'~' {
            (None, p + 1)
        } else {
            let (rs, p) = rule_id(l, p)?;
            if at(l, p) != b'~' {
                return None;
            }
            (Some(rs), p + 1)
        };
        let (n, p) = name_id(l, p)?;
        let (e, p) = opt_errs(l, p);
        Some((Kind::Compl(rs, n, e), p))
    }

    fn unit(&mut self, l: &[u8], mut p: usize) -> Option<(Unit, usize)> {
        while let Some(p1) = self.rule_set(l, p) {
            p = ws(l, p1);
        }
        let slot = self.pup;
        self.pup += 1;
        let off = p;
        let mut def = None;
        if let Some((n, p1)) = name_assgn(l, p) {
            if self.names[n as usize].is_some() {
                // the engine does not release the slot
                self.dups.push((n, off));
                self.fail(off, "this name is already defined");
                return None;
            }
            self.names[n as usize] = Some(slot);
            def = Some(n);
            p = p1;
        }
        p = ws(l, p);
        let pep = self.pep;
        let r = if at(l, p) == b'^' {
            Some((Kind::Start, p + 1))
        } else if at(l, p) == b'$' {
            Some((Kind::End, p + 1))
        } else if let Some((a, b, q)) = range_pat(l, p) {
            Some((Kind::Range(a, b), q))
        } else {
            any_pat(l, p, pep)
                .or_else(|| llim_pat(l, p))
                .or_else(|| {
                    if at(l, p) != b'<' {
                        return None;
                    }
                    let (n, q) = name_id(l, p + 1)?;
                    let (e, q) = if pep {
                        (Err3::default(), q)
                    } else {
                        opt_errs(l, q)
                    };
                    Some((Kind::Inv(n, e), q))
                })
                .or_else(|| {
                    let (n, q) = name_id(l, p)?;
                    let (e, q) = opt_errs(l, q);
                    Some((Kind::Repeat(n, e), q))
                })
                .or_else(|| {
                    let (w, q) = self.word(l, p)?;
                    let (e, q) = opt_errs(l, q);
                    Some((Kind::Word(w, e), q))
                })
                .or_else(|| self.compl_pat(l, p))
                .or_else(|| wt_pat(l, p))
        };
        match r {
            Some((kind, end)) => Some((
                Unit {
                    kind,
                    def,
                    slot,
                    off,
                    end,
                },
                end,
            )),
            None => {
                self.pup = slot;
                if !matches!(at(l, p), 0 | b'|' | b')') {
                    self.fail(p, "no pattern unit can be read here");
                }
                None
            }
        }
    }

    fn or_pat(&mut self, l: &[u8], p: usize) -> Option<(Unit, usize)> {
        let slot = self.pup;
        self.pup += 1;
        let r = self.or_body(l, p);
        match r {
            Some((a, b, bar, close)) => Some((
                Unit {
                    kind: Kind::Or(a, b, bar, close),
                    def: None,
                    slot,
                    off: p,
                    end: close + 1,
                },
                close + 1,
            )),
            None => {
                self.pup = slot;
                None
            }
        }
    }

    fn or_body(&mut self, l: &[u8], p: usize) -> Option<(Vec<Unit>, Vec<Unit>, usize, usize)> {
        if at(l, p) != b'(' {
            return None;
        }
        let (a, p) = self.cons_list(l, ws(l, p + 1))?;
        let bar = ws(l, p);
        if at(l, bar) != b'|' {
            self.fail(bar, "`|` is missing here (an alternative is `( A | B )`)");
            return None;
        }
        let (b, p) = self.cons_list(l, ws(l, bar + 1))?;
        let close = ws(l, p);
        if at(l, close) != b')' {
            self.fail(close, "`)` is missing here (an alternative is `( A | B )`)");
            return None;
        }
        Some((a, b, bar, close))
    }

    fn list(&mut self, l: &[u8], p: usize) -> Option<(Vec<Unit>, usize)> {
        let mut p = ws(l, p);
        let mut units = Vec::new();
        while let Some((u, p1)) = self.unit(l, p) {
            units.push(u);
            p = ws(l, p1);
        }
        (!units.is_empty()).then_some((units, p))
    }

    fn cons(&mut self, l: &[u8], mut p: usize) -> Option<(Vec<Unit>, usize)> {
        while let Some(p1) = self.rule_set(l, p) {
            p = ws(l, p1);
        }
        if let Some((u, p1)) = self.or_pat(l, p) {
            return Some((vec![u], p1));
        }
        self.list(l, p)
    }

    fn cons_list(&mut self, l: &[u8], p: usize) -> Option<(Vec<Unit>, usize)> {
        let mut p = ws(l, p);
        let mut units = Vec::new();
        while let Some((us, p1)) = self.cons(l, p) {
            units.extend(us);
            p = ws(l, p1);
        }
        (!units.is_empty()).then_some((units, p))
    }
}

// ---------------------------------------------------------------------
// what the engine computes from the parsed pattern

/// A parsed pattern with its units in slot order.
struct Pat<'a> {
    l: &'a [u8],
    ps: &'a Parsed,
    /// every unit, pre-order (alternatives included)
    flat: Vec<&'a Unit>,
}

fn flatten<'a>(units: &'a [Unit], out: &mut Vec<&'a Unit>) {
    for u in units {
        out.push(u);
        if let Kind::Or(a, b, ..) = &u.kind {
            flatten(a, out);
            flatten(b, out);
        }
    }
}

/// Uses of names (unit, name) that may read a match of an earlier attempt.
type Early<'a> = Vec<(&'a Unit, i32)>;

impl<'a> Pat<'a> {
    fn new(l: &'a [u8], ps: &'a Parsed) -> Pat<'a> {
        let mut flat = Vec::new();
        flatten(&ps.units, &mut flat);
        Pat { l, ps, flat }
    }

    fn slot_unit(&self, s: u32) -> Option<&'a Unit> {
        match self.flat.get(s as usize) {
            Some(u) if u.slot == s => Some(u),
            _ => self.flat.iter().find(|u| u.slot == s).copied(),
        }
    }

    /// The unit that defines name `n`.
    fn named(&self, n: i32) -> Option<&'a Unit> {
        let s = (*self.ps.names.get(usize::try_from(n).ok()?)?)?;
        self.slot_unit(s)
    }

    /// The name chain from `n` ends at a unit that is not a reference
    /// (false: an undefined name or a loop).
    fn name_ok(&self, n: i32) -> bool {
        let mut seen = [false; N_NAMES];
        let mut n = n;
        loop {
            if seen[n as usize] {
                return false;
            }
            seen[n as usize] = true;
            match self.named(n).map(|u| &u.kind) {
                None => return false,
                Some(Kind::Compl(_, m, _) | Kind::Repeat(m, _) | Kind::Inv(m, _)) => n = *m,
                Some(_) => return true,
            }
        }
    }

    /// The name chain from `n` comes back to a name seen before.
    fn loops(&self, n: i32) -> bool {
        let mut seen = [false; N_NAMES];
        let mut n = n;
        loop {
            if seen[n as usize] {
                return true;
            }
            seen[n as usize] = true;
            match self.named(n).map(|u| &u.kind) {
                Some(Kind::Compl(_, m, _) | Kind::Repeat(m, _) | Kind::Inv(m, _)) => n = *m,
                _ => return false,
            }
        }
    }

    /// `refs_ok`: every name and rule set used is defined, no loops.
    fn refs_ok(&self) -> bool {
        self.flat.iter().all(|u| match &u.kind {
            Kind::Compl(rs, n, _) => {
                self.name_ok(*n) && rs.is_none_or(|r| self.ps.rule_def[r as usize])
            }
            Kind::Repeat(n, _) | Kind::Inv(n, _) => self.name_ok(*n),
            Kind::Length(ns, _) => ns
                .iter()
                .all(|n| n.is_some_and(|n| self.named(n).is_some())),
            _ => true,
        })
    }

    /// `max_mat` (wrapping, as the engine; an alternative counts only the
    /// first unit of each branch).
    fn max_mat(&self, u: &Unit, depth: u32) -> i32 {
        match &u.kind {
            // the engine adds the 2nd number (deletions), as the C code
            Kind::Word(w, e) => (w.len() as i32).wrapping_add(e.del),
            Kind::Range(a, b) => a.wrapping_add(b.wrapping_sub(*a)),
            Kind::Any(..) => 1,
            Kind::Weight(_, rows, _) => rows.len() as i32,
            Kind::Compl(_, n, e) | Kind::Repeat(n, e) | Kind::Inv(n, e) => match self.named(*n) {
                Some(d) if depth < 10_000 => self.max_mat(d, depth + 1).wrapping_add(e.del),
                _ => 0,
            },
            Kind::Or(a, b, ..) => {
                let x = a.first().map_or(0, |u| self.max_mat(u, depth + 1));
                let y = b.first().map_or(0, |u| self.max_mat(u, depth + 1));
                x.max(y)
            }
            _ => 0,
        }
    }

    /// The value `parse_cmd` returns (0: "failed to parse pattern").
    fn max_mats(&self) -> i32 {
        self.ps
            .units
            .iter()
            .fold(0i32, |s, u| s.wrapping_add(self.max_mat(u, 0)))
    }

    /// The engine accepts the pattern.
    fn accepted(&self) -> bool {
        self.ps.ok && self.refs_ok() && self.max_mats() != 0
    }

    /// `walk_assigned`: names certainly matched after `units`; uses that
    /// are not are added to `out`.
    fn walk(&self, units: &'a [Unit], mut assigned: u64, out: &mut Early<'a>) -> u64 {
        for u in units {
            let mut need = |n: i32| {
                if !(0..=MAX_NAMES).contains(&n) || assigned & (1u64 << n) == 0 {
                    out.push((u, n));
                }
            };
            match &u.kind {
                Kind::Compl(_, n, _) | Kind::Repeat(n, _) | Kind::Inv(n, _) => need(*n),
                Kind::Length(ns, _) => ns.iter().for_each(|n| need(n.unwrap_or(-1))),
                Kind::Or(a, b, ..) => {
                    let x = self.walk(a, assigned, out);
                    let y = self.walk(b, assigned, out);
                    assigned = x & y;
                }
                _ => {}
            }
            if let Some(n) = u.def {
                assigned |= 1u64 << n;
            }
        }
        assigned
    }

    fn early(&self) -> Early<'a> {
        let mut out = Vec::new();
        self.walk(&self.ps.units, 0, &mut out);
        out
    }

    /// Why the pattern cannot be split for threads (None: it can).
    fn no_split(&self, early: bool) -> Option<&'static str> {
        if early {
            Some("a unit uses a match from an earlier attempt")
        } else if matches!(self.ps.units.first().map(|u| &u.kind), Some(Kind::Or(..))) {
            Some("the first unit is an alternative")
        } else if self.flat.iter().any(|u| matches!(u.kind, Kind::Start)) {
            Some("the pattern has `^`")
        } else {
            None
        }
    }

    /// (min, max) length of a match of `u` (a reversed range matches only
    /// its first length).
    fn len_range(&self, u: &Unit, depth: u32) -> (i64, i64) {
        match &u.kind {
            Kind::Range(a, b) => (*a as i64, (*a).max(*b) as i64),
            Kind::Word(w, e) => (
                (w.len() as i64 - e.del as i64).max(0),
                w.len() as i64 + e.ins as i64,
            ),
            Kind::Compl(_, n, e) | Kind::Repeat(n, e) | Kind::Inv(n, e) => match self.named(*n) {
                Some(d) if depth < 10_000 => {
                    let (a, b) = self.len_range(d, depth + 1);
                    ((a - e.del as i64).max(0), b + e.ins as i64)
                }
                _ => (0, 0),
            },
            Kind::Weight(_, rows, _) => (rows.len() as i64, rows.len() as i64),
            Kind::Any(..) => (1, 1),
            Kind::Or(a, b, ..) => {
                let (x, y) = (self.list_len(a, depth + 1), self.list_len(b, depth + 1));
                (x.0.min(y.0), x.1.max(y.1))
            }
            _ => (0, 0),
        }
    }

    fn list_len(&self, units: &[Unit], depth: u32) -> (i64, i64) {
        units.iter().fold((0, 0), |s, u| {
            let (a, b) = self.len_range(u, depth);
            (s.0 + a, s.1 + b)
        })
    }

    /// Can the engine skip gap lengths of `range` when `n` follows it?
    fn skippable(&self, n: &Unit, range: &Unit) -> bool {
        let own = |k: i32| self.named(k).is_some_and(|d| d.slot == range.slot);
        match &n.kind {
            Kind::Word(_, e) => e.exact(),
            Kind::Compl(None, k, e) => !own(*k) && e.ins == 0 && e.del == 0,
            Kind::Repeat(k, e) => !own(*k) && e.exact(),
            _ => false,
        }
    }

    /// Wide ranges (> 1000 lengths) whose gap lengths are tried one by
    /// one, with the unit after them.
    fn slow(&self, units: &'a [Unit], after: Option<&'a Unit>, skip: bool, out: &mut Early2<'a>) {
        for (i, u) in units.iter().enumerate() {
            let next = units.get(i + 1).or(after);
            match &u.kind {
                Kind::Range(a, b) if *b as i64 - *a as i64 > 1000 => {
                    if let Some(n) = next
                        && !(skip && self.skippable(n, u))
                    {
                        out.push((u, n));
                    }
                }
                Kind::Or(a, b, ..) => {
                    self.slow(a, next, skip, out);
                    self.slow(b, next, skip, out);
                }
                _ => {}
            }
        }
    }

    /// The pattern reads the same on both strands.
    fn symmetric(&self) -> bool {
        let us = &self.ps.units;
        let n = us.len();
        if self.ps.pep || n == 0 {
            return false;
        }
        (0..n.div_ceil(2)).all(|i| self.pair_sym(&us[i], &us[n - 1 - i], i == n - 1 - i))
    }

    fn pair_sym(&self, a: &Unit, b: &Unit, mid: bool) -> bool {
        let self_sym = |u: &Unit| match &u.kind {
            Kind::Range(..) => true,
            Kind::Word(w, _) => revcomp(w) == norm(w),
            _ => false,
        };
        if mid {
            return self_sym(a);
        }
        match (&a.kind, &b.kind) {
            (Kind::Range(a1, b1), Kind::Range(a2, b2)) => a1 == a2 && b1 == b2,
            (Kind::Word(w1, e1), Kind::Word(w2, e2)) => e1 == e2 && revcomp(w1) == norm(w2),
            (_, Kind::Compl(None, n, _)) => a.def == Some(*n),
            (_, Kind::Repeat(n, _)) => a.def == Some(*n) && self_sym(a),
            _ => false,
        }
    }

    fn line_of(&self, src: &Source, off: usize) -> u32 {
        src.src_line
            .get(off)
            .or(src.src_line.last())
            .copied()
            .unwrap_or(0)
    }

    /// The text of a unit, spaces made single, cut when long.
    fn text(&self, u: &Unit) -> String {
        let end = match u.kind {
            Kind::Or(..) => u.off + 1,
            _ => u.end,
        };
        snip(&self.l[u.off.min(self.l.len())..end.min(self.l.len())], 40)
    }
}

/// (wide range, unit after it)
type Early2<'a> = Vec<(&'a Unit, &'a Unit)>;

/// Upper case, U = T.
fn norm(w: &[u8]) -> Vec<u8> {
    w.iter()
        .map(|c| match c.to_ascii_uppercase() {
            b'U' => b'T',
            c => c,
        })
        .collect()
}

/// The reverse complement (IUPAC codes too), upper case.
fn revcomp(w: &[u8]) -> Vec<u8> {
    norm(w).iter().rev().map(|&c| compl(c)).collect()
}

/// `s` with runs of white space made one space, cut to `max` characters.
fn snip(s: &[u8], max: usize) -> String {
    let mut t = String::new();
    for w in s.split(|&c| isspace(c)).filter(|w| !w.is_empty()) {
        if !t.is_empty() {
            t.push(' ');
        }
        t.push_str(&String::from_utf8_lossy(w));
    }
    if t.chars().count() > max {
        t = t.chars().take(max.saturating_sub(3)).collect::<String>() + "...";
    }
    t
}

fn plural(n: i64, w: &str) -> String {
    let s = if n == 1 {
        ""
    } else if w.ends_with('h') {
        "es"
    } else {
        "s"
    };
    format!("{n} {w}{s}")
}

fn errs_words(e: Err3) -> String {
    if e.exact() {
        "exactly".into()
    } else {
        format!(
            "with up to {}, {}, {}",
            plural(e.mis as i64, "mismatch"),
            plural(e.del as i64, "deletion"),
            plural(e.ins as i64, "insertion")
        )
    }
}

fn iupac(c: u8) -> Option<&'static str> {
    Some(match c {
        b'R' => "R = A or G",
        b'Y' => "Y = C or T",
        b'K' => "K = G or T",
        b'M' => "M = A or C",
        b'S' => "S = C or G",
        b'W' => "W = A or T",
        b'B' => "B = C, G or T",
        b'D' => "D = A, G or T",
        b'H' => "H = A, C or T",
        b'V' => "V = A, C or G",
        b'N' => "N = any letter",
        b'U' => "U = T",
        _ => return None,
    })
}

/// "a, b or c" (`last` = "or") or "a, b and c".
fn join_words(items: &[String], last: &str) -> String {
    match items.split_last() {
        Some((z, r)) if !r.is_empty() => format!("{} {last} {z}", r.join(", ")),
        _ => items.concat(),
    }
}

/// A hint for the text where reading stopped.
fn hint(l: &[u8], p: usize, pep: bool) -> Option<String> {
    let rest = &l[p.min(l.len())..];
    let tok: Vec<u8> = rest.iter().take_while(|&&c| !isspace(c)).copied().collect();
    let t = String::from_utf8_lossy(&tok).into_owned();
    let c = at(l, p);
    let digits_then = |sep: &[u8]| {
        let k = tok.iter().take_while(|c| c.is_ascii_digit()).count();
        k > 0 && tok[k..].starts_with(sep)
    };
    Some(match c {
        0 => "the pattern ends too early".into(),
        b'[' => "write [mismatches,deletions,insertions] right after a word or a name, \
                 with three numbers and no spaces: ACGT[1,0,0]"
            .into(),
        b'{' => "a weight unit is {(a,c,g,t),(a,c,g,t),...} > cutoff: rows of 4 numbers \
                 (20 or 21 with -p), no space before `,` or `)` in a row"
            .into(),
        b'~' if pep => "`~` (reverse complement) cannot be used with -p".into(),
        _ if tok.contains(&b'|') => "put a space before and after `|`".into(),
        _ if (digits_then(b"-") || digits_then(b"..")) && !digits_then(b"...") => {
            "a range is written min...max with three dots, for example 50...30000".into()
        }
        b'p' if num(l, p + 1).is_some_and(|(n, _)| !(0..=MAX_NAMES).contains(&n)) => {
            "names are p0 to p50".into()
        }
        b'r' if num(l, p + 1).is_some() => {
            if pep {
                "rule sets cannot be used with -p".into()
            } else {
                "a rule set is written r1={au,ua,gc,cg} (no spaces inside { }) and used as r1~p2"
                    .into()
            }
        }
        _ if starts(l, p, b"length") => "write length(p1+p2) < N".into(),
        _ if !pep && tok.iter().all(|c| c.is_ascii_alphabetic()) && !tok.is_empty() => {
            let bad: Vec<char> = tok
                .iter()
                .filter(|&&c| !dna_letter(c))
                .map(|&c| c as char)
                .collect();
            if bad.is_empty() {
                return None;
            }
            format!(
                "`{t}`: `{}` is not a DNA letter (use A C G T U or the IUPAC codes R Y K M S W B D H V N)",
                bad[0]
            )
        }
        _ => return None,
    })
}

// ---------------------------------------------------------------------
// the checks

#[derive(Clone, Copy, PartialEq, PartialOrd)]
enum Level {
    Error,
    Warning,
    Note,
}

struct Msg {
    level: Level,
    /// source line (0: not about one unit)
    line: u32,
    text: String,
}

struct Checker<'s, 'a> {
    src: &'s Source<'s>,
    facts: &'s Facts,
    o: &'s RunOpts,
    pat: &'a Pat<'a>,
    msgs: Vec<Msg>,
}

impl Checker<'_, '_> {
    /// Add a message about the unit at `off` (None: about the file).
    fn add(&mut self, level: Level, off: Option<usize>, text: impl Into<String>) {
        let line = off.map_or(0, |o| self.pat.line_of(self.src, o));
        let text = text.into();
        self.msgs.push(Msg { level, line, text });
    }

    /// The errors: the engine would stop with "failed to parse pattern".
    fn errors(&mut self) {
        let (pat, ps) = (self.pat, self.pat.ps);
        if !ps.ok {
            if ps.dups.is_empty() {
                let (p, t) = stop_text(pat.l, ps);
                self.add(
                    Level::Error,
                    Some(p),
                    format!("the pattern cannot be read {t}"),
                );
            }
            for &(n, off) in &ps.dups {
                let first = ps.names[n as usize].and_then(|s| pat.slot_unit(s));
                let first =
                    first.map(|u| format!(" (first at line {})", pat.line_of(self.src, u.off)));
                let t = format!("p{n} is defined a second time{}", first.unwrap_or_default());
                self.add(
                    Level::Error,
                    Some(off),
                    t + "; a name can be defined only once",
                );
            }
            return;
        }
        for u in pat.flat.iter().copied() {
            let mut names: Vec<i32> = Vec::new();
            let e = |t: String| (Level::Error, Some(u.off), t);
            let mut out = Vec::new();
            match &u.kind {
                Kind::Compl(rs, n, _) => {
                    if let Some(r) = rs.filter(|&r| !ps.rule_def[r as usize]) {
                        out.push(e(format!(
                            "rule set r{r} is used but not defined (write r{r}={{au,ua,gc,cg}} before it)"
                        )));
                    }
                    names.push(*n);
                }
                Kind::Repeat(n, _) | Kind::Inv(n, _) => names.push(*n),
                Kind::Length(ns, _) => {
                    if length_bad(ns) {
                        out.push(e("length() can have at most 50 names".into()));
                    }
                    names.extend(ns.iter().flatten());
                }
                _ => {}
            }
            for n in names {
                if pat.named(n).is_none() {
                    out.push(e(format!(
                        "p{n} is used but not defined (write p{n}=... before it)"
                    )));
                } else if !matches!(u.kind, Kind::Length(..)) && pat.loops(n) {
                    out.push(e(format!("p{n} refers to itself (as in p1=~p1)")));
                }
            }
            for (l, o, t) in out {
                self.add(l, o, t);
            }
        }
        if pat.refs_ok() && pat.max_mats() == 0 {
            let t = "the pattern can match at most 0 letters (for example `^` alone or `0...0`); \
                     the program then says \"failed to parse pattern\"";
            self.add(Level::Error, None, t);
        }
        if !self.has_error() && !(self.facts.parsed && pat.accepted()) {
            self.add(
                Level::Error,
                None,
                "the program could not parse the pattern",
            );
        }
    }

    /// Label errors (`%@` comments).
    fn labels(&mut self) {
        let mine: Vec<(u32, usize)> = (self.pat.flat.iter())
            .filter(|u| !matches!(u.kind, Kind::Or(..)))
            .map(|u| (u.slot, u.off))
            .collect();
        let offs = if self.facts.parsed {
            &self.facts.units
        } else {
            &mine
        };
        let units: Vec<(u32, u32)> = offs
            .iter()
            .map(|&(s, o)| (s, self.pat.line_of(self.src, o)))
            .collect();
        let r = fmt::parse_labels(self.src.comments).and_then(|lb| lb.assign(&units).map(|_| ()));
        if let Err(e) = r {
            let (line, text) = split_line(&e);
            let text = format!("label: {text}");
            self.msgs.push(Msg {
                level: Level::Error,
                line,
                text,
            });
        }
    }

    /// Checks of a pattern that parses.
    fn checks(&mut self) {
        let pat = self.pat;
        let (warn, note) = (Level::Warning, Level::Note);
        for u in pat.flat.iter().copied() {
            let (t, at) = (pat.text(u), Some(u.off));
            match &u.kind {
                Kind::Range(a, b) if a > b => self.add(warn, at, format!(
                    "`{t}`: the first number is larger, so only length {a} is matched; write {b}...{a}"
                )),
                Kind::Range(a, b) if a < b && u.def.is_some() => self.add(note, at, format!(
                    "`{t}`: the shortest length that works is used, not the longest"
                )),
                Kind::Inv(n, _) => self.add(note, at, format!(
                    "`{t}` is p{n} in reverse order, not the reverse complement; \
                     for an inverted repeat use ~p{n}"
                )),
                Kind::Or(..) => self.add(note, at,
                    "alternative `( A | B )`: branch 1 is searched over the whole sequence before branch 2",
                ),
                Kind::Weight(..) => self.add(note, at,
                    "weight unit: N and IUPAC letters in the data get average scores, \
                     but other units never match N",
                ),
                _ => {}
            }
            // a range catches N from the data; p1 / <p1 then match any letter there
            if let Kind::Repeat(n, _) | Kind::Inv(n, _) = &u.kind
                && !self.o.strict_n
                && !pat.ps.pep
                && pat
                    .named(*n)
                    .is_some_and(|d| matches!(d.kind, Kind::Range(..)))
            {
                self.add(
                    note,
                    at,
                    format!(
                        "`{t}`: an N that p{n} caught from the data matches any letter here \
                     (false hits at assembly gaps); add --strict-n for genomes with gaps"
                    ),
                );
            }
            let (e, fixed) = match &u.kind {
                Kind::Word(w, e) => (*e, Some(w.len() as i64)),
                Kind::Compl(_, n, e) | Kind::Repeat(n, e) | Kind::Inv(n, e) => {
                    let r = pat.named(*n).map(|d| pat.len_range(d, 0));
                    (*e, r.filter(|r| r.0 == r.1).map(|r| r.0))
                }
                _ => continue,
            };
            if let Some(len) = fixed.filter(|&l| l > 0 && e.mis as i64 >= l) {
                let m = plural(e.mis as i64, "mismatch");
                self.add(
                    warn,
                    at,
                    format!("`{t}` allows {m} in {len} letters, so it matches any {len} letters"),
                );
            }
            if e.ins > 0 && e.del > 0 {
                self.add(
                    note,
                    at,
                    format!(
                        "`{t}` has deletions and insertions: the matcher takes one greedy alignment \
                     and can miss some matches"
                    ),
                );
            }
        }
        let early = pat.early();
        for (u, n) in &early {
            let t = pat.text(u);
            let later = pat.named(*n).is_some_and(|d| d.off > u.off);
            let what = match (&u.kind, later) {
                (Kind::Length(..), _) => format!("`{t}` reads p{n} before p{n} is matched"),
                (_, true) => format!("`{t}` uses p{n} before it is defined"),
                _ => format!("`{t}`: p{n} is defined only in one branch of an alternative"),
            };
            let t = what + "; it then uses a match from an earlier attempt or sequence";
            self.add(warn, Some(u.off), t);
        }
        let parsed = self.facts.parsed;
        let is_early = if parsed {
            self.facts.uses_earlier_state
        } else {
            !early.is_empty()
        };
        if is_early && early.is_empty() {
            self.add(
                warn,
                None,
                "a unit uses a match from an earlier attempt or sequence",
            );
        }
        // wide gaps that are slow
        let skip = !pat.ps.pep && !is_early;
        let mut slow = Vec::new();
        pat.slow(&pat.ps.units, None, skip, &mut slow);
        for (r, n) in slow {
            let why = if skip {
                format!(
                    "`{}` after it cannot be used to skip gap lengths",
                    pat.text(n)
                )
            } else if pat.ps.pep {
                "gap lengths are not skipped with -p".into()
            } else {
                "gap lengths are not skipped when a unit uses an earlier match".into()
            };
            let t = format!(
                "slow: each gap length of `{}` is tried ({why})",
                pat.text(r)
            );
            self.add(note, Some(r.off), t);
        }
        // names never used
        let mut used = [false; N_NAMES];
        for u in &pat.flat {
            match &u.kind {
                Kind::Compl(_, n, _) | Kind::Repeat(n, _) | Kind::Inv(n, _) => {
                    used[*n as usize] = true
                }
                Kind::Length(ns, _) => ns.iter().flatten().for_each(|&n| used[n as usize] = true),
                _ => {}
            }
        }
        for u in &pat.flat {
            if let Some(n) = u.def.filter(|&n| !used[n as usize]) {
                self.add(note, Some(u.off), format!("p{n} is defined but not used"));
            }
        }
        let why = pat.no_split(is_early);
        if !(if parsed {
            self.facts.can_split
        } else {
            why.is_none()
        }) {
            let why = why.unwrap_or("engine");
            self.add(
                note,
                None,
                format!(
                    "the pattern cannot be split for threads ({why}): \
                 with -t, each long sequence uses one thread"
                ),
            );
        }
        if pat.symmetric() && !(self.o.complements && self.o.dedup) {
            let t =
                "the pattern reads the same on both strands: with -c, every element is found twice";
            if !self.o.complements {
                self.add(note, None, format!("{t}; use --dedup"));
            } else {
                let fix = if self.o.format { "" } else { " with --format" };
                self.add(warn, None, format!("{t}; use --dedup{fix}"));
            }
        }
    }

    /// Checks of the file and the options.
    fn file_checks(&mut self) {
        if self.src.truncated {
            let t = "the pattern file is longer than 31 999 bytes; the text after that is not read";
            self.add(Level::Warning, None, t);
        }
        if let Some(v) = self.o.o_value.as_ref().filter(|v| v.first() == Some(&b'-')) {
            let v = String::from_utf8_lossy(v);
            let more = if v == "-c" {
                ", so the reverse strand is not searched"
            } else {
                ""
            };
            self.add(Level::Warning, None, format!(
                "-o takes a value: `{v}` is read as the value of -o, not as an option{more}; write -o 1"
            ));
        }
    }

    fn run(&mut self, all: bool) {
        self.file_checks();
        self.errors();
        if self.pat.ps.ok {
            self.labels();
        }
        if all && !self.has_error() {
            self.checks();
        }
        self.msgs
            .sort_by(|a, b| (a.line, a.level).partial_cmp(&(b.line, b.level)).unwrap());
    }

    fn has_error(&self) -> bool {
        self.msgs.iter().any(|m| m.level == Level::Error)
    }
}

/// Where reading stopped: (offset, "at `text` (byte N): why; hint").
fn stop_text(l: &[u8], ps: &Parsed) -> (usize, String) {
    let (p, why) = ps.stop.unwrap_or((0, "the pattern cannot be read"));
    let here = if p >= l.len() {
        "at the end of the pattern".to_string()
    } else {
        format!("at `{}`", snip(&l[p..], 30))
    };
    let h = hint(l, p, ps.pep)
        .map(|h| format!("; {h}"))
        .unwrap_or_default();
    (p, format!("{here} (byte {p}): {why}{h}"))
}

/// "pattern line N: text" -> (N, text)
fn split_line(e: &str) -> (u32, String) {
    if let Some(r) = e.strip_prefix("pattern line ")
        && let Some((n, t)) = r.split_once(": ")
        && let Ok(n) = n.parse()
    {
        return (n, t.to_string());
    }
    (0, e.to_string())
}

fn level_word(l: Level) -> &'static str {
    match l {
        Level::Error => "error",
        Level::Warning => "warning",
        Level::Note => "note",
    }
}

/// `--lint`: the report (for stdout) and whether it has an error.
pub fn lint(src: &Source, facts: &Facts, o: &RunOpts) -> (String, bool) {
    let ps = Parsed::new(src.line, o.protein);
    let pat = Pat::new(src.line, &ps);
    let mut c = Checker {
        src,
        facts,
        o,
        pat: &pat,
        msgs: Vec::new(),
    };
    c.run(true);
    let file = String::from_utf8_lossy(src.file);
    let mut out = String::new();
    let mut count = [0usize; 3];
    for m in &c.msgs {
        count[m.level as usize] += 1;
        if m.line > 0 {
            let _ = write!(out, "{file}:{}: ", m.line);
        } else {
            let _ = write!(out, "{file}: ");
        }
        let _ = writeln!(out, "{}: {}", level_word(m.level), m.text);
    }
    let _ = writeln!(
        out,
        "{} errors, {} warnings, {} notes",
        count[0], count[1], count[2]
    );
    (out, c.has_error())
}

// ---------------------------------------------------------------------
// explain

/// One row of the unit table: (unit number, offset, text, meaning).
type Row = (String, usize, String, String);

/// The first word of the `%@ TYPE` label on source line `line`.
fn label_on(src: &Source, line: u32) -> Option<String> {
    src.comments.iter().find_map(|(l, t)| {
        let r = t
            .strip_prefix(b"@")
            .filter(|r| r.first().is_some_and(u8::is_ascii_whitespace))?;
        let w = r.split(u8::is_ascii_whitespace).find(|w| !w.is_empty());
        w.filter(|_| *l == line)
            .map(|w| String::from_utf8_lossy(w).into_owned())
    })
}

impl Pat<'_> {
    /// The rows of the unit table; `n` counts the units.
    fn rows(&self, src: &Source, units: &[Unit], depth: usize, n: &mut usize, rows: &mut Vec<Row>) {
        let pad = "  ".repeat(depth);
        for u in units {
            if let Kind::Or(a, b, bar, close) = &u.kind {
                let s = String::new;
                rows.push((
                    s(),
                    u.off,
                    format!("{pad}("),
                    "one of two branches (branch 1 is searched first):".into(),
                ));
                self.rows(src, a, depth + 1, n, rows);
                rows.push((s(), *bar, format!("{pad}|"), "or branch 2:".into()));
                self.rows(src, b, depth + 1, n, rows);
                rows.push((s(), *close, format!("{pad})"), "end of the branches".into()));
                continue;
            }
            *n += 1;
            let mut m = self.meaning(u);
            if let Some(lb) = label_on(src, self.line_of(src, u.off)) {
                let _ = write!(m, "   [label: {lb}]");
            }
            rows.push((n.to_string(), u.off, format!("{pad}{}", self.text(u)), m));
        }
    }

    /// The unit in words.
    fn meaning(&self, u: &Unit) -> String {
        let call = u
            .def
            .map(|n| format!("; call them p{n}"))
            .unwrap_or_default();
        let m = match &u.kind {
            Kind::Start => "the start of the sequence".into(),
            Kind::End => "the end of the sequence".into(),
            Kind::Range(a, b) => {
                let (what, s) = if u.def.is_some() {
                    ("any", "length")
                } else {
                    ("a gap of", "gap")
                };
                let one = plural(*a as i64, "letter");
                match a.cmp(b) {
                    Ordering::Equal => format!("{what} {one}"),
                    Ordering::Less => {
                        format!("{what} {a} to {b} letters (the shortest {s} that works is used)")
                    }
                    Ordering::Greater => {
                        format!("{what} {one} (min > max: only length {a} is matched)")
                    }
                }
            }
            Kind::Word(w, e) => {
                let mut codes: Vec<&str> = Vec::new();
                for s in w.iter().filter_map(|&c| iupac(c)).filter(|_| !self.ps.pep) {
                    if !codes.contains(&s) {
                        codes.push(s);
                    }
                }
                let codes = if codes.is_empty() {
                    String::new()
                } else {
                    format!(" ({})", codes.join(", "))
                };
                format!(
                    "the letters {}, {}{codes}",
                    String::from_utf8_lossy(w),
                    errs_words(*e)
                )
            }
            Kind::Repeat(n, e) => format!("the same letters as p{n}, {}", errs_words(*e)),
            Kind::Compl(None, n, e) => {
                format!("the reverse complement of p{n}, {}", errs_words(*e))
            }
            Kind::Compl(Some(r), n, e) => {
                format!(
                    "p{n} read backwards, each letter paired with it by rule set r{r}, {}",
                    errs_words(*e)
                )
            }
            Kind::Inv(n, e) => {
                format!(
                    "the letters of p{n} in reverse order (not complemented), {}",
                    errs_words(*e)
                )
            }
            Kind::Weight(max, rows, cut) => {
                let cols = match rows[0].len() {
                    4 => "4 scores: A, C, G, T".to_string(),
                    k => format!("{k} scores, one per amino acid"),
                };
                let lim = max
                    .map(|m| format!(" and less than {m}"))
                    .unwrap_or_default();
                let k = rows.len() as i64;
                format!(
                    "{} scored with a weight matrix ({} of {cols}); the total score must be more than \
                     {cut}{lim}; N and IUPAC codes in the data get average scores",
                    plural(k, "letter"),
                    plural(k, "row")
                )
            }
            Kind::Length(ns, bound) => {
                let names: Vec<String> = ns
                    .iter()
                    .map(|n| n.map_or("?".into(), |n| format!("p{n}")))
                    .collect();
                let all = if names.len() == 1 {
                    "the length of"
                } else {
                    "the total length of"
                };
                format!(
                    "no letters; checks that {all} {} is less than {bound}",
                    join_words(&names, "and")
                )
            }
            Kind::Any(set, not) => {
                let items: Vec<String> = set.iter().map(|&c| (c as char).to_string()).collect();
                let what = if *not {
                    "one letter that is not"
                } else {
                    "one letter:"
                };
                format!("{what} {}", join_words(&items, "or"))
            }
            Kind::Or(..) => String::new(),
        };
        m + &call
    }
}

/// `--explain`: the description (for stdout) and whether it has an error.
pub fn explain(src: &Source, facts: &Facts, o: &RunOpts) -> (String, bool) {
    let ps = Parsed::new(src.line, o.protein);
    let pat = Pat::new(src.line, &ps);
    let mut out = String::new();
    let mode = if o.protein { "protein, -p" } else { "DNA" };
    let _ = writeln!(
        out,
        "pattern file: {} ({mode})",
        String::from_utf8_lossy(src.file)
    );
    let _ = writeln!(out, "joined pattern: {}", snip(src.line, 300));
    let mut rows = Vec::new();
    pat.rows(src, &ps.units, 0, &mut 0, &mut rows);
    let line = |off| match pat.line_of(src, off) {
        0 => String::new(),
        l => l.to_string(),
    };
    let wn = rows.iter().map(|r| r.0.len()).max().unwrap_or(0).max(4);
    let wl = rows
        .iter()
        .map(|r| line(r.1).len())
        .max()
        .unwrap_or(0)
        .max(4);
    let wt = rows
        .iter()
        .map(|r| r.2.chars().count())
        .max()
        .unwrap_or(0)
        .max(4);
    let _ = writeln!(
        out,
        "{:wn$}  {:wl$}  {:wt$}  meaning",
        "unit", "line", "text"
    );
    for (num, off, text, meaning) in &rows {
        let _ = writeln!(out, "{num:wn$}  {:wl$}  {text:wt$}  {meaning}", line(*off));
    }
    for rs in &ps.rules {
        let pairs: Vec<String> = rs
            .pairs
            .iter()
            .map(|p| format!("{}-{}", p[0] as char, p[1] as char))
            .collect();
        let l = pat.line_of(src, rs.off);
        let _ = writeln!(
            out,
            "rule set r{} (line {l}): letter pairs {}",
            rs.id,
            pairs.join(", ")
        );
    }
    if let Some((_, t)) = src
        .comments
        .iter()
        .find(|(_, t)| t.starts_with(b"@element"))
    {
        let _ = writeln!(out, "labels: {}", snip(&t[1..], 200));
    }
    let mut c = Checker {
        src,
        facts,
        o,
        pat: &pat,
        msgs: Vec::new(),
    };
    c.run(false);
    if c.has_error() {
        for m in c.msgs.iter().filter(|m| m.level == Level::Error) {
            let at = if m.line > 0 {
                format!(" (line {})", m.line)
            } else {
                String::new()
            };
            let _ = writeln!(out, "error{at}: {}", m.text);
        }
        return (out, true);
    }
    let (lo, hi) = pat.list_len(&ps.units, 0);
    let len = if lo == hi {
        plural(lo, "letter")
    } else {
        format!("{lo} to {hi} letters")
    };
    let _ = writeln!(out, "hit length: {len}");
    let start = if matches!(ps.units.first().map(|u| &u.kind), Some(Kind::Start)) {
        "the hit must begin at the start of the sequence (`^`)"
    } else {
        "the first unit is tried at each position, from left to right"
    };
    let _ = writeln!(
        out,
        "search: {start}; after a hit, the search goes on after its end (with -o: at the next position)"
    );
    if !o.protein {
        let strand = if o.complements {
            "both strands are searched (-c)"
        } else {
            "only the given strand is searched (-c: both strands)"
        };
        let _ = writeln!(out, "strands: {strand}");
    }
    let early = if facts.parsed {
        facts.uses_earlier_state
    } else {
        !pat.early().is_empty()
    };
    let why = pat.no_split(early);
    let _ = if facts.parsed && facts.can_split || !facts.parsed && why.is_none() {
        writeln!(
            out,
            "threads: with -t N, pieces of long sequences are searched at the same time"
        )
    } else {
        let why = why.unwrap_or("engine");
        writeln!(
            out,
            "threads: with -t N, each sequence is searched by one thread ({why})"
        )
    };
    (out, false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{DNA, Engine, PEPTIDE};

    /// The pattern file read as main.rs reads it: (line, source lines,
    /// comments, truncated).
    type Read = (Vec<u8>, Vec<u32>, Vec<(u32, Vec<u8>)>, bool);

    fn read_pat(b: &[u8]) -> Read {
        let (mut line, mut src, mut com) = (Vec::new(), Vec::new(), Vec::new());
        let mut n = 1u32;
        let mut k = 0;
        while k < b.len() && line.len() < 31_999 {
            let c = b[k];
            if c == b'%' {
                let e = b[k..]
                    .iter()
                    .position(|&c| c == b'\n')
                    .map_or(b.len(), |e| k + e);
                let mut t = b[k + 1..e].to_vec();
                if t.last() == Some(&b'\r') {
                    t.pop();
                }
                com.push((n, t));
                k = e;
                if k < b.len() {
                    n += 1;
                    k += 1;
                    line.push(b' ');
                    src.push(n);
                }
            } else if c == b'\n' || c == b'\r' {
                n += (c == b'\n') as u32;
                line.push(b' ');
                src.push(n);
                k += 1;
            } else {
                line.push(c);
                src.push(n);
                k += 1;
            }
        }
        if let Some(z) = line.iter().position(|&c| c == 0) {
            line.truncate(z);
        }
        (line, src, com, k < b.len())
    }

    fn opts(c: bool) -> RunOpts {
        RunOpts {
            protein: false,
            complements: c,
            dedup: false,
            format: false,
            o_value: None,
            strict_n: false,
        }
    }

    fn facts(line: &[u8], pep: bool) -> Facts {
        let mut e = Engine::new();
        let ok = e.parse_cmd(line, if pep { PEPTIDE } else { DNA }) != 0;
        Facts {
            parsed: ok,
            uses_earlier_state: ok && e.uses_earlier_state(),
            can_split: ok && e.can_split(),
            units: if ok { e.unit_offsets() } else { Vec::new() },
        }
    }

    fn run(
        text: &str,
        o: &RunOpts,
        f: fn(&Source, &Facts, &RunOpts) -> (String, bool),
    ) -> (String, bool) {
        let (line, sl, com, tr) = read_pat(text.as_bytes());
        let src = Source {
            file: b"x.pat",
            line: &line,
            src_line: &sl,
            comments: &com,
            truncated: tr,
        };
        f(&src, &facts(&line, o.protein), o)
    }

    fn lint_o(text: &str, o: &RunOpts) -> (String, bool) {
        run(text, o, lint)
    }

    fn lint_s(text: &str) -> String {
        lint_o(text, &opts(false)).0
    }

    fn tir_files() -> Vec<(String, Vec<u8>)> {
        let mut v = Vec::new();
        for d in ["tir_scan_patterns", "tir_scan_patterns/labelled"] {
            let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(d);
            for e in std::fs::read_dir(dir).unwrap() {
                let p = e.unwrap().path();
                if p.extension().is_some_and(|x| x == "pat") {
                    v.push((p.display().to_string(), std::fs::read(&p).unwrap()));
                }
            }
        }
        v.sort();
        assert_eq!(v.len(), 44);
        v
    }

    /// A small random generator (xorshift).
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn int(&mut self, a: i64, b: i64) -> i64 {
            a + (self.next() % (b - a + 1) as u64) as i64
        }
        fn p(&mut self, x: f64) -> bool {
            (self.next() % 1_000_000) as f64 / 1e6 < x
        }
        fn pick<'a>(&mut self, s: &[&'a str]) -> &'a str {
            s[self.next() as usize % s.len()]
        }
        fn ch(&mut self, s: &str) -> char {
            let b = s.as_bytes();
            b[self.next() as usize % b.len()] as char
        }
    }

    /// Random patterns, after `Pat` in tests/fuzz_compare.py.
    struct Gen<'r> {
        r: &'r mut Rng,
        pep: bool,
        defined: Vec<i64>,
        next_name: i64,
        rules: Vec<i64>,
    }

    const AA: &str = "ACDEFGHIKLMNPQRSTVWY";

    impl Gen<'_> {
        fn errs(&mut self) -> String {
            if self.r.p(0.6) {
                return String::new();
            }
            let (m, i, d) = (self.r.int(0, 2), self.r.int(0, 2), self.r.int(0, 2));
            format!("[{m},{i},{d}]")
        }
        fn word(&mut self, n: i64) -> String {
            let set = if self.pep {
                "ACDEFGHIKLMNPQRSTVWYXacd"
            } else if self.r.p(0.3) {
                "acgtuACGTUnNrRyYmMkKsSwWbBdDhHvV"
            } else {
                "acgtACGT"
            };
            (0..n).map(|_| self.r.ch(set)).collect()
        }
        fn unit(&mut self, depth: u32) -> String {
            let mut kinds = vec![
                "range", "range", "exact", "sim", "repeat", "inv", "weight", "llim", "or", "start",
                "end",
            ];
            if self.pep {
                kinds.extend(["any", "notany", "any"]);
            } else {
                kinds.extend(["compl", "compl", "rcompl"]);
            }
            let k = self.r.pick(&kinds);
            let mut name = String::new();
            if matches!(k, "range" | "exact" | "sim" | "weight" | "any" | "notany")
                && self.r.p(0.5)
                && self.next_name < 8
            {
                name = format!("p{}=", self.next_name);
                self.defined.push(self.next_name);
                self.next_name += 1;
            }
            let r = &mut *self.r;
            match k {
                "range" => {
                    let a = r.int(0, 6);
                    let b = if r.p(0.1) {
                        r.int(0, 6)
                    } else {
                        a + r.int(0, 6)
                    };
                    let sp = if r.p(0.1) { " " } else { "" };
                    format!("{name}{a}{sp}...{sp}{b}")
                }
                "exact" => {
                    let n = r.int(1, 4);
                    name + &self.word(n)
                }
                "sim" => {
                    let n = r.int(1, 6);
                    let e = format!("[{},{},{}]", r.int(0, 2), r.int(0, 2), r.int(0, 2));
                    name + &self.word(n) + &e
                }
                "repeat" | "inv" | "compl" | "rcompl" | "llim" => {
                    if self.defined.is_empty() || self.r.p(0.05) {
                        let a = self.r.int(1, 5);
                        let b = a + self.r.int(0, 3);
                        let s = format!("p{}={a}...{b}", self.next_name);
                        self.defined.push(self.next_name);
                        self.next_name += 1;
                        return s;
                    }
                    let n = self.defined[self.r.next() as usize % self.defined.len()];
                    let n = if self.r.p(0.05) { self.r.int(0, 9) } else { n };
                    let e = self.errs();
                    match k {
                        "repeat" => format!("p{n}{e}"),
                        "inv" => format!("<p{n}{e}"),
                        "compl" => format!("~p{n}{e}"),
                        "rcompl" => {
                            let rn = self.r.int(1, 3);
                            if !self.rules.contains(&rn) && self.r.p(0.95) {
                                self.rules.push(rn);
                            }
                            format!("r{rn}~p{n}{e}")
                        }
                        _ => {
                            let c = self.r.int(1, 3);
                            let names: Vec<String> = (0..c)
                                .map(|_| {
                                    format!(
                                        "p{}",
                                        self.defined[self.r.next() as usize % self.defined.len()]
                                    )
                                })
                                .collect();
                            format!("length({}) < {}", names.join("+"), self.r.int(1, 15))
                        }
                    }
                }
                "weight" => {
                    let tup = if !self.pep || r.p(0.3) { 4 } else { 20 };
                    let rows = r.int(1, 4);
                    let vecs: Vec<String> = (0..rows)
                        .map(|_| {
                            let v: Vec<String> =
                                (0..tup).map(|_| r.int(0, 100).to_string()).collect();
                            format!("({})", v.join(","))
                        })
                        .collect();
                    let body = format!("{{{}}}", vecs.join(","));
                    let cut = r.int(0, 60 * rows);
                    if r.p(0.3) {
                        format!("{name}{} > {body} > {cut}", cut + r.int(10, 200))
                    } else {
                        format!("{name}{body} > {cut}")
                    }
                }
                "or" if depth <= 1 => {
                    let a = self.list(depth + 1, 1, 2);
                    let b = self.list(depth + 1, 1, 2);
                    format!("({a} | {b})")
                }
                "start" => "^".into(),
                "end" => "$".into(),
                "any" | "notany" => {
                    let n = r.int(1, 4);
                    let s: String = (0..n).map(|_| r.ch(AA)).collect();
                    format!("{name}{k}({s})")
                }
                _ => "1...2".into(),
            }
        }
        fn list(&mut self, depth: u32, lo: i64, hi: i64) -> String {
            let n = self.r.int(lo, hi);
            (0..n)
                .map(|_| self.unit(depth))
                .collect::<Vec<_>>()
                .join(" ")
        }
        fn build(&mut self) -> String {
            let body = self.list(0, 1, 5);
            let mut txt = String::new();
            for rn in self.rules.clone() {
                let pairs = ["au", "ua", "gc", "cg", "gu", "ug", "ga", "ag", "AT", "TA"];
                let k = self.r.int(1, 4);
                let ps: Vec<&str> = (0..k).map(|_| self.r.pick(&pairs)).collect();
                txt += &format!("r{rn}={{{}}} ", ps.join(","));
            }
            txt += &body;
            if self.r.p(0.3) {
                // random damage
                let marks = [
                    "(", ")", "[", "]", "|", ",", "~", "x", "p9", "..", "{", "p1=", " ", "-",
                    "p51", "+",
                ];
                for _ in 0..self.r.int(1, 2) {
                    let i = self.r.int(0, txt.len() as i64) as usize;
                    let m = self.r.pick(&marks);
                    txt.insert_str(i, m);
                }
            }
            txt
        }
    }

    /// Random strings of pattern pieces.
    fn soup(r: &mut Rng, pep: bool) -> String {
        let toks = [
            "(",
            ")",
            "|",
            "p1=",
            "p2=",
            "~p1",
            "r1~p2",
            "r1={au,ua}",
            "r2={gc}",
            "1...5",
            "5...1",
            "0...0",
            "ACGT",
            "ac[1,0,0]",
            "length(p1+p2) < 5",
            "length(p1+) < 3",
            "{(1,2,3,4)} > 2",
            "9 > {(1,2,3,4),(4,3,2,1)} > 2",
            "^",
            "$",
            "<p1",
            "<p2[1,0,0]",
            "any(HK)",
            "notany(D)",
            "p1",
            "p2[0,1,1]",
            "[",
            ",",
            "...",
            "p51=",
            "p-0",
            "nnn",
            "MKV",
            "x",
            "(A | C)",
        ];
        let n = r.int(1, 8);
        let mut s = String::new();
        for _ in 0..n {
            s += r.pick(&toks);
            if !r.p(0.15) {
                s.push(if r.p(0.1) { '\n' } else { ' ' });
            }
        }
        let _ = pep;
        s
    }

    #[derive(Default)]
    struct Stats {
        n: usize,
        accepted: usize,
        skipped: usize,
    }

    fn check(line: &[u8], pep: bool, shared: &mut Engine, st: &mut Stats) {
        let ps = Parsed::new(line, pep);
        let pat = Pat::new(line, &ps);
        if ps.ok
            && pat
                .flat
                .iter()
                .any(|u| matches!(&u.kind, Kind::Length(ns, _) if length_bad(ns)))
        {
            st.skipped += 1;
            return;
        }
        // rule sets stay in an engine: use a new one when the text may have one
        let fresh = line
            .windows(2)
            .any(|w| w[0] == b'r' && (w[1].is_ascii_digit() || w[1] == b'-'));
        let mut own;
        let e = if fresh {
            own = Engine::new();
            &mut own
        } else {
            shared
        };
        let rc = e.parse_cmd(line, if pep { PEPTIDE } else { DNA });
        let show = String::from_utf8_lossy(line);
        st.n += 1;
        assert_eq!(
            rc != 0,
            pat.accepted(),
            "accept differs (pep {pep}): {show}"
        );
        if rc == 0 {
            return;
        }
        st.accepted += 1;
        assert_eq!(rc, pat.max_mats(), "max length: {show}");
        let mine: Vec<(u32, usize)> = pat
            .flat
            .iter()
            .filter(|u| !matches!(u.kind, Kind::Or(..)))
            .map(|u| (u.slot, u.off))
            .collect();
        assert_eq!(e.unit_offsets(), mine, "offsets: {show}");
        let early = !pat.early().is_empty();
        assert_eq!(e.uses_earlier_state(), early, "earlier state: {show}");
        assert_eq!(
            e.can_split(),
            pat.no_split(early).is_none(),
            "can_split: {show}"
        );
    }

    const HAND: &[&str] = &[
        "p1=3...3 CACTA[0,0,0] p2=7...7 500...15000 ~p2 TAGTG[0,0,0] p1",
        "p1=4...7 3...8 ~p1",
        "r1={au,ua,gc,cg,gu,ug,ga,ag} r2={au,ua,gc,cg} p1=2...3 0...4 p2=2...5 1...5 r1~p2 0...4 r2~p1",
        "p1=2...3 p3=0...4 p2=2...5 p4=1...5 ~p2 p5=0...4 length(p3+p4+p5) < 7 ~p1",
        "p1=175 > {(50,0,50,0),(0,0,0,75),(50,50,0,0)} > 100 0...10 p1",
        "p1=4...4 (GGCC | CCGG) ~p1[0,0,1]",
        "( p1=3...3 3...8 ~p1 | p1=5...5 4...4 ~p1 GGG )",
        "(GAGA | (GCGCA | TTCGA))",
        "(GAGA|GCGCA)",
        "^ ACGT 0...4 $",
        "^",
        "0...0",
        "p1=~p1",
        "p1=AC p1=GT",
        "~p3",
        "r4~p1",
        "p1=3...3 r4~p1",
        "p1=3...3 <p1[1,0,0]",
        "p1 = 3 ... 5 p1",
        "p1=3...3 1...2 length(p1+) < 4",
        "p0=AC p-0",
        "CACTA[0,0]",
        "CACTA [0,0,0]",
        "500-15000",
        "{(1,2,3,4) ,(1,2,3,4)} > 2",
        "{(1, 2,3,4)} > 2",
        "{( 1,2,3 ,4)} > 2",
        "{(1,2,3,4)}>2",
        "10>{(1,2,3,4)}>2",
        "AC r1={au} ( GT | CA ) r1~p1",
        "p1=AC r1={au} ( GT | CA ) r1~p1",
        "p1=4...4 ( AC | p2=GT ) p2",
        "p1=4...4 ( AC p2=GT | TT ) length(p1+p2) < 9",
        "p51=AC",
        "p50=AC p50",
        "AC xx",
        "AC\x0bGT",
        "ACu[0,0,0]",
    ];

    const HAND_PEP: &[&str] = &[
        "ACG 0...4 any(HQD) 1...3 notany(HK)",
        "p1=0...4 any(HQD) 1...3 notany(HK) p1",
        "CGXXYWG[1,0,0]",
        "ACp1",
        "p1=ACD <p1[1,0,0]",
        "p1=ACD <p1",
        "anything",
        "lengthy",
        "any(HK",
        "{(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)} > 3",
        "~p1",
        "r1={au} AC",
    ];

    fn self_check() {
        let mut shared = Engine::new();
        let mut st = Stats::default();
        for p in HAND {
            check(p.as_bytes(), false, &mut shared, &mut st);
        }
        for p in HAND_PEP {
            check(p.as_bytes(), true, &mut shared, &mut st);
        }
        for (_, b) in tir_files() {
            check(&read_pat(&b).0, false, &mut shared, &mut st);
        }
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
        for (f, pep) in [("test_dna_patterns", false), ("test_prot_patterns", true)] {
            for l in std::fs::read(root.join(f)).unwrap().split(|&c| c == b'\n') {
                check(&read_pat(l).0, pep, &mut shared, &mut st);
            }
        }
        let unhex = |s: &str| -> Vec<u8> {
            (0..s.len() / 2)
                .map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap())
                .collect()
        };
        for f in ["golden.tsv", "fixed.tsv"] {
            let t = std::fs::read_to_string(root.join(f)).unwrap();
            for l in t.lines() {
                let mut c = l.split('\t');
                let (pat, args) = (unhex(c.next().unwrap()), unhex(c.next().unwrap_or("")));
                let pep = args.split(|&c| c == 0).any(|a| a == b"-p");
                check(&read_pat(&pat).0, pep, &mut shared, &mut st);
            }
        }
        let fixed = st.n;
        let mut r = Rng(0x9e37_79b9_7f4a_7c15);
        for i in 0..4000 {
            let pep = r.p(0.25);
            let t = if i % 3 == 2 {
                soup(&mut r, pep)
            } else {
                Gen {
                    r: &mut r,
                    pep,
                    defined: Vec::new(),
                    next_name: 1,
                    rules: Vec::new(),
                }
                .build()
            };
            check(&read_pat(t.as_bytes()).0, pep, &mut shared, &mut st);
        }
        eprintln!(
            "self check: {} patterns ({} fixed, {} random), {} accepted, {} skipped (bad length())",
            st.n,
            fixed,
            st.n - fixed,
            st.accepted,
            st.skipped
        );
        assert!(st.accepted > 1000 && st.n - st.accepted > 500);
    }

    #[test]
    fn same_as_engine() {
        std::thread::Builder::new()
            .stack_size(1 << 28)
            .spawn(self_check)
            .unwrap()
            .join()
            .unwrap();
    }

    fn count(s: &str, what: &str) -> usize {
        s.lines().filter(|l| l.contains(what)).count()
    }

    #[test]
    fn tir_patterns_are_clean() {
        for (f, b) in tir_files() {
            let t = String::from_utf8(b).unwrap();
            let (s, err) = lint_o(&t, &opts(false));
            assert!(!err, "{f}: {s}");
            assert_eq!(
                count(&s, ": error:") + count(&s, ": warning:"),
                0,
                "{f}: {s}"
            );
            assert!(s.contains("reads the same on both strands"), "{f}: {s}");
            assert!(!s.contains("slow"), "{f}: {s}");
            let (s, _) = lint_o(&t, &opts(true));
            assert_eq!(count(&s, ": warning:"), 1, "{f}: {s}");
            assert!(s.contains("use --dedup with --format"), "{f}: {s}");
            let mut o = opts(true);
            o.dedup = true;
            assert_eq!(count(&lint_o(&t, &o).0, ": warning:"), 0, "{f}");
            let (e, err) = run(&t, &opts(false), explain);
            assert!(!err && e.contains("hit length:"), "{f}: {e}");
        }
    }

    /// `pat` gives a message with `level` that contains `what`.
    fn fires(pat: &str, level: &str, what: &str) {
        let s = lint_s(pat);
        assert!(
            s.lines()
                .any(|l| l.contains(&format!(": {level}: ")) && l.contains(what)),
            "{pat:?} should give {level} `{what}`:\n{s}"
        );
    }

    #[test]
    fn errors() {
        fires(
            "p1=3...3 CACTA[0,0] p1",
            "error",
            "cannot be read at `[0,0] p1`",
        );
        fires(
            "p1=3...3 CACTA[0,0] p1",
            "error",
            "[mismatches,deletions,insertions]",
        );
        fires("AC 500-15000 GT", "error", "three dots");
        fires("(GAGA|GCGCA)", "error", "space before and after `|`");
        fires("( AC | GT", "error", "`)` is missing");
        fires("AC XZ", "error", "not a DNA letter");
        fires("p1=3...3 ~p2", "error", "p2 is used but not defined");
        fires(
            "p1=3...3 r1~p1",
            "error",
            "rule set r1 is used but not defined",
        );
        fires("p1=~p1", "error", "p1 refers to itself");
        fires("p1=3...3 p2=~p3 p3=~p2", "error", "refers to itself");
        fires("^", "error", "at most 0 letters");
        fires("0...0", "error", "at most 0 letters");
        fires(
            "p1=3...3 p1=4...4",
            "error",
            "p1 is defined a second time (first at line 1)",
        );
        fires(
            &format!("p1=3...3 length(p1{}) < 4", "+p1".repeat(50)),
            "error",
            "length()",
        );
        fires(
            "p1=3...3 %@ tsd\nAC %@ tsd x\n",
            "error",
            "label: `x` is not key=value",
        );
        fires(
            "p1=3...3 AC\n%@ tsd\n",
            "error",
            "label: label on a line where no pattern unit starts",
        );
        let (s, err) = lint_o("p1=3...3 r1~p1", &opts(false));
        assert!(err && s.ends_with("1 errors, 0 warnings, 0 notes\n"), "{s}");
        assert!(s.starts_with("x.pat:1: error:"), "{s}");
        // the line of the unit
        fires("p1=3...3\nAC\n~p4\n", "error", "p4");
        assert!(lint_s("p1=3...3\nAC\n~p4\n").starts_with("x.pat:3: error"));
    }

    #[test]
    fn warnings() {
        fires("AC 5...2 GT", "warning", "only length 5 is matched");
        fires("~p1 p1=3...3", "warning", "uses p1 before it is defined");
        fires(
            "( p1=3...3 | AC ) ~p1",
            "warning",
            "defined only in one branch",
        );
        fires(
            "length(p1) < 4 p1=3...3",
            "warning",
            "reads p1 before p1 is matched",
        );
        fires("ACG[3,0,0] 0...5 GT", "warning", "matches any 3 letters");
        fires(
            "p1=3...3 0...5 ~p1[3,0,0]",
            "warning",
            "matches any 3 letters",
        );
        let mut o = opts(false);
        o.o_value = Some(b"-c".to_vec());
        let s = lint_o("AC 0...4 GT", &o).0;
        assert!(
            s.contains("warning: -o takes a value") && s.contains("reverse strand"),
            "{s}"
        );
        let big = "A ".repeat(16_100);
        let s = lint_s(&big);
        assert!(
            s.contains("warning: the pattern file is longer than 31 999 bytes"),
            "{s}"
        );
        // symmetric with -c
        let s = lint_o("p1=4...4 CA 10...20 TG p1", &opts(true)).0;
        assert!(s.contains("warning: the pattern reads the same"), "{s}");
        let s = lint_o("p1=4...4 CA 10...20 CA p1", &opts(true)).0;
        assert!(!s.contains("reads the same"), "{s}");
        let s = lint_o("ARC 10...20 GYT", &opts(true)).0;
        assert!(s.contains("reads the same"), "{s}");
    }

    #[test]
    fn notes() {
        fires(
            "p1=10...10 50...30000 p1[1,0,0]",
            "note",
            "slow: each gap length of `50...30000`",
        );
        fires("50...30000 AC[1,0,0]", "note", "slow");
        fires("r1={au} p1=5...5 50...30000 r1~p1", "note", "slow");
        fires("p1=5...5 50...30000 ~p1[0,1,0]", "note", "slow");
        fires("p1=0...5000 ~p1", "note", "slow");
        fires("p1=3...3 0...5 <p1", "note", "not the reverse complement");
        fires("p1=2...10 0...5 p1", "note", "shortest length that works");
        fires(
            "(AC | GT) 0...4 TT",
            "note",
            "branch 1 is searched over the whole sequence",
        );
        fires(
            "(AC | GT) 0...4 TT",
            "note",
            "cannot be split for threads (the first unit is an alternative)",
        );
        fires("^ AC", "note", "cannot be split");
        fires("ACGT[1,1,1]", "note", "greedy alignment");
        fires("{(1,2,3,4),(4,3,2,1)} > 3", "note", "average scores");
        fires(
            "p1=3...3 p2=AC 0...5 p1",
            "note",
            "p2 is defined but not used",
        );
        fires("AC 0...5 GT", "note", "reads the same");
        fires("p1=3...3 0...5 p1", "note", "add --strict-n");
        fires("p1=3...3 0...5 <p1", "note", "add --strict-n");
        let mut o = opts(false);
        o.strict_n = true;
        assert!(!lint_o("p1=3...3 0...5 p1", &o).0.contains("--strict-n"));
        // ~p1 never matches a caught N; a word cannot catch N
        assert!(!lint_s("p1=3...3 0...5 ~p1").contains("--strict-n"));
        assert!(!lint_s("p1=ACG 0...5 p1").contains("--strict-n"));
        // no false notes
        for p in [
            "p1=10...10 50...30000 ~p1[1,0,0] AC",
            "p1=3...3 50...30000 p1",
            "ACGT[1,1,0] 50...30000 GT 0...2",
            "AC 50...30000",
        ] {
            let s = lint_s(p);
            assert!(
                !s.contains("slow") && !s.contains("greedy") && !s.contains(": error:"),
                "{p}: {s}"
            );
        }
    }

    #[test]
    fn protein() {
        let mut o = opts(false);
        o.protein = true;
        let (s, err) = lint_o("p1=0...4 any(HQD) 1...3 notany(HK) p1", &o);
        assert!(!err && s.contains("0 errors, 0 warnings"), "{s}");
        let (s, err) = lint_o("AC ~p1", &o);
        assert!(err && s.contains("cannot be used with -p"), "{s}");
        let (e, err) = run("p1=0...4 any(HQD) 1...3 notany(HK) p1", &o, explain);
        assert!(
            !err && e.contains("one letter: H, Q or D") && e.contains("not H or K"),
            "{e}"
        );
        assert!(!e.contains("strands"), "{e}");
    }

    #[test]
    fn explain_words() {
        let o = opts(false);
        let (e, err) = run(
            "r1={au,ua} p1=3...3 RCN[1,0,0] p2=2...5 ( AC | p3=GT ) 9 > {(1,2,3,4)} > 2 \
             length(p1+p2) < 9 <p1 r1~p2 ~p2[0,1,1] p1 5...2 $",
            &o,
            explain,
        );
        assert!(!err, "{e}");
        for w in [
            "any 3 letters; call them p1",
            "the letters RCN, with up to 1 mismatch, 0 deletions, 0 insertions (R = A or G, N = any letter)",
            "any 2 to 5 letters (the shortest length that works is used); call them p2",
            "one of two branches (branch 1 is searched first):",
            "or branch 2:",
            "the letters GT, exactly; call them p3",
            "1 letter scored with a weight matrix (1 row of 4 scores: A, C, G, T); the total score must be more than 2 and less than 9",
            "the total length of p1 and p2 is less than 9",
            "the letters of p1 in reverse order (not complemented), exactly",
            "p2 read backwards, each letter paired with it by rule set r1",
            "the reverse complement of p2, with up to 0 mismatches, 1 deletion, 1 insertion",
            "the same letters as p1, exactly",
            "a gap of 5 letters (min > max: only length 5 is matched)",
            "the end of the sequence",
            "rule set r1 (line 1): letter pairs a-t, t-a",
            "threads: with -t N, pieces of long sequences",
        ] {
            assert!(e.contains(w), "missing `{w}`:\n{e}");
        }
        let (e, err) = run("p1=3...3\nCACTA[0,0] p1\n", &o, explain);
        assert!(err, "{e}");
        assert!(
            e.contains("error (line 2): the pattern cannot be read at `[0,0] p1` (byte 14)"),
            "{e}"
        );
        assert!(e.contains("1     1     p1=3...3"), "{e}");
        let (e, err) = run("p1=3...3 ~p2", &o, explain);
        assert!(
            err && e.contains("error (line 1): p2 is used but not defined"),
            "{e}"
        );
    }

    #[test]
    fn explain_labels() {
        let t = "%@element hAT_TIR_transposon Name=DTA\np1=8...8 %@ target_site_duplication\n\
                 50...300\np1 %@ target_site_duplication\n";
        let (e, err) = run(t, &opts(true), explain);
        assert!(!err, "{e}");
        assert!(
            e.contains("labels: element hAT_TIR_transposon Name=DTA"),
            "{e}"
        );
        assert_eq!(count(&e, "[label: target_site_duplication]"), 2, "{e}");
        assert!(e.contains("strands: both strands are searched (-c)"), "{e}");
    }
}
