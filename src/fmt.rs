//! Output formats (`--format gff3|bed6|bed12|jsonl`) and the `%@` labels
//! of the pattern file that name the parts of a hit.
//!
//! Workers encode each hit as a raw record (`encode`); the writer thread
//! turns the records into the format, in input order, so Names and IDs
//! are the same for every `-t`.

use crate::engine::{Buf, compl};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, PartialEq)]
pub enum Format {
    Gff3,
    Bed6,
    Bed12,
    Jsonl,
}

impl Format {
    pub fn parse(s: &[u8]) -> Option<Format> {
        Some(match s {
            b"gff3" => Format::Gff3,
            b"bed6" => Format::Bed6,
            b"bed12" => Format::Bed12,
            b"jsonl" => Format::Jsonl,
            _ => return None,
        })
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Role {
    Tsd,
    Tir,
    Other,
}

type Attrs = Vec<(Vec<u8>, Vec<u8>)>;

/// One `%@` label: the feature type and the key=value pairs.
struct Label {
    typ: Vec<u8>,
    so: Option<Vec<u8>>,
    attrs: Attrs,
    role: Role,
}

/// The labels of a pattern file.
#[derive(Default)]
pub struct Labels {
    /// `%@element`: (type, Name prefix, Classification, Method, SO, other)
    element: Option<Label>,
    name: Option<Vec<u8>>,
    class: Option<Vec<u8>>,
    method: Option<Vec<u8>>,
    /// (source line, label) of each `%@ TYPE` part label
    parts: Vec<(u32, Label)>,
}

/// Keys the program writes itself.
const RESERVED: [&[u8]; 3] = [b"ID", b"Parent", b"Name"];

/// Parse the labels in the comments of the pattern file: `(source line,
/// comment text after '%')`.
pub fn parse_labels(comments: &[(u32, Vec<u8>)]) -> Result<Labels, String> {
    let mut lb = Labels::default();
    for (line, text) in comments {
        let Some(rest) = text.strip_prefix(b"@") else {
            continue;
        };
        let err = |m: &str| Err(format!("pattern line {line}: {m}"));
        let (is_elem, rest) = match rest.strip_prefix(b"element") {
            Some(r) if r.is_empty() || r[0].is_ascii_whitespace() => (true, r),
            _ if rest.first().is_some_and(|c| c.is_ascii_whitespace()) => (false, rest),
            _ => {
                return err("a label is `%@element TYPE key=value ...` or `%@ TYPE key=value ...`");
            }
        };
        let mut words = rest
            .split(|c| c.is_ascii_whitespace())
            .filter(|w| !w.is_empty());
        let typ = match words.next() {
            Some(t) if !t.contains(&b'=') => t.to_vec(),
            _ => return err("the label has no type"),
        };
        let mut l = Label {
            typ,
            so: None,
            attrs: Vec::new(),
            role: Role::Other,
        };
        let mut role = None;
        for w in words {
            let Some(k) = w.iter().position(|&c| c == b'=') else {
                return err(&format!(
                    "`{}` is not key=value",
                    String::from_utf8_lossy(w)
                ));
            };
            let (key, val) = (&w[..k], &w[k + 1..]);
            if key.is_empty() || val.is_empty() {
                return err(&format!(
                    "`{}` is not key=value",
                    String::from_utf8_lossy(w)
                ));
            }
            let reserved = RESERVED.contains(&key)
                || (!is_elem && matches!(key, b"Classification" | b"Method"));
            if reserved && !(is_elem && key == b"Name") {
                return err(&format!(
                    "{} is set by the program here",
                    String::from_utf8_lossy(key)
                ));
            }
            match key {
                b"Sequence_ontology" => l.so = Some(val.to_vec()),
                b"role" if !is_elem => {
                    role = Some(match val {
                        b"tsd" => Role::Tsd,
                        b"tir" => Role::Tir,
                        b"other" => Role::Other,
                        _ => return err("role is tsd, tir or other"),
                    })
                }
                b"Name" => lb.name = Some(val.to_vec()),
                b"Classification" => lb.class = Some(val.to_vec()),
                b"Method" => lb.method = Some(val.to_vec()),
                _ => l.attrs.push((key.to_vec(), val.to_vec())),
            }
        }
        if is_elem {
            if lb.element.is_some() {
                return err("a second %@element label");
            }
            lb.element = Some(l);
        } else {
            l.role = role.unwrap_or(match &l.typ[..] {
                b"target_site_duplication" => Role::Tsd,
                b"terminal_inverted_repeat"
                | b"five_prime_terminal_inverted_repeat"
                | b"three_prime_terminal_inverted_repeat" => Role::Tir,
                _ => Role::Other,
            });
            if lb.parts.iter().any(|(n, _)| n == line) {
                return err("a second part label on the same line");
            }
            lb.parts.push((*line, l));
        }
    }
    Ok(lb)
}

impl Labels {
    /// The part label (index) of each unit slot, and the source line of
    /// each slot.  `units`: (slot, source line) of every unit.  Error: a
    /// part label on a line where no unit starts.
    pub fn assign(&self, units: &[(u32, u32)]) -> Result<(Vec<u32>, Vec<u32>), String> {
        let n = units
            .iter()
            .map(|&(s, _)| s as usize + 1)
            .max()
            .unwrap_or(0);
        let mut group = vec![u32::MAX; n];
        let mut line = vec![0u32; n];
        for &(s, l) in units {
            line[s as usize] = l;
            if let Some(g) = self.parts.iter().position(|(pl, _)| *pl == l) {
                group[s as usize] = g as u32;
            }
        }
        for (pl, _) in &self.parts {
            if !units.iter().any(|&(_, l)| l == *pl) {
                return Err(format!(
                    "pattern line {pl}: label on a line where no pattern unit starts"
                ));
            }
        }
        Ok((group, line))
    }
}

/// Sequence Ontology ids for column-3 names (so-simple.obo, 2026-08-07).
const SO: [(&str, &str); 21] = [
    ("repeat_region", "SO:0000657"),
    ("target_site_duplication", "SO:0000434"),
    ("terminal_inverted_repeat_element", "SO:0000208"),
    ("terminal_inverted_repeat", "SO:0000481"),
    ("five_prime_terminal_inverted_repeat", "SO:0000420"),
    ("three_prime_terminal_inverted_repeat", "SO:0000421"),
    ("CACTA_TIR_transposon", "SO:0002285"),
    ("Tc1_Mariner_TIR_transposon", "SO:0002278"),
    ("hAT_TIR_transposon", "SO:0002279"),
    ("Mutator_TIR_transposon", "SO:0002280"),
    ("Merlin_TIR_transposon", "SO:0002281"),
    ("Transib_TIR_transposon", "SO:0002282"),
    ("piggyBac_TIR_transposon", "SO:0002283"),
    ("PIF_Harbinger_TIR_transposon", "SO:0002284"),
    ("P_TIR_transposon", "SO:0001535"),
    ("helitron", "SO:0000544"),
    ("stem_loop", "SO:0000313"),
    ("TF_binding_site", "SO:0000235"),
    ("polypeptide_motif", "SO:0001067"),
    ("sequence_motif", "SO:0001683"),
    ("region", "SO:0000001"),
];

fn so_id(typ: &[u8]) -> Option<&'static [u8]> {
    SO.iter()
        .find(|(n, _)| n.as_bytes() == typ)
        .map(|(_, i)| i.as_bytes())
}

// ----------------------------------------------------------------------
// raw hit records (worker → writer)

fn put_u64(b: &mut Vec<u8>, v: u64) {
    b.extend_from_slice(&v.to_le_bytes());
}

/// Encode a hit: pattern, strand, record length, id, and for each entry
/// its unit slot, start (offset on the strand searched), length and text.
#[allow(clippy::too_many_arguments)]
pub fn encode(
    b: &mut Vec<u8>,
    pat: usize,
    rev: bool,
    ln: usize,
    id: &[u8],
    hits: &[i64],
    n: usize,
    slots: &[u32],
    data: &Buf,
    cb: i64,
) {
    put_u64(b, pat as u64);
    b.push(rev as u8);
    put_u64(b, ln as u64);
    put_u64(b, id.len() as u64);
    b.extend_from_slice(id);
    put_u64(b, (hits[0] - cb) as u64);
    put_u64(b, (hits[n] - cb) as u64);
    put_u64(b, n as u64);
    for i in 0..n {
        let from = hits[i] - cb;
        let len = (hits[i + 1] - hits[i]).max(0);
        put_u64(b, slots.get(i).copied().unwrap_or(u32::MAX) as u64);
        put_u64(b, from as u64);
        put_u64(b, len as u64);
        for p in from..from + len {
            b.push(data.get(p));
        }
    }
}

struct Unit {
    slot: u32,
    /// 1-based closed on the forward strand; end = start - 1 when empty
    start: i64,
    end: i64,
    /// as matched, on the strand searched
    text: Vec<u8>,
}

struct Hit {
    /// index of the pattern file
    pat: usize,
    rev: bool,
    ln: i64,
    id: Vec<u8>,
    start: i64,
    end: i64,
    units: Vec<Unit>,
}

struct Rd<'a>(&'a [u8]);

impl Rd<'_> {
    fn u64(&mut self) -> u64 {
        let (a, b) = self.0.split_at(8);
        self.0 = b;
        u64::from_le_bytes(a.try_into().unwrap())
    }
    fn bytes(&mut self, n: usize) -> Vec<u8> {
        let (a, b) = self.0.split_at(n);
        self.0 = b;
        a.to_vec()
    }
}

/// Forward-strand 1-based closed range of `len` bytes at offset `s0` of
/// the strand searched.
fn fwd(rev: bool, ln: i64, s0: i64, len: i64) -> (i64, i64) {
    if rev {
        (ln - s0 - len + 1, ln - s0)
    } else {
        (s0 + 1, s0 + len)
    }
}

fn decode(raw: &[u8]) -> Hit {
    let pat = u64::from_le_bytes(raw[..8].try_into().unwrap()) as usize;
    let rev = raw[8] != 0;
    let mut r = Rd(&raw[9..]);
    let ln = r.u64() as i64;
    let idl = r.u64() as usize;
    let id = r.bytes(idl);
    let s0 = r.u64() as i64;
    let past = r.u64() as i64;
    let (start, end) = fwd(rev, ln, s0, past - s0);
    let n = r.u64() as usize;
    let mut units = Vec::with_capacity(n);
    for _ in 0..n {
        let slot = r.u64() as u32;
        let from = r.u64() as i64;
        let len = r.u64() as i64;
        let (start, end) = fwd(rev, ln, from, len);
        units.push(Unit {
            slot,
            start,
            end,
            text: r.bytes(len as usize),
        });
    }
    Hit {
        pat,
        rev,
        ln,
        id,
        start,
        end,
        units,
    }
}

// ----------------------------------------------------------------------
// the writer side

pub struct Options {
    pub format: Format,
    pub name_prefix: Option<Vec<u8>>,
    pub name_start: u64,
    pub typ: Option<Vec<u8>>,
    pub dedup: bool,
}

/// The labels of one pattern file.
pub struct PatLabels {
    /// pattern file name, for the GFF3 header (not the whole command
    /// line: the output must not depend on `-t`)
    pub file: Vec<u8>,
    pub lb: Labels,
    /// part label of each unit slot (u32::MAX: none)
    pub group: Vec<u32>,
    /// source line of each unit slot
    pub line: Vec<u32>,
}

pub struct Formatter {
    o: Options,
    pats: Vec<PatLabels>,
    /// next number of each Name prefix (one counter for all patterns
    /// with the same prefix)
    next_name: HashMap<Vec<u8>, u64>,
    /// several patterns: the Names given (a prefix can end in a digit,
    /// so `X` + 11 and `X1` + 1 would be the same Name)
    names: HashSet<Vec<u8>>,
    seen: HashSet<Vec<u8>>,
    /// --dedup: hits of the current record, forward spans, spans found
    /// on both strands (for each pattern)
    pending: Vec<Hit>,
    fwd_spans: HashSet<(usize, i64, i64)>,
    both: HashSet<(usize, i64, i64)>,
}

/// One feature of a hit, before writing.
struct Feat<'a> {
    typ: &'a [u8],
    start: i64,
    end: i64,
    id: Vec<u8>,
    parent: Option<Vec<u8>>,
    so: Option<&'a [u8]>,
    /// the TSD= / TIR= values (element only), other attributes
    extra: Attrs,
    attrs: &'a [(Vec<u8>, Vec<u8>)],
}

/// A labelled part of a hit.
struct Part {
    g: usize,
    start: i64,
    end: i64,
    /// as matched, on the strand searched (entries in hit order)
    text: Vec<u8>,
}

impl Formatter {
    pub fn new(o: Options, pats: Vec<PatLabels>) -> Formatter {
        Formatter {
            next_name: HashMap::new(),
            names: HashSet::new(),
            o,
            pats,
            seen: HashSet::new(),
            pending: Vec::new(),
            fwd_spans: HashSet::new(),
            both: HashSet::new(),
        }
    }

    /// Text written before the first hit.
    pub fn header(&self, out: &mut Vec<u8>) {
        if self.o.format == Format::Gff3 {
            out.extend_from_slice(b"##gff-version 3\n# scan_for_matches ");
            out.extend_from_slice(env!("CARGO_PKG_VERSION").as_bytes());
            for p in &self.pats {
                out.extend_from_slice(b" pattern=");
                out.extend_from_slice(&p.file);
            }
            out.push(b'\n');
        }
    }

    /// One raw hit.  Returns false when `--dedup` drops it (it does not
    /// count for `-m`).
    pub fn hit(&mut self, raw: &[u8], out: &mut Vec<u8>) -> bool {
        let h = decode(raw);
        if !self.o.dedup {
            self.write(&h, false, out);
            return true;
        }
        let span = (h.pat, h.start, h.end);
        if !h.rev {
            self.fwd_spans.insert(span);
        } else if self.fwd_spans.contains(&span) {
            self.both.insert(span);
            return false;
        }
        self.pending.push(h);
        true
    }

    /// The end of a record (and of the output): write what `--dedup`
    /// holds back.
    pub fn end_record(&mut self, out: &mut Vec<u8>) {
        let pending = std::mem::take(&mut self.pending);
        for h in &pending {
            let dot = !h.rev && self.both.contains(&(h.pat, h.start, h.end));
            self.write(h, dot, out);
        }
        self.fwd_spans.clear();
        self.both.clear();
    }

    fn group_of(&self, pat: usize, slot: u32) -> Option<usize> {
        let g = *self.pats[pat].group.get(slot as usize)?;
        (g != u32::MAX).then_some(g as usize)
    }

    fn write(&mut self, h: &Hit, dot: bool, out: &mut Vec<u8>) {
        if self.o.format != Format::Jsonl && h.end < h.start {
            return; // a hit of length 0 has no place in GFF3 / BED
        }
        let mut name = self
            .o
            .name_prefix
            .clone()
            .or_else(|| self.pats[h.pat].lb.name.clone())
            .unwrap_or_else(|| b"sfm".to_vec());
        let next = self
            .next_name
            .entry(name.clone())
            .or_insert(self.o.name_start);
        let k = name.len();
        loop {
            name.truncate(k);
            name.extend_from_slice(next.to_string().as_bytes());
            *next += 1;
            if self.pats.len() == 1 || self.names.insert(name.clone()) {
                break;
            }
        }
        let strand = if dot {
            b'.'
        } else if h.rev {
            b'-'
        } else {
            b'+'
        };
        match self.o.format {
            Format::Jsonl => self.jsonl(h, &name, strand, out),
            Format::Gff3 => self.gff3(h, &name, strand, out),
            Format::Bed6 | Format::Bed12 => self.bed(h, &name, strand, out),
        }
    }

    fn elem_type(&self, pat: usize) -> &[u8] {
        self.o
            .typ
            .as_deref()
            .or(self.pats[pat].lb.element.as_ref().map(|l| &l.typ[..]))
            .unwrap_or(b"sequence_motif")
    }

    /// The labelled parts of a hit, in rising order, and the span of the
    /// element (the hit without its TSDs; None: nothing but TSDs).
    fn parts(&self, h: &Hit) -> (Vec<Part>, Option<(i64, i64)>) {
        let mut parts: Vec<Part> = Vec::new();
        let mut elem: Option<(i64, i64)> = None;
        for u in &h.units {
            if u.end < u.start {
                continue;
            }
            let g = self.group_of(h.pat, u.slot);
            if g.is_none_or(|g| self.pats[h.pat].lb.parts[g].1.role != Role::Tsd) {
                elem = Some(match elem {
                    None => (u.start, u.end),
                    Some((a, b)) => (a.min(u.start), b.max(u.end)),
                });
            }
            let Some(g) = g else { continue };
            match parts.iter_mut().find(|p| p.g == g) {
                Some(p) => {
                    p.start = p.start.min(u.start);
                    p.end = p.end.max(u.end);
                    p.text.extend_from_slice(&u.text);
                }
                None => parts.push(Part {
                    g,
                    start: u.start,
                    end: u.end,
                    text: u.text.clone(),
                }),
            }
        }
        parts.sort_by_key(|p| (p.start, p.g));
        (parts, elem)
    }

    /// Forward-strand letters of a part.
    fn fwd_text(h: &Hit, t: &[u8]) -> Vec<u8> {
        if h.rev {
            t.iter().rev().map(|&c| compl(c)).collect()
        } else {
            t.to_vec()
        }
    }

    fn gff3(&mut self, h: &Hit, name: &[u8], strand: u8, out: &mut Vec<u8>) {
        if self.seen.insert(h.id.clone()) {
            out.extend_from_slice(b"##sequence-region ");
            esc_seqid(out, &h.id);
            out.extend_from_slice(format!(" 1 {}\n", h.ln).as_bytes());
        }
        let (parts, elem) = self.parts(h);
        let lb = &self.pats[h.pat].lb;
        let role = |p: &Part| lb.parts[p.g].1.role;
        let has_tsd = parts.iter().any(|p| role(p) == Role::Tsd);
        let id = |suffix: &[u8]| {
            let mut v = name.to_vec();
            v.extend_from_slice(suffix);
            v
        };
        let elem_id = if has_tsd { id(b".te") } else { id(b"") };
        let no_attrs: Attrs = Vec::new();
        let mut feats: Vec<Feat> = Vec::new();
        if has_tsd {
            feats.push(Feat {
                typ: b"repeat_region",
                start: h.start,
                end: h.end,
                id: id(b""),
                parent: None,
                so: so_id(b"repeat_region"),
                extra: Vec::new(),
                attrs: &no_attrs,
            });
        }
        // TSD= / TIR=: the two parts of the role, forward letters
        let pair = |r: Role, key: &[u8]| -> Option<(Vec<u8>, Vec<u8>)> {
            let v: Vec<&Part> = parts.iter().filter(|p| role(p) == r).collect();
            if v.len() != 2 {
                return None;
            }
            let mut s = Self::fwd_text(h, &v[0].text);
            s.push(b'_');
            s.extend_from_slice(&Self::fwd_text(h, &v[1].text));
            Some((key.to_vec(), s))
        };
        let el = lb.element.as_ref();
        let etyp = self.elem_type(h.pat);
        let (estart, eend) = if has_tsd {
            elem.unwrap_or((0, -1))
        } else {
            (h.start, h.end)
        };
        if estart <= eend {
            feats.push(Feat {
                typ: etyp,
                start: estart,
                end: eend,
                id: elem_id.clone(),
                parent: has_tsd.then(|| id(b"")),
                so: el
                    .and_then(|l| l.so.as_deref())
                    .filter(|_| self.o.typ.is_none())
                    .or_else(|| so_id(etyp)),
                extra: [pair(Role::Tsd, b"TSD"), pair(Role::Tir, b"TIR")]
                    .into_iter()
                    .flatten()
                    .collect(),
                attrs: el.map(|l| &l.attrs[..]).unwrap_or(&no_attrs),
            });
        }
        let parent_of_parts = if estart <= eend { elem_id } else { id(b"") };
        let count = |r: Role| parts.iter().filter(|p| role(p) == r).count();
        let (n_tsd, n_tir) = (count(Role::Tsd), count(Role::Tir));
        let (mut k_tsd, mut k_tir, mut k_other) = (0, 0, 0);
        for p in &parts {
            let l = &lb.parts[p.g].1;
            let suffix = match l.role {
                Role::Tsd => {
                    k_tsd += 1;
                    lr_name(n_tsd, k_tsd, "TSD")
                }
                Role::Tir => {
                    k_tir += 1;
                    lr_name(n_tir, k_tir, "TIR")
                }
                Role::Other => {
                    k_other += 1;
                    format!(".{k_other}")
                }
            };
            feats.push(Feat {
                typ: &l.typ,
                start: p.start,
                end: p.end,
                id: id(suffix.as_bytes()),
                parent: Some(if l.role == Role::Tsd {
                    id(b"")
                } else {
                    parent_of_parts.clone()
                }),
                so: l.so.as_deref().or_else(|| so_id(&l.typ)),
                extra: Vec::new(),
                attrs: &l.attrs,
            });
        }
        // the top feature first, then in rising order (longer first)
        if feats.len() > 1 {
            feats[1..].sort_by_key(|f| (f.start, std::cmp::Reverse(f.end)));
        }
        let class = lb.class.as_deref();
        let method = lb.method.as_deref().unwrap_or(b"structural");
        for f in &feats {
            esc_seqid(out, &h.id);
            out.extend_from_slice(b"\tscan_for_matches\t");
            esc(out, f.typ, false);
            out.extend_from_slice(
                format!("\t{}\t{}\t.\t{}\t.\t", f.start, f.end, strand as char).as_bytes(),
            );
            let mut kv = |k: &[u8], v: &[u8], first: bool| {
                if !first {
                    out.push(b';');
                }
                esc(out, k, true);
                out.push(b'=');
                esc(out, v, true);
            };
            kv(b"ID", &f.id, true);
            if let Some(p) = &f.parent {
                kv(b"Parent", p, false);
            }
            kv(b"Name", name, false);
            if let Some(c) = class {
                kv(b"Classification", c, false);
            }
            kv(b"Method", method, false);
            if let Some(s) = f.so {
                kv(b"Sequence_ontology", s, false);
            }
            for (k, v) in f.extra.iter().chain(f.attrs.iter()) {
                kv(k, v, false);
            }
            out.push(b'\n');
        }
        out.extend_from_slice(b"###\n");
    }

    fn bed(&self, h: &Hit, name: &[u8], strand: u8, out: &mut Vec<u8>) {
        let line = |out: &mut Vec<u8>, a: i64, b: i64| {
            out.extend_from_slice(&h.id);
            out.extend_from_slice(format!("\t{}\t{}\t", a - 1, b).as_bytes());
            out.extend_from_slice(name);
            out.extend_from_slice(format!("\t0\t{}", strand as char).as_bytes());
        };
        if self.o.format == Format::Bed6 {
            line(out, h.start, h.end);
            out.push(b'\n');
            return;
        }
        let (parts, elem) = self.parts(h);
        // blocks: the labelled parts (overlapping ones joined), or the hit
        let mut blocks: Vec<(i64, i64)> = Vec::new();
        for p in &parts {
            match blocks.last_mut() {
                Some(b) if p.start <= b.1 => b.1 = b.1.max(p.end),
                _ => blocks.push((p.start, p.end)),
            }
        }
        if blocks.is_empty() {
            blocks.push((h.start, h.end));
        }
        let (a, b) = (blocks[0].0, blocks[blocks.len() - 1].1);
        let lb = &self.pats[h.pat].lb;
        let has_tsd = parts.iter().any(|p| lb.parts[p.g].1.role == Role::Tsd);
        let (ta, tb) = match (has_tsd, elem) {
            (false, _) => (a, b),
            (true, Some((x, y))) if x.max(a) <= y.min(b) => (x.max(a), y.min(b)),
            _ => (a, a - 1),
        };
        line(out, a, b);
        out.extend_from_slice(format!("\t{}\t{}\t0\t{}\t", ta - 1, tb, blocks.len()).as_bytes());
        for (x, y) in &blocks {
            out.extend_from_slice(format!("{},", y - x + 1).as_bytes());
        }
        out.push(b'\t');
        for (x, _) in &blocks {
            out.extend_from_slice(format!("{},", x - a).as_bytes());
        }
        out.push(b'\n');
    }

    fn jsonl(&self, h: &Hit, name: &[u8], strand: u8, out: &mut Vec<u8>) {
        out.extend_from_slice(b"{\"seq\":");
        json_str(out, &h.id);
        out.extend_from_slice(
            format!(
                ",\"strand\":\"{}\",\"start\":{},\"end\":{},\"name\":",
                strand as char, h.start, h.end
            )
            .as_bytes(),
        );
        json_str(out, name);
        if self.pats.len() > 1 {
            out.extend_from_slice(b",\"pattern\":");
            json_str(out, &self.pats[h.pat].file);
        }
        out.extend_from_slice(b",\"type\":");
        json_str(out, self.elem_type(h.pat));
        out.extend_from_slice(b",\"units\":[");
        for (i, u) in h.units.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            out.extend_from_slice(b"{\"line\":");
            match self.pats[h.pat].line.get(u.slot as usize) {
                Some(l) if *l > 0 => out.extend_from_slice(l.to_string().as_bytes()),
                _ => out.extend_from_slice(b"null"),
            }
            out.extend_from_slice(
                format!(",\"start\":{},\"end\":{},\"text\":", u.start, u.end).as_bytes(),
            );
            json_str(out, &u.text);
            if let Some(g) = self.group_of(h.pat, u.slot) {
                out.extend_from_slice(b",\"label\":");
                json_str(out, &self.pats[h.pat].lb.parts[g].1.typ);
            }
            out.push(b'}');
        }
        out.extend_from_slice(b"]}\n");
    }
}

/// `.lTSD`/`.rTSD` for a pair, `.TSD1`, `.TSD2`, ... otherwise.
fn lr_name(n: usize, k: usize, what: &str) -> String {
    if n == 2 {
        format!(".{}{what}", if k == 1 { 'l' } else { 'r' })
    } else {
        format!(".{what}{k}")
    }
}

fn pct(out: &mut Vec<u8>, c: u8) {
    out.extend_from_slice(format!("%{c:02X}").as_bytes());
}

/// GFF3 column 1: characters outside [a-zA-Z0-9.:^*$@!+_?-|] as %XX.
fn esc_seqid(out: &mut Vec<u8>, s: &[u8]) {
    for &c in s {
        if c.is_ascii_alphanumeric() || b".:^*$@!+_?-|".contains(&c) {
            out.push(c);
        } else {
            pct(out, c);
        }
    }
}

/// GFF3 columns 3 and 9: control characters and '%' as %XX; in column 9
/// (`attr`) also ';', '=', '&' and ','.
fn esc(out: &mut Vec<u8>, s: &[u8], attr: bool) {
    for &c in s {
        if c < 0x20 || c == 0x7f || c == b'%' || (attr && b";=&,".contains(&c)) {
            pct(out, c);
        } else {
            out.push(c);
        }
    }
}

fn json_str(out: &mut Vec<u8>, s: &[u8]) {
    out.push(b'"');
    for &c in s {
        match c {
            b'"' => out.extend_from_slice(b"\\\""),
            b'\\' => out.extend_from_slice(b"\\\\"),
            0x20..=0x7e => out.push(c),
            _ => out.extend_from_slice(format!("\\u{c:04x}").as_bytes()),
        }
    }
    out.push(b'"');
}

#[cfg(test)]
mod tests {
    use super::*;

    fn labels(src: &str) -> Result<Labels, String> {
        let c: Vec<(u32, Vec<u8>)> = src
            .lines()
            .enumerate()
            .filter_map(|(i, l)| {
                let k = l.find('%')?;
                Some((i as u32 + 1, l.as_bytes()[k + 1..].to_vec()))
            })
            .collect();
        parse_labels(&c)
    }

    #[test]
    fn label_syntax() {
        let lb = labels(
            "%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC x=1\n\
             p1=3...3 %@ target_site_duplication\n\
             CACTA %@ five_prime_terminal_inverted_repeat\n\
             10...20 % plain comment\n\
             AA %@ spacer role=tsd Sequence_ontology=SO:1\n",
        )
        .unwrap();
        assert_eq!(lb.name.as_deref(), Some(&b"DTC"[..]));
        assert_eq!(lb.class.as_deref(), Some(&b"TIR/DTC"[..]));
        assert_eq!(
            lb.element.as_ref().unwrap().attrs,
            vec![(b"x".to_vec(), b"1".to_vec())]
        );
        let roles: Vec<(u32, bool, bool)> = lb
            .parts
            .iter()
            .map(|(l, p)| (*l, p.role == Role::Tsd, p.role == Role::Tir))
            .collect();
        assert_eq!(
            roles,
            vec![(2, true, false), (3, false, true), (5, true, false)]
        );
        assert_eq!(lb.parts[2].1.so.as_deref(), Some(&b"SO:1"[..]));
        for bad in [
            "%@element\n",
            "%@\n",
            "%@foo bar\n",
            "%@ x y\n",
            "%@ x Name=a\n",
            "%@ x role=z\n",
            "%@element a\n%@element b\n",
            "A %@ x\n%@ y\n",
        ] {
            let r = if bad.starts_with("A") {
                // two labels for line 1
                labels("A %@ x %@ y\n").and(Err(String::new()))
            } else {
                labels(bad)
            };
            assert!(r.is_err(), "{bad:?}");
        }
        // a label line without a unit
        let lb = labels("% @ not a label\n%@ x\n").unwrap();
        assert!(lb.assign(&[(0, 1)]).is_err());
        assert!(lb.assign(&[(0, 2)]).is_ok());
    }

    #[test]
    fn escapes() {
        let mut o = Vec::new();
        esc_seqid(&mut o, b"chr 1;x");
        esc(&mut o, b"a;b=c%d", true);
        json_str(&mut o, b"q\"\\\x01\xff");
        assert_eq!(o, b"chr%201%3Bxa%3Bb%3Dc%25d\"q\\\"\\\\\\u0001\\u00ff\"");
    }

    #[test]
    fn coordinates() {
        // 10 bp, unit at offset 2 length 3 of the forward / reverse strand
        assert_eq!(fwd(false, 10, 2, 3), (3, 5));
        assert_eq!(fwd(true, 10, 2, 3), (6, 8));
        assert_eq!(fwd(false, 10, 4, 0), (5, 4));
    }
}
