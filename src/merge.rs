//! `--merge`: join the GFF3 output of several runs (patterns), remove
//! elements found twice and give genome-wide unique Names.

use std::collections::{BTreeSet, HashMap};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use crate::gz::{self, GzReader};

const USAGE: &str = "\
usage: scan_for_matches --merge [--overlap F] [--name-prefix P] [--name-start N]
                        [--output FILE] FILE.gff3 [FILE2.gff3 ...]

Join GFF3 files (priority = order; '-' is stdin, gzip is detected), remove
duplicate elements and give new unique Names.
  --overlap F       also drop an element that overlaps a kept one by at
                    least F (0 < F <= 1) of the length of each
  --name-prefix P   Name prefix for all elements (default: old Name
                    without trailing digits)
  --name-start N    first number of each prefix (default 1)
  --output FILE     write to FILE instead of stdout
";

/// The command line options.
#[derive(Debug, PartialEq)]
struct Opts {
    overlap: Option<f64>,
    prefix: Option<Vec<u8>>,
    start: u64,
    output: Option<String>,
    files: Vec<String>,
}

/// What the argument parser decided.
#[derive(Debug, PartialEq)]
enum Parsed {
    Run(Opts),
    Help,
}

/// Split `--opt=value` or take the next argument as the value.
fn opt_value(
    args: &[Vec<u8>],
    i: &mut usize,
    name: &str,
    inline: Option<&[u8]>,
) -> Result<Vec<u8>, String> {
    if let Some(v) = inline {
        return Ok(v.to_vec());
    }
    *i += 1;
    args.get(*i)
        .cloned()
        .ok_or_else(|| format!("option {name} needs a value"))
}

/// Parse the arguments after `--merge`.
fn parse_args(args: &[Vec<u8>]) -> Result<Parsed, String> {
    let mut o = Opts {
        overlap: None,
        prefix: None,
        start: 1,
        output: None,
        files: Vec::new(),
    };
    let mut i = 0;
    let mut only_files = false;
    while i < args.len() {
        let a = &args[i];
        if only_files || a == b"-" || !a.starts_with(b"-") {
            o.files.push(String::from_utf8_lossy(a).into_owned());
        } else if a == b"--" {
            only_files = true;
        } else if a == b"-h" || a == b"--help" {
            return Ok(Parsed::Help);
        } else {
            let (name, inline) = match a.iter().position(|&c| c == b'=') {
                Some(p) if a.starts_with(b"--") => (&a[..p], Some(&a[p + 1..])),
                _ => (&a[..], None),
            };
            let name = String::from_utf8_lossy(name).into_owned();
            let val = match name.as_str() {
                "--overlap" | "--name-prefix" | "--name-start" | "--output" => {
                    opt_value(args, &mut i, &name, inline)?
                }
                _ => return Err(format!("unknown option {}", String::from_utf8_lossy(a))),
            };
            let text = String::from_utf8_lossy(&val).into_owned();
            match name.as_str() {
                "--overlap" => match text.parse::<f64>() {
                    Ok(f) if f.is_finite() && f > 0.0 && f <= 1.0 => o.overlap = Some(f),
                    _ => return Err(format!("--overlap must be a number in (0, 1]: {text}")),
                },
                "--name-start" => {
                    o.start = text.parse().map_err(|_| {
                        format!("--name-start must be a non-negative integer: {text}")
                    })?
                }
                "--output" => o.output = Some(text),
                _ => o.prefix = Some(val),
            }
        }
        i += 1;
    }
    if o.files.is_empty() {
        return Err("no input file".to_string());
    }
    Ok(Parsed::Run(o))
}

/// A parsed GFF3 line (borrowed from the file buffer).
struct Feat<'a> {
    range: (usize, usize),
    lineno: usize,
    id: Option<&'a [u8]>,
    parents: Vec<&'a [u8]>,
    seqid: &'a [u8],
    start: i64,
    end: i64,
    /// number of `###` lines before it
    block: usize,
}

/// One top feature with all its descendants.
#[derive(Debug)]
struct Elem {
    file: usize,
    /// byte ranges of the lines in the file buffer, input order
    lines: Vec<(usize, usize)>,
    /// index of the top feature in `lines`
    top: usize,
    seqid: u32,
    start: i64,
    end: i64,
}

/// Byte index where the ninth column (attributes) starts.
fn attr_start(line: &[u8]) -> usize {
    let mut tabs = 0;
    for (p, &c) in line.iter().enumerate() {
        if c == b'\t' {
            tabs += 1;
            if tabs == 8 {
                return p + 1;
            }
        }
    }
    line.len()
}

/// The ninth column of a line, empty when absent or `.`.
fn attr_col(line: &[u8]) -> &[u8] {
    let a = &line[attr_start(line)..];
    if a == b"." { &[] } else { a }
}

/// The attribute pieces of an attribute column.
fn pieces(attrs: &[u8]) -> impl Iterator<Item = &[u8]> {
    attrs.split(|&c| c == b';').filter(|p| !p.is_empty())
}

/// Key and value of one attribute piece.
fn key_value(piece: &[u8]) -> (&[u8], &[u8]) {
    match piece.iter().position(|&c| c == b'=') {
        Some(p) => (&piece[..p], &piece[p + 1..]),
        None => (piece, &[]),
    }
}

/// Value of the attribute `key` in a line.
fn attr_of<'a>(line: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    pieces(attr_col(line))
        .map(key_value)
        .find(|(k, _)| *k == key)
        .map(|(_, v)| v)
}

/// Parse one data line.
fn parse_feat(data: &[u8], range: (usize, usize), lineno: usize) -> Result<Feat<'_>, String> {
    let line = &data[range.0..range.1];
    let cols: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
    if cols.len() != 9 {
        return Err(format!(
            "expected 9 tab-separated columns, found {}",
            cols.len()
        ));
    }
    let num = |c: &[u8], what: &str| -> Result<i64, String> {
        std::str::from_utf8(c)
            .ok()
            .and_then(|s| s.parse::<i64>().ok())
            .filter(|&v| v >= 1)
            .ok_or_else(|| format!("bad {what} '{}'", String::from_utf8_lossy(c)))
    };
    let start = num(cols[3], "start")?;
    let end = num(cols[4], "end")?;
    if start > end {
        return Err(format!("start {start} is greater than end {end}"));
    }
    let id = attr_of(line, b"ID");
    let parents = match attr_of(line, b"Parent") {
        Some(v) => v.split(|&c| c == b',').filter(|p| !p.is_empty()).collect(),
        None => Vec::new(),
    };
    Ok(Feat {
        range,
        lineno,
        id,
        parents,
        seqid: cols[0],
        start,
        end,
        block: 0,
    })
}

/// Index of the element root of every feature (follows first parents).
fn find_roots(feats: &[Feat], parent: &[usize], name: &str) -> Result<Vec<usize>, String> {
    let mut root = vec![usize::MAX; feats.len()];
    let mut busy = vec![false; feats.len()];
    for s in 0..feats.len() {
        let mut path = Vec::new();
        let mut cur = s;
        while root[cur] == usize::MAX && parent[cur] != usize::MAX {
            if busy[cur] {
                return Err(format!("{name}:{}: Parent loop", feats[cur].lineno));
            }
            busy[cur] = true;
            path.push(cur);
            cur = parent[cur];
        }
        let r = if root[cur] == usize::MAX {
            cur
        } else {
            root[cur]
        };
        root[cur] = r;
        for p in path {
            root[p] = r;
        }
    }
    Ok(root)
}

/// Interned seqids and the first `##sequence-region` line of each.
#[derive(Default)]
struct Seqs {
    ids: HashMap<Vec<u8>, u32>,
    names: Vec<Vec<u8>>,
    regions: HashMap<u32, Vec<u8>>,
}

impl Seqs {
    fn id(&mut self, s: &[u8]) -> u32 {
        if let Some(&i) = self.ids.get(s) {
            return i;
        }
        let i = self.names.len() as u32;
        self.ids.insert(s.to_vec(), i);
        self.names.push(s.to_vec());
        i
    }
}

/// Read the lines of one file into elements.
fn read_elems(file: usize, name: &str, data: &[u8], seqs: &mut Seqs) -> Result<Vec<Elem>, String> {
    let mut feats = Vec::new();
    let (mut pos, mut lineno, mut block) = (0, 0, 0);
    while pos < data.len() {
        let mut end = data[pos..]
            .iter()
            .position(|&c| c == b'\n')
            .map_or(data.len(), |p| pos + p);
        let next = end + 1;
        if end > pos && data[end - 1] == b'\r' {
            end -= 1;
        }
        let line = &data[pos..end];
        lineno += 1;
        if line == b"##FASTA" {
            break;
        }
        if line == b"###" {
            block += 1;
        }
        if line.starts_with(b"##sequence-region") {
            let sid = line[17..]
                .split(|c| c.is_ascii_whitespace())
                .find(|w| !w.is_empty());
            if let Some(sid) = sid {
                let id = seqs.id(sid);
                seqs.regions.entry(id).or_insert_with(|| line.to_vec());
            }
        } else if !line.is_empty() && line[0] != b'#' {
            let mut f = parse_feat(data, (pos, end), lineno)
                .map_err(|m| format!("{name}:{lineno}: {m}"))?;
            f.block = block;
            feats.push(f);
        }
        pos = next;
    }
    group(file, name, &feats, seqs)
}

/// Attach every feature to the element of its (first) parent.
fn group(file: usize, name: &str, feats: &[Feat], seqs: &mut Seqs) -> Result<Vec<Elem>, String> {
    let mut by_id: HashMap<&[u8], usize> = HashMap::new();
    for (i, f) in feats.iter().enumerate() {
        if let Some(id) = f.id {
            let first = *by_id.entry(id).or_insert(i);
            // lines of one feature (same ID) have the same parents and no
            // `###` between them; two elements with one ID (files joined
            // with cat) are an error
            if feats[first].parents != f.parents || feats[first].block != f.block {
                let msg = format!(
                    "ID '{}' is used twice (line {}); IDs must be unique in one file",
                    String::from_utf8_lossy(id),
                    feats[first].lineno
                );
                return Err(format!("{name}:{}: {msg}", f.lineno));
            }
        }
    }
    let mut parent = vec![usize::MAX; feats.len()];
    let mut all: Vec<Vec<usize>> = vec![Vec::new(); feats.len()];
    for (i, f) in feats.iter().enumerate() {
        for p in &f.parents {
            let Some(&j) = by_id.get(p) else {
                let msg = format!(
                    "Parent '{}' is not an ID in this file",
                    String::from_utf8_lossy(p)
                );
                return Err(format!("{name}:{}: {msg}", f.lineno));
            };
            all[i].push(j);
        }
        if let Some(&j) = all[i].first() {
            parent[i] = j;
        }
    }
    let root = find_roots(feats, &parent, name)?;
    for (i, ps) in all.iter().enumerate() {
        if ps.iter().any(|&j| root[j] != root[i]) {
            let msg = "Parents are in different elements";
            return Err(format!("{name}:{}: {msg}", feats[i].lineno));
        }
    }
    let mut elem_of = vec![usize::MAX; feats.len()];
    let mut elems: Vec<Elem> = Vec::new();
    for (i, f) in feats.iter().enumerate() {
        if root[i] == i {
            elem_of[i] = elems.len();
            let seqid = seqs.id(f.seqid);
            elems.push(Elem {
                file,
                lines: Vec::new(),
                top: 0,
                seqid,
                start: f.start,
                end: f.end,
            });
        }
    }
    for (i, f) in feats.iter().enumerate() {
        let e = &mut elems[elem_of[root[i]]];
        if root[i] == i {
            e.top = e.lines.len();
        }
        e.lines.push(f.range);
    }
    Ok(elems)
}

/// Do two spans match: equal, or (with `f`) a reciprocal overlap of `f`.
fn same_span(a: (i64, i64), b: (i64, i64), f: Option<f64>) -> bool {
    if a == b {
        return true;
    }
    let Some(f) = f else { return false };
    let ov = (a.1.min(b.1) - a.0.max(b.0) + 1) as f64;
    let need = |s: (i64, i64)| f * (s.1 - s.0 + 1) as f64 - 1e-9;
    ov > 0.0 && ov >= need(a) && ov >= need(b)
}

/// For every element: the index of the kept element it duplicates.
/// Elements are in priority order.
fn find_dups(elems: &[Elem], nseq: usize, f: Option<f64>) -> Vec<Option<usize>> {
    let mut kept: Vec<BTreeSet<(i64, i64, usize)>> = vec![BTreeSet::new(); nseq];
    let mut dup = vec![None; elems.len()];
    for (i, e) in elems.iter().enumerate() {
        let set = &mut kept[e.seqid as usize];
        let lo = match f {
            Some(f) => e.start - ((e.end - e.start + 1) as f64 / f) as i64 - 2,
            None => e.start,
        };
        let hi = if f.is_some() { e.end } else { e.start };
        let found = set
            .range((lo, i64::MIN, 0)..=(hi, i64::MAX, usize::MAX))
            .filter(|&&(s, t, _)| same_span((s, t), (e.start, e.end), f))
            .map(|&(_, _, k)| k)
            .min();
        match found {
            Some(k) => dup[i] = Some(k),
            None => {
                set.insert((e.start, e.end, i));
            }
        }
    }
    dup
}

/// Compare seqids: digit runs as numbers, then plain bytes.
fn natural_cmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].is_ascii_digit() && b[j].is_ascii_digit() {
            fn run(s: &[u8], p: usize) -> (usize, &[u8]) {
                let n = s[p..].iter().take_while(|c| c.is_ascii_digit()).count();
                let z = s[p..p + n].iter().take_while(|&&c| c == b'0').count();
                (p + n, &s[p + z..p + n])
            }
            let (ni, da) = run(a, i);
            let (nj, db) = run(b, j);
            let c = da.len().cmp(&db.len()).then_with(|| da.cmp(db));
            if c.is_ne() {
                return c;
            }
            (i, j) = (ni, nj);
        } else {
            if a[i] != b[j] {
                return a[i].cmp(&b[j]);
            }
            i += 1;
            j += 1;
        }
    }
    (a.len() - i).cmp(&(b.len() - j)).then_with(|| a.cmp(b))
}

/// GFF3 escape of an attribute value.
fn escape(v: &[u8]) -> Vec<u8> {
    let mut o = Vec::with_capacity(v.len());
    for &c in v {
        if c < 0x20 || c == 0x7f || b",;=&%".contains(&c) {
            o.extend_from_slice(format!("%{c:02X}").as_bytes());
        } else {
            o.push(c);
        }
    }
    o
}

/// Undo `%XX` escapes.
fn unescape(v: &[u8]) -> Vec<u8> {
    let hex = |c: u8| (c as char).to_digit(16);
    let mut o = Vec::with_capacity(v.len());
    let mut i = 0;
    while i < v.len() {
        if v[i] == b'%' && i + 2 < v.len() && hex(v[i + 1]).is_some() && hex(v[i + 2]).is_some() {
            o.push((hex(v[i + 1]).unwrap() * 16 + hex(v[i + 2]).unwrap()) as u8);
            i += 3;
        } else {
            o.push(v[i]);
            i += 1;
        }
    }
    o
}

/// The old and new names of one element.
struct Rename<'a> {
    old_top: Option<&'a [u8]>,
    new: &'a [u8],
}

impl Rename<'_> {
    /// New value of an ID or Parent value.
    fn map(&self, v: &[u8], out: &mut Vec<u8>) {
        out.extend_from_slice(self.new);
        match self.old_top {
            Some(t) if v == t => {}
            Some(t) if v.len() > t.len() && v.starts_with(t) && v[t.len()] == b'.' => {
                out.extend_from_slice(&v[t.len()..])
            }
            _ => {
                out.push(b'.');
                out.extend_from_slice(v);
            }
        }
    }
}

/// Write one line with renamed ID, Parent and Name; the top line also
/// gets `Name=` (if missing) and the `Merged=` values.
fn rewrite(line: &[u8], top: bool, r: &Rename, merged: &[u8], out: &mut Vec<u8>) {
    let attrs = attr_col(line);
    out.extend_from_slice(&line[..attr_start(line)]);
    let has = |k: &[u8]| pieces(attrs).any(|p| key_value(p).0 == k);
    let add_name = top && !has(b"Name");
    let mut add_merged = top && !merged.is_empty();
    let mut first = true;
    let mut sep = |out: &mut Vec<u8>| {
        if !std::mem::take(&mut first) {
            out.push(b';');
        }
    };
    if add_name && !has(b"ID") {
        sep(out);
        out.extend_from_slice(b"Name=");
        out.extend_from_slice(r.new);
    }
    for p in pieces(attrs) {
        sep(out);
        let (k, v) = key_value(p);
        match k {
            b"Name" => {
                out.extend_from_slice(b"Name=");
                out.extend_from_slice(r.new);
            }
            b"ID" => {
                out.extend_from_slice(b"ID=");
                r.map(v, out);
                if add_name {
                    out.extend_from_slice(b";Name=");
                    out.extend_from_slice(r.new);
                }
            }
            b"Parent" => {
                out.extend_from_slice(b"Parent=");
                for (n, x) in v.split(|&c| c == b',').enumerate() {
                    if n > 0 {
                        out.push(b',');
                    }
                    r.map(x, out);
                }
            }
            b"Merged" if add_merged => {
                add_merged = false;
                out.extend_from_slice(p);
                out.push(b',');
                out.extend_from_slice(merged);
            }
            _ => out.extend_from_slice(p),
        }
    }
    if add_merged {
        sep(out);
        out.extend_from_slice(b"Merged=");
        out.extend_from_slice(merged);
    }
    out.push(b'\n');
}

/// Name prefix of an old name: no trailing digits, else `sfm`.
fn default_prefix(old: &[u8]) -> Vec<u8> {
    let n = old
        .iter()
        .rposition(|c| !c.is_ascii_digit())
        .map_or(0, |p| p + 1);
    if n == 0 {
        b"sfm".to_vec()
    } else {
        old[..n].to_vec()
    }
}

/// Open an input (`-` is stdin) and read it all; gzip is unpacked.
fn read_input(name: &str) -> std::io::Result<Vec<u8>> {
    let inner: Box<dyn Read> = if name == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(std::fs::File::open(name)?)
    };
    let mut r = BufReader::new(inner);
    let mut data = Vec::new();
    if gz::is_gzip(r.fill_buf()?) {
        GzReader::new(r)?.read_to_end(&mut data)?;
    } else {
        r.read_to_end(&mut data)?;
    }
    Ok(data)
}

/// File name without directories (for `Merged=`).
fn base_name(name: &str) -> &str {
    if name == "-" {
        "stdin"
    } else {
        name.rsplit('/').next().unwrap_or(name)
    }
}

/// The old Name of an element (top Name, else ID).
fn old_name<'a>(data: &'a [u8], e: &Elem) -> &'a [u8] {
    let l = e.lines[e.top];
    let line = &data[l.0..l.1];
    attr_of(line, b"Name")
        .or_else(|| attr_of(line, b"ID"))
        .unwrap_or(&[])
}

/// The merge itself: read all files, write the result; returns a message
/// on a data error.
fn merge(o: &Opts, w: &mut dyn Write) -> Result<(), String> {
    let mut bufs = Vec::new();
    let mut seqs = Seqs::default();
    let mut elems: Vec<Elem> = Vec::new();
    for (fi, name) in o.files.iter().enumerate() {
        let data = read_input(name).map_err(|e| format!("{name}: {e}"))?;
        elems.extend(read_elems(fi, name, &data, &mut seqs)?);
        bufs.push(data);
    }
    let dup = find_dups(&elems, seqs.names.len(), o.overlap);
    let mut merged: Vec<Vec<u8>> = vec![Vec::new(); elems.len()];
    let mut removed = 0;
    for (i, d) in dup.iter().enumerate() {
        if let Some(k) = *d {
            let e = &elems[i];
            let mut v = base_name(&o.files[e.file]).as_bytes().to_vec();
            v.push(b':');
            v.extend(unescape(old_name(&bufs[e.file], e)));
            if !merged[k].is_empty() {
                merged[k].push(b',');
            }
            merged[k].extend(escape(&v));
            removed += 1;
        }
    }
    let mut order: Vec<usize> = (0..elems.len()).filter(|&i| dup[i].is_none()).collect();
    order.sort_by(|&a, &b| {
        let (x, y) = (&elems[a], &elems[b]);
        natural_cmp(&seqs.names[x.seqid as usize], &seqs.names[y.seqid as usize])
            .then(x.start.cmp(&y.start))
            .then(y.end.cmp(&x.end))
            .then(a.cmp(&b))
    });
    let mut out = Vec::new();
    let names: Vec<&str> = o.files.iter().map(|s| s.as_str()).collect();
    out.extend_from_slice(
        format!(
            "##gff-version 3\n# scan_for_matches {} merge {}\n",
            env!("CARGO_PKG_VERSION"),
            names.join(" ")
        )
        .as_bytes(),
    );
    let mut seen = std::collections::HashSet::new();
    for &i in &order {
        let id = elems[i].seqid;
        if let Some(l) = seqs.regions.get(&id).filter(|_| seen.insert(id)) {
            out.extend_from_slice(l);
            out.push(b'\n');
        }
    }
    let mut counters: HashMap<Vec<u8>, u64> = HashMap::new();
    for &i in &order {
        let e = &elems[i];
        let data = &bufs[e.file];
        let old = old_name(data, e);
        let prefix = match &o.prefix {
            Some(p) => escape(p),
            None => default_prefix(old),
        };
        let c = counters.entry(prefix.clone()).or_insert(o.start);
        let mut new = prefix;
        new.extend_from_slice(c.to_string().as_bytes());
        *c += 1;
        let tl = e.lines[e.top];
        let r = Rename {
            old_top: attr_of(&data[tl.0..tl.1], b"ID"),
            new: &new,
        };
        for (n, &(a, b)) in e.lines.iter().enumerate() {
            rewrite(&data[a..b], n == e.top, &r, &merged[i], &mut out);
        }
        out.extend_from_slice(b"###\n");
        if out.len() > 1 << 20 {
            w.write_all(&out).map_err(|e| e.to_string())?;
            out.clear();
        }
    }
    w.write_all(&out)
        .and_then(|_| w.flush())
        .map_err(|e| e.to_string())?;
    eprintln!(
        "merge: {} elements read, {removed} duplicates removed, {} written",
        elems.len(),
        order.len()
    );
    Ok(())
}

/// `scan_for_matches --merge ARGS...`; `args` are the arguments after
/// `--merge`.  Returns the exit status.
pub fn main(args: &[Vec<u8>]) -> i32 {
    let o = match parse_args(args) {
        Ok(Parsed::Help) => {
            print!("{USAGE}");
            return 0;
        }
        Ok(Parsed::Run(o)) => o,
        Err(m) => {
            eprintln!("merge: {m}");
            return 2;
        }
    };
    for f in o.files.iter().filter(|f| *f != "-") {
        if let Err(e) = std::fs::File::open(f) {
            eprintln!("merge: cannot open {f}: {e}");
            return 2;
        }
    }
    let mut w: Box<dyn Write> = match &o.output {
        Some(p) => match std::fs::File::create(p) {
            Ok(f) => Box::new(BufWriter::new(f)),
            Err(e) => {
                eprintln!("merge: cannot write {p}: {e}");
                return 2;
            }
        },
        None => Box::new(BufWriter::new(std::io::stdout())),
    };
    match merge(&o, &mut *w) {
        Ok(()) => 0,
        Err(m) => {
            eprintln!("merge: {m}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;

    fn args(a: &[&str]) -> Vec<Vec<u8>> {
        a.iter().map(|s| s.as_bytes().to_vec()).collect()
    }

    #[test]
    fn natural_order() {
        assert_eq!(natural_cmp(b"chr2", b"chr10"), Ordering::Less);
        assert_eq!(natural_cmp(b"chr10", b"chr2"), Ordering::Greater);
        assert_eq!(natural_cmp(b"a01", b"a1"), Ordering::Less);
        assert_eq!(natural_cmp(b"chr1", b"chr1"), Ordering::Equal);
        assert_eq!(natural_cmp(b"chr1", b"chr1a"), Ordering::Less);
    }

    #[test]
    fn escapes() {
        assert_eq!(
            escape(b"a,b;c=d&e%f\tg"),
            b"a%2Cb%3Bc%3Dd%26e%25f%09g".to_vec()
        );
        assert_eq!(unescape(b"a%2Cb%zz%4"), b"a,b%zz%4".to_vec());
    }

    #[test]
    fn id_mapping() {
        let r = Rename {
            old_top: Some(b"DTC5"),
            new: b"X2",
        };
        let m = |v: &[u8]| {
            let mut o = Vec::new();
            r.map(v, &mut o);
            String::from_utf8(o).unwrap()
        };
        assert_eq!(m(b"DTC5"), "X2");
        assert_eq!(m(b"DTC5.te"), "X2.te");
        assert_eq!(m(b"DTC55"), "X2.DTC55");
        assert_eq!(m(b"other"), "X2.other");
    }

    #[test]
    fn rewrite_lines() {
        let r = Rename {
            old_top: Some(b"A1"),
            new: b"B7",
        };
        let mut o = Vec::new();
        let top = b"c\ts\tt\t1\t9\t.\t+\t.\tID=A1;K=v";
        rewrite(top, true, &r, b"f:A9", &mut o);
        assert_eq!(
            o,
            b"c\ts\tt\t1\t9\t.\t+\t.\tID=B7;Name=B7;K=v;Merged=f:A9\n"
        );
        o.clear();
        let sub = b"c\ts\tt\t1\t9\t.\t+\t.\tID=A1.te;Parent=A1,x;Name=A1;Merged=q";
        rewrite(sub, false, &r, b"f:A9", &mut o);
        assert_eq!(
            o,
            b"c\ts\tt\t1\t9\t.\t+\t.\tID=B7.te;Parent=B7,B7.x;Name=B7;Merged=q\n"
        );
        o.clear();
        rewrite(b"c\ts\tt\t1\t9\t.\t+\t.\t.", true, &r, b"", &mut o);
        assert_eq!(o, b"c\ts\tt\t1\t9\t.\t+\t.\tName=B7\n");
    }

    #[test]
    fn spans() {
        assert!(same_span((1, 10), (1, 10), None));
        assert!(!same_span((1, 10), (2, 10), None));
        assert!(same_span((1, 10), (2, 10), Some(0.9)));
        assert!(!same_span((1, 100), (1, 10), Some(0.9)));
        assert!(!same_span((1, 10), (20, 30), Some(0.1)));
    }

    #[test]
    fn prefixes() {
        assert_eq!(default_prefix(b"DTC5"), b"DTC");
        assert_eq!(default_prefix(b"DTH_98"), b"DTH_");
        assert_eq!(default_prefix(b"123"), b"sfm");
        assert_eq!(default_prefix(b""), b"sfm");
    }

    #[test]
    fn options() {
        let Ok(Parsed::Run(o)) =
            parse_args(&args(&["--overlap=0.9", "--name-start", "5", "-", "a"]))
        else {
            panic!()
        };
        assert_eq!((o.overlap, o.start, o.files.len()), (Some(0.9), 5, 2));
        assert_eq!(parse_args(&args(&["-h"])), Ok(Parsed::Help));
        for bad in [
            &["a", "--overlap", "0"][..],
            &["a", "--overlap=1.5"],
            &["a", "--name-start", "-1"],
            &["a", "--bogus"],
            &["a", "--output"],
            &[],
        ] {
            assert!(parse_args(&args(bad)).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn forward_reference_and_loop() {
        let ok = b"c\ts\tt\t2\t3\t.\t.\t.\tID=k;Parent=t\nc\ts\tt\t1\t9\t.\t.\t.\tID=t\n";
        let mut seqs = Seqs::default();
        let e = read_elems(0, "f", ok, &mut seqs).unwrap();
        assert_eq!((e.len(), e[0].lines.len(), e[0].top), (1, 2, 1));
        let lp = b"c\ts\tt\t2\t3\t.\t.\t.\tID=a;Parent=b\nc\ts\tt\t1\t9\t.\t.\t.\tID=b;Parent=a\n";
        assert!(
            read_elems(0, "f", lp, &mut seqs)
                .unwrap_err()
                .contains("loop")
        );
    }
}
