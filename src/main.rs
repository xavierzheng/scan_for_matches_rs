//! scan_for_matches -- Rust port of the program by Ross Overbeek (ANL).
//!
//!   scan_for_matches [-c] [-p] [-o N] [-i ignore_file] [-t threads] pattern_file < fasta_input > hits
//!   scan_for_matches [options] --input fasta_input --output hits pattern_file
//!
//! This file is the port of scan_for_matches.c (option handling, FASTA
//! reading and hit printing).  The pattern language lives in engine.rs.

mod engine;
mod fmt;
mod gz;
mod lint;
mod merge;
mod sys;

use engine::{Buf, DNA, Engine, PEPTIDE, SharedIndex, compl};
use std::io::Read;

const MAX_PAT_LINE_LN: usize = 32000;
const EOF: i32 = -1;

#[inline(always)]
fn isspace(c: i32) -> bool {
    matches!(c, 0x20 | 0x09 | 0x0a | 0x0b | 0x0c | 0x0d)
}

/// Characters left out of sequences.
#[inline(always)]
fn is_seq_space(b: u8) -> bool {
    b == b' ' || b == b'\n' || b == b'\r'
}

/// Byte reader with C `getc`/`ungetc` semantics.
struct Input<R: Read> {
    r: R,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
    back: Option<i32>,
    /// a damaged-input error (gzip); the input ends there
    err: Option<String>,
}

impl<R: Read> Input<R> {
    fn new(r: R) -> Input<R> {
        Input {
            r,
            buf: vec![0u8; 1 << 16],
            pos: 0,
            len: 0,
            back: None,
            err: None,
        }
    }

    fn fill(&mut self) -> bool {
        loop {
            match self.r.read(&mut self.buf) {
                Ok(0) => return false,
                Ok(n) => {
                    self.pos = 0;
                    self.len = n;
                    return true;
                }
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(e) if e.kind() == std::io::ErrorKind::InvalidData => {
                    // damaged gzip input: reported in input order
                    self.err = Some(e.to_string());
                    return false;
                }
                // like getc(): a read error ends the input
                Err(_) => return false,
            }
        }
    }

    #[inline(always)]
    fn getc(&mut self) -> i32 {
        if let Some(c) = self.back.take() {
            return c;
        }
        if self.pos >= self.len && !self.fill() {
            return EOF;
        }
        let c = self.buf[self.pos];
        self.pos += 1;
        c as i32
    }

    fn ungetc(&mut self, c: i32) {
        if c != EOF {
            self.back = Some(c);
        }
    }

    /// fscanf(fp, "%s", s): returns 1 on success, 0/EOF otherwise
    fn scan_s(&mut self, out: &mut Vec<u8>) -> i32 {
        out.clear();
        let mut c;
        loop {
            c = self.getc();
            if c == EOF {
                return EOF;
            }
            if !isspace(c) {
                break;
            }
        }
        loop {
            out.push(c as u8);
            c = self.getc();
            if c == EOF {
                break;
            }
            if isspace(c) {
                self.ungetc(c);
                break;
            }
        }
        1
    }

    /// fscanf(fp, ">%s", s)
    fn scan_gt_s(&mut self, out: &mut Vec<u8>) -> i32 {
        let c = self.getc();
        if c == EOF {
            return EOF;
        }
        if c != b'>' as i32 {
            self.ungetc(c);
            return 0;
        }
        self.scan_s(out)
    }

    /// Skip anything before the first line that starts with '>' (blank
    /// lines, comments).  The '>' is left unread.
    fn skip_to_first_header(&mut self) {
        let mut line_start = true;
        loop {
            let c = self.getc();
            if c == EOF {
                return;
            }
            if c == b'>' as i32 && line_start {
                self.ungetc(c);
                return;
            }
            line_start = c == b'\n' as i32 || c == b'\r' as i32;
        }
    }

    /// Sequence body: everything up to EOF or '>', dropping ' ', '\n' and
    /// '\r' (Windows line ends).
    /// Returns the char that stopped the scan.
    fn read_body(&mut self, out: &mut Vec<u8>) -> i32 {
        out.clear();
        if let Some(c) = self.back.take() {
            if c == b'>' as i32 {
                return c;
            }
            if !is_seq_space(c as u8) {
                out.push(c as u8);
            }
        }
        loop {
            if self.pos >= self.len && !self.fill() {
                return EOF;
            }
            let chunk = &self.buf[self.pos..self.len];
            match chunk.iter().position(|&b| b == b'>') {
                Some(k) => {
                    out.extend(chunk[..k].iter().copied().filter(|&b| !is_seq_space(b)));
                    self.pos += k + 1;
                    return b'>' as i32;
                }
                None => {
                    out.extend(chunk.iter().copied().filter(|&b| !is_seq_space(b)));
                    self.pos = self.len;
                }
            }
        }
    }
}

/// Where formatted hits go: stdout, or a buffer (threads).
trait Sink {
    fn write(&mut self, b: &[u8]);
    /// called after each complete hit; false: the hit was dropped
    /// (`--dedup`) and does not count for `-m`
    fn end_hit(&mut self) -> bool {
        true
    }
}

impl Sink for sys::Out {
    fn write(&mut self, b: &[u8]) {
        sys::Out::write(self, b)
    }
}

/// The hits of one record, formatted; `ends[k]` is where hit k ends.
#[derive(Default)]
struct HitBuf {
    buf: Vec<u8>,
    ends: Vec<usize>,
}

impl Sink for HitBuf {
    fn write(&mut self, b: &[u8]) {
        self.buf.extend_from_slice(b);
    }
    fn end_hit(&mut self) -> bool {
        self.ends.push(self.buf.len());
        true
    }
}

/// The final output: the hits as the C program prints them, or in one of
/// the formats of `--format` (each hit arrives as a raw record).
enum Writer {
    Plain(sys::Out),
    Fmt {
        out: sys::Out,
        f: Box<fmt::Formatter>,
        raw: Vec<u8>,
        text: Vec<u8>,
    },
}

impl Writer {
    fn new(f: Option<fmt::Formatter>) -> Writer {
        let mut out = sys::Out::new();
        match f {
            None => Writer::Plain(out),
            Some(f) => {
                let mut text = Vec::new();
                f.header(&mut text);
                out.write(&text);
                text.clear();
                Writer::Fmt {
                    out,
                    f: Box::new(f),
                    raw: Vec::new(),
                    text,
                }
            }
        }
    }

    /// After the last hit of a record (and before the program ends).
    fn end_record(&mut self) {
        if let Writer::Fmt { out, f, text, .. } = self {
            f.end_record(text);
            out.write(text);
            text.clear();
        }
    }
}

impl Sink for Writer {
    fn write(&mut self, b: &[u8]) {
        match self {
            Writer::Plain(out) => out.write(b),
            Writer::Fmt { raw, .. } => raw.extend_from_slice(b),
        }
    }
    fn end_hit(&mut self) -> bool {
        match self {
            Writer::Plain(_) => true,
            Writer::Fmt { out, f, raw, text } => {
                let counted = f.hit(raw, text);
                raw.clear();
                out.write(text);
                text.clear();
                counted
            }
        }
    }
}

struct Printer<S: Sink> {
    out: S,
    line: Vec<u8>,
    /// raw records for `--format` instead of the C output
    raw: bool,
    /// index of the pattern file (several patterns)
    pat: usize,
}

impl<S: Sink> Printer<S> {
    /// Print one hit the way the C code does: `>id:[a,b]` on one line,
    /// then every matched piece followed by a space, then a newline.
    #[allow(clippy::too_many_arguments)]
    fn hit(
        &mut self,
        id: &[u8],
        a: i64,
        b: i64,
        hits: &[i64],
        n: usize,
        data: &Buf,
        cdata: i64,
    ) -> bool {
        self.line.clear();
        self.line.push(b'>');
        self.line.extend_from_slice(id);
        self.line
            .extend_from_slice(format!(":[{},{}]\n", a, b).as_bytes());
        self.out.write(&self.line);
        self.line.clear();
        for i1 in 0..n {
            let from = hits[i1] - cdata;
            let to = from + (hits[i1 + 1] - hits[i1]).max(0);
            for p in from..to {
                self.line.push(data.get(p));
                if self.line.len() >= 1 << 16 {
                    self.out.write(&self.line);
                    self.line.clear();
                }
            }
            self.line.push(b' ');
        }
        self.line.push(b'\n');
        self.out.write(&self.line);
        self.out.end_hit()
    }
}

/// One step of reading the FASTA input.
enum ReadItem {
    Record {
        id: Vec<u8>,
        body: Vec<u8>,
    },
    /// the input cannot be read on: message for stderr, exit status 1
    Error(Vec<u8>),
    End,
}

/// Reads records the way the C program does (see scan_for_matches.c).
struct FastaReader {
    inp: Input<Box<dyn Read>>,
    got_gt: bool,
}

impl FastaReader {
    fn new(r: Box<dyn Read>) -> FastaReader {
        let mut inp = Input::new(r);
        inp.skip_to_first_header();
        FastaReader { inp, got_gt: false }
    }

    fn next(&mut self) -> ReadItem {
        let mut id = Vec::new();
        let r = if !self.got_gt {
            self.inp.scan_gt_s(&mut id)
        } else {
            self.inp.scan_s(&mut id)
        };
        if r != 1 {
            return self.end_or_error();
        }
        // skip the rest of the header line
        loop {
            let c = self.inp.getc();
            if c == b'\n' as i32 || c == EOF {
                break;
            }
        }
        let mut body = Vec::new();
        let stop = self.inp.read_body(&mut body);
        self.got_gt = stop == b'>' as i32;
        if let Some(e) = self.inp.err.take() {
            return ReadItem::Error(format!("{e}\n").into_bytes());
        }
        if body.len() > engine::MAX_SEQ_LEN {
            let mut msg = b"sequence ".to_vec();
            msg.extend_from_slice(&id);
            msg.extend_from_slice(
                format!(
                    " is too long (more than {} characters)\n",
                    engine::MAX_SEQ_LEN
                )
                .as_bytes(),
            );
            return ReadItem::Error(msg);
        }
        ReadItem::Record { id, body }
    }

    fn end_or_error(&mut self) -> ReadItem {
        match self.inp.err.take() {
            Some(e) => ReadItem::Error(format!("{e}\n").into_bytes()),
            None => ReadItem::End,
        }
    }
}

#[derive(Clone, Copy)]
struct Opts {
    protein: bool,
    complements: bool,
    show_overlaps: bool,
    /// `--format`: raw records, unit of each hit entry tracked
    raw: bool,
    /// `-t`: threads that decompress bgzip input (1: no extra threads)
    threads: usize,
}

/// Turn the record in `data` (length `ln`) into its reverse complement.
fn reverse_complement(data: &mut Buf, ln: usize) {
    let s = &mut data.v[..ln];
    s.reverse();
    for c in s {
        *c = compl(*c);
    }
}

/// Print the hit in `hits` (`n` units); on the reverse strand (`rev`) the
/// positions are counted on the forward strand, as in the C code.
/// Returns whether the hit counts for `-m`.
#[allow(clippy::too_many_arguments)]
fn print_hit<S: Sink>(
    pr: &mut Printer<S>,
    eng: &Engine,
    id: &[u8],
    rev: bool,
    ln: usize,
    hits: &[i64],
    n: usize,
    data: &Buf,
) -> bool {
    let cb = eng.cdata_base();
    if pr.raw {
        pr.line.clear();
        fmt::encode(
            &mut pr.line,
            pr.pat,
            rev,
            ln,
            id,
            hits,
            n,
            eng.hit_slots(),
            data,
            cb,
        );
        pr.out.write(&pr.line);
        return pr.out.end_hit();
    }
    if !rev {
        pr.hit(
            id,
            1 + hits[0] - cb,
            1 + (hits[n] - 1 - cb),
            hits,
            n,
            data,
            cb,
        )
    } else {
        let l = ln as i64;
        pr.hit(
            id,
            1 + (l - 1) - (hits[0] - cb),
            1 + (l - 1) - (hits[n] - 1 - cb),
            hits,
            n,
            data,
            cb,
        )
    }
}

/// Search one record (forward strand, then the reverse complement with
/// `-c`) and print its hits while `max_hits` allows.  Returns whether a
/// hit was printed.  This is the body of the C main loop.
#[allow(clippy::too_many_arguments)]
fn scan_record<S: Sink>(
    eng: &mut Engine,
    data: &mut Buf,
    hits: &mut Vec<i64>,
    pr: &mut Printer<S>,
    id: &[u8],
    body: &[u8],
    opts: Opts,
    max_hits: &mut i32,
) -> bool {
    // the persistent data buffer, NUL terminated
    data.store(body);
    let ln = body.iter().position(|&b| b == 0).unwrap_or(body.len());
    if !opts.protein {
        eng.comp_data(data);
    } else {
        eng.copy_data(data);
    }

    let mut hit_in_line = false;
    let mut i = eng.first_match(ln as i32, hits);
    while *max_hits > 0 && i > 0 {
        hit_in_line = true;
        if print_hit(pr, eng, id, false, ln, hits, i as usize, data) {
            *max_hits -= 1;
        }
        i = if !opts.show_overlaps {
            eng.cont_match(hits)
        } else {
            eng.next_match(hits)
        };
    }

    if opts.complements {
        reverse_complement(data, ln);
        eng.comp_data(data);

        let mut i = eng.first_match(ln as i32, hits);
        while *max_hits > 0 && i > 0 {
            hit_in_line = true;
            if print_hit(pr, eng, id, true, ln, hits, i as usize, data) {
                *max_hits -= 1;
            }
            i = eng.cont_match(hits);
        }
    }
    hit_in_line
}

/// Count a record without hits for `-n`.
fn missed(stop_after: &mut i32) {
    *stop_after = stop_after.wrapping_sub(1);
    if *stop_after == 0 {
        eprintln!("exceeded limit of lines failing to match");
        std::process::exit(1);
    }
}

fn read_error(msg: &[u8]) -> ! {
    use std::io::Write;
    let _ = std::io::stderr().write_all(msg);
    std::process::exit(1);
}

/// Several patterns: after searching a record at least this long, an
/// engine gives the memory of the record back.
const RELEASE_LEN: usize = 1 << 22;

/// One thread: records are read, searched (by each pattern in turn) and
/// printed one after the other.
fn run_sequential(
    mut engs: Vec<Engine>,
    out: Writer,
    input: Option<std::fs::File>,
    opts: Opts,
    ignore: &std::collections::HashSet<Vec<u8>>,
    mut max_hits: i32,
    mut stop_after: i32,
) {
    let mut data = Buf::new();
    let mut hits: Vec<i64> = vec![0; 2000];
    let mut pr = Printer {
        out,
        line: Vec::new(),
        raw: opts.raw,
        pat: 0,
    };
    for eng in &mut engs {
        eng.set_track_units(opts.raw);
    }
    let several = engs.len() > 1;
    let mut rd = FastaReader::new(fasta_source(input, opts.threads));
    while max_hits > 0 {
        match rd.next() {
            ReadItem::End => break,
            ReadItem::Error(msg) => read_error(&msg),
            ReadItem::Record { id, body } => {
                if ignore.contains(&id) {
                    continue;
                }
                let mut hit = false;
                for (k, eng) in engs.iter_mut().enumerate() {
                    if max_hits <= 0 {
                        break;
                    }
                    pr.pat = k;
                    hit |= scan_record(
                        eng,
                        &mut data,
                        &mut hits,
                        &mut pr,
                        &id,
                        &body,
                        opts,
                        &mut max_hits,
                    );
                    // (a pattern that reads matches of earlier records
                    // keeps them)
                    if several && body.len() >= RELEASE_LEN && !eng.uses_earlier_state() {
                        eng.release_data();
                    }
                }
                pr.out.end_record();
                if !hit {
                    missed(&mut stop_after);
                }
            }
        }
    }
}

/// Records of at least two pieces are searched in pieces of this many
/// start positions by several threads (when the pattern allows it).
/// Tests set smaller pieces with the environment variable SFM_PIECE.
fn piece_size() -> usize {
    std::env::var("SFM_PIECE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n >= 1)
        .unwrap_or(1 << 20)
}

/// Warm-up before each piece (see `PieceEngine::search`).  The writer
/// turns it on when it has to search pieces again itself too often
/// (long, dense hits: the chain of a worker meets the one-thread chain
/// only after many hits); then the reader also makes larger pieces, so
/// the warm-up costs about a quarter more work.  The environment
/// variable SFM_WARM=N (tests) fixes it to N positions.
struct Warm {
    /// positions searched before each piece
    len: std::sync::atomic::AtomicUsize,
    fixed: bool,
}

/// Largest warm-up.
const MAX_WARM: usize = 2 << 20;

impl Warm {
    fn new() -> Warm {
        let fixed = std::env::var("SFM_WARM")
            .ok()
            .and_then(|v| v.parse::<usize>().ok());
        Warm {
            len: std::sync::atomic::AtomicUsize::new(fixed.unwrap_or(0)),
            fixed: fixed.is_some(),
        }
    }

    fn get(&self) -> usize {
        self.len.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// Piece size for the next long record.
    fn piece(&self, base: usize) -> usize {
        if self.fixed {
            base
        } else {
            base.max(4 * self.get())
        }
    }
}

/// What the writer counts to set the warm-up.
#[derive(Default)]
struct WarmStats {
    /// positions in pieces, positions searched again, largest distance
    /// searched again in one piece
    seen: u64,
    again: u64,
    max_again: usize,
}

impl WarmStats {
    /// The writer joined a piece of `len` positions and searched `again`
    /// of them itself.  When more than 1/20 of all positions so far were
    /// searched again, the warm-up becomes twice the largest distance
    /// searched again in one piece (at most MAX_WARM; it never gets
    /// smaller).
    fn joined(&mut self, warm: &Warm, len: usize, again: usize) {
        if warm.fixed {
            return;
        }
        self.seen += len as u64;
        self.again += again as u64;
        self.max_again = self.max_again.max(again);
        if self.again * 20 > self.seen && self.seen >= 8 << 20 {
            let w = (2 * self.max_again).min(MAX_WARM).max(warm.get());
            warm.len.store(w, std::sync::atomic::Ordering::Relaxed);
        }
    }
}

/// A record.  A long record is searched in pieces; each strand (letters
/// and coded sequence) and its k-mer index are made once and shared by
/// all engines (and patterns).
struct LongRec {
    /// unique for each record
    serial: usize,
    id: Vec<u8>,
    body: Vec<u8>,
    ln: usize,
    index: [SharedIndex; 2],
    strands: [std::sync::OnceLock<std::sync::Arc<Strand>>; 2],
}

/// One strand of a long record: the letters (for the output) and the
/// coded sequence the engines search.
struct Strand {
    data: Buf,
    seq: engine::SharedSeq,
}

/// Work for one worker.  The jobs of a record come one after the other:
/// pattern by pattern, and for each pattern strand by strand.
enum Job {
    /// a whole record, searched by the patterns `pats` in turn
    Whole(WholeJob),
    Piece(PieceJob),
}

struct WholeJob {
    rec: std::sync::Arc<LongRec>,
    pats: Vec<usize>,
    /// the first / last job of the record
    first_of_record: bool,
    end_of_record: bool,
}

/// Hits of pattern `pat` starting at `from..=last` of one strand of a
/// long record.
struct PieceJob {
    rec: std::sync::Arc<LongRec>,
    pat: usize,
    rev: bool,
    from: usize,
    last: usize,
    /// the first piece of its strand (for this pattern)
    first: bool,
    /// the first / last job of the record
    first_of_record: bool,
    end_of_record: bool,
}

/// The hits a worker found in a piece, formatted, with their start and
/// end (offsets on the strand, end exclusive).  Without `-o` they are
/// the hits the one-thread search finds when it reaches the piece with
/// nothing to skip.
struct PieceHits {
    out: HitBuf,
    starts: Vec<i64>,
    pasts: Vec<i64>,
    /// the worker found `starts[0]` searching from here: no hit starts
    /// in `e0..starts[0]` (`from` without warm-up; i64::MAX: unknown)
    e0: i64,
}

/// What a worker or the reader sends to the writer, in input order.
enum Done {
    /// the hits of a whole job; first and last job of the record
    Hits(HitBuf, bool, bool),
    Piece(PieceJob, PieceHits),
    Error(Vec<u8>),
    End,
}

/// An engine with the strand of a long record it holds (loaded once for
/// all its pieces).
struct PieceEngine {
    eng: Engine,
    /// the strand loaded (shared with the other engines)
    strand: Option<std::sync::Arc<Strand>>,
    hits: Vec<i64>,
    /// serial number and strand of the record loaded
    loaded: Option<(usize, bool)>,
    /// the last hit `join_piece` found with this engine was not empty
    last_hit_nonempty: bool,
    /// `--format`: raw records
    raw: bool,
    /// index of the pattern file
    pat: usize,
}

impl PieceEngine {
    fn new(mut eng: Engine, raw: bool, threads: usize, pat: usize) -> PieceEngine {
        eng.set_track_units(raw);
        eng.set_index_threads(threads);
        PieceEngine {
            raw,
            pat,
            eng,
            strand: None,
            hits: vec![0; 2000],
            loaded: None,
            last_hit_nonempty: true,
        }
    }

    fn load(&mut self, rec: &LongRec, rev: bool) {
        let key = (rec.serial, rev);
        if self.loaded == Some(key) {
            return;
        }
        let eng = &self.eng;
        let st = rec.strands[rev as usize].get_or_init(|| {
            let mut data = Buf::new();
            data.store(&rec.body);
            if rev {
                reverse_complement(&mut data, rec.ln);
            }
            let seq = eng.shared_seq(&data);
            std::sync::Arc::new(Strand { data, seq })
        });
        self.eng.set_shared_seq(Some(&st.seq));
        self.eng
            .set_shared_index(Some(rec.index[rev as usize].clone()));
        self.strand = Some(std::sync::Arc::clone(st));
        self.loaded = Some(key);
    }

    /// The letters of the strand loaded.
    fn data(&self) -> &Buf {
        &self.strand.as_ref().expect("a strand is loaded").data
    }

    /// Before searching a whole record with this engine.
    fn unload(&mut self) {
        if self.loaded.is_some() {
            self.eng.set_shared_index(None);
            self.eng.set_shared_seq(None);
            self.strand = None;
            self.loaded = None;
        }
    }

    /// Give the memory of the record back (several patterns: each thread
    /// holds a long record in one engine at a time).
    fn release(&mut self) {
        self.unload();
        self.eng.release_data();
    }

    /// The hits of a piece: all of them (`all`, forward strand with -o),
    /// or the chain of non-overlapping hits from its first start.  With
    /// `warm` > 0 the chain starts `warm` positions before the piece (hits
    /// that start before it are not kept), so that it has most likely
    /// met the one-thread chain when the piece starts.
    fn search(&mut self, job: &PieceJob, all: bool, warm: usize) -> PieceHits {
        let rec = &*job.rec;
        self.load(rec, job.rev);
        let mut pr = Printer {
            out: HitBuf::default(),
            line: Vec::new(),
            raw: self.raw,
            pat: self.pat,
        };
        let (mut starts, mut pasts) = (Vec::new(), Vec::new());
        let cb = self.eng.cdata_base();
        let from = job.from as i64;
        let warm = if all || job.first { 0 } else { warm };
        let mut e0 = from;
        let mut i = self.eng.first_match_in(
            rec.ln as i32,
            job.from.saturating_sub(warm) as i64,
            job.last as i64,
            &mut self.hits,
        );
        while i > 0 {
            let n = i as usize;
            let (s, p) = (self.hits[0] - cb, self.hits[n] - cb);
            if s < from {
                // warm-up hit: the search goes on from its end (from its
                // start when it is empty: then e0 is not known)
                e0 = if p > s { p.max(from) } else { i64::MAX };
                i = self.eng.cont_match(&mut self.hits);
                continue;
            }
            print_hit(
                &mut pr,
                &self.eng,
                &rec.id,
                job.rev,
                rec.ln,
                &self.hits,
                n,
                self.data(),
            );
            starts.push(self.hits[0] - cb);
            pasts.push(self.hits[n] - cb);
            i = if all {
                self.eng.next_match(&mut self.hits)
            } else {
                self.eng.cont_match(&mut self.hits)
            };
        }
        PieceHits {
            out: pr.out,
            starts,
            pasts,
            e0,
        }
    }
}

/// Several threads: a reader thread, `threads` workers that each search
/// whole records, or pieces of long records, with their own engines (one
/// for each pattern), and this thread, which prints the results in input
/// order.  Hit limit (`-m`) and miss limit (`-n`) are applied in input
/// order, so the output equals the one-thread output.
#[allow(clippy::too_many_arguments)]
fn run_threads(
    mut out: Writer,
    input: Option<std::fs::File>,
    threads: usize,
    lines: Vec<Vec<u8>>,
    seq_type: i32,
    opts: Opts,
    split: Vec<bool>,
    ignore: std::collections::HashSet<Vec<u8>>,
    mut max_hits: i32,
    mut stop_after: i32,
) {
    use std::sync::mpsc;
    use std::sync::{Arc, Mutex};
    let n_pats = lines.len();
    let (job_tx, job_rx) = mpsc::sync_channel::<(usize, Job)>(2 * threads);
    // at most this many jobs between the reader and the writer: when an
    // early piece is slow, the finished later ones wait in memory
    let (ticket_tx, ticket_rx) = mpsc::sync_channel::<()>(4 * threads);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (done_tx, done_rx) = mpsc::channel::<(usize, Done)>();
    let split: Vec<bool> = split.iter().map(|&s| s && !opts.protein).collect();
    let base_piece = piece_size();
    let warm: Arc<Vec<Warm>> = Arc::new((0..n_pats).map(|_| Warm::new()).collect());
    let mut wstats: Vec<WarmStats> = (0..n_pats).map(|_| WarmStats::default()).collect();

    let tx = done_tx.clone();
    let rwarm = Arc::clone(&warm);
    std::thread::spawn(move || {
        let mut rd = FastaReader::new(fasta_source(input, opts.threads));
        let mut k = 0usize;
        let mut serial = 0usize;
        loop {
            match rd.next() {
                ReadItem::Record { id, body } => {
                    if ignore.contains(&id) {
                        continue;
                    }
                    let ln = body.iter().position(|&b| b == 0).unwrap_or(body.len());
                    let rec = Arc::new(LongRec {
                        serial,
                        id,
                        body,
                        ln,
                        index: Default::default(),
                        strands: Default::default(),
                    });
                    serial += 1;
                    let jobs = record_jobs(&rec, &split, &rwarm, base_piece, opts.complements);
                    for job in jobs {
                        if ticket_tx.send(()).is_err() || job_tx.send((k, job)).is_err() {
                            return;
                        }
                        k += 1;
                    }
                }
                ReadItem::Error(msg) => {
                    let _ = tx.send((k, Done::Error(msg)));
                    return;
                }
                ReadItem::End => {
                    let _ = tx.send((k, Done::End));
                    return;
                }
            }
        }
    });

    let lines = Arc::new(lines);
    let new_engines = move |lines: &[Vec<u8>]| -> Vec<PieceEngine> {
        lines
            .iter()
            .enumerate()
            .map(|(k, line)| {
                let mut eng = Engine::new();
                eng.parse_cmd(line, seq_type);
                PieceEngine::new(eng, opts.raw, opts.threads, k)
            })
            .collect()
    };
    for _ in 0..threads {
        let rx = Arc::clone(&job_rx);
        let tx = done_tx.clone();
        let lines = Arc::clone(&lines);
        let warm = Arc::clone(&warm);
        std::thread::Builder::new()
            .stack_size(1 << 30) // the parser recurses as deeply as the C code
            .spawn(move || {
                let mut pes = new_engines(&lines);
                // the engine that holds a long record (several patterns)
                let mut held: Option<usize> = None;
                let mut data = Buf::new();
                let mut hits: Vec<i64> = vec![0; 2000];
                loop {
                    let job = rx.lock().unwrap().recv();
                    let Ok((k, job)) = job else { return };
                    let done = match job {
                        Job::Whole(w) => {
                            let mut pr = Printer {
                                out: HitBuf::default(),
                                line: Vec::new(),
                                raw: opts.raw,
                                pat: 0,
                            };
                            for &p in &w.pats {
                                let pe = &mut pes[p];
                                pe.unload();
                                pr.pat = p;
                                let mut unlimited = i32::MAX;
                                scan_record(
                                    &mut pe.eng,
                                    &mut data,
                                    &mut hits,
                                    &mut pr,
                                    &w.rec.id,
                                    &w.rec.body,
                                    opts,
                                    &mut unlimited,
                                );
                                if n_pats > 1 && w.rec.ln >= RELEASE_LEN {
                                    pe.release();
                                }
                            }
                            Done::Hits(pr.out, w.first_of_record, w.end_of_record)
                        }
                        Job::Piece(job) => {
                            if let Some(h) = held.filter(|&h| h != job.pat) {
                                pes[h].release();
                            }
                            if n_pats > 1 {
                                held = Some(job.pat);
                            }
                            let all = opts.show_overlaps && !job.rev;
                            let ph = pes[job.pat].search(&job, all, warm[job.pat].get());
                            Done::Piece(job, ph)
                        }
                    };
                    if tx.send((k, done)).is_err() {
                        return;
                    }
                }
            })
            .expect("thread");
    }
    drop(done_tx);

    let mut waiting = std::collections::BTreeMap::new();
    let mut next = 0usize;
    // a long record: end of the last hit printed on this strand, and
    // whether the record has a hit; engines for searching again (one for
    // each pattern, made when needed)
    let mut past = 0i64;
    let mut rec_hit = false;
    let mut weng: Vec<Option<PieceEngine>> = (0..n_pats).map(|_| None).collect();
    let mut wheld: Option<usize> = None;
    let mut wpr = Printer {
        out: HitBuf::default(),
        line: Vec::new(),
        raw: opts.raw,
        pat: 0,
    };
    while max_hits > 0 {
        let item = match waiting.remove(&next) {
            Some(d) => d,
            None => match done_rx.recv() {
                Ok((k, d)) => {
                    waiting.insert(k, d);
                    continue;
                }
                Err(_) => break,
            },
        };
        next += 1;
        if matches!(item, Done::Hits(..) | Done::Piece(..)) {
            let _ = ticket_rx.try_recv();
        }
        let end_of_record = match item {
            Done::End => break,
            Done::Error(msg) => {
                out.end_record();
                read_error(&msg)
            }
            Done::Hits(h, first, end) => {
                if first {
                    rec_hit = false;
                }
                rec_hit |= print_hits(&mut out, &h, 0, &mut max_hits) > 0;
                end
            }
            Done::Piece(job, ph) => {
                if job.first_of_record {
                    rec_hit = false;
                }
                if job.first {
                    past = 0;
                }
                if opts.show_overlaps && !job.rev {
                    // every hit, in search order
                    rec_hit |= print_hits(&mut out, &ph.out, 0, &mut max_hits) > 0;
                } else {
                    if let Some(h) = wheld.filter(|&h| h != job.pat)
                        && let Some(we) = weng[h].as_mut()
                    {
                        we.release();
                    }
                    if n_pats > 1 {
                        wheld = Some(job.pat);
                    }
                    let we = weng[job.pat].get_or_insert_with(|| {
                        let mut e = Engine::new();
                        e.parse_cmd(&lines[job.pat], seq_type);
                        PieceEngine::new(e, opts.raw, opts.threads, job.pat)
                    });
                    wpr.pat = job.pat;
                    let mut again = 0;
                    rec_hit |= join_piece(
                        &mut out,
                        &job,
                        &ph,
                        we,
                        &mut wpr,
                        &mut past,
                        &mut max_hits,
                        &mut again,
                    );
                    wstats[job.pat].joined(&warm[job.pat], job.last + 1 - job.from, again);
                }
                job.end_of_record
            }
        };
        if end_of_record {
            out.end_record();
            if !rec_hit && max_hits > 0 {
                missed(&mut stop_after);
            }
        }
    }
    out.end_record();
}

/// The jobs of a record, in output order: for each pattern, the whole
/// record (patterns next to each other that do not split share one job)
/// or its pieces, strand by strand.
fn record_jobs(
    rec: &std::sync::Arc<LongRec>,
    split: &[bool],
    warm: &[Warm],
    base_piece: usize,
    complements: bool,
) -> Vec<Job> {
    let ln = rec.ln;
    let mut jobs: Vec<Job> = Vec::new();
    let mut whole: Vec<usize> = Vec::new();
    let flush = |whole: &mut Vec<usize>, jobs: &mut Vec<Job>| {
        if !whole.is_empty() {
            jobs.push(Job::Whole(WholeJob {
                rec: std::sync::Arc::clone(rec),
                pats: std::mem::take(whole),
                first_of_record: false,
                end_of_record: false,
            }));
        }
    };
    for (pat, &sp) in split.iter().enumerate() {
        let piece = warm[pat].piece(base_piece);
        if !sp || ln < 2 * piece {
            whole.push(pat);
            continue;
        }
        flush(&mut whole, &mut jobs);
        let strands: &[bool] = if complements {
            &[false, true]
        } else {
            &[false]
        };
        for &rev in strands {
            let mut from = 0;
            while from <= ln {
                // a hit of length 0 can start just past the end
                let mut last = from + piece - 1;
                if last + 1 >= ln {
                    last = ln;
                }
                jobs.push(Job::Piece(PieceJob {
                    rec: std::sync::Arc::clone(rec),
                    pat,
                    rev,
                    from,
                    last,
                    first: from == 0,
                    first_of_record: false,
                    end_of_record: false,
                }));
                from = last + 1;
            }
        }
    }
    flush(&mut whole, &mut jobs);
    let n = jobs.len();
    for (i, j) in jobs.iter_mut().enumerate() {
        let (f, e) = match j {
            Job::Whole(w) => (&mut w.first_of_record, &mut w.end_of_record),
            Job::Piece(p) => (&mut p.first_of_record, &mut p.end_of_record),
        };
        *f = i == 0;
        *e = i + 1 == n;
    }
    jobs
}

/// Print the hits of `h` from hit `from` on while `max_hits` allows;
/// returns how many were printed (or dropped by `--dedup`).
fn print_hits(out: &mut Writer, h: &HitBuf, from: usize, max_hits: &mut i32) -> usize {
    let mut start = if from == 0 { 0 } else { h.ends[from - 1] };
    let mut printed = 0;
    for &end in &h.ends[from..] {
        if *max_hits <= 0 {
            break;
        }
        out.write(&h.buf[start..end]);
        start = end;
        if out.end_hit() {
            *max_hits -= 1;
        }
        printed += 1;
    }
    printed
}

/// Continue the one-thread search ("next hit that starts at or after
/// the end of the last one") through a piece.  `past` is the end of the
/// last hit printed on the strand.  The worker's hits are right from the
/// first one that the true search reaches; before that the piece is
/// searched again from `past` with `we`.  Returns whether a hit was
/// printed; `searched` is set to the number of positions of the piece
/// searched again.
#[allow(clippy::too_many_arguments)]
fn join_piece(
    out: &mut Writer,
    job: &PieceJob,
    ph: &PieceHits,
    we: &mut PieceEngine,
    wpr: &mut Printer<HitBuf>,
    past: &mut i64,
    max_hits: &mut i32,
    searched: &mut usize,
) -> bool {
    let last = job.last as i64;
    let from = job.from as i64;
    let mut printed = false;
    let mut again = false; // `we` holds a search of this piece
    let mut again_from = 0i64;
    loop {
        if again {
            *searched = ((*past).min(last + 1) - again_from).max(0) as usize;
        }
        if *past > last || *max_hits <= 0 {
            return printed;
        }
        // the worker's first hit at or after `past`
        let j = ph.starts.partition_point(|&s| s < *past);
        // the worker searched on from the end of hit j-1 (from `e0` for
        // its first hit), at or before `past`: from here on its hits are
        // the true ones.  (When the last hit has length 0 the search goes
        // on at its own start, so only an ending after the start is
        // used.)  No hit starts between `past` and `from` (the earlier
        // pieces were joined up to there).
        let prev = if j == 0 { ph.e0 } else { ph.pasts[j - 1] };
        let synced = prev <= (*past).max(from);
        if synced && (!again || we.last_hit_nonempty) {
            if j < ph.starts.len() {
                let n = print_hits(out, &ph.out, j, max_hits);
                if n > 0 {
                    printed = true;
                    *past = ph.pasts[j + n - 1];
                }
            }
            return printed;
        }
        // search the piece again from `past` (from its start when `past`
        // is before it: hits that start earlier belong to the earlier
        // pieces, also a hit of length 0 at `past`)
        let i = if !again {
            again = true;
            again_from = (*past).max(from);
            we.load(&job.rec, job.rev);
            we.eng
                .first_match_in(job.rec.ln as i32, again_from, last, &mut we.hits)
        } else {
            we.eng.cont_match(&mut we.hits)
        };
        if i <= 0 {
            // searched to the end of the piece
            *searched = (last + 1 - again_from).max(0) as usize;
            return printed;
        }
        let n = i as usize;
        let cb = we.eng.cdata_base();
        wpr.out.buf.clear();
        wpr.out.ends.clear();
        print_hit(
            wpr,
            &we.eng,
            &job.rec.id,
            job.rev,
            job.rec.ln,
            &we.hits,
            n,
            we.data(),
        );
        out.write(&wpr.out.buf);
        if out.end_hit() {
            *max_hits -= 1;
        }
        printed = true;
        let start = we.hits[0] - cb;
        *past = we.hits[n] - cb;
        we.last_hit_nonempty = *past > start;
    }
}

/// An error in the long options or the labels of the pattern file.
const HELP: &str = "\
scan_for_matches VERSION: search DNA or protein sequences for a pattern

usage:
  scan_for_matches [options] pattern_file < fasta_input > hits
  scan_for_matches [options] --input fasta_input --output hits pattern_file
  scan_for_matches --format F [options] pattern_file pattern_file ... < fasta_input
  scan_for_matches --lint [-p] [-c] [--dedup] pattern_file ...
  scan_for_matches --explain [-p] pattern_file ...
  scan_for_matches --merge [merge options] run1.gff3 run2.gff3 ... (--merge --help)

The FASTA input can be plain text, gzip or bgzip (detected).

options:
  -c                 also search the opposite strand (reverse complement);
                     not with -p
  -p                 protein sequences
  -n N               stop (exit status 1) after N sequences without a hit
  -m N               report at most N hits
  -o N               show overlapping hits (the value is not used, but it
                     must be there)
  -i FILE            file of sequence ids to skip
  -t N               use N threads (default 1); the output is the same as
                     with one thread; bgzip input is decompressed by N threads
  --input FILE       read the FASTA input from FILE (default: stdin; -: stdin)
  --output FILE      write the hits to FILE (default: stdout; -: stdout)
  -h, --help         show this help and exit

output formats (default: the original output):
  --format F         gff3, bed6, bed12 or jsonl
  --name-prefix P    Name of each hit = P + number (default: the Name= of the
                     %@element label, else sfm)
  --name-start N     first number (default 1)
  --type T           column 3 of the element (default: the %@element type,
                     else sequence_motif)
  --dedup            with -c: a reverse-strand hit with the same span as a
                     forward hit is dropped; the element gets strand .

several pattern files (needs --format): the input is read once; each record
is searched by each pattern in turn.  Names with the same prefix are numbered
by one counter, so they stay unique.

checking patterns (no input is read):
  --lint             report traps in the pattern (errors, warnings, notes);
                     exit status 1 when there is an error
  --explain          describe each unit of the pattern in plain words

Long options can be written --format gff3 or --format=gff3.
Note: -i and -o are the options of the original program (ids to skip,
overlapping hits), not input and output.

example:
  scan_for_matches -t 8 -c --dedup --format gff3 --name-prefix DTC \\
      --input genome.fna.gz --output DTC.gff3 DTC.pat

Pattern language and labels: see README.md.
";

fn print_help() -> ! {
    print!("{}", HELP.replace("VERSION", env!("CARGO_PKG_VERSION")));
    std::process::exit(0);
}

fn long_error(msg: &str) -> ! {
    eprintln!("scan_for_matches: {msg}");
    std::process::exit(2);
}

fn usage(errflag: i32, optind: i32, argc: i32) -> ! {
    eprintln!("errflag={} optind={} argc={}", errflag, optind, argc);
    eprintln!(
        "usage: scan_for_matches -c [for complementary strand] -p [protein][-n stop_after_n_misses] [-m max_hits] [-i file_of_ids_to_ignore] -o [overlapping hits] [-t threads] pattern_file < fasta_input > hits"
    );
    std::process::exit(2);
}

/// The FASTA input: the `--input` file, else stdin.
fn fasta_source(input: Option<std::fs::File>, threads: usize) -> Box<dyn Read> {
    match input {
        Some(f) => open_fasta_input(std::io::BufReader::with_capacity(1 << 20, f), threads),
        None => open_fasta_input(std::io::stdin().lock(), threads),
    }
}

/// The FASTA input: plain text, or gzip / bgzip compressed (detected by
/// the first two bytes).  With `threads` > 1, bgzip blocks are
/// decompressed by that many threads.
fn open_fasta_input<R: Read + 'static>(mut r: R, threads: usize) -> Box<dyn Read> {
    let mut head = Vec::with_capacity(2);
    while head.len() < 2 {
        let mut b = [0u8; 1];
        match r.read(&mut b) {
            Ok(0) => break,
            Ok(_) => head.push(b[0]),
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(_) => break,
        }
    }
    let r = std::io::Cursor::new(head.clone()).chain(r);
    if gz::is_gzip(&head) {
        let g: std::io::Result<Box<dyn Read>> = if threads > 1 {
            gz::BgzfReader::new(r, threads).map(|g| Box::new(g) as Box<dyn Read>)
        } else {
            gz::GzReader::new(r).map(|g| Box::new(g) as Box<dyn Read>)
        };
        match g {
            Ok(g) => g,
            Err(e) => {
                eprintln!("{e}");
                std::process::exit(1);
            }
        }
    } else {
        Box::new(r)
    }
}

fn open_path(p: &[u8]) -> Option<std::fs::File> {
    use std::os::unix::ffi::OsStrExt;
    std::fs::File::open(std::ffi::OsStr::from_bytes(p)).ok()
}

fn main() {
    // The parser and max_mat() recurse as deeply as the C code does; give
    // them a large stack so only a truly endless recursion (which crashes
    // the C program too) ends in SIGSEGV.
    let t = std::thread::Builder::new()
        .stack_size(1 << 30)
        .spawn(real_main)
        .expect("thread");
    let _ = t.join();
}

/// A pattern file as read: the joined line, the source line of each of
/// its bytes, and the comments (for the `%@` labels).
struct PatText {
    file: Vec<u8>,
    line: Vec<u8>,
    src_line: Vec<u32>,
    comments: Vec<(u32, Vec<u8>)>,
    /// text after byte 31 999 was lost
    truncated: bool,
}

/// Read a pattern: '%' starts a comment, newlines become spaces.
fn read_pattern(file: Vec<u8>, f: std::fs::File) -> PatText {
    let mut line: Vec<u8> = Vec::new();
    let mut src_line: Vec<u32> = Vec::new();
    let mut comments: Vec<(u32, Vec<u8>)> = Vec::new();
    let mut fp = Input::new(f);
    let mut n_line = 1u32;
    let mut i = fp.getc();
    while i != EOF && line.len() < MAX_PAT_LINE_LN - 1 {
        if i == b'%' as i32 {
            let mut text = Vec::new();
            loop {
                i = fp.getc();
                if i == b'\n' as i32 || i == EOF {
                    break;
                }
                text.push(i as u8);
            }
            if text.last() == Some(&b'\r') {
                text.pop();
            }
            comments.push((n_line, text));
            if i != EOF {
                n_line += 1;
                i = b' ' as i32;
            }
        } else if i == b'\n' as i32 || i == b'\r' as i32 {
            if i == b'\n' as i32 {
                n_line += 1;
            }
            i = b' ' as i32;
        } else {
            line.push(i as u8);
            src_line.push(n_line);
            i = fp.getc();
        }
    }
    // text after the cut (other than spaces and comments) is lost
    let mut truncated = false;
    let mut in_comment = false;
    while i != EOF {
        if i == b'\n' as i32 {
            in_comment = false;
        } else if i == b'%' as i32 {
            in_comment = true;
        } else if !in_comment && !isspace(i) {
            truncated = true;
            break;
        }
        i = fp.getc();
    }
    // the C string ends at the first NUL
    if let Some(k) = line.iter().position(|&b| b == 0) {
        line.truncate(k);
    }
    PatText {
        file,
        line,
        src_line,
        comments,
        truncated,
    }
}

/// `--lint` / `--explain` for each pattern file; exits.
fn lint_patterns(pats: &[PatText], lint: bool, ro: &lint::RunOpts) -> ! {
    let mut error = false;
    for (k, p) in pats.iter().enumerate() {
        let mut eng = Engine::new();
        let parsed = eng.parse_cmd(&p.line, if ro.protein { PEPTIDE } else { DNA }) != 0;
        let facts = lint::Facts {
            parsed,
            uses_earlier_state: parsed && eng.uses_earlier_state(),
            can_split: parsed && eng.can_split(),
            units: if parsed {
                eng.unit_offsets()
            } else {
                Vec::new()
            },
        };
        let src = lint::Source {
            file: &p.file,
            line: &p.line,
            src_line: &p.src_line,
            comments: &p.comments,
            truncated: p.truncated,
        };
        let (text, err) = if lint {
            lint::lint(&src, &facts, ro)
        } else {
            lint::explain(&src, &facts, ro)
        };
        if k > 0 {
            println!();
        }
        print!("{text}");
        error |= err;
    }
    use std::io::Write;
    let _ = std::io::stdout().flush();
    std::process::exit(error as i32);
}

fn real_main() {
    let mut args = sys::Args::from_env();
    if args.argc() > 1 && args.get(1) == b"--merge" {
        let rest: Vec<Vec<u8>> = (2..args.argc() as usize).map(|i| args.get(i)).collect();
        std::process::exit(merge::main(&rest));
    }
    // long options (output formats); the C options are read by getopt
    let long = args
        .take_long(
            &[
                "format",
                "name-prefix",
                "name-start",
                "type",
                "input",
                "output",
            ],
            &["dedup", "help", "lint", "explain"],
        )
        .unwrap_or_else(|e| long_error(&e));
    // -h only as an option of its own, not as the value of an option
    // (the C program ignored -h)
    let short_h = (1..args.argc() as usize).any(|i| {
        let before = if i > 1 { args.get(i - 1) } else { Vec::new() };
        let takes_value = matches!(before.as_slice(), b"-n" | b"-m" | b"-o" | b"-i" | b"-t");
        args.get(i) == b"-h" && !takes_value
    });
    if short_h || long.iter().any(|(n, _)| n == "help") {
        print_help();
    }
    let argc = args.argc();
    let mut fopts = fmt::Options {
        format: fmt::Format::Gff3,
        name_prefix: None,
        name_start: 1,
        typ: None,
        dedup: false,
    };
    let mut format = None;
    let mut input_path: Option<Vec<u8>> = None;
    let mut output_path: Option<Vec<u8>> = None;
    let mut fmt_opts = 0;
    let (mut do_lint, mut do_explain) = (false, false);
    for (name, v) in &long {
        if !matches!(name.as_str(), "input" | "output" | "lint" | "explain") {
            fmt_opts += 1;
        }
        match name.as_str() {
            "input" => input_path = Some(v.clone()).filter(|v| v != b"-"),
            "output" => output_path = Some(v.clone()).filter(|v| v != b"-"),
            "format" => {
                format = Some(
                    fmt::Format::parse(v)
                        .unwrap_or_else(|| long_error("--format is gff3, bed6, bed12 or jsonl")),
                )
            }
            "name-prefix" => fopts.name_prefix = Some(v.clone()),
            "name-start" => {
                fopts.name_start = std::str::from_utf8(v)
                    .ok()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or_else(|| long_error("--name-start needs a number"))
            }
            "type" => fopts.typ = Some(v.clone()),
            "lint" => do_lint = true,
            "explain" => do_explain = true,
            _ => fopts.dedup = true,
        }
    }
    if do_lint && do_explain {
        long_error("use --lint or --explain, not both");
    }
    let checking = do_lint || do_explain;
    if format.is_none() && fmt_opts > 0 && !checking {
        long_error("--dedup, --name-prefix, --name-start and --type need --format");
    }

    let mut stop_after: i32 = 100_000_000;
    let mut max_hits: i32 = 100_000_000;
    let mut show_overlaps = false;
    let mut o_value: Option<Vec<u8>> = None;
    let mut complements = false;
    let mut protein = false;
    let mut errflag: i32 = 0;
    let mut ig_fp: Option<std::fs::File> = None;
    let mut threads: usize = 1;

    loop {
        let c = args.getopt(b"pcn:m:o:i:t:\0");
        if c == -1 {
            break;
        }
        match c as u8 {
            b'o' => {
                show_overlaps = true;
                o_value = sys::optarg_bytes();
            }
            b'c' => complements = true,
            b'p' => protein = true,
            b'n' => {
                if sys::sscanf_optarg_int(&mut stop_after) != 1 {
                    errflag = 1;
                    eprintln!("invalid value on -n option (make it a positive integer)");
                }
            }
            b'm' => {
                if sys::sscanf_optarg_int(&mut max_hits) != 1 {
                    errflag = 1;
                    eprintln!("invalid value on -m option (make it a positive integer)");
                }
            }
            b't' => {
                match sys::optarg_bytes()
                    .and_then(|v| std::str::from_utf8(&v).ok()?.trim().parse::<usize>().ok())
                {
                    Some(n) if n >= 1 => threads = n,
                    _ => {
                        errflag = 1;
                        eprintln!("invalid value on -t option (make it a positive integer)");
                    }
                }
            }
            b'i' => {
                ig_fp = sys::optarg_bytes().and_then(|p| open_path(&p));
                if ig_fp.is_none() {
                    errflag = 1;
                    eprintln!("invalid file name for ids to ignore");
                }
            }
            _ => {}
        }
    }

    let optind = sys::optind_get();
    if errflag != 0 || optind >= argc {
        usage(errflag, optind, argc);
    }
    if complements && protein {
        eprintln!("-c (complementary strand) cannot be used with -p (protein sequences)");
        usage(errflag, optind, argc);
    }
    // one or more pattern files
    let pat_files: Vec<(Vec<u8>, std::fs::File)> = (optind as usize..argc as usize)
        .map(|i| {
            let name = args.get(i);
            match open_path(&name) {
                Some(f) => (name, f),
                None if optind as usize + 1 == argc as usize => usage(errflag, optind, argc),
                None => long_error(&format!(
                    "cannot open pattern file {}",
                    String::from_utf8_lossy(&name)
                )),
            }
        })
        .collect();
    let several = pat_files.len() > 1;
    if several && format.is_none() && !checking {
        long_error("several pattern files need --format");
    }
    let pats: Vec<PatText> = pat_files
        .into_iter()
        .map(|(name, f)| read_pattern(name, f))
        .collect();

    if checking {
        let ro = lint::RunOpts {
            protein,
            complements,
            dedup: fopts.dedup,
            format: format.is_some(),
            o_value,
        };
        lint_patterns(&pats, do_lint, &ro);
    }

    // ids to ignore
    let mut ignore: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
    let mut n_ignore = 0usize;
    if let Some(f) = ig_fp {
        let mut inp = Input::new(f);
        let mut id = Vec::new();
        while inp.scan_s(&mut id) == 1 {
            ignore.insert(id.clone());
            n_ignore += 1;
        }
        if n_ignore > 0 {
            eprintln!("ignoring {} id(s)", n_ignore);
        }
    }

    let seq_type = if protein { PEPTIDE } else { DNA };
    let mut engs: Vec<Engine> = Vec::new();
    for p in &pats {
        let mut eng = Engine::new();
        let rc = eng.parse_cmd(&p.line, seq_type);
        if rc == 0 {
            let mut msg = b"failed to parse pattern: ".to_vec();
            if several {
                msg.extend_from_slice(&p.file);
                msg.extend_from_slice(b": ");
            }
            msg.extend_from_slice(&p.line);
            msg.push(b'\n');
            use std::io::Write;
            let _ = std::io::stderr().write_all(&msg);
            std::process::exit(1);
        }
        engs.push(eng);
    }

    let formatter = format.map(|f| {
        fopts.format = f;
        let labels = pats
            .iter()
            .zip(&engs)
            .map(|(p, eng)| {
                let err = |e: String| -> ! {
                    if several {
                        long_error(&format!("{}: {e}", String::from_utf8_lossy(&p.file)))
                    } else {
                        long_error(&e)
                    }
                };
                let lb = fmt::parse_labels(&p.comments).unwrap_or_else(|e| err(e));
                let units: Vec<(u32, u32)> = eng
                    .unit_offsets()
                    .into_iter()
                    .map(|(slot, off)| (slot, p.src_line.get(off).copied().unwrap_or(0)))
                    .collect();
                let (group, line) = lb.assign(&units).unwrap_or_else(|e| err(e));
                fmt::PatLabels {
                    file: p.file.clone(),
                    lb,
                    group,
                    line,
                }
            })
            .collect();
        fmt::Formatter::new(fopts, labels)
    });
    if let Some(f) = &formatter {
        fmt::set_text_slots(f.text_slots());
    }
    let opts = Opts {
        protein,
        complements,
        show_overlaps,
        raw: formatter.is_some(),
        threads,
    };
    let input = input_path.map(|p| {
        open_path(&p).unwrap_or_else(|| {
            long_error(&format!(
                "cannot open input file {}",
                String::from_utf8_lossy(&p)
            ))
        })
    });
    if let Some(p) = output_path {
        use std::os::unix::ffi::OsStrExt;
        let name = String::from_utf8_lossy(&p).into_owned();
        let f = std::fs::File::create(std::ffi::OsStr::from_bytes(&p))
            .unwrap_or_else(|e| long_error(&format!("cannot write output file {name}: {e}")));
        if !sys::redirect_stdout(f) {
            long_error(&format!("cannot write output file {name}"));
        }
    }
    let out = Writer::new(formatter);
    // patterns that use matches of earlier sequences stay on one thread
    if threads > 1 && !engs.iter().any(|e| e.uses_earlier_state()) && max_hits > 0 {
        let split = engs.iter().map(|e| e.can_split()).collect();
        run_threads(
            out,
            input,
            threads,
            pats.into_iter().map(|p| p.line).collect(),
            seq_type,
            opts,
            split,
            ignore,
            max_hits,
            stop_after,
        );
    } else {
        run_sequential(engs, out, input, opts, &ignore, max_hits, stop_after);
    }
}
