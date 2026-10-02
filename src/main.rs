//! scan_for_matches -- Rust port of the program by Ross Overbeek (ANL).
//!
//!   scan_for_matches [-c] [-p] [-o N] [-i ignore_file] [-t threads] pattern_file < fasta_input > hits
//!
//! This file is the port of scan_for_matches.c (option handling, FASTA
//! reading and hit printing).  The pattern language lives in engine.rs.

mod engine;
mod gz;
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
    /// called after each complete hit
    fn end_hit(&mut self) {}
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
    fn end_hit(&mut self) {
        self.ends.push(self.buf.len());
    }
}

struct Printer<S: Sink> {
    out: S,
    line: Vec<u8>,
}

impl<S: Sink> Printer<S> {
    /// Print one hit the way the C code does: `>id:[a,b]` on one line,
    /// then every matched piece followed by a space, then a newline.
    #[allow(clippy::too_many_arguments)]
    fn hit(&mut self, id: &[u8], a: i64, b: i64, hits: &[i64], n: usize, data: &Buf, cdata: i64) {
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
        self.out.end_hit();
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
}

/// Turn the record in `data` (length `ln`) into its reverse complement.
fn reverse_complement(data: &mut Buf, ln: usize) {
    if ln > 0 {
        let (mut a, mut b) = (0usize, ln - 1);
        while a <= b {
            let tmp = compl(data.v[a]);
            data.v[a] = compl(data.v[b]);
            data.v[b] = tmp;
            a += 1;
            if b == 0 {
                break;
            }
            b -= 1;
        }
    }
}

/// Print the hit in `hits` (`n` units); on the reverse strand (`rev`) the
/// positions are counted on the forward strand, as in the C code.
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
) {
    let cb = eng.cdata_base();
    if !rev {
        pr.hit(
            id,
            1 + hits[0] - cb,
            1 + (hits[n] - 1 - cb),
            hits,
            n,
            data,
            cb,
        );
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
        );
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
        *max_hits -= 1;
        print_hit(pr, eng, id, false, ln, hits, i as usize, data);
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
            *max_hits -= 1;
            print_hit(pr, eng, id, true, ln, hits, i as usize, data);
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

/// One thread: records are read, searched and printed one after the other.
fn run_sequential(
    mut eng: Engine,
    opts: Opts,
    ignore: &std::collections::HashSet<Vec<u8>>,
    mut max_hits: i32,
    mut stop_after: i32,
) {
    let mut data = Buf::new();
    let mut hits: Vec<i64> = vec![0; 2000];
    let mut pr = Printer {
        out: sys::Out::new(),
        line: Vec::new(),
    };
    let mut rd = FastaReader::new(open_fasta_input(std::io::stdin().lock()));
    while max_hits > 0 {
        match rd.next() {
            ReadItem::End => break,
            ReadItem::Error(msg) => read_error(&msg),
            ReadItem::Record { id, body } => {
                if ignore.contains(&id) {
                    continue;
                }
                if !scan_record(
                    &mut eng,
                    &mut data,
                    &mut hits,
                    &mut pr,
                    &id,
                    &body,
                    opts,
                    &mut max_hits,
                ) {
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

/// A long record searched in pieces; the k-mer index of each strand is
/// built once and shared by all engines.
struct LongRec {
    /// unique for each long record
    serial: usize,
    id: Vec<u8>,
    body: Vec<u8>,
    ln: usize,
    index: [SharedIndex; 2],
}

/// Work for one worker.
enum Job {
    /// a whole record: id, body
    Whole(Vec<u8>, Vec<u8>),
    Piece(PieceJob),
}

/// Hits starting at `from..=last` of one strand of a long record.
struct PieceJob {
    rec: std::sync::Arc<LongRec>,
    rev: bool,
    from: usize,
    last: usize,
    /// the first piece of its strand
    first: bool,
    /// the last piece of the record (both strands)
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
}

/// What a worker or the reader sends to the writer, in input order.
enum Done {
    Hits(HitBuf),
    Piece(PieceJob, PieceHits),
    Error(Vec<u8>),
    End,
}

/// An engine with the strand of a long record it holds (loaded once for
/// all its pieces).
struct PieceEngine {
    eng: Engine,
    data: Buf,
    hits: Vec<i64>,
    /// serial number and strand of the record loaded
    loaded: Option<(usize, bool)>,
    /// the last hit `join_piece` found with this engine was not empty
    last_hit_nonempty: bool,
}

impl PieceEngine {
    fn new(eng: Engine) -> PieceEngine {
        PieceEngine {
            eng,
            data: Buf::new(),
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
        self.data.store(&rec.body);
        if rev {
            reverse_complement(&mut self.data, rec.ln);
        }
        self.eng
            .set_shared_index(Some(rec.index[rev as usize].clone()));
        self.eng.comp_data(&self.data);
        self.loaded = Some(key);
    }

    /// Before searching a whole record with this engine.
    fn unload(&mut self) {
        if self.loaded.is_some() {
            self.eng.set_shared_index(None);
            self.loaded = None;
        }
    }

    /// The hits of a piece: all of them (`all`, forward strand with -o),
    /// or the chain of non-overlapping hits from its first start.
    fn search(&mut self, job: &PieceJob, all: bool) -> PieceHits {
        let rec = &*job.rec;
        self.load(rec, job.rev);
        let mut pr = Printer {
            out: HitBuf::default(),
            line: Vec::new(),
        };
        let (mut starts, mut pasts) = (Vec::new(), Vec::new());
        let cb = self.eng.cdata_base();
        let mut i = self.eng.first_match_in(
            rec.ln as i32,
            job.from as i64,
            job.last as i64,
            &mut self.hits,
        );
        while i > 0 {
            let n = i as usize;
            print_hit(
                &mut pr, &self.eng, &rec.id, job.rev, rec.ln, &self.hits, n, &self.data,
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
        }
    }
}

/// Several threads: a reader thread, `threads` workers that each search
/// whole records, or pieces of long records, with their own engine, and
/// this thread, which prints the results in input order.  Hit limit
/// (`-m`) and miss limit (`-n`) are applied in input order, so the output
/// equals the one-thread output.
#[allow(clippy::too_many_arguments)]
fn run_threads(
    threads: usize,
    line: Vec<u8>,
    seq_type: i32,
    opts: Opts,
    split: bool,
    ignore: std::collections::HashSet<Vec<u8>>,
    mut max_hits: i32,
    mut stop_after: i32,
) {
    use std::sync::mpsc;
    use std::sync::{Arc, Mutex};
    let (job_tx, job_rx) = mpsc::sync_channel::<(usize, Job)>(2 * threads);
    let job_rx = Arc::new(Mutex::new(job_rx));
    let (done_tx, done_rx) = mpsc::channel::<(usize, Done)>();
    let split = split && !opts.protein;
    let piece = piece_size();

    let tx = done_tx.clone();
    std::thread::spawn(move || {
        let mut rd = FastaReader::new(open_fasta_input(std::io::stdin().lock()));
        let mut k = 0usize;
        loop {
            match rd.next() {
                ReadItem::Record { id, body } => {
                    if ignore.contains(&id) {
                        continue;
                    }
                    let ln = body.iter().position(|&b| b == 0).unwrap_or(body.len());
                    if !split || ln < 2 * piece {
                        if job_tx.send((k, Job::Whole(id, body))).is_err() {
                            return;
                        }
                        k += 1;
                        continue;
                    }
                    let rec = Arc::new(LongRec {
                        serial: k,
                        id,
                        body,
                        ln,
                        index: Default::default(),
                    });
                    let strands: &[bool] = if opts.complements {
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
                            let job = PieceJob {
                                rec: Arc::clone(&rec),
                                rev,
                                from,
                                last,
                                first: from == 0,
                                end_of_record: last == ln && rev == opts.complements,
                            };
                            if job_tx.send((k, Job::Piece(job))).is_err() {
                                return;
                            }
                            k += 1;
                            from = last + 1;
                        }
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

    let line = Arc::new(line);
    for _ in 0..threads {
        let rx = Arc::clone(&job_rx);
        let tx = done_tx.clone();
        let line = Arc::clone(&line);
        std::thread::Builder::new()
            .stack_size(1 << 30) // the parser recurses as deeply as the C code
            .spawn(move || {
                let mut eng = Engine::new();
                eng.parse_cmd(&line, seq_type);
                let mut pe = PieceEngine::new(eng);
                let mut data = Buf::new();
                let mut hits: Vec<i64> = vec![0; 2000];
                loop {
                    let job = rx.lock().unwrap().recv();
                    let Ok((k, job)) = job else { return };
                    let done = match job {
                        Job::Whole(id, body) => {
                            pe.unload();
                            let mut pr = Printer {
                                out: HitBuf::default(),
                                line: Vec::new(),
                            };
                            let mut unlimited = i32::MAX;
                            scan_record(
                                &mut pe.eng,
                                &mut data,
                                &mut hits,
                                &mut pr,
                                &id,
                                &body,
                                opts,
                                &mut unlimited,
                            );
                            Done::Hits(pr.out)
                        }
                        Job::Piece(job) => {
                            let ph = pe.search(&job, opts.show_overlaps && !job.rev);
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

    let mut out = sys::Out::new();
    let mut waiting = std::collections::BTreeMap::new();
    let mut next = 0usize;
    // a long record: end of the last hit printed on this strand, and
    // whether the record has a hit; an engine for searching again
    let mut past = 0i64;
    let mut rec_hit = false;
    let mut weng: Option<PieceEngine> = None;
    let mut wpr = Printer {
        out: HitBuf::default(),
        line: Vec::new(),
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
        match item {
            Done::End => break,
            Done::Error(msg) => read_error(&msg),
            Done::Hits(h) => {
                let printed = print_hits(&mut out, &h, 0, &mut max_hits);
                if printed == 0 {
                    missed(&mut stop_after);
                }
            }
            Done::Piece(job, ph) => {
                if job.first {
                    past = 0;
                    if !job.rev {
                        rec_hit = false;
                    }
                }
                if opts.show_overlaps && !job.rev {
                    // every hit, in search order
                    rec_hit |= print_hits(&mut out, &ph.out, 0, &mut max_hits) > 0;
                } else {
                    let we = weng.get_or_insert_with(|| {
                        let mut e = Engine::new();
                        e.parse_cmd(&line, seq_type);
                        PieceEngine::new(e)
                    });
                    rec_hit |=
                        join_piece(&mut out, &job, &ph, we, &mut wpr, &mut past, &mut max_hits);
                }
                if job.end_of_record && !rec_hit && max_hits > 0 {
                    missed(&mut stop_after);
                }
            }
        }
    }
}

/// Print the hits of `h` from hit `from` on while `max_hits` allows;
/// returns how many were printed.
fn print_hits(out: &mut sys::Out, h: &HitBuf, from: usize, max_hits: &mut i32) -> usize {
    let mut start = if from == 0 { 0 } else { h.ends[from - 1] };
    let mut printed = 0;
    for &end in &h.ends[from..] {
        if *max_hits <= 0 {
            break;
        }
        out.write(&h.buf[start..end]);
        start = end;
        *max_hits -= 1;
        printed += 1;
    }
    printed
}

/// Continue the one-thread search ("next hit that starts at or after
/// the end of the last one") through a piece.  `past` is the end of the
/// last hit printed on the strand.  The worker's hits are right from the
/// first one that the true search reaches; before that the piece is
/// searched again from `past` with `we`.  Returns whether a hit was
/// printed.
#[allow(clippy::too_many_arguments)]
fn join_piece(
    out: &mut sys::Out,
    job: &PieceJob,
    ph: &PieceHits,
    we: &mut PieceEngine,
    wpr: &mut Printer<HitBuf>,
    past: &mut i64,
    max_hits: &mut i32,
) -> bool {
    let last = job.last as i64;
    let mut printed = false;
    let mut again = false; // `we` holds a search of this piece
    loop {
        if *past > last || *max_hits <= 0 {
            return printed;
        }
        // the worker's first hit at or after `past`
        let j = ph.starts.partition_point(|&s| s < *past);
        // the worker searched on from the end of hit j-1, at or before
        // `past`: from here on its hits are the true ones.  (When the
        // last hit has length 0 the search goes on at its own start,
        // so only an ending after the start is used.)
        let synced = j == 0 || ph.pasts[j - 1] <= *past;
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
        // search the piece again from `past`
        let i = if !again {
            again = true;
            we.load(&job.rec, job.rev);
            we.eng
                .first_match_in(job.rec.ln as i32, *past, last, &mut we.hits)
        } else {
            we.eng.cont_match(&mut we.hits)
        };
        if i <= 0 {
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
            &we.data,
        );
        out.write(&wpr.out.buf);
        *max_hits -= 1;
        printed = true;
        let start = we.hits[0] - cb;
        *past = we.hits[n] - cb;
        we.last_hit_nonempty = *past > start;
    }
}

fn usage(errflag: i32, optind: i32, argc: i32) -> ! {
    eprintln!("errflag={} optind={} argc={}", errflag, optind, argc);
    eprintln!(
        "usage: scan_for_matches -c [for complementary strand] -p [protein][-n stop_after_n_misses] [-m max_hits] [-i file_of_ids_to_ignore] -o [overlapping hits] [-t threads] pattern_file < fasta_input > hits"
    );
    std::process::exit(2);
}

/// The FASTA input: plain text, or gzip / bgzip compressed (detected by
/// the first two bytes).
fn open_fasta_input<R: Read + 'static>(mut r: R) -> Box<dyn Read> {
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
        match gz::GzReader::new(r) {
            Ok(g) => Box::new(g),
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

fn real_main() {
    let mut args = sys::Args::from_env();
    let argc = args.argc();

    let mut stop_after: i32 = 100_000_000;
    let mut max_hits: i32 = 100_000_000;
    let mut show_overlaps = false;
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
            b'o' => show_overlaps = true,
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
    let pat_file = match open_path(&args.get(optind as usize)) {
        Some(f) => f,
        None => usage(errflag, optind, argc),
    };

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

    // read the pattern: '%' starts a comment, newlines become spaces
    let mut line: Vec<u8> = Vec::new();
    {
        let mut fp = Input::new(pat_file);
        let mut i = fp.getc();
        while i != EOF && line.len() < MAX_PAT_LINE_LN - 1 {
            if i == b'%' as i32 {
                loop {
                    i = fp.getc();
                    if i == b'\n' as i32 || i == EOF {
                        break;
                    }
                }
                if i != EOF {
                    i = b' ' as i32;
                }
            } else if i == b'\n' as i32 || i == b'\r' as i32 {
                i = b' ' as i32;
            } else {
                line.push(i as u8);
                i = fp.getc();
            }
        }
    }
    // the C string ends at the first NUL
    if let Some(k) = line.iter().position(|&b| b == 0) {
        line.truncate(k);
    }

    let mut eng = Engine::new();
    let rc = eng.parse_cmd(&line, if protein { PEPTIDE } else { DNA });
    if rc == 0 {
        let mut msg = b"failed to parse pattern: ".to_vec();
        msg.extend_from_slice(&line);
        msg.push(b'\n');
        use std::io::Write;
        let _ = std::io::stderr().write_all(&msg);
        std::process::exit(1);
    }

    let opts = Opts {
        protein,
        complements,
        show_overlaps,
    };
    // patterns that use matches of earlier sequences stay on one thread
    if threads > 1 && !eng.uses_earlier_state() && max_hits > 0 {
        let split = eng.can_split();
        run_threads(
            threads,
            line,
            if protein { PEPTIDE } else { DNA },
            opts,
            split,
            ignore,
            max_hits,
            stop_after,
        );
    } else {
        run_sequential(eng, opts, &ignore, max_hits, stop_after);
    }
}
