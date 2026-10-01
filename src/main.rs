//! scan_for_matches -- Rust port of the program by Ross Overbeek (ANL).
//!
//!   scan_for_matches [-c] [-p] [-o N] [-i ignore_file] pattern_file < fasta_input > hits
//!
//! This file is the port of scan_for_matches.c (option handling, FASTA
//! reading and hit printing).  The pattern language lives in engine.rs.

mod engine;
mod sys;

use engine::{compl, Buf, Engine, CDATA_BASE, DNA, PEPTIDE};
use std::io::Read;

const MAX_PAT_LINE_LN: usize = 32000;
const EOF: i32 = -1;

#[inline(always)]
fn isspace(c: i32) -> bool {
    matches!(c, 0x20 | 0x09 | 0x0a | 0x0b | 0x0c | 0x0d)
}

/// Byte reader with C `getc`/`ungetc` semantics.
struct Input<R: Read> {
    r: R,
    buf: Vec<u8>,
    pos: usize,
    len: usize,
    back: Option<i32>,
}

impl<R: Read> Input<R> {
    fn new(r: R) -> Input<R> {
        Input { r, buf: vec![0u8; 1 << 16], pos: 0, len: 0, back: None }
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

    /// Sequence body: everything up to EOF or '>', dropping ' ' and '\n'.
    /// Returns the char that stopped the scan.
    fn read_body(&mut self, out: &mut Vec<u8>) -> i32 {
        out.clear();
        if let Some(c) = self.back.take() {
            if c == b'>' as i32 {
                return c;
            }
            if c != b' ' as i32 && c != b'\n' as i32 {
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
                    out.extend(chunk[..k].iter().copied().filter(|&b| b != b' ' && b != b'\n'));
                    self.pos += k + 1;
                    return b'>' as i32;
                }
                None => {
                    out.extend(chunk.iter().copied().filter(|&b| b != b' ' && b != b'\n'));
                    self.pos = self.len;
                }
            }
        }
    }
}

struct Printer {
    out: sys::Out,
    line: Vec<u8>,
}

impl Printer {
    /// Print one hit the way the C code does:
    ///   >id:[a,b]
    ///   seg1 seg2 ... segN \n
    fn hit(&mut self, id: &[u8], a: i64, b: i64, hits: &[i64], n: usize, data: &Buf) {
        self.line.clear();
        self.line.push(b'>');
        self.line.extend_from_slice(id);
        self.line.extend_from_slice(format!(":[{},{}]\n", a, b).as_bytes());
        self.out.write(&self.line);
        self.line.clear();
        for i1 in 0..n {
            let j = (hits[i1 + 1] - hits[i1]) as i32;
            // `for (...; j; j--)` with an int counter
            let count = j as u32 as u64;
            let mut p = hits[i1] - CDATA_BASE;
            for _ in 0..count {
                match data.try_get(p) {
                    Some(b) => self.line.push(b),
                    None => {
                        // C faults here after printf() has buffered every
                        // character before this one
                        self.out.write(&self.line);
                        sys::segv();
                    }
                }
                p += 1;
                if self.line.len() >= 1 << 16 {
                    self.out.write(&self.line);
                    self.line.clear();
                }
            }
            self.line.push(b' ');
        }
        self.line.push(b'\n');
        self.out.write(&self.line);
    }
}

fn usage(errflag: i32, optind: i32, argc: i32) -> ! {
    eprintln!("errflag={} optind={} argc={}", errflag, optind, argc);
    eprintln!(
        "usage: scan_for_matches -c [for complementary strand] -p [protein][-n stop_after_n_misses] [-m max_hits] [-i file_of_ids_to_ignore] -o [overlapping hits] pattern_file < fasta_input > hits"
    );
    std::process::exit(2);
}

fn open_path(p: &[u8]) -> Option<std::fs::File> {
    use std::os::unix::ffi::OsStrExt;
    std::fs::File::open(std::ffi::OsStr::from_bytes(p)).ok()
}

fn main() {
    // The parser and max_mat() recurse as deeply as the C code does; give
    // them a large stack so only a truly endless recursion (which crashes
    // the C program too) ends in SIGSEGV.
    let t = std::thread::Builder::new().stack_size(1 << 30).spawn(real_main).expect("thread");
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

    loop {
        let c = args.getopt(b"pcnmo:i:\0");
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
    let pat_file = match open_path(&args.get(optind as usize)) {
        Some(f) => f,
        None => usage(errflag, optind, argc),
    };

    // ids to ignore
    let mut ignore: Vec<Vec<u8>> = Vec::new();
    if let Some(f) = ig_fp {
        let mut inp = Input::new(f);
        let mut id = Vec::new();
        while inp.scan_s(&mut id) == 1 {
            ignore.push(id.clone());
        }
        if !ignore.is_empty() {
            eprintln!("ignoring {} id(s)", ignore.len());
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
            } else if i == b'\n' as i32 {
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

    let mut data = Buf::new();
    // `char *hits[2000]` is not initialised in C.  When the first search
    // finds nothing, `past_last = hits[0]` copies the leftover value, which
    // is this constant in the reference build.
    let mut hits: Vec<i64> = vec![0; 2000];
    hits[0] = 0x0f00_7fff_ffff_fff8;
    let mut pr = Printer { out: sys::Out::new(), line: Vec::new() };
    let mut inp = Input::new(std::io::stdin().lock());
    let mut id: Vec<u8> = Vec::new();
    let mut body: Vec<u8> = Vec::new();
    let mut got_gt = false;

    loop {
        if max_hits <= 0 {
            break;
        }
        let r = if !got_gt { inp.scan_gt_s(&mut id) } else { inp.scan_s(&mut id) };
        if r != 1 {
            break;
        }
        // skip the rest of the header line
        loop {
            let c = inp.getc();
            if c == b'\n' as i32 {
                break;
            }
            if c == EOF {
                sys::hang();
            }
        }
        let stop = inp.read_body(&mut body);
        got_gt = stop == b'>' as i32;

        // copy into the persistent data buffer, NUL terminated; past the
        // malloc'd block the C program faults while reading the input
        if body.len() as i64 >= engine::ALLOC_LEN {
            sys::segv();
        }
        data.v[..body.len()].copy_from_slice(&body);
        data.v[body.len()] = 0;
        let ln = body.iter().position(|&b| b == 0).unwrap_or(body.len());

        if ignore.iter().any(|x| *x == id) {
            continue;
        }

        if !protein {
            eng.comp_data(&data);
        } else {
            eng.copy_data(&data);
        }

        let mut hit_in_line = false;
        let mut i = eng.first_match(ln as i32, &mut hits);
        while max_hits > 0 && i > 0 {
            hit_in_line = true;
            max_hits -= 1;
            let n = i as usize;
            pr.hit(&id, 1 + hits[0] - CDATA_BASE, 1 + (hits[n] - 1 - CDATA_BASE), &hits, n, &data);
            i = if !show_overlaps { eng.cont_match(&mut hits) } else { eng.next_match(&mut hits) };
        }

        if complements {
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
            eng.comp_data(&data);

            let mut i = eng.first_match(ln as i32, &mut hits);
            while max_hits > 0 && i > 0 {
                hit_in_line = true;
                max_hits -= 1;
                let n = i as usize;
                let l = ln as i64;
                pr.hit(&id, 1 + (l - 1) - (hits[0] - CDATA_BASE), 1 + (l - 1) - (hits[n] - 1 - CDATA_BASE), &hits, n, &data);
                i = eng.cont_match(&mut hits);
            }
        }

        if !hit_in_line {
            stop_after = stop_after.wrapping_sub(1);
            if stop_after == 0 {
                eprintln!("exceeded limit of lines failing to match");
                std::process::exit(1);
            }
        }
    }
}
