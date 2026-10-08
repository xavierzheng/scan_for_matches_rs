//! gzip / bgzip input.  bgzip files (as made by htslib) are a series of
//! gzip members; all of them are read, one after the other.  Decompression
//! uses the system zlib library.  `BgzfReader` decompresses the blocks of
//! a bgzip file in several threads; a member that is not a bgzip block
//! (plain gzip) and everything after it is read by `GzReader`.

use std::io::{self, Read};
use std::os::raw::{c_char, c_int, c_uint, c_ulong, c_void};

#[repr(C)]
struct ZStream {
    next_in: *const u8,
    avail_in: c_uint,
    total_in: c_ulong,
    next_out: *mut u8,
    avail_out: c_uint,
    total_out: c_ulong,
    msg: *const c_char,
    state: *mut c_void,
    zalloc: *mut c_void,
    zfree: *mut c_void,
    opaque: *mut c_void,
    data_type: c_int,
    adler: c_ulong,
    reserved: c_ulong,
}

#[link(name = "z")]
unsafe extern "C" {
    fn zlibVersion() -> *const c_char;
    fn inflateInit2_(
        strm: *mut ZStream,
        window_bits: c_int,
        version: *const c_char,
        stream_size: c_int,
    ) -> c_int;
    fn inflate(strm: *mut ZStream, flush: c_int) -> c_int;
    fn inflateReset(strm: *mut ZStream) -> c_int;
    fn inflateEnd(strm: *mut ZStream) -> c_int;
    fn crc32(crc: c_ulong, buf: *const u8, len: c_uint) -> c_ulong;
}

const Z_OK: c_int = 0;
const Z_STREAM_END: c_int = 1;
const Z_BUF_ERROR: c_int = -5;
const Z_NO_FLUSH: c_int = 0;
const Z_FINISH: c_int = 4;
/// raw deflate data (the payload of a bgzip block)
const RAW_WINDOW: c_int = -15;
/// 15-bit window, gzip header expected
const GZIP_WINDOW: c_int = 15 + 16;

/// Does the input start with the gzip magic bytes?
pub fn is_gzip(head: &[u8]) -> bool {
    head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
}

pub struct GzReader<R: Read> {
    inner: R,
    inbuf: Vec<u8>,
    in_pos: usize,
    in_len: usize,
    in_eof: bool,
    /// zlib keeps a pointer to the stream, so it must not move
    strm: Box<ZStream>,
    /// inside a member (false after a member's end)
    in_member: bool,
}

fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("gzip input: {msg}"))
}

/// A zlib inflate stream (`window`: GZIP_WINDOW or RAW_WINDOW).
fn new_stream(window: c_int) -> io::Result<Box<ZStream>> {
    let mut strm = Box::new(ZStream {
        next_in: std::ptr::null(),
        avail_in: 0,
        total_in: 0,
        next_out: std::ptr::null_mut(),
        avail_out: 0,
        total_out: 0,
        msg: std::ptr::null(),
        state: std::ptr::null_mut(),
        zalloc: std::ptr::null_mut(),
        zfree: std::ptr::null_mut(),
        opaque: std::ptr::null_mut(),
        data_type: 0,
        adler: 0,
        reserved: 0,
    });
    let rc = unsafe {
        inflateInit2_(
            &mut *strm,
            window,
            zlibVersion(),
            std::mem::size_of::<ZStream>() as c_int,
        )
    };
    if rc != Z_OK {
        return Err(bad("cannot start zlib"));
    }
    Ok(strm)
}

impl<R: Read> GzReader<R> {
    pub fn new(inner: R) -> io::Result<GzReader<R>> {
        let strm = new_stream(GZIP_WINDOW)?;
        Ok(GzReader {
            inner,
            inbuf: vec![0; 1 << 17],
            in_pos: 0,
            in_len: 0,
            in_eof: false,
            strm,
            in_member: false,
        })
    }

    /// Refill the input buffer completely (or up to the end of the input):
    /// zlib then always gets the same pieces, whatever sizes the reads of
    /// a pipe return, so damaged input gives the same output every time.
    fn fill(&mut self) -> io::Result<()> {
        if self.in_pos < self.in_len || self.in_eof {
            return Ok(());
        }
        self.in_pos = 0;
        self.in_len = 0;
        while self.in_len < self.inbuf.len() {
            match self.inner.read(&mut self.inbuf[self.in_len..]) {
                Ok(0) => {
                    self.in_eof = true;
                    break;
                }
                Ok(n) => self.in_len += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
}

impl<R: Read> Read for GzReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            self.fill()?;
            let avail = self.in_len - self.in_pos;
            if avail == 0 && self.in_eof {
                if self.in_member {
                    return Err(bad("unexpected end of compressed data"));
                }
                return Ok(0);
            }
            let s = &mut *self.strm;
            s.next_in = self.inbuf[self.in_pos..].as_ptr();
            s.avail_in = avail as c_uint;
            s.next_out = out.as_mut_ptr();
            s.avail_out = out.len().min(c_uint::MAX as usize) as c_uint;
            let want = s.avail_out as usize;
            let rc = unsafe { inflate(s, Z_NO_FLUSH) };
            let used = avail - s.avail_in as usize;
            let made = want - s.avail_out as usize;
            self.in_pos += used;
            if used > 0 || made > 0 {
                self.in_member = true;
            }
            match rc {
                Z_OK => {}
                Z_STREAM_END => {
                    // next member (bgzip block), if any
                    unsafe { inflateReset(s) };
                    self.in_member = false;
                }
                Z_BUF_ERROR if made == 0 && used == 0 => {
                    // needs more input; loop reads it (or reports the end)
                    if self.in_eof {
                        return Err(bad("unexpected end of compressed data"));
                    }
                }
                _ => return Err(bad("damaged compressed data")),
            }
            if made > 0 {
                return Ok(made);
            }
        }
    }
}

impl<R: Read> Drop for GzReader<R> {
    fn drop(&mut self) {
        unsafe { inflateEnd(&mut *self.strm) };
    }
}

/// The size of the bgzip block at the start of `h`: `Some(Some(n))` for a
/// block of `n` bytes, `Some(None)` when `h` does not start with a bgzip
/// block header, `None` when more bytes are needed to tell.
fn bgzf_block_size(h: &[u8]) -> Option<Option<usize>> {
    const MAGIC: [u8; 4] = [0x1f, 0x8b, 8, 4]; // gzip, deflate, FEXTRA only
    let n = h.len().min(4);
    if h[..n] != MAGIC[..n] {
        return Some(None);
    }
    if h.len() < 12 {
        return None;
    }
    let xlen = u16::from_le_bytes([h[10], h[11]]) as usize;
    if h.len() < 12 + xlen {
        return None;
    }
    let mut x = &h[12..12 + xlen];
    while x.len() >= 4 {
        let slen = u16::from_le_bytes([x[2], x[3]]) as usize;
        if x[0] == b'B' && x[1] == b'C' && slen == 2 && x.len() >= 6 {
            let size = u16::from_le_bytes([x[4], x[5]]) as usize + 1;
            // header, payload, CRC32 and ISIZE
            return Some((size >= 12 + xlen + 8).then_some(size));
        }
        if x.len() < 4 + slen {
            break;
        }
        x = &x[4 + slen..];
    }
    Some(None)
}

/// Decompress one bgzip block (`block`: the whole block) and append its
/// data to `out`.
fn inflate_block(strm: &mut ZStream, block: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    let xlen = u16::from_le_bytes([block[10], block[11]]) as usize;
    let n = block.len();
    let payload = &block[12 + xlen..n - 8];
    let crc = u32::from_le_bytes(block[n - 8..n - 4].try_into().unwrap());
    let isize = u32::from_le_bytes(block[n - 4..].try_into().unwrap()) as usize;
    // deflate cannot expand data more than about 1032 times
    if isize > 1100 * payload.len() + 64 {
        return Err(bad("damaged compressed data"));
    }
    let start = out.len();
    out.resize(start + isize + 1, 0);
    unsafe { inflateReset(strm) };
    strm.next_in = payload.as_ptr();
    strm.avail_in = payload.len() as c_uint;
    strm.next_out = out[start..].as_mut_ptr();
    strm.avail_out = (isize + 1) as c_uint;
    let rc = unsafe { inflate(strm, Z_FINISH) };
    let made = isize + 1 - strm.avail_out as usize;
    if rc != Z_STREAM_END || strm.avail_in != 0 || made != isize {
        return Err(bad("damaged compressed data"));
    }
    out.truncate(start + isize);
    let got = unsafe { crc32(0, out[start..].as_ptr(), isize as c_uint) } as u32;
    if got != crc {
        return Err(bad("damaged compressed data"));
    }
    Ok(())
}

/// Some bgzip blocks: their bytes, one after the other, and their sizes.
struct Batch {
    bytes: Vec<u8>,
    sizes: Vec<usize>,
}

/// The data of a batch; with an error, the data of the blocks before the
/// bad one (so the output does not depend on how blocks were batched).
type Done = (Vec<u8>, Option<io::Error>);
type Job = (Batch, mpsc::Sender<Done>);

use std::collections::VecDeque;
use std::sync::{Arc, Mutex, mpsc};

/// Blocks per batch (about 1 MB of data).
const BATCH_BLOCKS: usize = 16;

/// gzip / bgzip input decompressed by `threads` worker threads.  The
/// reader cuts the input into batches of whole bgzip blocks; the workers
/// decompress them; `read` gives the data back in input order.  Output
/// is the same as `GzReader` for any good input.
pub struct BgzfReader<R: Read> {
    /// the compressed input (None once handed to `serial`)
    inner: Option<R>,
    /// input bytes read but not yet put in a batch
    pending: Vec<u8>,
    in_eof: bool,
    /// from the first member that is not a bgzip block to the end
    serial: Option<GzReader<io::Chain<io::Cursor<Vec<u8>>, R>>>,
    jobs: Option<mpsc::Sender<Job>>,
    /// results of the batches given to the workers, in input order
    queue: VecDeque<mpsc::Receiver<Done>>,
    max_queue: usize,
    cur: Vec<u8>,
    pos: usize,
    /// returned once `cur` is used up
    err: Option<io::Error>,
}

impl<R: Read> BgzfReader<R> {
    pub fn new(inner: R, threads: usize) -> io::Result<BgzfReader<R>> {
        let threads = threads.max(1);
        let (tx, rx) = mpsc::channel::<Job>();
        let rx = Arc::new(Mutex::new(rx));
        for _ in 0..threads {
            let rx = Arc::clone(&rx);
            std::thread::Builder::new()
                .name("bgzf".into())
                .spawn(move || {
                    // zlib state is not Send: each worker makes its own
                    let Ok(mut strm) = new_stream(RAW_WINDOW) else {
                        return;
                    };
                    loop {
                        let job = rx.lock().unwrap().recv();
                        let Ok((batch, done)) = job else { break };
                        let mut out = Vec::with_capacity(batch.sizes.len() << 16);
                        let mut at = 0;
                        let mut err = None;
                        for &n in &batch.sizes {
                            let ok = out.len();
                            if let Err(e) =
                                inflate_block(&mut strm, &batch.bytes[at..at + n], &mut out)
                            {
                                out.truncate(ok);
                                err = Some(e);
                                break;
                            }
                            at += n;
                        }
                        let _ = done.send((out, err));
                    }
                    unsafe { inflateEnd(&mut *strm) };
                })
                .map_err(|_| bad("cannot start threads"))?;
        }
        Ok(BgzfReader {
            inner: Some(inner),
            pending: Vec::new(),
            in_eof: false,
            serial: None,
            jobs: Some(tx),
            queue: VecDeque::new(),
            max_queue: 2 * threads,
            cur: Vec::new(),
            pos: 0,
            err: None,
        })
    }

    /// Read more input into `pending`; false at the end of the input.
    fn more(&mut self) -> io::Result<bool> {
        let Some(r) = self.inner.as_mut() else {
            return Ok(false);
        };
        if self.in_eof {
            return Ok(false);
        }
        let start = self.pending.len();
        self.pending.resize(start + (1 << 20), 0);
        loop {
            match r.read(&mut self.pending[start..]) {
                Ok(0) => {
                    self.pending.truncate(start);
                    self.in_eof = true;
                    return Ok(false);
                }
                Ok(n) => {
                    self.pending.truncate(start + n);
                    return Ok(true);
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => {
                    self.pending.truncate(start);
                    return Err(e);
                }
            }
        }
    }

    /// Give the next batch of blocks to the workers.  At the first input
    /// that is not a whole bgzip block, the rest goes to `serial`.
    /// False when there is nothing more to give.
    fn send_batch(&mut self) -> io::Result<bool> {
        if self.inner.is_none() {
            return Ok(false);
        }
        let mut sizes = Vec::new();
        let mut used = 0;
        let mut not_bgzf = false;
        let mut at_end = false;
        while sizes.len() < BATCH_BLOCKS {
            let h = &self.pending[used..];
            match bgzf_block_size(h) {
                Some(Some(n)) if h.len() >= n => {
                    sizes.push(n);
                    used += n;
                }
                Some(None) => {
                    not_bgzf = true;
                    break;
                }
                // more bytes needed
                _ => {
                    if !sizes.is_empty() {
                        break; // send the blocks that are complete
                    }
                    if !self.more()? {
                        at_end = true;
                        break;
                    }
                }
            }
        }
        if !sizes.is_empty() {
            let bytes = self.pending[..used].to_vec();
            self.pending.drain(..used);
            let (tx, rx) = mpsc::channel();
            let jobs = self.jobs.as_ref().unwrap();
            if jobs.send((Batch { bytes, sizes }, tx)).is_err() {
                return Err(bad("threads stopped"));
            }
            self.queue.push_back(rx);
            return Ok(true);
        }
        // a member that is not a bgzip block, or the input ends inside a
        // block: GzReader reads (or reports) the rest
        if not_bgzf || (at_end && !self.pending.is_empty()) {
            let rest = std::mem::take(&mut self.pending);
            let inner = self.inner.take().unwrap();
            self.serial = Some(GzReader::new(io::Cursor::new(rest).chain(inner))?);
        } else {
            self.inner = None; // end of input
        }
        Ok(false)
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            if self.pos < self.cur.len() {
                let n = out.len().min(self.cur.len() - self.pos);
                out[..n].copy_from_slice(&self.cur[self.pos..self.pos + n]);
                self.pos += n;
                return Ok(n);
            }
            if let Some(e) = self.err.take() {
                // stop here: nothing after a damaged block is read
                self.inner = None;
                self.serial = None;
                self.queue.clear();
                return Err(e);
            }
            while self.queue.len() < self.max_queue && self.send_batch()? {}
            match self.queue.pop_front() {
                Some(rx) => {
                    let (data, err) = rx.recv().map_err(|_| bad("threads stopped"))?;
                    self.cur = data;
                    self.pos = 0;
                    self.err = err;
                }
                None => {
                    return match self.serial.as_mut() {
                        Some(g) => g.read(out),
                        None => Ok(0),
                    };
                }
            }
        }
    }
}

impl<R: Read> Drop for BgzfReader<R> {
    fn drop(&mut self) {
        // the workers stop when the job channel closes
        self.jobs = None;
    }
}
