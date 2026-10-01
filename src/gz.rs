//! gzip / bgzip input.  bgzip files (as made by htslib) are a series of
//! gzip members; all of them are read, one after the other.  Decompression
//! uses the system zlib library.

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
}

const Z_OK: c_int = 0;
const Z_STREAM_END: c_int = 1;
const Z_BUF_ERROR: c_int = -5;
const Z_NO_FLUSH: c_int = 0;
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

impl<R: Read> GzReader<R> {
    pub fn new(inner: R) -> io::Result<GzReader<R>> {
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
                GZIP_WINDOW,
                zlibVersion(),
                std::mem::size_of::<ZStream>() as c_int,
            )
        };
        if rc != Z_OK {
            return Err(bad("cannot start zlib"));
        }
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

    fn fill(&mut self) -> io::Result<()> {
        if self.in_pos < self.in_len || self.in_eof {
            return Ok(());
        }
        loop {
            match self.inner.read(&mut self.inbuf) {
                Ok(0) => {
                    self.in_eof = true;
                    return Ok(());
                }
                Ok(n) => {
                    self.in_pos = 0;
                    self.in_len = n;
                    return Ok(());
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                Err(e) => return Err(e),
            }
        }
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
