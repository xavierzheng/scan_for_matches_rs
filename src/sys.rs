//! Thin bindings to the C library that the original program relied on.
//!
//! The original uses the platform `getopt`, `sscanf` and stdio for its
//! command line and its output.  We call the same libc functions so that
//! option handling, diagnostics and stdout buffering behave identically
//! (libc is always linked by std, so no extra crate is needed).

use std::ffi::CString;
use std::os::raw::{c_char, c_int, c_void};

unsafe extern "C" {
    fn getopt(argc: c_int, argv: *const *mut c_char, optstring: *const c_char) -> c_int;
    static mut optarg: *mut c_char;
    static mut optind: c_int;
    fn sscanf(s: *const c_char, fmt: *const c_char, ...) -> c_int;
    fn fdopen(fd: c_int, mode: *const c_char) -> *mut c_void;
    fn fwrite(p: *const c_void, size: usize, n: usize, f: *mut c_void) -> usize;
    fn signal(sig: c_int, handler: usize) -> usize;
    fn raise(sig: c_int) -> c_int;
    fn dup2(old: c_int, new: c_int) -> c_int;
    fn mmap(
        addr: *mut c_void,
        len: usize,
        prot: c_int,
        flags: c_int,
        fd: c_int,
        off: i64,
    ) -> *mut c_void;
    fn madvise(addr: *mut c_void, len: usize, advice: c_int) -> c_int;
    fn mprotect(addr: *mut c_void, len: usize, prot: c_int) -> c_int;
    fn munmap(addr: *mut c_void, len: usize) -> c_int;
}

const SIGSEGV: c_int = 11;
#[cfg(target_os = "macos")]
const SIGBUS: c_int = 10;
#[cfg(not(target_os = "macos"))]
const SIGBUS: c_int = 7;
const SIG_DFL: usize = 0;

/// Terminate the way the C program does when it dereferences a bad pointer.
pub fn segv() -> ! {
    unsafe {
        signal(SIGSEGV, SIG_DFL);
        raise(SIGSEGV);
    }
    std::process::abort()
}

/// A protection fault (write into a read-only segment).
pub fn sigbus() -> ! {
    unsafe {
        signal(SIGBUS, SIG_DFL);
        raise(SIGBUS);
    }
    std::process::abort()
}

/// argv kept alive for the whole run (getopt may permute it on glibc).
pub struct Args {
    _owned: Vec<CString>,
    ptrs: Vec<*mut c_char>,
}

impl Args {
    pub fn from_env() -> Args {
        use std::os::unix::ffi::OsStringExt;
        let owned: Vec<CString> = std::env::args_os()
            .map(|a| CString::new(a.into_vec()).unwrap_or_default())
            .collect();
        let mut ptrs: Vec<*mut c_char> = owned.iter().map(|c| c.as_ptr() as *mut c_char).collect();
        ptrs.push(std::ptr::null_mut());
        Args {
            _owned: owned,
            ptrs,
        }
    }

    pub fn argc(&self) -> c_int {
        (self.ptrs.len() - 1) as c_int
    }

    /// Bytes of argv[i] (after any permutation done by getopt).
    pub fn get(&self, i: usize) -> Vec<u8> {
        unsafe { std::ffi::CStr::from_ptr(self.ptrs[i]).to_bytes().to_vec() }
    }

    /// Take the long options (`--name value`, `--name=value`, `--flag`)
    /// out of argv, before `getopt` sees them; `--` ends the options.
    /// `with_value`: the names that take a value.
    #[allow(clippy::type_complexity)]
    pub fn take_long(
        &mut self,
        with_value: &[&str],
        flags: &[&str],
    ) -> Result<Vec<(String, Vec<u8>)>, String> {
        let mut out = Vec::new();
        let mut i = 1;
        while i < self.ptrs.len() - 1 {
            let a = self.get(i);
            if a == b"--" {
                break;
            }
            let Some(rest) = a.strip_prefix(b"--") else {
                i += 1;
                continue;
            };
            let (name, val) = match rest.iter().position(|&c| c == b'=') {
                Some(k) => (&rest[..k], Some(rest[k + 1..].to_vec())),
                None => (rest, None),
            };
            let name = String::from_utf8_lossy(name).into_owned();
            let mut n = 1;
            let val = if with_value.contains(&name.as_str()) {
                match val {
                    Some(v) => v,
                    None if i + 1 < self.ptrs.len() - 1 => {
                        n = 2;
                        self.get(i + 1)
                    }
                    None => return Err(format!("option --{name} needs a value")),
                }
            } else if flags.contains(&name.as_str()) && val.is_none() {
                Vec::new()
            } else {
                return Err(format!("unknown option --{name}"));
            };
            out.push((name, val));
            self.ptrs.drain(i..i + n);
        }
        Ok(out)
    }

    pub fn getopt(&mut self, optstring: &[u8]) -> c_int {
        unsafe {
            getopt(
                self.argc(),
                self.ptrs.as_ptr(),
                optstring.as_ptr() as *const c_char,
            )
        }
    }
}

pub fn optind_get() -> c_int {
    unsafe { optind }
}

pub fn optarg_bytes() -> Option<Vec<u8>> {
    unsafe {
        let p = optarg;
        if p.is_null() {
            None
        } else {
            Some(std::ffi::CStr::from_ptr(p).to_bytes().to_vec())
        }
    }
}

/// `sscanf(optarg,"%d",&v)` exactly as the C code does it (including the
/// crash when optarg is NULL).
pub fn sscanf_optarg_int(v: &mut i32) -> c_int {
    unsafe {
        let mut x: c_int = *v;
        let r = sscanf(optarg, c"%d".as_ptr(), &mut x as *mut c_int);
        *v = x;
        r
    }
}

/// Make `f` the process's stdout (file descriptor 1).  Call before the
/// first `Out::new`.
pub fn redirect_stdout(f: std::fs::File) -> bool {
    use std::os::unix::io::AsRawFd;
    let ok = unsafe { dup2(f.as_raw_fd(), 1) } == 1;
    drop(f); // fd 1 keeps the file open
    ok
}

/// stdout as a C stdio stream, so buffering matches the original.
pub struct Out {
    fp: *mut c_void,
}

impl Out {
    pub fn new() -> Out {
        let fp = unsafe { fdopen(1, c"w".as_ptr()) };
        Out { fp }
    }

    pub fn write(&mut self, b: &[u8]) {
        if !self.fp.is_null() && !b.is_empty() {
            unsafe {
                fwrite(b.as_ptr() as *const c_void, 1, b.len(), self.fp);
            }
        }
    }
}

/// Give the pages of `len` bytes at `p` (page aligned) back to the
/// system; they read as zero bytes when touched again.
pub fn release_pages(p: *mut u8, len: usize) {
    const MADV_DONTNEED: c_int = 4;
    if len > 0 {
        unsafe { madvise(p as *mut c_void, len, MADV_DONTNEED) };
    }
}

/// Make `len` bytes at `p` (page aligned) read only: a write is a
/// SIGSEGV.
pub fn protect_read_only(p: *mut u8, len: usize) {
    const PROT_READ: c_int = 0x1;
    unsafe { mprotect(p as *mut c_void, len, PROT_READ) };
}

/// Unmap memory from `map_zeroed`.
pub fn unmap(p: *mut u8, len: usize) {
    unsafe { munmap(p as *mut c_void, len) };
}

/// Zeroed read/write memory (address space only; pages are used when
/// touched), placed at `hint` if that range is free.  Never unmapped.
pub fn map_zeroed(hint: usize, len: usize) -> *mut u8 {
    #[cfg(target_os = "macos")]
    const MAP_ANON: c_int = 0x1000;
    #[cfg(not(target_os = "macos"))]
    const MAP_ANON: c_int = 0x20;
    #[cfg(target_os = "macos")]
    const MAP_NORESERVE: c_int = 0;
    #[cfg(not(target_os = "macos"))]
    const MAP_NORESERVE: c_int = 0x4000;
    const MAP_PRIVATE: c_int = 0x2;
    const PROT_RW: c_int = 0x1 | 0x2;
    let p = unsafe {
        mmap(
            hint as *mut c_void,
            len,
            PROT_RW,
            MAP_PRIVATE | MAP_ANON | MAP_NORESERVE,
            -1,
            0,
        )
    };
    if p as isize == -1 {
        std::process::abort();
    }
    p as *mut u8
}
