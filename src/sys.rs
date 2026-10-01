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
}

const SIGSEGV: c_int = 11;
const SIG_DFL: usize = 0;

/// Terminate the way the C program does when it dereferences a bad pointer.
pub fn segv() -> ! {
    unsafe {
        signal(SIGSEGV, SIG_DFL);
        raise(SIGSEGV);
    }
    std::process::abort()
}

/// The C program spins forever in `while (getc(stdin) != '\n');` at EOF.
/// Buffered output is never flushed, exactly as in the original.
pub fn hang() -> ! {
    loop {
        std::thread::sleep(std::time::Duration::from_secs(3600));
    }
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
        Args { _owned: owned, ptrs }
    }

    pub fn argc(&self) -> c_int {
        (self.ptrs.len() - 1) as c_int
    }

    /// Bytes of argv[i] (after any permutation done by getopt).
    pub fn get(&self, i: usize) -> Vec<u8> {
        unsafe { std::ffi::CStr::from_ptr(self.ptrs[i]).to_bytes().to_vec() }
    }

    pub fn getopt(&mut self, optstring: &[u8]) -> c_int {
        unsafe { getopt(self.argc(), self.ptrs.as_ptr(), optstring.as_ptr() as *const c_char) }
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
