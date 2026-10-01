//! The test suite shipped with the C program (`run_tests` / `testit`):
//! every pattern line is run on its own and the combined output must equal
//! `test_output` byte for byte.

use std::io::Write;
use std::process::{Command, Stdio};

fn run(args: &[&str], pattern: &[u8], input: &[u8]) -> Vec<u8> {
    let dir = std::env::temp_dir().join(format!("sfm_suite_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let pat = dir.join("tmp.pat");
    std::fs::write(&pat, pattern).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"))
        .args(args)
        .arg(&pat)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap().stdout
}

#[test]
fn original_test_suite() {
    let data = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/");
    let read = |f: &str| std::fs::read(format!("{data}{f}")).unwrap();
    let mut out = Vec::new();
    let dna_in = read("test_dna_input");
    for line in read("test_dna_patterns").split_inclusive(|&b| b == b'\n') {
        out.extend_from_slice(b"=====================\n");
        out.extend_from_slice(line);
        out.extend_from_slice(b"\n\n");
        out.extend(run(&[], line, &dna_in));
    }
    out.extend_from_slice(b"==========<<<< Protein Matches >>>>==========\n");
    let prot_in = read("test_prot_input");
    for line in read("test_prot_patterns").split_inclusive(|&b| b == b'\n') {
        out.extend_from_slice(b"=====================\n");
        out.extend_from_slice(line);
        out.extend_from_slice(b"\n\n");
        out.extend(run(&["-p"], line, &prot_in));
    }
    let expected = read("test_output");
    assert!(out == expected, "output differs from test_output");
}
