//! `--lint` and `--explain` on the command line: no input is read; exit
//! status 1 when there is an error.

use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_scan_for_matches");
const DTC: &str = "tir_scan_patterns/labelled/00_DTC_published_2025.pat";

fn run(args: &[&str]) -> (Option<i32>, String, String) {
    let out = Command::new(BIN)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(args)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    (
        out.status.code(),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(out.stderr).unwrap(),
    )
}

fn pattern(name: &str, text: &str) -> String {
    let d = std::env::temp_dir().join(format!("sfm_lint_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    let p = d.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_str().unwrap().to_string()
}

#[test]
fn lint_tir_pattern() {
    let (st, out, err) = run(&["--lint", DTC]);
    assert_eq!(st, Some(0), "{err}");
    assert!(out.contains("note: the pattern reads the same on both strands"));
    assert!(out.contains("add --strict-n"), "{out}");
    assert!(out.ends_with("0 errors, 0 warnings, 2 notes\n"), "{out}");
    // -c without --dedup: a warning; with --dedup (no --format needed): none
    let (st, out, _) = run(&["--lint", "-c", DTC]);
    assert_eq!(st, Some(0));
    assert!(out.ends_with("0 errors, 1 warnings, 1 notes\n"), "{out}");
    let (st, out, err) = run(&["--lint", "-c", "--dedup", "--strict-n", DTC]);
    assert_eq!(st, Some(0), "{err}");
    assert!(out.ends_with("0 errors, 0 warnings, 0 notes\n"), "{out}");
}

#[test]
fn lint_errors_and_several_files() {
    let bad = pattern("bad.pat", "p1=3...3 ~p9\n");
    let rev = pattern("rev.pat", "p1=10...3 2...5 ~p1\n");
    let (st, out, _) = run(&["--lint", &bad, &rev]);
    assert_eq!(st, Some(1));
    assert!(out.contains("error: p9 is used but not defined"), "{out}");
    assert!(out.contains("warning: `p1=10...3`"), "{out}");
    // a pattern that fails to parse
    let junk = pattern("junk.pat", "p1=3...3 GAG@ p1\n");
    let (st, out, _) = run(&["--explain", &junk]);
    assert_eq!(st, Some(1));
    assert!(out.contains("error"), "{out}");
    // -o -c: -c is the value of -o
    let (_, out, _) = run(&["--lint", "-o", "-c", DTC]);
    assert!(out.contains("-o takes a value"), "{out}");
}

#[test]
fn explain_and_option_errors() {
    let (st, out, _) = run(&["--explain", DTC]);
    assert_eq!(st, Some(0));
    assert!(out.contains("~p2           the reverse complement of p2, exactly"));
    assert!(out.contains("hit length: 530 to 15030 letters"));
    let (st, _, err) = run(&["--lint", "--explain", DTC]);
    assert_eq!(st, Some(2));
    assert!(err.contains("not both"));
}

/// `[m,d,i]` as in README.original: the 2nd number lets the data miss a
/// letter (deletion), the 3rd lets it have an extra letter (insertion).
#[test]
fn explain_deletions_and_insertions() {
    let p = pattern("di.pat", "ACGT[1,1,2]\n");
    let (st, out, _) = run(&["--explain", &p]);
    assert_eq!(st, Some(0));
    assert!(
        out.contains("ACGT, with up to 1 mismatch, 1 deletion, 2 insertions"),
        "{out}"
    );
    assert!(out.contains("hit length: 3 to 6 letters"), "{out}");
}
