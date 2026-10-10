//! --help, --input and --output.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_scan_for_matches");
const PAT: &str = "p1=4...6 2...6 ~p1";

fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("sfm_cli_{}_{name}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fasta() -> Vec<u8> {
    let mut x: u64 = 7;
    let mut out = String::new();
    for r in 0..20 {
        out.push_str(&format!(">rec{r}\n"));
        for _ in 0..30 {
            let line: String = (0..60)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"acgt"[(x % 4) as usize] as char
                })
                .collect();
            out.push_str(&line);
            out.push('\n');
        }
    }
    out.into_bytes()
}

/// Run with `stdin` as input; (exit status, stdout, stderr).
fn run(args: &[&str], stdin: &[u8]) -> (Option<i32>, Vec<u8>, String) {
    let mut child = Command::new(BIN)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut i = child.stdin.take().unwrap();
    let data = stdin.to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = i.write_all(&data);
    });
    let out = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    (
        out.status.code(),
        out.stdout,
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn help_lists_the_options() {
    for a in [&["--help"][..], &["-h"], &["-c", "-h", "x.pat"]] {
        let (code, out, err) = run(a, b"");
        let out = String::from_utf8(out).unwrap();
        assert_eq!(code, Some(0), "{a:?}");
        assert!(err.is_empty());
        for opt in [
            "-c ", "-p ", "-n N", "-m N", "-o N", "-i FILE", "-t N", "--input", "--output",
            "--format", "--dedup", "--help",
        ] {
            assert!(out.contains(opt), "{opt} in help");
        }
        assert!(out.contains(env!("CARGO_PKG_VERSION")));
    }
}

#[test]
fn h_as_value_of_an_option_is_not_help() {
    // -o takes a value (not used): "-o -h" is not a request for help
    let d = dir("hval");
    let pat = d.join("p");
    std::fs::write(&pat, format!("{PAT}\n")).unwrap();
    let (code, out, _) = run(&["-o", "-h", pat.to_str().unwrap()], &fasta());
    assert_eq!(code, Some(0));
    assert!(!String::from_utf8_lossy(&out).contains("usage"));
}

#[test]
fn input_and_output_files_give_the_same_hits() {
    let d = dir("io");
    let pat = d.join("p");
    std::fs::write(&pat, format!("{PAT}\n")).unwrap();
    let fa = d.join("in.fa");
    std::fs::write(&fa, fasta()).unwrap();
    let p = pat.to_str().unwrap();
    for extra in [&["-c"][..], &["-c", "-t", "3"], &["-c", "--format", "gff3"]] {
        let mut a = extra.to_vec();
        a.push(p);
        let (code, want, _) = run(&a, &fasta());
        assert_eq!(code, Some(0));
        assert!(!want.is_empty());

        // --input
        let mut b = extra.to_vec();
        b.extend(["--input", fa.to_str().unwrap(), p]);
        assert_eq!(run(&b, b""), (Some(0), want.clone(), String::new()));

        // --output (with --input=, and - for stdin)
        let outf = d.join("out.txt");
        let o = format!("--output={}", outf.to_str().unwrap());
        let mut c = extra.to_vec();
        c.extend(["--input", "-", &o, p]);
        let (code, so, se) = run(&c, &fasta());
        assert_eq!((code, so.len(), se.as_str()), (Some(0), 0, ""));
        assert_eq!(std::fs::read(&outf).unwrap(), want);
    }
}

#[test]
fn missing_input_file_is_an_error() {
    let d = dir("miss");
    let pat = d.join("p");
    std::fs::write(&pat, format!("{PAT}\n")).unwrap();
    let (code, out, err) = run(&["--input", "/no/such/file.fa", pat.to_str().unwrap()], b"");
    assert_eq!(code, Some(2));
    assert!(out.is_empty());
    assert_eq!(
        err,
        "scan_for_matches: cannot open input file /no/such/file.fa\n"
    );
    let (code, _, err) = run(&[pat.to_str().unwrap(), "--input"], b"");
    assert_eq!(code, Some(2));
    assert_eq!(err, "scan_for_matches: option --input needs a value\n");
}

#[test]
fn progress_goes_to_stderr_only() {
    let d = dir("prog");
    let pat = d.join("p");
    std::fs::write(&pat, format!("{PAT}\n")).unwrap();
    let p = pat.to_str().unwrap();
    let (_, plain, _) = run(&["-c", p], &fasta());
    // every hit counts (no --dedup)
    let hits = plain.iter().filter(|&&b| b == b'>').count();
    assert!(hits > 0);
    for extra in [&["-c"][..], &["-c", "-t", "3"], &["-c", "--format", "gff3"]] {
        let mut a = extra.to_vec();
        a.push(p);
        let (_, want, _) = run(&a, &fasta());
        let mut b = vec!["--progress"];
        b.extend(&a);
        // small pieces: the long-record path of -t is used too
        let mut child = Command::new(BIN)
            .args(&b)
            .env("SFM_PIECE", "300")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(&fasta()).unwrap();
        let out = child.wait_with_output().unwrap();
        assert_eq!(out.status.code(), Some(0));
        assert_eq!(out.stdout, want, "{extra:?}");
        let err = String::from_utf8(out.stderr).unwrap();
        let lines: Vec<&str> = err.lines().collect();
        assert_eq!(lines.len(), 2, "{err}");
        assert!(lines[0].ends_with("0 records (0.0 Mb) done, 0 hits; started"), "{err}");
        let fin = format!("20 records (0.0 Mb) done, {hits} hits; finished");
        assert!(lines[1].ends_with(&fin), "{err} / {fin}");
    }
}
