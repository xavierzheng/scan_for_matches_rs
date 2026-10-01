//! Golden cases recorded from the C reference program with
//! `tests/make_golden.py` and split with `tests/split_golden.py`: pattern,
//! options and input, with the expected stdout and exit status (exit code
//! or killing signal).
//!
//! * `golden.tsv` – cases where this version must equal the C program.
//! * `fixed.tsv` – cases changed by the fixes of version 0.1.0 (expected
//!   results are this version's).

use std::io::Write;
use std::os::unix::process::ExitStatusExt;
use std::process::{Command, Stdio};

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect()
}

struct Case {
    line: usize,
    pattern: Vec<u8>,
    args: Vec<String>,
    input: Vec<u8>,
    status: String,
    stdout: Vec<u8>,
}

fn run(c: &Case, dir: &std::path::Path) -> Result<(), String> {
    let pat = dir.join(format!("case{}.pat", c.line));
    std::fs::write(&pat, &c.pattern).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"))
        .args(&c.args)
        .arg(&pat)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let input = c.input.clone();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&input);
    });
    let out = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    let status = match out.status.code() {
        Some(code) => code.to_string(),
        None => format!("sig{}", out.status.signal().unwrap_or(0)),
    };
    if status != c.status || out.stdout != c.stdout {
        return Err(format!(
            "golden line {}: pattern {:?} args {:?}: status {} (want {}), stdout {} bytes (want {}){}",
            c.line,
            String::from_utf8_lossy(&c.pattern),
            c.args,
            status,
            c.status,
            out.stdout.len(),
            c.stdout.len(),
            if out.stdout != c.stdout {
                ", stdout differs"
            } else {
                ""
            }
        ));
    }
    Ok(())
}

#[test]
fn golden_same_as_c() {
    check("golden.tsv", 500);
}

#[test]
fn golden_fixed() {
    check("fixed.tsv", 200);
}

fn check(file: &str, min_cases: usize) {
    let text = std::fs::read_to_string(format!("{}/tests/data/{file}", env!("CARGO_MANIFEST_DIR")))
        .unwrap();
    let cases: Vec<Case> = text
        .lines()
        .enumerate()
        .map(|(i, l)| {
            let f: Vec<&str> = l.split('\t').collect();
            let args = unhex(f[1]);
            let args = if args.is_empty() {
                Vec::new()
            } else {
                args.split(|&b| b == 0)
                    .map(|a| String::from_utf8(a.to_vec()).unwrap())
                    .collect()
            };
            Case {
                line: i + 1,
                pattern: unhex(f[0]),
                args,
                input: unhex(f[2]),
                status: f[3].to_string(),
                stdout: unhex(f[4]),
            }
        })
        .collect();
    assert!(cases.len() >= min_cases);
    let dir = std::env::temp_dir().join(format!("sfm_golden_{}_{file}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    // at most 8 workers (leave cores for other work)
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(8);
    let chunk = cases.len().div_ceil(threads);
    let failures: Vec<String> = std::thread::scope(|s| {
        let hs: Vec<_> = cases
            .chunks(chunk)
            .map(|part| {
                let dir = dir.clone();
                s.spawn(move || {
                    part.iter()
                        .filter_map(|c| run(c, &dir).err())
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        hs.into_iter().flat_map(|h| h.join().unwrap()).collect()
    });
    let _ = std::fs::remove_dir_all(&dir);
    assert!(
        failures.is_empty(),
        "{file}: {} of {} cases differ:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}
