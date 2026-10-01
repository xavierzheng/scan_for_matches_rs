//! Golden cases recorded from the C reference program with
//! `tests/make_golden.py`: pattern, options and input, with the expected
//! stdout and exit status (exit code or killing signal).

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
fn golden_cases() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/golden.tsv"
    ))
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
    assert!(cases.len() > 500);
    let dir = std::env::temp_dir().join(format!("sfm_golden_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
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
        "{} of {} golden cases differ:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}
