//! FASTA input compressed with gzip or bgzip (a series of gzip members).

use std::io::Write;
use std::process::{Command, Stdio};

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut child = Command::new("gzip")
        .arg("-c")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("gzip command");
    let mut stdin = child.stdin.take().unwrap();
    let data = data.to_vec();
    let feeder = std::thread::spawn(move || stdin.write_all(&data).unwrap());
    let out = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    out.stdout
}

fn scan(args: &[&str], pattern: &str, input: &[u8]) -> (Option<i32>, Vec<u8>, String) {
    scan_in_pieces(args, pattern, input, 0)
}

/// As `scan`, but the input is written in pieces of `piece` bytes with a
/// short pause after each (0: all at once).
fn scan_in_pieces(
    args: &[&str],
    pattern: &str,
    input: &[u8],
    piece: usize,
) -> (Option<i32>, Vec<u8>, String) {
    let dir = std::env::temp_dir().join(format!(
        "sfm_gz_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let pat = dir.join("pattern");
    std::fs::write(&pat, format!("{pattern}\n")).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"))
        .args(args)
        .arg(&pat)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    let input = input.to_vec();
    let feeder = std::thread::spawn(move || {
        if piece == 0 {
            let _ = stdin.write_all(&input);
            return;
        }
        for p in input.chunks(piece) {
            if stdin.write_all(p).and_then(|_| stdin.flush()).is_err() {
                return;
            }
            std::thread::sleep(std::time::Duration::from_micros(200));
        }
    });
    let out = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    (
        out.status.code(),
        out.stdout,
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// Some records of random DNA with hairpins.
fn fasta() -> Vec<u8> {
    let mut x: u64 = 99;
    let mut out = String::new();
    for r in 0..40 {
        out.push_str(&format!(">rec{r} test record\n"));
        let s: String = (0..3000)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"acgtACGT"[(x % 8) as usize] as char
            })
            .collect();
        for line in s.as_bytes().chunks(70) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
    }
    out.into_bytes()
}

const PAT: &str = "p1=4...6 2...6 ~p1";

#[test]
fn gzip_gives_same_output() {
    let plain = fasta();
    let a = scan(&["-c"], PAT, &plain);
    let b = scan(&["-c"], PAT, &gzip(&plain));
    assert_eq!(a.0, Some(0));
    assert!(!a.1.is_empty());
    assert_eq!(a, b);
}

#[test]
fn several_gzip_members_like_bgzip() {
    let plain = fasta();
    let cut = plain.len() / 3;
    let cut2 = 2 * plain.len() / 3;
    let mut multi = gzip(&plain[..cut]);
    multi.extend(gzip(&plain[cut..cut2]));
    multi.extend(gzip(&plain[cut2..]));
    multi.extend(gzip(b"")); // bgzip ends with an empty block
    assert_eq!(scan(&["-c"], PAT, &plain), scan(&["-c"], PAT, &multi));
}

#[test]
fn truncated_gzip_is_an_error() {
    let gz = gzip(&fasta());
    let (code, _, err) = scan(&[], PAT, &gz[..gz.len() / 2]);
    assert_eq!(code, Some(1));
    assert_eq!(err, "gzip input: unexpected end of compressed data\n");
}

#[test]
fn damaged_gzip_is_an_error() {
    let mut gz = gzip(&fasta());
    for b in &mut gz[100..200] {
        *b = 0;
    }
    let (code, _, err) = scan(&[], PAT, &gz);
    assert_eq!(code, Some(1));
    assert_eq!(err, "gzip input: damaged compressed data\n");
}

#[test]
fn damaged_gzip_gives_the_same_output_however_it_arrives() {
    // the output before the error must not depend on the sizes of the
    // reads from the pipe
    let mut gz = gzip(&fasta());
    let n = gz.len();
    for b in &mut gz[n / 2..n / 2 + 50] {
        *b = 0x55;
    }
    let whole = scan(&["-c"], PAT, &gz);
    assert_eq!(whole.0, Some(1));
    assert!(!whole.1.is_empty(), "hits before the damage are printed");
    for piece in [1000, 4093] {
        assert_eq!(
            whole,
            scan_in_pieces(&["-c"], PAT, &gz, piece),
            "pieces of {piece}"
        );
    }
}
