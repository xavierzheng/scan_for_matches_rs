//! --strict-n: a letter other than A, C, G, T that a name caught is one
//! mismatch when `p1` / `<p1` uses the name again.

use std::io::Write;
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_scan_for_matches");

/// Hits (`id:[a,b]`) of `pat` on `fasta`; exit status and stderr.
fn hits(pat: &str, args: &[&str], fasta: &str, env: &[(&str, &str)]) -> (Vec<String>, i32, String) {
    let d = std::env::temp_dir().join(format!("sfm_strict_{}_{:x}", std::process::id(), {
        use std::hash::{BuildHasher, RandomState};
        RandomState::new().hash_one((pat, args))
    }));
    std::fs::create_dir_all(&d).unwrap();
    let pp = d.join("p.pat");
    std::fs::write(&pp, format!("{pat}\n")).unwrap();
    let mut child = Command::new(BIN)
        .args(args)
        .arg(&pp)
        .envs(env.iter().copied())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(fasta.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    std::fs::remove_file(&pp).unwrap();
    std::fs::remove_dir(&d).unwrap();
    let h = String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter_map(|l| l.strip_prefix('>').map(str::to_string))
        .collect();
    (
        h,
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stderr).into(),
    )
}

/// The ids with a hit, without and with --strict-n.
fn ids(pat: &str, fasta: &str) -> (String, String) {
    let f = |args: &[&str]| {
        let (h, _, _) = hits(pat, args, fasta, &[]);
        h.iter()
            .map(|h| h.split(':').next().unwrap())
            .collect::<Vec<_>>()
            .join(" ")
    };
    (f(&[]), f(&["--strict-n"]))
}

#[test]
fn the_rule() {
    // name p1 = 4 letters, 6 letters between
    let fa = ">ok\nGCTACCCCCCGCTA\n>n1\nGNTACCCCCCGCTA\n>n2\nGNNACCCCCCGCTA\n>r1\nGRTACCCCCCGCTA\n";
    // exact TSD: 1 N is too many with --strict-n
    assert_eq!(
        ids("^ p1=4...4 6...6 p1 $", fa),
        ("ok n1 n2".into(), "ok".into())
    );
    // 1 mismatch allowed: 1 N or 1 R uses it, 2 N are too many
    assert_eq!(
        ids("^ p1=4...4 6...6 p1[1,0,0] $", fa),
        ("ok n1 n2 r1".into(), "ok n1 r1".into())
    );
    // <p1 the same
    let fa = ">ok\nGCTACCCCCCATCG\n>n1\nGNTACCCCCCATCG\n";
    assert_eq!(
        ids("^ p1=4...4 6...6 <p1 $", fa),
        ("ok n1".into(), "ok".into())
    );
    assert_eq!(
        ids("^ p1=4...4 6...6 <p1[1,0,0] $", fa),
        ("ok n1".into(), "ok n1".into())
    );
    // ~p1 never matched a caught N (also without the option)
    let fa = ">ok\nGCTACCCCCCTAGC\n>n1\nGNTACCCCCCTAGC\n";
    assert_eq!(
        ids("^ p1=4...4 6...6 ~p1[1,0,0] $", fa),
        ("ok".into(), "ok".into())
    );
    // N in the gap: no effect
    let fa = ">gap\nGCTACCNNNNGCTA\n";
    assert_eq!(
        ids("^ p1=4...4 6...6 p1 $", fa),
        ("gap".into(), "gap".into())
    );
}

/// An assembly gap before a real element: with a wide gap range the
/// search skips gap lengths; the false TSD must not hide the real one.
#[test]
fn gap_skipping_and_threads() {
    let mut s = String::from("ACGT");
    s += &"N".repeat(300);
    let mut x: u64 = 3;
    let mut rnd = |n: usize| -> String {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize] as char
            })
            .collect()
    };
    s += &rnd(500);
    s += "GATTACAGATCC";
    s += &rnd(200);
    s += "GATTACAGATCC";
    s += &rnd(3000);
    let fa = format!(">chr\n{s}\n");
    let pat = "p1=12...12 150...250 p1";
    let start = |h: &String| -> usize { h[5..h.find(',').unwrap()].parse().unwrap() };
    let (plain, _, _) = hits(pat, &[], &fa, &[]);
    let (strict, _, _) = hits(pat, &["--strict-n"], &fa, &[]);
    // without the option the N run is a TSD; with it only real TSDs
    assert!(plain.first().is_some_and(|h| start(h) <= 304), "{plain:?}");
    assert!(strict.iter().any(|h| h == "chr:[805,1028]"), "{strict:?}");
    assert!(strict.iter().all(|h| start(h) > 304), "{strict:?}");
    // threads and small pieces give the same hits
    for piece in ["50", "300", "1000"] {
        let (t, _, _) = hits(
            pat,
            &["--strict-n", "-t", "4"],
            &fa,
            &[("SFM_PIECE", piece)],
        );
        assert_eq!(t, strict, "SFM_PIECE={piece}");
    }
}

#[test]
fn not_with_protein() {
    let (_, code, err) = hits(
        "p1=2...2 1...3 p1",
        &["--strict-n", "-p"],
        ">x\nMKMK\n",
        &[],
    );
    assert_eq!(code, 2);
    assert!(err.contains("--strict-n is for DNA"), "{err}");
}
