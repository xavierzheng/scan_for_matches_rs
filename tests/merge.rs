//! `--merge`: join GFF3 files, remove duplicates, new unique Names.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_scan_for_matches");

fn dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("sfm_merge_{}_{name}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Write `text` to a file in the test directory; returns the path.
fn file(d: &std::path::Path, name: &str, text: &[u8]) -> String {
    let p = d.join(name);
    std::fs::write(&p, text).unwrap();
    p.to_str().unwrap().to_string()
}

/// Run `scan_for_matches --merge ARGS` with `stdin`; (status, stdout, stderr).
fn merge(args: &[&str], stdin: &[u8]) -> (Option<i32>, String, String) {
    let mut child = Command::new(BIN)
        .arg("--merge")
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(stdin).unwrap();
    let o = child.wait_with_output().unwrap();
    (
        o.status.code(),
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    )
}

/// One element: a top feature `name` with a child `name.te`.
fn elem(seq: &str, s: u32, e: u32, strand: &str, name: &str) -> String {
    format!(
        "{seq}\tx\trepeat_region\t{s}\t{e}\t.\t{strand}\t.\tID={name};Name={name};K=v\n\
         {seq}\tx\tTE\t{}\t{}\t.\t{strand}\t.\tID={name}.te;Parent={name};Name={name}\n###\n",
        s + 1,
        e - 1
    )
}

/// The Name= values of the top features (`repeat_region` lines).
fn tops(out: &str) -> Vec<String> {
    out.lines()
        .filter(|l| l.split('\t').nth(2) == Some("repeat_region"))
        .map(|l| l.split('\t').take(5).collect::<Vec<_>>().join(" "))
        .collect()
}

#[test]
fn same_span_different_strand_priority_and_merged() {
    let d = dir("dup");
    let a = file(&d, "a.gff3", elem("c1", 10, 50, "+", "AA1").as_bytes());
    let b = file(&d, "b.gff3", elem("c1", 10, 50, ".", "BB7").as_bytes());
    let (st, out, err) = merge(&[&a, &b], b"");
    assert_eq!(st, Some(0), "{err}");
    assert!(
        err.contains("2 elements read, 1 duplicates removed, 1 written"),
        "{err}"
    );
    assert_eq!(tops(&out), ["c1 x repeat_region 10 50"]);
    assert!(
        out.contains("\t+\t.\tID=AA1;Name=AA1;K=v;Merged=b.gff3:BB7\n"),
        "{out}"
    );
    assert!(!out.contains("BB7."), "{out}");
    // reversed priority: the other one wins
    let (_, out, _) = merge(&[&b, &a], b"");
    assert!(
        out.contains("\t.\t.\tID=BB1;Name=BB1;K=v;Merged=a.gff3:AA1\n"),
        "{out}"
    );
}

#[test]
fn overlap_is_reciprocal() {
    let mut text = elem("c1", 1, 100, "+", "L1");
    text += &elem("c1", 10, 20, "+", "S1");
    let d = dir("ov");
    let a = file(&d, "a.gff3", text.as_bytes());
    let (_, out, err) = merge(&["--overlap", "0.9", &a], b"");
    assert!(err.contains("0 duplicates"), "{err}");
    assert_eq!(tops(&out).len(), 2);
    let (_, _, err) = merge(&["--overlap=0.1", &a], b"");
    assert!(err.contains("1 duplicates"), "{err}");
    // both ways: 95 of 100 and 95 of 100 merge; shift of 20 does not
    let x = file(&d, "x.gff3", elem("c1", 1, 100, "+", "X1").as_bytes());
    let y = file(&d, "y.gff3", elem("c1", 6, 105, "+", "Y1").as_bytes());
    let z = file(&d, "z.gff3", elem("c1", 21, 120, "+", "Z1").as_bytes());
    assert!(
        merge(&["--overlap", "0.9", &x, &y], b"")
            .2
            .contains("1 duplicates")
    );
    assert!(
        merge(&["--overlap", "0.9", &x, &z], b"")
            .2
            .contains("0 duplicates")
    );
    assert!(
        merge(&["--overlap", "0.9", &z, &x], b"")
            .2
            .contains("0 duplicates")
    );
    // without --overlap only the exact span counts
    assert!(merge(&[&x, &y], b"").2.contains("0 duplicates"));
}

#[test]
fn natural_seqid_order_and_regions() {
    let mut text = String::from("##sequence-region chr10 1 999\n##sequence-region chr2 1 500\n");
    text += &elem("chr10", 5, 20, "+", "T1");
    text += &elem("chr2", 30, 40, "+", "T2");
    text += &elem("chr2", 30, 60, "+", "T3");
    text += &elem("chr2", 30, 60, "+", "T4");
    text += &elem("chr2", 5, 9, "+", "T5");
    let d = dir("nat");
    let a = file(&d, "a.gff3", text.as_bytes());
    let (st, out, _) = merge(&[&a], b"");
    assert_eq!(st, Some(0));
    assert_eq!(
        tops(&out),
        [
            "chr2 x repeat_region 5 9",
            "chr2 x repeat_region 30 60",
            "chr2 x repeat_region 30 40",
            "chr10 x repeat_region 5 20"
        ]
    );
    let regions: Vec<&str> = out.lines().filter(|l| l.starts_with("##seq")).collect();
    assert_eq!(
        regions,
        [
            "##sequence-region chr2 1 500",
            "##sequence-region chr10 1 999"
        ]
    );
}

#[test]
fn header_and_renaming() {
    let mut text = String::from("##gff-version 3\n# comment\n\n");
    text += "c1\tx\tTSD\t10\t12\t.\t+\t.\tID=D5.lTSD;Parent=D5;Name=D5;Note=a b\n";
    text += "c1\tx\tTIR\t13\t15\t.\t+\t.\tID=D5.l;Parent=D5.te\n";
    text += "c1\tx\tTE\t13\t40\t.\t+\t.\tID=D5.te;Parent=D5\n";
    text += "c1\tx\trepeat_region\t10\t43\t.\t+\t.\tID=D5\n###\n##FASTA\n>c1\nACGT\n";
    let d = dir("ren");
    let a = file(&d, "a.gff3", text.as_bytes());
    let (st, out, err) = merge(&["--name-start", "7", &a], b"");
    assert_eq!(st, Some(0), "{err}");
    let want = "##gff-version 3\n".to_string()
        + &format!(
            "# scan_for_matches {} merge {a}\n",
            env!("CARGO_PKG_VERSION")
        )
        + "c1\tx\tTSD\t10\t12\t.\t+\t.\tID=D7.lTSD;Parent=D7;Name=D7;Note=a b\n"
        + "c1\tx\tTIR\t13\t15\t.\t+\t.\tID=D7.l;Parent=D7.te\n"
        + "c1\tx\tTE\t13\t40\t.\t+\t.\tID=D7.te;Parent=D7\n"
        + "c1\tx\trepeat_region\t10\t43\t.\t+\t.\tID=D7;Name=D7\n###\n";
    assert_eq!(out, want);
}

#[test]
fn multi_level_and_foreign_ids() {
    let text = "c\tx\tg\t5\t9\t.\t+\t.\tID=m1;Parent=t1,t2\n\
                c\tx\tt\t1\t20\t.\t+\t.\tID=t2;Parent=g1\n\
                c\tx\tt\t1\t20\t.\t+\t.\tID=t1;Parent=g1\n\
                c\tx\tgene\t1\t20\t.\t+\t.\tID=gene_9;Name=gene_9\n\
                c\tx\tgene\t1\t20\t.\t+\t.\tID=g1\n";
    // 'gene_9' and 'g1' are two tops with the same span: the second is a duplicate.
    let (st, out, err) = merge(&["--name-prefix", "N_", "-"], text.as_bytes());
    assert_eq!(st, Some(0), "{err}");
    assert!(
        err.contains("2 elements read, 1 duplicates removed, 1 written"),
        "{err}"
    );
    assert!(out.contains("ID=N_1;Name=N_1;Merged=stdin:g1\n"), "{out}");
    // a single element with a forward reference
    let text = text.replace(
        "\t1\t20\t.\t+\t.\tID=gene_9;Name=gene_9",
        "\t100\t200\t.\t+\t.\tID=zz",
    );
    let (_, out, _) = merge(
        &["--name-prefix=N_", "--name-start=3", "-"],
        text.as_bytes(),
    );
    assert!(out.contains("\tID=N_3.m1;Parent=N_3.t1,N_3.t2\n"), "{out}");
    assert!(out.contains("\tID=N_3.t1;Parent=N_3\n"), "{out}");
    assert!(out.ends_with("###\n"));
}

#[test]
fn counters_per_prefix() {
    let mut text = elem("c1", 1, 9, "+", "DTC5");
    text += &elem("c1", 20, 29, "+", "DTH_98");
    text += &elem("c1", 30, 39, "+", "DTC2");
    text += &elem("c1", 40, 49, "+", "77");
    let d = dir("cnt");
    let a = file(&d, "a.gff3", text.as_bytes());
    let (_, out, _) = merge(&[&a], b"");
    let names: Vec<&str> = out
        .lines()
        .filter_map(|l| l.split("\tID=").nth(1))
        .filter(|l| !l.contains("Parent"))
        .collect();
    let ids: Vec<&str> = names.iter().map(|l| l.split(';').next().unwrap()).collect();
    assert_eq!(ids, ["DTC1", "DTH_1", "DTC2", "sfm1"]);
}

/// A gzip file with one stored block (no compression), no crates needed.
fn gzip(data: &[u8]) -> Vec<u8> {
    let mut crc = !0u32;
    for &b in data {
        crc ^= b as u32;
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xEDB8_8320 & (crc & 1).wrapping_neg());
        }
    }
    let mut o = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255, 1];
    o.extend((data.len() as u16).to_le_bytes());
    o.extend((!(data.len() as u16)).to_le_bytes());
    o.extend(data);
    o.extend((!crc).to_le_bytes());
    o.extend((data.len() as u32).to_le_bytes());
    o
}

#[test]
fn gzip_input_file_and_stdin() {
    let d = dir("gz");
    let a = file(
        &d,
        "a.gff3.gz",
        &gzip(elem("c1", 10, 50, "+", "AA1").as_bytes()),
    );
    let (st, out, err) = merge(&[&a], b"");
    assert_eq!(st, Some(0), "{err}");
    assert_eq!(tops(&out), ["c1 x repeat_region 10 50"]);
    let (_, out2, _) = merge(&["-"], &gzip(elem("c1", 10, 50, "+", "AA1").as_bytes()));
    assert_eq!(tops(&out2), tops(&out));
}

#[test]
fn output_file() {
    let d = dir("outf");
    let a = file(&d, "a.gff3", elem("c1", 10, 50, "+", "AA1").as_bytes());
    let o = d.join("o.gff3");
    let (st, out, _) = merge(&["--output", o.to_str().unwrap(), &a], b"");
    assert_eq!((st, out.as_str()), (Some(0), ""));
    assert!(
        std::fs::read_to_string(o)
            .unwrap()
            .contains("ID=AA1;Name=AA1")
    );
}

#[test]
fn data_errors_exit_1() {
    let d = dir("err");
    let cases: [(&str, &str); 4] = [
        ("c\tx\tt\t1\t9\t.\t+\t.\n", "a.gff3:1: "),
        ("c\tx\tt\tzz\t9\t.\t+\t.\tID=a\n", "a.gff3:1: bad start"),
        (
            "# x\nc\tx\tt\t1\t9\t.\t+\t.\tID=a;Parent=nope\n",
            "a.gff3:2: ",
        ),
        ("c\tx\tt\t5\t1\t.\t+\t.\tID=a\n", "a.gff3:1: "),
    ];
    for (text, want) in cases {
        let a = file(&d, "a.gff3", text.as_bytes());
        let (st, out, err) = merge(&[&a], b"");
        assert_eq!(st, Some(1), "{text}");
        assert!(out.is_empty());
        assert!(err.starts_with("merge: ") && err.contains(want), "{err}");
    }
    // parents in two elements
    let t = "c\tx\tt\t1\t9\t.\t+\t.\tID=a\nc\tx\tt\t1\t9\t.\t+\t.\tID=b\nc\tx\tt\t2\t3\t.\t+\t.\tID=k;Parent=a,b\n";
    let a = file(&d, "a.gff3", t.as_bytes());
    let (st, _, err) = merge(&[&a], b"");
    assert_eq!(st, Some(1));
    assert!(err.contains("a.gff3:3: "), "{err}");
    // two runs joined with cat: the same ID on two elements
    let t = format!(
        "{}{}",
        elem("c", 1, 9, "+", "X1"),
        elem("c", 20, 29, "+", "X1")
    );
    let a = file(&d, "a.gff3", t.as_bytes());
    let (st, _, err) = merge(&[&a], b"");
    assert_eq!(st, Some(1));
    assert!(err.contains("a.gff3:4: ID 'X1' is used twice"), "{err}");
    // one feature on two lines (same ID, same parent) is fine
    let t = "c\tx\tt\t1\t30\t.\t+\t.\tID=g\nc\tx\tCDS\t2\t5\t.\t+\t.\tID=cds;Parent=g\nc\tx\tCDS\t9\t12\t.\t+\t.\tID=cds;Parent=g\n";
    let a = file(&d, "a.gff3", t.as_bytes());
    let (st, _, err) = merge(&[&a], b"");
    assert_eq!(st, Some(0), "{err}");
}

#[test]
fn option_errors_exit_2_and_help() {
    let d = dir("opt");
    let a = file(&d, "a.gff3", elem("c1", 10, 50, "+", "AA1").as_bytes());
    let cases: [&[&str]; 8] = [
        &[],
        &["--bogus", &a],
        &["--overlap", &a],
        &["--overlap", "0", &a],
        &["--overlap=2", &a],
        &["--name-start", "x", &a],
        &["--name-start=-1", &a],
        &["/nonexistent/none.gff3"],
    ];
    for c in cases {
        let (st, out, err) = merge(c, b"");
        assert_eq!(st, Some(2), "{c:?}: {err}");
        assert!(out.is_empty() && err.starts_with("merge: "), "{c:?}: {err}");
    }
    for h in ["-h", "--help"] {
        let (st, out, _) = merge(&[h], b"");
        assert_eq!(st, Some(0));
        assert!(out.starts_with("usage: scan_for_matches --merge"), "{out}");
    }
}

/// Run the scanner with a labelled pattern on `fa`; the GFF3 text.
fn scan(d: &std::path::Path, pat: &str, fa: &str) -> String {
    let p = file(d, "p.pat", pat.as_bytes());
    let mut c = Command::new(BIN)
        .args(["--format", "gff3", &p])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    c.stdin.take().unwrap().write_all(fa.as_bytes()).unwrap();
    String::from_utf8(c.wait_with_output().unwrap().stdout).unwrap()
}

#[test]
fn round_trip_two_patterns() {
    let d = dir("rt");
    let filler = "TTTTGGGAAA".repeat(8);
    let el = |tsd: &str| format!("{tsd}CACTAGCTCAA{filler}TTGAGCTAGTG{tsd}");
    let fa = format!(
        ">chr10\n{}{}{}\n>chr2\n{}{}{}\n",
        filler,
        el("AGT"),
        filler,
        filler,
        el("GGA"),
        filler
    );
    let pat = |n: u32| {
        format!(
            "%@element CACTA_TIR_transposon Name=DTC\n\
             p1=3...3 %@ target_site_duplication\n\
             CACTA p2={n}...{n} %@ five_prime_terminal_inverted_repeat\n\
             20...200\n\
             ~p2 TAGTG %@ three_prime_terminal_inverted_repeat\n\
             p1 %@ target_site_duplication\n"
        )
    };
    let a = file(&d, "a.gff3", scan(&d, &pat(4), &fa).as_bytes());
    let b = file(&d, "b.gff3", scan(&d, &pat(6), &fa).as_bytes());
    let (st, out, err) = merge(&[&a, &b], b"");
    assert_eq!(st, Some(0), "{err}");
    assert!(
        err.contains("4 elements read, 2 duplicates removed, 2 written"),
        "{err}"
    );
    let t = tops(&out);
    assert_eq!(t.len(), 2);
    assert!(
        t[0].starts_with("chr2 ") && t[1].starts_with("chr10 "),
        "{t:?}"
    );
    assert!(
        out.contains("ID=DTC1;Name=DTC1;") && out.contains("Merged=b.gff3:DTC"),
        "{out}"
    );
    assert!(
        out.contains("ID=DTC2.te;Parent=DTC2;") && out.contains("ID=DTC2.rTSD;"),
        "{out}"
    );
    assert_eq!(out.matches("\n###\n").count(), 2);
}
