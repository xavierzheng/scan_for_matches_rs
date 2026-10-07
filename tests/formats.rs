//! `--format gff3|bed6|bed12|jsonl` and the `%@` labels.

use std::io::Write;
use std::process::{Command, Stdio};

type Out = (Option<i32>, String, String);

fn scan(args: &[&str], pattern: &str, input: &str) -> Out {
    let dir = std::env::temp_dir().join(format!(
        "sfm_fmt_{}_{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("pat"), pattern).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_scan_for_matches"))
        .current_dir(&dir)
        .args(args)
        .arg("pat")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let o = child.wait_with_output().unwrap();
    (
        o.status.code(),
        String::from_utf8(o.stdout).unwrap(),
        String::from_utf8(o.stderr).unwrap(),
    )
}

const CACTA: &str = "\
%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC
p1=3...3               %@ target_site_duplication
CACTA p2=7...7         %@ five_prime_terminal_inverted_repeat
10...50
~p2 TAGTG              %@ three_prime_terminal_inverted_repeat
p1                     %@ target_site_duplication
";

/// 20 bp, element (TSD GAT, TIRs CACTAACGTTGC / GCAACGTTAGTG), 15 bp
const SEQ: &str = ">chr1 x\nCCTTGGAACCTTGGAACCTT\
GATCACTAACGTTGCAAAAATTTTTCCCCCGGGGGAAAAATTTTTGCAACGTTAGTGGAT\
CCTTGGAACCTTGGA\n";

#[test]
fn gff3_cacta() {
    let (st, out, _) = scan(&["--format", "gff3"], CACTA, SEQ);
    assert_eq!(st, Some(0));
    let cols = |l: &str| l.split('\t').map(String::from).collect::<Vec<_>>();
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(lines[0], "##gff-version 3");
    assert_eq!(
        lines[1],
        format!(
            "# scan_for_matches {} pattern=pat",
            env!("CARGO_PKG_VERSION")
        )
    );
    assert_eq!(lines[2], "##sequence-region chr1 1 95");
    let want = [
        (
            "repeat_region",
            21,
            80,
            "ID=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657",
        ),
        (
            "target_site_duplication",
            21,
            23,
            "ID=DTC1.lTSD;Parent=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000434",
        ),
        (
            "CACTA_TIR_transposon",
            24,
            77,
            "ID=DTC1.te;Parent=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0002285;TSD=GAT_GAT;TIR=CACTAACGTTGC_GCAACGTTAGTG",
        ),
        (
            "five_prime_terminal_inverted_repeat",
            24,
            35,
            "ID=DTC1.lTIR;Parent=DTC1.te;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000420",
        ),
        (
            "three_prime_terminal_inverted_repeat",
            66,
            77,
            "ID=DTC1.rTIR;Parent=DTC1.te;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000421",
        ),
        (
            "target_site_duplication",
            78,
            80,
            "ID=DTC1.rTSD;Parent=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000434",
        ),
    ];
    for (l, (t, a, b, attrs)) in lines[3..9].iter().zip(want) {
        let c = cols(l);
        assert_eq!(
            c,
            vec![
                "chr1",
                "scan_for_matches",
                t,
                &a.to_string(),
                &b.to_string(),
                ".",
                "+",
                ".",
                attrs
            ]
        );
    }
    assert_eq!(lines[9], "###");
    assert_eq!(lines.len(), 10);
}

#[test]
fn dedup_and_names() {
    // the element reads the same on both strands: -c finds it twice
    let (_, out, _) = scan(&["--format", "bed6", "-c"], CACTA, SEQ);
    assert_eq!(out, "chr1\t20\t80\tDTC1\t0\t+\nchr1\t20\t80\tDTC2\t0\t-\n");
    let (_, out, _) = scan(
        &[
            "--format",
            "bed6",
            "-c",
            "--dedup",
            "--name-prefix",
            "X",
            "--name-start=7",
        ],
        CACTA,
        SEQ,
    );
    assert_eq!(out, "chr1\t20\t80\tX7\t0\t.\n");
    // -m counts the elements written, after --dedup; the search stops
    // at the limit, so chr2 is not searched on the reverse strand
    let two = format!("{SEQ}{}", SEQ.replace("chr1", "chr2"));
    let (_, out, _) = scan(
        &["--format", "bed6", "-c", "--dedup", "-m", "2"],
        CACTA,
        &two,
    );
    assert_eq!(out, "chr1\t20\t80\tDTC1\t0\t.\nchr2\t20\t80\tDTC2\t0\t+\n");
    // BED12: blocks are the labelled parts, thick part without the TSDs
    let (_, out, _) = scan(&["--format", "bed12"], CACTA, SEQ);
    assert_eq!(
        out,
        "chr1\t20\t80\tDTC1\t0\t+\t23\t77\t0\t4\t3,12,12,3,\t0,3,45,57,\n"
    );
}

/// Reverse complement.
fn rc(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            x => x,
        })
        .collect()
}

/// A tiny parser for the JSON lines written by the program.
fn field<'a>(obj: &'a str, key: &str) -> &'a str {
    let k = format!("\"{key}\":");
    let s = &obj[obj.find(&k).unwrap() + k.len()..];
    if let Some(r) = s.strip_prefix('"') {
        &r[..r.find('"').unwrap()]
    } else {
        &s[..s.find([',', '}']).unwrap()]
    }
}

#[test]
fn jsonl_units_read_back() {
    // alternatives, zero-width units, a label in each branch, both strands
    let pat = "\
^ % start
p1=2...4 ( AC[1,0,0] %@ left
| GG %@ right
) 0...6
~p1 %@ back
";
    let mut x = 7u64;
    let mut seq = String::new();
    for _ in 0..3000 {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seq.push(b"ACGT"[(x >> 60) as usize & 3] as char);
    }
    let input: String = (0..30)
        .map(|k| format!(">s{k}\n{}\n", &seq[k * 100..k * 100 + 100]))
        .collect();
    let (st, out, err) = scan(&["--format", "jsonl", "-c", "-o", "1"], pat, &input);
    assert_eq!(st, Some(0), "{err}");
    let mut seen = std::collections::HashSet::new();
    for l in out.lines() {
        let id: usize = field(l, "seq")[1..].parse().unwrap();
        let rec = &seq[id * 100..id * 100 + 100];
        let rev = field(l, "strand") == "-";
        let units = &l[l.find("\"units\":[").unwrap()..];
        for u in units.split("{\"line\"").skip(1) {
            let u = format!("{{\"line\"{u}");
            let (a, b): (usize, usize) = (
                field(&u, "start").parse().unwrap(),
                field(&u, "end").parse().unwrap(),
            );
            let text = field(&u, "text");
            let fwd = &rec[a - 1..b];
            assert_eq!(text, if rev { rc(fwd) } else { fwd.to_string() }, "{l}");
            let line = field(&u, "line");
            let label = u.contains("\"label\"").then(|| field(&u, "label"));
            seen.insert((line.to_string(), label.map(String::from)));
        }
    }
    // each branch's unit is mapped to its own line and label
    for want in [
        ("1", None),
        ("2", Some("left")),
        ("3", Some("right")),
        ("4", None),
        ("5", Some("back")),
    ] {
        assert!(
            seen.contains(&(want.0.to_string(), want.1.map(String::from))),
            "{want:?} not in {seen:?}"
        );
    }
}

#[test]
fn label_errors_and_plain_output() {
    let bad = "AAA %@ x y\n";
    let (st, _, err) = scan(&["--format", "gff3"], bad, ">a\nAAA\n");
    assert_eq!(st, Some(2));
    assert!(err.contains("pattern line 1"), "{err}");
    let (st, _, err) = scan(&["--format", "gff3"], "% @x\n%@ x\nAAA\n", ">a\nAAA\n");
    assert_eq!(st, Some(2));
    assert!(err.contains("no pattern unit starts"), "{err}");
    // without --format the labels are plain comments, as in 0.2.0
    let (st, out, _) = scan(&[], bad, ">a\nAAA\n");
    assert_eq!((st, out.as_str()), (Some(0), ">a:[1,3]\nAAA \n"));
    let (st, _, err) = scan(&["--dedup"], "AAA\n", ">a\nAAA\n");
    assert_eq!(st, Some(2));
    assert!(err.contains("need --format"), "{err}");
    let (st, _, err) = scan(&["--format", "xml"], "AAA\n", "");
    assert_eq!(st, Some(2));
    assert!(err.contains("gff3, bed6"), "{err}");
}

#[test]
fn no_labels() {
    let (_, out, _) = scan(
        &["--format", "gff3", "--type", "stem_loop"],
        "p1=3...3 2...4 ~p1\n",
        ">s 1\nAACGTTTACGTT\n",
    );
    let body: Vec<&str> = out.lines().skip(3).collect();
    assert_eq!(
        body,
        vec![
            "s\tscan_for_matches\tstem_loop\t2\t11\t.\t+\t.\tID=sfm1;Name=sfm1;Method=structural;Sequence_ontology=SO:0000313",
            "###"
        ]
    );
}
