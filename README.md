# scan_for_matches (Rust port)

[![CI](https://github.com/xavierzheng/scan_for_matches_rs/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/xavierzheng/scan_for_matches_rs/actions/workflows/ci.yml)

A Rust port of **scan_for_matches** by Ross Overbeek (Argonne National
Laboratory): scan nucleotide or protein sequences in FASTA format for
patterns (ranges, hairpins/reverse complements with pairing rules, repeats,
inexact matches, weight matrices, alternatives, length limits, `^`/`$`,
`any()`/`notany()` for proteins).

The pattern language is described in [`README.original`](README.original)
(the original documentation).

## Versions

* **0.0.0** – exact reproduction of the original C program, including its
  bugs, crashes and fixed limits.
* **0.1.0** – the bugs and input problems of the original
  are fixed (see [`CHANGELOG.md`](CHANGELOG.md)).  For all other input the
  results are byte-identical to the original C program: same hits, same
  output format, same exit status, same messages.
* **0.2.0** – faster (gap skipping, `-t N` threads, gzip input); same
  output as 0.1.0.
* **0.3.0** (this version) – output as GFF3, BED or JSON lines
  (`--format`), with the parts of a hit named by `%@` labels in the
  pattern file.  Without `--format` the output is the same as 0.2.0.

## Build

```sh
cargo build --release
# binary: target/release/scan_for_matches
```

No third-party crates are used.

## Usage

Same command line as the original:

```sh
scan_for_matches [-c] [-p] [-n N] [-m N] [-o N] [-i ids_to_ignore] [-t N] pattern_file < fasta_input > hits
```

| option | meaning |
|---|---|
| `-c` | also search the opposite strand (reverse complement); not with `-p` |
| `-p` | protein sequences |
| `-n N` | stop (exit status 1) after N sequences without a hit |
| `-m N` | report at most N hits |
| `-o N` | show overlapping hits (the value is not used) |
| `-i file` | file of sequence ids to skip |
| `-t N` | use N threads: several records at the same time, and long records in pieces; default 1.  The output is the same as with one thread |

The FASTA input on stdin can be plain text or compressed with gzip or
bgzip (for example NCBI `*.fna.gz` files); compression is detected
automatically:

```sh
scan_for_matches -c pat_file < genome.fna.gz
```

### Output formats (`--format`)

| option | meaning |
|---|---|
| `--format F` | `gff3`, `bed6`, `bed12` or `jsonl` instead of the original output |
| `--name-prefix P` | Name of each hit = P + number (default: the `Name=` of the label, else `sfm`) |
| `--name-start N` | first number (default 1); numbers have no fixed width: `DTC1` … `DTC123456` |
| `--type T` | column 3 of the element (default: the `%@element` type, else `sequence_motif`) |
| `--dedup` | with `-c`: a reverse-strand hit with the same span as a forward hit is dropped; the element gets strand `.`.  `-m` counts the elements written |

The output is the same for every `-t N`.  Long options can be written
`--format gff3` or `--format=gff3`.

Labels are `%` comments, so the C program and older versions read the
same pattern file (they ignore the labels):

```
%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC
p1=3...3               %@ target_site_duplication
CACTA[0,0,0] p2=7...7  %@ five_prime_terminal_inverted_repeat
500...15000
~p2 TAGTG[0,0,0]       %@ three_prime_terminal_inverted_repeat
p1                     %@ target_site_duplication
```

* `%@element TYPE key=value ...` (once, anywhere): the element.  Keys:
  `Name` (Name prefix), `Classification`, `Method` (default
  `structural`), `Sequence_ontology`; other keys are copied to column 9.
* `%@ TYPE key=value ...` at the end of a line: all units that start on
  that line make one feature (here the 5' TIR is `CACTA[0,0,0] p2`).
  `role=tsd|tir|other` sets the role; by default the type
  `target_site_duplication` is a TSD and the three
  `*terminal_inverted_repeat` types are TIRs.
* Lines without a label (the spacer above) give no feature.
* Types are free text.  `Sequence_ontology=` is written when the label
  gives it or when the type is one of 21 SO names the program knows
  (TE, TSD, TIR, `repeat_region`, `stem_loop`, `TF_binding_site`,
  `sequence_motif`, `region`, ...).
* A pattern without labels also works: one feature per hit.

GFF3 of one hit (`>chr1:[21,80]`):

```
##gff-version 3
# scan_for_matches 0.3.0 pattern=cacta.pat
##sequence-region chr1 1 95
chr1  scan_for_matches  repeat_region                         21  80  .  +  .  ID=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657
chr1  scan_for_matches  target_site_duplication               21  23  .  +  .  ID=DTC1.lTSD;Parent=DTC1;Name=DTC1;...
chr1  scan_for_matches  CACTA_TIR_transposon                  24  77  .  +  .  ID=DTC1.te;Parent=DTC1;Name=DTC1;...;Sequence_ontology=SO:0002285;TSD=GAT_GAT;TIR=CACTAACGTTGC_GCAACGTTAGTG
chr1  scan_for_matches  five_prime_terminal_inverted_repeat   24  35  .  +  .  ID=DTC1.lTIR;Parent=DTC1.te;Name=DTC1;...
chr1  scan_for_matches  three_prime_terminal_inverted_repeat  66  77  .  +  .  ID=DTC1.rTIR;Parent=DTC1.te;Name=DTC1;...
chr1  scan_for_matches  target_site_duplication               78  80  .  +  .  ID=DTC1.rTSD;Parent=DTC1;Name=DTC1;...
###
```

* Coordinates are 1-based and closed on the forward strand, also for
  `-c` hits (start <= end).
* The top line carries the bare Name as ID: `repeat_region` (TSD to TSD)
  when the pattern has TSD labels, else the element.  The element is the
  hit without its TSDs.  Parts: `.lTSD`/`.rTSD`, `.lTIR`/`.rTIR` (left /
  right on the forward strand), other parts `.1`, `.2`, ...
* `TSD=` and `TIR=` (on the element, as in EDTA): the two TSDs / TIRs,
  forward-strand letters, left_right.
* Every line of a hit has the same `Name`, `Classification` and `Method`.
* Names are unique within one run.  For several patterns, give each run
  its own `--name-prefix`, or `--name-start` to go on numbering.

BED6: `chrom start-1 end Name 0 strand`, one line per hit.  BED12: the
blocks are the labelled parts (the line spans them), the thick part is
the element without TSDs; without labels, one block = the hit.

JSON lines: one object per hit, for Python/pandas
(`pd.read_json(f, lines=True)`):

```
{"seq":"chr1","strand":"+","start":21,"end":80,"name":"DTC1","type":"CACTA_TIR_transposon","units":[{"line":2,"start":21,"end":23,"text":"GAT","label":"target_site_duplication"},...]}
```

`units` are all entries of the hit in pattern order: `line` = line of
the pattern file, `text` = as matched on the strand searched (for a `-`
hit: the reverse complement of the forward letters at start..end); an
empty unit has end = start - 1.

Example (from the original README):

```sh
$ echo 'p1=4...7 3...8 ~p1' > pat_file
$ scan_for_matches pat_file < tests/data/test_dna_input
>tst1:[6,27]
cguaacc ggttaacc gguuacg 
>tst2:[6,27]
CGUAACC GGTTAACC GGUUACG 
```

## Differences from the original C program

Only these (details in `CHANGELOG.md`):

* no fixed limits: any number of pattern units, pattern codes and weights
  that fit on a pattern line, long hairpin stems and many mismatches,
  sequence ids of any length, any number of ids to ignore, sequences up to
  2 147 483 645 characters;
* `-n` and `-m` work (the original always crashed);
* undefined names (`~p3` without `p3`), names that refer to themselves
  (`p1=~p1`), undefined rule sets, negative lengths and negative
  mismatch/insert/delete counts are pattern errors ("failed to parse
  pattern", exit status 1) instead of crashes or reads outside memory;
* `p50` and `r50` work like the other names and rule sets;
* inexact matches with mismatches, inserts and deletes together try every
  choice (the original lost the "delete" choice after a mismatch, so it
  missed some hits);
* inexact words with more inserts than letters do not read past the end of
  the sequence;
* a header line at the end of the input without a newline ends the input
  (the original looped for ever);
* Windows line ends (`\r\n`) in sequences and patterns, and text or blank
  lines before the first `>`, are accepted;
* protein data is matched in upper case (lowercase data did not match);
  the output still shows the data as it is;
* `-c` with `-p` is an error (exit status 2);
* with `-c`, a lowercase `s` is shown as `s` (the original showed `S`);
* bytes 0x80 and above are unknown characters, like `x`;
* (0.2.0) gzip / bgzip compressed FASTA input is accepted;
* (0.2.0) `-t N` searches records in parallel.  A pattern that reads
  a name before it is certainly matched (defined later, or only in one
  branch of `( | )`) uses matches of earlier sequences, so it is always
  searched on one thread;
* (0.2.0) gap lengths that cannot lead to a hit are skipped (see
  Speed); the output does not change.

Unchanged on purpose (see `TODO.md`, group C): a pattern whose longest
match is 0 characters is a pattern error, the value of `-o` is not used, a
`length()` unit can see the last match of a name that is not matched yet.

## Tests

```sh
cargo test --release
```

* `tests/original_suite.rs` – the test suite shipped with the C program
  (`run_tests`, `test_output`).
* `tests/golden.rs` – cases recorded from the C reference program
  (`tests/make_golden.py`, split with `tests/split_golden.py`):
  `golden.tsv` (550 cases where this version equals C) and `fixed.tsv`
  (234 cases changed by the fixes).
* `tests/fixes.rs` – one test per fix, with expected results worked out
  from the pattern language.
* `tests/gzip_input.rs` – gzip and multi-member (bgzip) input, damaged
  input.
* `tests/threads.rs` – `-t 2/4/8` give the same output as `-t 1`
  (options, hit and miss limits, ignore list, damaged input, patterns
  that use earlier sequences); `tests/fuzz_threads.py` does the same with
  random patterns and input (`FUZZ_FORMAT=1`: with random `--format`
  and `--dedup`).
* `tests/formats.rs` – `--format`: GFF3 of a CACTA hit, BED, `--dedup`
  with `-m`, label errors, and JSON lines whose unit coordinates, read
  back from the FASTA on both strands, give the printed text (with
  alternatives `( | )` and zero-width units).
* `tests/fuzz_skip.py OLD NEW N SEED` – compares a build without gap
  skipping (for example 0.1.0) with a new one on random patterns with
  wide ranges followed by words, reverse complements and repeats.
* unit test in `src/engine.rs` – the inexact matcher against a separate
  recursive implementation of its search order (200 000 random cases).

Differential testing against the C program:

```sh
tests/build_reference.sh ../scan_for_matches_original target/reference/scan_for_matches
python3 tests/fuzz_compare.py --compat target/reference/scan_for_matches target/release/scan_for_matches 10000 1
```

With `--compat`, `fuzz_compare.py` generates only input that no fix
touches, so stdout, stderr and exit status must equal the C program's
(40 000 cases: no difference).  Without `--compat` it also generates
over-limit patterns, odd bytes and so on, as used for version 0.0.0.

The reference program is the original built on macOS arm64 with Apple
clang, `cc -std=gnu89 -O2` (the original does not compile with the default
flags of current compilers).

## Speed

Typical patterns run at about the speed of the C program (from 0.8× to
1.3× of its time, depending on the pattern).  Inexact matches with
mismatches, inserts and deletes together can take longer than in the
original, because the choices it skipped are now tried.

(0.2.0) A range followed by an exact word, an exact reverse complement,
a reverse complement with mismatches only (`~p1[1,0,0]`) or an exact
repeat, for example `p1=8...12 50...30000 ~p1`, is much faster: only the
gap lengths where the next unit can match are tried (same order, same
output).  For long gaps a 5-mer position index of the record is used; it
takes about 4 bytes per base (a 74 Mb chromosome: about 300 MB).  A range
followed by anything else (an inexact word, a weight matrix, another
range) is searched as before.

(0.2.0) With `-t N`, records are searched at the same time, and a long
record (2 Mb or more) is cut into pieces of 1 Mb of start positions that
the threads share; the output is the same as with one thread.  Each
thread keeps a copy of the record it searches (about 2 bytes per base);
the index of a record is shared.

22 TIR patterns of B. napus (`tir_scan_patterns`, for example
`p1=8...8 p2=5...7 50...30000 ~p2 p1`), Apple M4, same hits for all
versions; times are the sum over the 22 patterns:

| input | C (original) | 0.1.0 | 0.2.0 |
|---|---|---|---|
| chromosome A1, 1 Mb, 22 patterns, 1 thread | 3162 s | 3963 s | 3.2 s |
| chromosome A1, 31 Mb, 22 patterns, 1 thread | – | – | 107 s |
| genome, 1.0 Gb, 22 patterns, `-t 8` | – | – | 9.6 min (max 5.6 GB) |

Threads, whole genome, one pattern (02_DTA_short5to7, 81 514 hits):

| `-t` | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| seconds | 349 | 184 | 125 | 110 |
| max memory | 0.6 GB | 1.3 GB | 1.7 GB | 4.2 GB |

The M4 has 4 fast and 6 slow cores, and the search is limited by memory
speed, so more than 4 threads gains little on this machine.

## Layout

* `src/main.rs` – port of `scan_for_matches.c` (options, FASTA, output)
* `src/engine.rs` – port of `ggpunit.c` (pattern parser and matcher)
* `src/fmt.rs` – `--format` output and the `%@` labels
* `src/gz.rs` – gzip / bgzip input (system zlib)
* `src/sys.rs` – the few C library calls used (`getopt`, `sscanf`, stdio
  output, signals, `mmap`)
