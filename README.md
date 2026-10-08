# scan_for_matches (Rust port)

[![CI](https://github.com/xavierzheng/scan_for_matches_rs/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/xavierzheng/scan_for_matches_rs/actions/workflows/ci.yml)

A Rust port of **scan_for_matches** by Ross Overbeek (Argonne National
Laboratory): scan nucleotide or protein sequences in FASTA format for
patterns (ranges, hairpins/reverse complements with pairing rules, repeats,
inexact matches, weight matrices, alternatives, length limits, `^`/`$`,
`any()`/`notany()` for proteins).

The pattern language is described in [`README.original`](README.original)
(the original documentation).

## Quick start

**Why use it**

* Same pattern language and same results as the original C program (tested
  against it), without its crashes and size limits.
* Fast: gap skipping (up to 1000× on TIR patterns), `-t N` threads, gzip /
  bgzip input.  22 TIR patterns on a 1.0 Gb genome: about 11 min on 20
  threads.
* Output as GFF3, BED or JSON lines, with named parts (TSD, TIR, ...) and
  unique Names; several patterns in one run; `--merge` for one
  genome-wide file.
* `--lint` and `--explain` check and describe a pattern before a long run.
* 22 ready-made TIR transposon patterns in
  [`tir_scan_patterns/`](tir_scan_patterns/README.md).

**Install**

* macOS (Apple silicon): download `scan_for_matches-<version>-macos-arm64.tar.gz`
  from [Releases](https://github.com/xavierzheng/scan_for_matches_rs/releases).
* Linux or other systems: install Rust (<https://rustup.rs>), then

  ```sh
  git clone https://github.com/xavierzheng/scan_for_matches_rs.git
  cd scan_for_matches_rs
  cargo build --release      # binary: target/release/scan_for_matches
  ```

  No other libraries are needed (the system zlib is used).

**Use**

```sh
# a hairpin: a 4-7 bp stem, a 3-8 bp loop, the reverse complement of the stem
echo 'p1=4...7 3...8 ~p1' > hairpin.pat
scan_for_matches hairpin.pat < seqs.fa > hits.txt

# check a pattern first (no input is read)
scan_for_matches --explain hairpin.pat
scan_for_matches --lint -c hairpin.pat

# TIR transposons in a genome: GFF3, 8 threads, all 22 patterns in one run,
# then one file with genome-wide unique Names
scan_for_matches -t 8 -c --dedup --format gff3 \
    --input genome.fa.gz --output tir_all.gff3 tir_scan_patterns/labelled/*.pat
scan_for_matches --merge --output tir.gff3 tir_all.gff3

scan_for_matches --help      # all options
```

## Versions

Newest first:

* **0.4.0** (this version) – several patterns in one run (the input is
  read once), `--lint` and `--explain` for patterns, `--merge` for the
  GFF3 of several runs (genome-wide unique Names), and
  `tools/jaspar2sfm.py` (JASPAR matrix → weight unit).  Without the new
  options the output is the same as 0.3.0 (only the version in the GFF3
  header changes; extra file arguments, which 0.3.0 ignored, are now
  pattern files).
* **0.3.0** – output as GFF3, BED or JSON lines
  (`--format`), with the parts of a hit named by `%@` labels in the
  pattern file; `--help`, `--input`, `--output`; faster threads (bgzip
  input, k-mer index, long dense hits); the 22 TIR patterns in
  [`tir_scan_patterns/`](tir_scan_patterns/README.md).  Without the new
  options the output is the same as 0.2.0.
* **0.2.0** – faster (gap skipping, `-t N` threads, gzip input); same
  output as 0.1.0.
* **0.1.0** – the bugs and input problems of the original
  are fixed (see [`CHANGELOG.md`](CHANGELOG.md)).  For all other input the
  results are byte-identical to the original C program: same hits, same
  output format, same exit status, same messages.
* **0.0.0** – exact reproduction of the original C program, including its
  bugs, crashes and fixed limits.

## Build

```sh
cargo build --release
# binary: target/release/scan_for_matches
```

No third-party crates are used.

## Usage

Same command line as the original, or with input and output files:

```sh
scan_for_matches [-c] [-p] [-n N] [-m N] [-o N] [-i ids_to_ignore] [-t N] pattern_file < fasta_input > hits
scan_for_matches [options] --input fasta_input --output hits pattern_file
scan_for_matches --format F [options] pattern_file pattern_file ... < fasta_input > hits
scan_for_matches --lint [-p] [-c] [--dedup] pattern_file ...
scan_for_matches --explain [-p] pattern_file ...
scan_for_matches --merge [merge options] run1.gff3 run2.gff3 ... > all.gff3
```

`scan_for_matches --help` shows all options.

| option | meaning |
|---|---|
| `-c` | also search the opposite strand (reverse complement); not with `-p` |
| `-p` | protein sequences |
| `-n N` | stop (exit status 1) after N sequences without a hit |
| `-m N` | report at most N hits |
| `-o N` | show overlapping hits (the value is not used) |
| `-i file` | file of sequence ids to skip |
| `-t N` | use N threads: several records at the same time, and long records in pieces; default 1.  The output is the same as with one thread |
| `--input FILE` | read the FASTA input from FILE instead of stdin (`-`: stdin) |
| `--output FILE` | write the hits to FILE instead of stdout (`-`: stdout) |
| `-h`, `--help` | show the options and exit |

`-i` and `-o` keep the meaning they have in the original program (ids
to skip, overlapping hits); input and output files are `--input` and
`--output` only.

The FASTA input on stdin can be plain text or compressed with gzip or
bgzip (for example NCBI `*.fna.gz` files); compression is detected
automatically.  With `-t N`, the blocks of a bgzip file are
decompressed by N threads (plain gzip uses one thread):

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
# scan_for_matches 0.4.0 pattern=cacta.pat
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
* Names are unique within one run (also with several pattern files in
  one run: one counter for each Name prefix).  To join separate runs,
  use `--merge` (below), or give each run its own `--name-prefix`.

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

### Several patterns in one run

```sh
scan_for_matches -t 20 -c --dedup --format gff3 \
    --input genome.fna.gz --output tir.gff3 tir_scan_patterns/labelled/*.pat
```

* Needs `--format` (the original output cannot show which pattern made a
  hit).  The input is read and decompressed once; each record is
  searched by each pattern in turn, and its hits are written pattern by
  pattern, in the order of the files.
* The hits of each pattern are the same as in a run with that pattern
  alone; `--dedup` works for each pattern.  Hits found by two patterns
  are both written (use `--merge` to remove them).
* Names: the prefix of each pattern (its `Name=` label, or
  `--name-prefix` for all), one counter for each prefix.
* GFF3 header: `pattern=` for each file; JSON lines get a `"pattern"`
  field (only with several files).
* `-m` counts the hits of all patterns; `-n` counts records where no
  pattern has a hit.
* B. napus genome, 22 TIR patterns, `-t 20`: 648 s in one run (max
  6.4 GB), 762 s as 22 runs (max 4.9 GB); same hits.

### Checking patterns: `--lint` and `--explain`

No input is read.  Several pattern files can be given.

* `--explain` describes each unit in plain words: source line, text,
  meaning (mismatches, IUPAC codes, names, gaps, labels), the hit length
  range, and how the search goes.
* `--lint` reports traps as `file:line: error|warning|note: text` and a
  count.  Exit status 1 when there is an error, else 0.  Give it the
  options of the real run (`-p`, `-c`, `--dedup`, `-o`), as some checks
  depend on them.
  - errors: the pattern does not parse (where, and why), a name defined
    twice, an undefined name or rule set, a name that refers to itself,
    a longest match of 0 letters, bad `%@` labels;
  - warnings: a reversed range `10...3` (only length 10 is matched), a
    name used before it is defined or defined in one branch of `( | )`,
    text cut at 31 999 bytes, mismatches >= word length, `-o -c` (`-c`
    becomes the value of `-o`), `-c` without `--dedup` for a pattern that
    reads the same on both strands (every element found twice);
  - notes: a wide gap that cannot skip gap lengths (slow), `<pN` (reverse,
    not reverse complement), `pN=a...b` uses the shortest length that
    works, alternative order, inserts and deletes together, weight units
    score `N` as an average, unused names, patterns that do not split for
    threads.

```
$ scan_for_matches --lint -c tir_scan_patterns/labelled/00_DTC_published_2025.pat
tir_scan_patterns/labelled/00_DTC_published_2025.pat: warning: the pattern reads the same on both strands: with -c, every element is found twice; use --dedup with --format
0 errors, 1 warnings, 0 notes
```

### Merging several runs: `--merge`

```
scan_for_matches --merge [--overlap F] [--name-prefix P] [--name-start N] \
                 [--output FILE] FILE.gff3 [FILE2.gff3 ...]
```

`--merge` joins the GFF3 files of several patterns or runs into one file with
genome-wide unique Names. The files are in priority order (`-` is stdin; gzip
and bgzip are detected). Any GFF3 with `ID`/`Parent` works.

- An *element* is a top feature (no `Parent`) and all its descendants. IDs are
  local to each input file, so two runs can both have `DTC1`.
- Two elements are the *same* when seqid, start and end of the top feature are
  equal (strand is ignored). With `--overlap F` (0 < F <= 1) they are also the
  same when they overlap by at least F of the length of *each* one.
- The element of the earlier file (then the earlier line) is kept. The kept top
  feature gets `Merged=<file>:<old Name>` for each dropped element.
- Output order: seqid in natural order (`chr2` before `chr10`), then start
  ascending, end descending, then input order.
- New Names are prefix + counter (`DTC1`, `DTC2`, ...), counted in output order.
  The prefix is `--name-prefix`, or the old Name without trailing digits
  (`DTH_98` gives `DTH_`; `sfm` if nothing is left). `--name-start` sets the
  first number (default 1). `ID`, `Parent` and `Name` are renamed in all lines
  of an element; all other columns and attributes are copied as they are.
- An ID must be unique in one file (lines of one feature with the same
  ID and Parent are allowed).  To merge several runs, give the files one
  by one; files joined with `cat` repeat IDs and give an error.
- A summary goes to stderr: `merge: N elements read, N duplicates removed, N written`.
- Exit status: 0 ok, 1 bad input (`merge: FILE:LINE: message`), 2 bad options.

```
scan_for_matches -c --dedup --format gff3 00_DTC_published_2025.pat < genome.fa > a.gff3
scan_for_matches -c --dedup --format gff3 09_DTC_CACTA_seed10.pat   < genome.fa > b.gff3
scan_for_matches --merge --overlap 0.9 a.gff3 b.gff3 --output all.gff3
```

The output of one run with several pattern files can be merged too (one
input file): elements found by two patterns are then removed.

### JASPAR matrix to weight unit (`tools/jaspar2sfm.py`)

`tools/jaspar2sfm.py` turns a JASPAR position frequency matrix into a weight
unit `{(a,c,g,t),...} > CUTOFF` for a pattern file.  Python 3 standard library only.

```
python3 -I tools/jaspar2sfm.py MA0549.1.jaspar --pvalue 1e-4 > bzr2.pat
python3 -I tools/jaspar2sfm.py MA0549.1.jaspar --revcomp >> bzr2.pat   # other strand
```

* Input: JASPAR text (`>ID NAME` and four rows `A [ ... ]`); `-` reads stdin.
  Several matrices in one file give one block each; `--name ID_OR_NAME` keeps one.
* Weights: `round(scale * log2(freq / background))`, `--scale 100`; the
  frequency uses `--pseudocount 0.25` (spread by the background);
  `--background A,C,G,T` defaults to equal.
* Cutoff: the smallest integer T with `P(score > T) <= --pvalue` under the
  background, by an exact DP over the integer scores (the engine tests `score > T`).
  The `%` comment lines show the exact p-value of T and the maximum score.
* `--revcomp` converts the matrix of the opposite strand (reverse columns, A<->T, C<->G).
  Use it instead of `-c` when the unit is part of a longer pattern.
* Limits: no score is printed, the first window above T is reported; N and IUPAC
  letters score as averages of the weights; the pattern file is cut at 31 999
  characters (a warning goes to stderr above 30 000).  Exit status 1 on bad input.
* Check: MA0549.1 at p = 1e-4 gave 101 random hits per 1 Mb (100 expected).
  Tests: `python3 -I -m unittest discover -s tools`.

### Example

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
* `tests/multi.rs` – several pattern files: the hits of each pattern
  equal the single runs, for `-t 1/4` and many piece sizes; Names,
  header, `-m`, errors.
* `tests/lint_cli.rs` and unit tests in `src/lint.rs` – `--lint` /
  `--explain`; the lint parser accepts exactly the patterns the engine
  accepts, with the same units (about 4 900 patterns, 4 000 of them
  random).
* `tests/merge.rs` and unit tests in `src/merge.rs` – `--merge`.
* `tools/test_jaspar2sfm.py` – `python3 -I -m unittest discover -s tools`.
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
flags of current compilers).  On Linux x86_64 (HPC, the C program built
with gcc 9.5), version 0.4.0 gives the same result in 20 000 `--compat`
cases.  `fuzz_skip.py` there uses version 0.1.0 (no gap skipping) built
from its git tag as the old binary.

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

22 TIR patterns of B. napus ([`tir_scan_patterns/`](tir_scan_patterns/README.md), for example
`p1=8...8 p2=5...7 50...30000 ~p2 p1`), whole genome (1.0 Gb, bgzip),
`-c --dedup --format gff3 -t 20`, HPC (Intel Xeon Gold 6238R), same hits
for all versions (newest first):

| version | time | max memory |
|---|---|---|
| 0.4.0, the 22 pattern files in one run | 648 s | 6.4 GB |
| 0.3.0, 22 runs (sum) | 762 s | 4.9 GB |
| 0.3.0 without the warm-up and bgzip threads, 22 runs | 1546 s | 10 GB |

Apple M4, same hits for all versions; times are the sum over the 22
patterns:

| input | 0.2.0 | 0.1.0 | C (original) |
|---|---|---|---|
| chromosome A1, 1 Mb, 22 patterns, 1 thread | 3.2 s | 3963 s | 3162 s |
| chromosome A1, 31 Mb, 22 patterns, 1 thread | 107 s | – | – |
| genome, 1.0 Gb, 22 patterns, `-t 8` | 9.6 min (max 5.6 GB) | – | – |

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
* `src/lint.rs` – `--lint` and `--explain`
* `src/merge.rs` – `--merge`
* `tools/jaspar2sfm.py` – JASPAR matrix → weight unit
* `src/gz.rs` – gzip / bgzip input (system zlib; bgzip blocks in threads)
* `src/sys.rs` – the few C library calls used (`getopt`, `sscanf`, stdio
  output, signals, `mmap`)
