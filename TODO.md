# TODO

## Later

- [x] `--lint` and `--explain` for patterns (0.4.0).
- [x] Merge tool: several patterns/runs → genome-wide unique Names
      (`--merge`, 0.4.0).
- [ ] `--format`: a hit of length 0 is written only in JSON lines (GFF3
      and BED have no place for it).
- [x] Several patterns in one run (read the genome once) (0.4.0).
- [x] Converter JASPAR matrix → integer weight unit with a p-value
      cutoff (`tools/jaspar2sfm.py`, 0.4.0).
- [ ] A pattern-writing skill for small models (JSON slot spec → pattern)
      and `--synth` test sequences.

## Platforms

- [x] Linux x86_64: decide if Linux needs its own reference.
      Decided 2026-10-09: no.  The results below (C program built on
      Linux, CI on Linux) are the same as on macOS, except the 2
      undefined-behaviour cases of 0.0.0.
      CI result (2026-10-01, GitHub `ubuntu-latest`, AMD EPYC 9V74 with
      AVX512; there is no SIMD code (0.2.0 compares 8 bases with plain
      64-bit integers), so CPU features do not matter):
      - builds and runs; the original test suite passes;
      - 786 of 788 golden cases (recorded on macOS arm64) match;
      - the 2 differences are undefined-behaviour cases:
        `p50` crash (Linux stdio writes 4096 bytes before the SIGSEGV,
        macOS none) and data byte 0xA5 (reads a byte of a heap pointer;
        Linux heap addresses differ);
      - version 0.1.0 (no undefined behaviour left): all tests pass on
        Linux x86_64 in CI; version 0.2.0 too (CI 2026-10-03);
      - 2026-10-08, HPC rossini: compared with the C program built on
        Linux (module scan_for_matches/20260930, gcc 9.5): `fuzz_compare
        --compat` 2 x 10 000 cases, all the same (version 0.4.0).
        `fuzz_skip` (0.1.0 built on the node vs 0.4.0): plain 800 and
        `--chain` 1000, no DIFF, no NEW-timeout.
      - 2026-10-08, HPC Intel Xeon Gold 6238R (Rust 1.98): `cargo test
        --release` passes; 22 TIR patterns on the B. napus genome give
        the same output (md5) as 0.3.0 with the new threading code.

## A. Bugs: crash, hang or silently wrong result (fixed in 0.1.0)

- [x] `-n` and `-m` always crash (SIGSEGV): the `getopt` string
      `"pcnmo:i:"` lacks `:` after `n` and `m`, so `optarg` is NULL.
- [x] A header line at end of file without a newline makes the program
      loop forever (`while (getc(stdin) != '\n')`).
- [x] Exactly 100 pattern units: the last, failed parse attempt uses slot
      100 (past `pu_s[100]`) and overwrites pattern codes, so hits are
      silently wrong.
- [x] More than 100 pattern units, more than 600 pattern code bytes, very
      large weight lists: writes past `pu_s`, `cv`, `iv` (wrong pattern,
      abort or crash).
- [x] Hit with more than 100 units: `revhits[100]` overflows (abort).
- [x] Reverse complement with mismatch/insert/delete longer than 100:
      `result[100]` in `loose_match` overflows (abort).
- [x] More than 101 open choices in `loose_match`: `stack[100]` overflows
      (abort).
- [x] Sequence ids of 1000 or more characters overflow `id[1000]` (random
      bytes in the output, or crash).
- [x] More than 20 000 ids in the `-i` file overflow `ignore[20000]`
      (crash).
- [x] Sequences longer than 250 000 000 characters overflow the data
      buffer (crash).
- [x] Undefined name (`~p3` without `p3`) or undefined rule set
      (`r7~p1` without `r7={...}`): crash instead of a pattern error.
- [x] `p1=~p1` (a name that refers to itself): endless recursion, crash.
- [x] `p50` and `r50` are accepted (`<= MAX_NAMES`) but the arrays have
      50 entries; they overwrite `past_last` and `start_srch`.
- [x] Negative ranges (`-3...5`) and similarity patterns with more inserts
      than letters read outside the sequence buffer.
- [x] Inexact match with both inserts and deletes: after `Pop` the
      "delete" choice is written to a released stack entry, so this search
      path is lost and some real hits are missed.

## B. Input problems (fixed in 0.1.0)

- [x] Windows line ends: `\r` is kept inside the sequence, so matches
      across lines fail (only spaces and `\n` are removed).
- [x] A blank line (or any text) before the first `>` gives no output and
      no error.
- [x] Protein mode: pattern letters are changed to uppercase but the data
      is not, so lowercase protein data never matches.
- [x] `-c` together with `-p` runs the DNA reverse complement on protein
      data; it should be rejected.
- [x] `compl()` maps lowercase `s` to uppercase `S` in `-c` output.
- [x] Bytes 0x80 and above in the data or pattern index `punit_to_code`
      with a negative (signed `char`) subscript and read other globals.

## C. Odd behaviour (kept on purpose, decided 2026-10-09)

Not to do: these are the behaviour of the original C program, and the
output stays the same as the C program.  A record, not a task list.
`--lint` warns about several of them.

- A pattern whose maximum match length is 0 (for example `^` alone or
      `0...0`) is reported as "failed to parse pattern".
- The value of `-o N` is not used; the usage text is not correct
      (`-n`/`-m` shown with values, `-o` without).
      (2026-10-08: `--help` now shows the correct options; the short
      usage text printed on an option error is still the old one.)
- `length(...)` can use the last match length of a name from an
      earlier attempt or an earlier sequence, when that name has not been
      matched yet in the current attempt.
- In an alternative `( | )`, the second alternative is not retried at
      later positions (dead branch for `alt == 2` in BACKTRACK); scanning
      relies on the first unit inside the alternative.
- `-o` (overlapping hits) works on the forward strand only: with
      `-c`, hits on the reverse strand never overlap (the C code always
      uses `cont_match` there).  The original README does not say this.
- `-n`: the original README says `-n 10` limits the output to 10 hits,
      but the C code (and its usage text) uses `-n N` as "stop after N
      sequences without a hit" and `-m N` as the hit limit.  Version 0.1.0
      follows the code.  Decide which meaning to keep, and correct the
      README or the code.
- The inexact matcher always takes a matching character and never
      tries to skip it with an insert/delete, so some valid alignments
      within the mismatch/insert/delete limits are not found (for example
      when an early match forces a later failure).  The original works the
      same way; a complete search would change results and speed.
- A reversed range `a...b` with `a > b` is accepted without a
      warning and matches only the length `a`.
- `-o` takes a value, so `-o -c pattern` uses `-c` as that value:
      the reverse strand is not searched and there is no error.
- A pattern that reads the same on both strands (all TIR patterns)
      reports each element twice with `-c` (`[a,b]` and `[b,a]`).
- The pattern file is joined into one line (newlines become spaces)
      and cut at 31 999 bytes without a warning.
- [x] A name that caught `N` from the data (a range takes any letter)
      works as a wildcard when it is used again as `p1` or `<p1` (not
      `~p1`: it never matches a caught `N`).  A TSD made of an assembly
      gap matches anything (`p1=8...8 20...50 p1` hits `NNNNNNNNNN...`
      and random letters).  Same in the C program; kept as the default.
      `--strict-n` (0.4.1) makes each such letter 1 mismatch.
- A weight-matrix unit scores IUPAC codes and `N` in the data as
      averages of the matching weights, while all other units never
      match `N`.

## `--strict-n` and release 0.4.1 (decided 2026-10-08, done 2026-10-09)

Problem (TODO C, PATTERNS.md mistake 12): a name that caught `N` (or
another IUPAC code) from the data works as a wildcard when it is used
again (`p1`, `~p1`, `<p1`), so TSDs/TIRs made of assembly gaps give
false hits.  Same in the C program.

- [x] Option `--strict-n` (OFF by default: without it the output stays
      byte-identical to the C program and 0.4.0).
      Rule: a letter caught by a name that is not A, C, G, T (U) counts
      as ONE MISMATCH at that place when the name is used again.  No
      other threshold: the mismatches the pattern allows decide.

      | pattern | letters the name caught | result |
      |---|---|---|
      | `p1=3...3 ... p1` (exact TSD) | `GNT` (1 N) | no match |
      | `p2[1,0,0]` (1 mismatch) | 1 N, the rest right | match (uses the 1 mismatch) |
      | `p2[1,0,0]` | 2 N | no match |
      | `~p2`, `~p2[1,0,0]` | any N | no match (as before: `~pN` never matched a caught N) |
      | any pattern | N inside a gap (`500...15000`) | no effect |
      | R, Y, ... (other IUPAC codes) | same as N: 1 code = 1 mismatch | |

      In the engine, not a filter after the search: a hit that is found
      and then dropped makes the search jump past its end, so a real
      element next to it could be lost.  Also the gap-skipping paths
      (masks from the named text, k-mer index) must follow the rule.
- [x] Tests: without the option, byte-identical (cargo test, fuzz_compare
      vs C, fuzz_threads, fuzz_skip, 22-pattern md5); with it, fuzz
      against a small separate reference implementation; maize and
      B. napus: how many hits go away, check that they are next to N.
- [x] `--lint`: a note for patterns that use a name again: "add
      --strict-n for genomes with gaps (N)".
- [x] Docs as clear as the table above (the user asked for this level):
      README options, PATTERNS.md (mistake 12 + the table), TIR examples
      and Quick start use `--strict-n`, CHANGELOG.
- [x] Release 0.4.1 (memory work, PATTERNS.md, license, `--explain`
      fix, `--strict-n`): version, CHANGELOG "Unreleased" → 0.4.1,
      README Versions (newest first), tag after CI passes.

## Documentation: pattern guide (planned 2026-10-08)

- [x] `PATTERNS.md` (repo root; `docs/` is not pushed): the pattern
      language for biology undergraduates.  Table of contents with links
      to each section; text diagrams (TSD, TIR, gap, hairpin); every
      example run with the binary and checked with `--explain`; the traps
      that `--lint` reports.  Written by a Sonnet agent.
- [x] (round 1: 10/10, 9 sure; round 2 after fixes: 10/10, 10 sure;
      grader and kit outside the repo) Check with a Haiku agent that gets only `PATTERNS.md`: 10 tasks
      ("write a pattern for ..."); its patterns are graded with the binary
      (`--lint` without errors, the planted sites in a test sequence are
      found); wrong answers show unclear parts; 1-2 rounds of fixes.  Also
      a list of words an undergraduate does not know.
- [x] Link `PATTERNS.md` from the README (intro and Quick start).

 (planned 2026-10-08, in this order)

Measured on maize (2.3 Gb, 10 chromosomes of about 300 Mb, bgzip),
`-t 20`, BED6, version 0.4.0 (`../bench_zma/`):

| pattern | without `-o` | with `-o` | hits (without → with `-o`) |
|---|---|---|---|
| 00_DTC | 16 s, 13 GB | 16 s, 13 GB | 2 884 → 3 021 |
| 08_DTH_seed12 | 55 s, 15 GB | 57 s, 16 GB | 101 618 → 1 303 532 |
| 02_DTA_short5to7 | 140 s, 16 GB | 222 s, 97 GB | 163 883 → 40 209 313 |

- [x] 1. `--format` with `-o`: memory (done: 97 GB → 3.9 GB with 2).
      Workers sent each hit as a raw
      record with the text of EVERY unit, also the 50...30000 gap
      (`fmt::encode`), and finished pieces wait in the writer while an
      earlier piece is slow (about 3 x threads pieces).  Fix: send only
      the text the format uses (BED6/BED12: none; GFF3: the TSD/TIR
      parts for `TSD=`/`TIR=`; JSON lines: all), and limit the pieces
      that wait.  The writer stays the only thread that writes (order,
      Names, `-m`, `--dedup` unchanged).  Temp files for each piece:
      only if this is not enough.  Output must not change.
- [x] 2. One copy of each strand of a long record, shared by all threads
      (done: maize 00_DTC 13.3 GB → 2.5 GB).
      Now every thread that searches pieces of a record holds its own
      copy (`PieceEngine.data` + the engine's coded sequence, about 2
      bytes per base: 300 Mb x 2 x 20 threads = 12 GB).  The engine reads
      the coded sequence only (writes into it happen only in the
      emulated C undefined-behaviour cases: then use a private copy).
      One thread per record is not the answer: few, unequal chromosomes
      leave cores idle, and it saves little memory.  Output must not
      change.

## Speed: known limits (0.2.0; a record, no plan to change, 2026-10-09)

`--lint` gives a note for the slow wide gaps.

- Gap skipping works only when the range is followed by an exact
      word, an exact `~pN`, `~pN[m,0,0]` or an exact `pN`.  A range
      followed by an inexact word, `pN` with errors, `~pN` with inserts
      or deletes, a weight matrix or another range is searched length by
      length (slow for wide ranges such as `50...30000`).
- With `-t N` on many long records and a fast pattern, splitting
      long records into pieces is slower than one thread per record
      (each thread copied each long record, until the shared strands;
      not measured again; whole genome, 00_DTC `-c`:
      5.8 s without pieces, 9.3 s with them).  Sharing one copy needs a
      change of the engine memory layout (see "Memory: next work", 2).
- [x] Without `-o`, with `-t N`, the writer thread searched a piece
      again by itself until its chain of non-overlapping hits met the
      worker's chain (`join_piece`); with long, dense hits that took
      about 300 kb (1 Mb pieces).  Fixed by an adaptive warm-up before
      each piece (02_DTA_short5to7 308 s → 106 s, 08_DTH_Tourist_seed8
      160 s → 35 s, B. napus, `-t 20`).
- Each engine maps 4 GB of address space (`CD_CAP`) and never unmaps
      it (fine for a command-line run).
