# TODO

## Platforms

- [ ] Linux x86_64: decide if Linux needs its own reference.
      CI result (2026-10-01, GitHub `ubuntu-latest`, AMD EPYC 9V74 with
      AVX512; there is no SIMD code, so CPU features do not matter):
      - builds and runs; the original test suite passes;
      - 786 of 788 golden cases (recorded on macOS arm64) match;
      - the 2 differences are undefined-behaviour cases:
        `p50` crash (Linux stdio writes 4096 bytes before the SIGSEGV,
        macOS none) and data byte 0xA5 (reads a byte of a heap pointer;
        Linux heap addresses differ);
      - version 0.1.0 (no undefined behaviour left): all tests pass on
        Linux x86_64 in CI;
      - not yet compared with the C program built on Linux (its memory
        layout, and so its undefined-behaviour results, differ from the
        macOS build).  Normal input is expected to match.

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

## C. Odd behaviour (can be kept)

- [ ] A pattern whose maximum match length is 0 (for example `^` alone or
      `0...0`) is reported as "failed to parse pattern".
- [ ] The value of `-o N` is not used; the usage text is not correct
      (`-n`/`-m` shown with values, `-o` without).
- [ ] `length(...)` can use the last match length of a name from an
      earlier attempt or an earlier sequence, when that name has not been
      matched yet in the current attempt.
- [ ] In an alternative `( | )`, the second alternative is not retried at
      later positions (dead branch for `alt == 2` in BACKTRACK); scanning
      relies on the first unit inside the alternative.
