---
name: sfm-patterns
description: Write, check and explain search patterns for scan_for_matches (Rust port, github.com/xavierzheng/scan_for_matches_rs) - the pattern language for DNA/RNA/protein motifs, hairpins, direct repeats, target site duplications (TSD), terminal inverted repeats (TIR) and transposon structures. Use when the user wants a .pat file, asks how to find a structure in a genome with scan_for_matches, or wants a pattern explained or fixed.
---

# Writing scan_for_matches patterns

You help a user write a pattern file (`.pat`) for `scan_for_matches`.
The program reads FASTA and reports every place where the pattern fits.
Full guide (if you can read files or the web): `PATTERNS.md` in the
repository.  This file has all you need to write correct patterns.

## Workflow

1. **Get the facts (the "slots").**  Ask only for what is missing:

   | slot | example |
   |---|---|
   | what to find | a CACTA (DTC) transposon, a hairpin, a motif |
   | fixed motifs, and where | `CACTA` at the 5' end |
   | TSD length | 3 bp (exact) |
   | TIR (stem) length and mismatches | 7 bp after the motif, 1 mismatch |
   | inside / loop / gap length | 500 to 15000 bp |
   | strands | both (`-c`) or one |
   | DNA or protein | DNA |

2. **Write the pattern** from the units and templates below.  Build it
   left to right, in the 5' → 3' order of the forward strand.
3. **Check it** (see "Checking"): `--explain` (is this what was meant?)
   and `--lint` (traps).  If you can run commands, run them.  If you
   cannot (a chat app), give the user the commands and ask for the output.
4. **Test it** on a tiny FASTA with the structure planted in it, then on
   the same FASTA with one letter broken (the hit must go away).
5. **Give the run command** (see "Running").
6. Explain the pattern to the user unit by unit, in plain words.

Never claim a pattern works without step 3 or 4.  Say which step was
not done.

## The units

Units are separated by spaces.  Upper and lower case are the same.

| unit | meaning |
|---|---|
| `GAGA` | these letters exactly; IUPAC codes allowed (`R`=A/G, `Y`=C/T, `S`=C/G, `W`=A/T, `K`=G/T, `M`=A/C, `B`,`D`,`H`,`V`, `N`=any); `U` = `T` |
| `GAGA[m,d,i]` | with at most `m` mismatches, `d` deletions (a pattern letter missing in the data), `i` insertions (an extra letter in the data) |
| `a...b` | any `a` to `b` letters (a gap); three dots; small number first; fixed length: `8...8` |
| `p1=UNIT` | give the letters that UNIT matched the name `p1` (`p0` to `p50`); `p1=8...8` catches whatever 8 letters are there |
| `p1` | the same letters as `p1` again (direct repeat, TSD) |
| `~p1` | the reverse complement of `p1` (hairpin stem, TIR) |
| `<p1` | `p1` reversed, NOT complemented (rarely wanted) |
| `p1[m,d,i]`, `~p1[m,d,i]` | the same with errors; prefer `[m,0,0]` |
| `r1={au,ua,gc,cg,gu,ug}` then `r1~p1` | reverse complement with your own pair rules (RNA G-U); pair `xy` = `x` in the name, `y` in the data; `u` = `t` here too |
| `length(p1+p2) < 20` | the named parts together are shorter than 20 (strictly less) |
| `^`, `$` | start / end of each FASTA record |
| `{(a,c,g,t),(a,c,g,t),...} > T` | weight matrix: integer scores for A,C,G,T per position; sum must be more than T |
| `%` | comment to the end of the line |

Alternatives: `( GAGA | GCGCA )` matches A or B.  Spaces around `|`.
Exactly two branches; nest for more: `( TAA | ( TAG | TGA ) )`.

Protein (`-p`): letters are amino acids (`N` = asparagine, `X` = any),
`any(HQD)` = one of H, Q, D; `notany(HK)` = one letter not H or K.  No
`~`, no `-c`, no `--strict-n` with `-p`.

## Templates (copy and change the numbers)

| you want | pattern |
|---|---|
| a motif with up to 1 mismatch | `GAGACC[1,0,0]` |
| a hairpin: stem 4-7, loop 3-8 | `p1=4...7 3...8 ~p1` |
| a hairpin with 1 bad pair | `p1=8...8 3...8 ~p1[1,0,0]` |
| an RNA hairpin with G-U pairs | `r1={au,ua,gc,cg,gu,ug} p1=6...6 3...8 r1~p1` |
| a direct repeat, gap 10-50 | `p1=6...6 10...50 p1` |
| a TIR element (no TSD) | `p1=10...10 50...3000 ~p1[1,0,0]` |
| a TIR transposon with TSD | `p1=3...3 CACTA p2=7...7 500...15000 ~p2 TAGTG p1` |
| errors in only part of a TIR (motif exact) | `p1=3...3 CACTG p2=8...8 300...8000 ~p2[1,0,0] CAGTG p1` |
| a TSD of fixed letters (`TA`) | `TA p1=12...12 100...1000 ~p1[2,0,0] TA` |

A whole record (start codon, 30-60 letters, a stop codon):
`^ ATG 30...60 ( TAA | ( TAG | TGA ) ) $`

The TIR transposon, unit by unit:

```
5'  TSD   TIR                                   TIR    TSD  3'
    GAT   CACTA ACGTTGC  ... inside ...  GCAACGT TAGTG  GAT
    p1    CACTA p2           500...15000    ~p2   TAGTG  p1
```

## Rules that cause most mistakes

0. **A TSD is the same letters on both sides**: name the left copy and
   repeat it (`p1=3...3 ... p1`).  A TIR pair is a name and its reverse
   complement (`p2=7...7 ... ~p2`).
1. **Write the other end yourself.**  The program does not guess it.  A
   motif `CACTA` at the 5' end needs its reverse complement `TAGTG` at
   the 3' end (reverse `CACTA` → `ATCAC`, complement → `TAGTG`).  Or name
   it: `p3=CACTA ... ~p3`.
2. **TIR / hairpin = `~p1`**, not `p1` (direct repeat) and not `<p1`
   (reversed only).  `p1`=`GAT` → `p1` matches `GAT`, `~p1` matches
   `ATC`, `<p1` matches `TAG`.
3. **Define a name before you use it, and only once** (also not in both
   branches of `( | )`).
4. **The gap is only the inside**, not the element length.  Element =
   TSD + TIR + inside + TIR + TSD.  `--explain` prints the hit length.
5. **Small number first**: `10...3` matches only length 10.
6. **A range with a name uses the shortest length that works**
   (`p1=4...7` takes 4 if 4 fits).  For a fixed stem use `p1=7...7`.
7. **One pattern per file.**  Lines of a file are joined with spaces.
8. **Mismatches ≥ length match anything**: `ACG[3,0,0]`.
9. **Speed after a wide gap** (a gap of a few hundred is fine either
   way).  Fast: the gap is followed by an exact word, an exact `p1`, an
   exact `~p1`, or `~p1[m,0,0]` (mismatches only) when `p1` has at least
   5 x (m + 1) letters (10 for `[1,0,0]`) or 5 exact letters follow it
   (`~p2[1,0,0] TAGTG`).  Slow: a word with errors (`TAGTG[1,0,0]`),
   `p1[1,0,0]`, any unit with deletions/insertions, or a short
   `~p2[m,0,0]` with few exact letters after it.  `--lint` says "slow"
   for each of these.  So put the mismatches on the `~p2[m,0,0]` right
   after the gap, and keep the words exact.
10. **Hits do not overlap**: after a hit the search goes on after its end.
    `-o 1` shows overlapping hits (forward strand only; large output).
11. **Branch order**: `( A | B )` searches branch A over the whole record
    first.  For every site of both, use two pattern files.
12. **Assembly gaps (N)**: a range can catch `N`; a later `p1` / `<p1`
    then accepts any letter there.  Run genomes with `--strict-n`: each
    caught letter that is not A/C/G/T then counts as 1 mismatch.
    (`~p1` never matches a caught `N`.)  `N` in the data never matches a
    word letter.
13. **Symmetric patterns** (all TIR patterns, hairpins) are found twice
    with `-c`; add `--dedup` (needs `--format`).

## Labels (for GFF3 output)

`%@` comments name the parts; old versions ignore them.

```
%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC
p1=3...3               %@ target_site_duplication
CACTA p2=7...7         %@ five_prime_terminal_inverted_repeat
500...15000
~p2 TAGTG              %@ three_prime_terminal_inverted_repeat
p1                     %@ target_site_duplication
```

`%@ TYPE` at the end of a line labels the units that start on that line.
Lines without a label (the gap) give no feature.  The two TSD and TIR
types above are recognised as TSD / TIR.

## Checking

```sh
scan_for_matches --explain my.pat     # each unit in words, hit length
scan_for_matches --lint -c --dedup --strict-n --format gff3 my.pat
```

Give `--lint` the options of the real run.  Exit status 1 = error (the
pattern does not parse: fix it).  Read every warning; notes are hints.
"failed to parse pattern" at run time = a syntax error: run `--lint`.

Test on a planted sequence:

```sh
printf '>t\nccccGATCACTAACGTTGC%sGCAACGTTAGTGGATcccc\n' "$(printf 'A%.0s' $(seq 600))" > t.fa
scan_for_matches my.pat < t.fa      # one hit, at the planted place
```

Then break one letter of a TIR or TSD and run again: the hit must go.
Use `^ ... $` in a pattern to test on a sequence that is exactly the
structure.

## Running

```sh
scan_for_matches -t 8 -c --dedup --strict-n --format gff3 --name-prefix DTC \
    --input genome.fa.gz --output DTC.gff3 DTC.pat
```

- `-c` both strands: use it for genome scans (an element or motif can
  sit on either strand).  `--dedup` (with `--format`) keeps one copy of
  an element that reads the same on both strands (TIR elements,
  hairpins).
- **`--strict-n` for every genome run** (no false hits at assembly gaps).
- `-t N` threads (same output); input may be plain, gzip or bgzip.
- `--format gff3|bed6|bed12|jsonl`; without it, the original text output
  (`>id:[start,end]` and the matched units).
- Several pattern files in one run (needs `--format`); then
  `scan_for_matches --merge --output all.gff3 run.gff3` for unique Names.
- `-o -c` is a trap: `-o` takes a value; write `-o 1 -c`.

## Answer format

Give: the pattern file (in a code block), one line per unit saying what
it matches, the check commands and what output to expect, and the run
command (for a genome: `-c`, `--strict-n`, and `--dedup --format gff3`
when the pattern is symmetric).  If you ran the checks, show their output.
