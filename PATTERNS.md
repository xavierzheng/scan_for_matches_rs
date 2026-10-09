# Pattern guide for scan_for_matches

This guide teaches the pattern language of `scan_for_matches`. It is for biology students who have never used the tool. You need no programming skills. You need a terminal and a small text editor.

All commands and outputs in this guide are real. They were run with version 0.4.1 of the Rust port. In the examples, `scan_for_matches` is the name of the program.

For the full original text (1990s), read [README.original](README.original). For all options and output formats, read [README.md](README.md).

## Table of contents

1. [Patterns you will use most](#patterns-you-will-use-most)
2. [What is a pattern](#what-is-a-pattern)
3. [Your first pattern in 2 minutes](#your-first-pattern-in-2-minutes)
4. [How the search works](#how-the-search-works)
5. [Letters and words](#letters-and-words)
6. [Mismatches, deletions and insertions](#mismatches-deletions-and-insertions)
7. [Ranges](#ranges)
8. [Names and repeats](#names-and-repeats)
9. [Reverse complement](#reverse-complement)
10. [Reverse](#reverse)
11. [Alternatives](#alternatives)
12. [Length limit](#length-limit)
13. [Start and end of the sequence](#start-and-end-of-the-sequence)
14. [Weight matrices](#weight-matrices)
15. [Protein patterns](#protein-patterns)
16. [Comments and pattern files](#comments-and-pattern-files)
17. [Diagrams](#diagrams)
18. [Worked example CACTA transposon](#worked-example-cacta-transposon)
19. [Common mistakes](#common-mistakes)
20. [Options used in this guide](#options-used-in-this-guide)
21. [Checking a pattern](#checking-a-pattern)
22. [Glossary](#glossary)

## Patterns you will use most

Copy one of these and change the numbers.

| you want | pattern | section |
|---|---|---|
| an exact motif, with IUPAC letters | `TATAWAWR` | [Letters and words](#letters-and-words) |
| a motif with up to 1 mismatch | `GAGACC[1,0,0]` | [Mismatches, deletions and insertions](#mismatches-deletions-and-insertions) |
| any letters, exactly 8 (a fixed length is written `8...8`) | `8...8` | [Ranges](#ranges) |
| a hairpin: stem 4 to 7, loop 3 to 8 | `p1=4...7 3...8 ~p1` | [Reverse complement](#reverse-complement) |
| a hairpin with 1 mismatched pair | `p1=8...8 3...8 ~p1[1,0,0]` | [Reverse complement](#reverse-complement) |
| a direct repeat with a gap of 10 to 50 | `p1=6...6 10...50 p1` | [Names and repeats](#names-and-repeats) |
| a TIR transposon with TSD | `p1=3...3 CACTA p2=7...7 500...15000 ~p2 TAGTG p1` | [Worked example CACTA transposon](#worked-example-cacta-transposon) |
| a whole record | `^ ATG 30...60 TAA $` | [Start and end of the sequence](#start-and-end-of-the-sequence) |

Each pattern goes in its own file. Check it first with `--explain` and `--lint` (see [Checking a pattern](#checking-a-pattern)).

## What is a pattern

A pattern is a short text that describes a piece of DNA, RNA or protein. The program reads a FASTA file and finds every piece that fits the pattern. Such a piece is a "hit".

A pattern is a list of "units". You write the units one after the other, with spaces between them. Each unit matches one part of the sequence. The parts follow each other from left to right.

For example, `GAGA 3...8 TTT` has three units. It means: the letters GAGA, then 3 to 8 any letters, then the letters TTT.

You can ask for more than fixed letters. You can ask for "the same piece again" (a repeat). You can ask for "the reverse complement of a piece seen before" (a hairpin). This is why the tool finds hairpins, transposon ends and target site duplications.

A pattern lives in a text file. A FASTA file holds the sequences to search.

## Your first pattern in 2 minutes

We search for a hairpin. A hairpin is a stem, a loop, and the same stem again, but reversed and complemented. (The terms are in the [glossary](#glossary).)

Make the pattern file:

```
$ echo 'p1=4...7 3...8 ~p1' > hairpin.pat
```

Make a small FASTA file `first.fa`:

```
>seq1
GGAGGCGTAACCGGTTAACCGGTTACGAGG
>seq2
AAGGAAGGAAGGAAGG
```

Run the program. The pattern file is the argument. The FASTA file goes in with `<`.

```
$ scan_for_matches hairpin.pat < first.fa
>seq1:[6,27]
CGTAACC GGTTAACC GGTTACG 
```

How to read the output:

* `>seq1:[6,27]` is the hit. It is in `seq1`, from letter 6 to letter 27. The first letter is number 1.
* The next line shows the matched letters, one block for each unit.
* `seq2` has no hit, so the program prints nothing for it.

The pattern has three units:

| unit | meaning | matched here |
|---|---|---|
| `p1=4...7` | take 4 to 7 any letters and call them `p1` | `CGTAACC` |
| `3...8` | then take 3 to 8 any letters | `GGTTAACC` |
| `~p1` | then the reverse complement of `p1` | `GGTTACG` |

Check: the reverse complement of `CGTAACC` is `GGTTACG`. Reverse it: `CCAATGC`. Then swap A with T and C with G: `GGTTACG`.

Ask the program to explain a pattern with `--explain`. It reads no FASTA file.

```
$ scan_for_matches --explain hairpin.pat
pattern file: hairpin.pat (DNA)
joined pattern: p1=4...7 3...8 ~p1
unit  line  text      meaning
1     1     p1=4...7  any 4 to 7 letters (the shortest length that works is used); call them p1
2     1     3...8     a gap of 3 to 8 letters (the shortest gap that works is used)
3     1     ~p1       the reverse complement of p1, exactly
hit length: 11 to 22 letters
search: the first unit is tried at each position, from left to right; after a hit, the search goes on after its end (with -o: at the next position)
strands: only the given strand is searched (-c: both strands)
threads: with -t N, pieces of long sequences are searched at the same time
```

The last line is about threads. A thread is one CPU core at work. With `-t N` the program uses N cores.

Add `-c` to search the other strand too:

```
$ scan_for_matches -c hairpin.pat < first.fa
>seq1:[6,27]
CGTAACC GGTTAACC GGTTACG 
>seq1:[27,6]
CGTAACC GGTTAACC GGTTACG 
```

The second hit is on the other strand. Its start (27) is larger than its end (6). The hairpin reads the same on both strands, so you get it twice. See [Diagrams](#diagrams) for how to remove the double hit.

## How the search works

The program works like this:

1. It starts at letter 1 of the sequence.
2. It tries to match unit 1 here. Then unit 2 right after it. Then unit 3, and so on.
3. If a unit fails, the program goes back to the unit before and tries another choice for it.
4. If nothing works, it moves one letter to the right and starts again.
5. If all units match, it reports a hit.
6. After a hit, the search goes on after the end of the hit.

Two rules matter a lot.

**A range tries the shortest length first.** `p1=2...6` tries 2 letters, then 3, then 4, and so on. The first length that lets the whole pattern match is the one you get. You do not get the longest one.

```
$ echo 'p1=2...6 0...3 p1' > l.pat
$ printf '>s\nACGTACGTAA\n' > d.fa
$ scan_for_matches l.pat < d.fa
>s:[1,6]
AC GT AC 
```

The data also holds `ACG T ACG` (a longer repeat). The program reports the short one, because it tried 2 letters first.

**After a hit, the search goes on after the end of the hit.** Hits that overlap the first hit are not shown. Use `-o 1` to see them. With `-o 1` the program also shows the other ways to match at the same place:

```
$ scan_for_matches -o 1 l.pat < d.fa
>s:[1,6]
AC GT AC 
>s:[1,7]
ACG T ACG 
>s:[1,8]
ACGT  ACGT 
>s:[2,7]
CG TA CG 
...
```

The `...` means that the list goes on. Write the value after `-o` (for example `-o 1`). The program does not use the value, but it must be there.

**Both strands.** Without `-c`, the program searches only the strand you give it. With `-c`, it also searches the reverse complement. A hit on the other strand has a start that is larger than its end.

```
$ echo ACGGT > p.pat
$ printf '>s\naaACCGTaa\n' > s.fa
$ scan_for_matches -c p.pat < s.fa
>s:[7,3]
ACGGT 
```

The data holds `ACCGT`. Its reverse complement is `ACGGT`, so the hit is on the other strand.

## Letters and words

A word is a run of letters, for example `GAGA`. It matches the same letters in the data.

```
$ echo GAGA > p.pat
$ printf '>s\nttgagaTTgaga\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[3,6]
gaga 
>s:[9,12]
gaga 
```

Rules for words:

* Upper case and lower case are the same.
* `U` is the same as `T`. You can use RNA or DNA letters in the pattern and in the data.
* A word ends at a space, a tab or `[` (see the next section).
* You can use IUPAC letters (IUPAC = International Union of Pure and Applied Chemistry; these are the standard one-letter codes). One letter can stand for several bases.

| code | bases | code | bases |
|---|---|---|---|
| A | A | K | G or T |
| C | C | M | A or C |
| G | G | B | C, G or T |
| T, U | T | D | A, G or T |
| R | A or G (purine) | H | A, C or T |
| Y | C or T (pyrimidine) | V | A, C or G |
| S | C or G | N | A, C, G or T |
| W | A or T | | |

Every letter in this table was tested. The letters `X`, `I` and `Z` are not valid in a DNA pattern.

```
$ echo RRRRYYYY > p.pat
$ printf '>s\nttGAAGCTCCtt\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[3,10]
GAAGCTCC 
```

`RRRR` matched `GAAG` (all purines). `YYYY` matched `CTCC` (all pyrimidines).

```
$ echo TATAA > p.pat
$ printf '>s\nggtauaagg\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[3,7]
tauaa 
```

Watch out: the IUPAC codes work in the pattern only. An `N` in the data never matches a letter in a normal word. `ACGT` does not match `ACNT`, and `ACGT[1,0,0]` does not match it either (both were tested). A pattern `ACNT` does match the data `ACGT`.

## Mismatches, deletions and insertions

Add `[m,d,i]` after a word to allow errors. The three numbers are the maximum counts of:

| position | name | what it allows |
|---|---|---|
| 1st | mismatches `m` (also called substitutions) | a letter of the data is different from the pattern letter |
| 2nd | deletions `d` | a letter of the pattern is missing in the data (the data is **shorter**) |
| 3rd | insertions `i` | an extra letter in the data (the data is **longer**) |

Test with the word `ACGT`. The tests use `^` and `$` (start and end of the sequence) so that the whole sequence must fit:

```
$ echo '^ ACGT[0,1,0] $' > p.pat
$ printf '>s\nAGT\n' > s.fa           # one letter less than ACGT
$ scan_for_matches p.pat < s.fa
>s:[1,3]
 AGT  
$ printf '>s\nACAGT\n' > s.fa         # one letter more than ACGT
$ scan_for_matches p.pat < s.fa
$
```

`[0,1,0]` accepts `AGT` (a letter of the pattern is missing in the data). It does not accept `ACAGT`. Now the third number:

```
$ echo '^ ACGT[0,0,1] $' > p.pat
$ printf '>s\nACAGT\n' > s.fa         # one letter more than ACGT
$ scan_for_matches p.pat < s.fa
>s:[1,5]
 ACAGT  
$ printf '>s\nAGT\n' > s.fa           # one letter less than ACGT
$ scan_for_matches p.pat < s.fa
$
```

So `[0,0,1]` accepts `ACAGT` but not `AGT`.

The summary:

* `[1,0,0]`: one different letter. `ACGT` matches `ACCT`.
* `[0,1,0]`: one letter missing in the data. `ACGT` matches `AGT`.
* `[0,0,1]`: one extra letter in the data. `ACGT` matches `ACAGT`.

A mismatch example:

```
$ echo 'ACGT[1,0,0]' > q3.pat
$ printf '>q\nCCACCTCC\n' > q3.fa
$ scan_for_matches q3.pat < q3.fa
>q:[3,6]
ACCT 
```

This is the same order as in [README.original](README.original): mismatches, then deletions, then insertions. The same definitions are used everywhere in this guide: a deletion is a pattern letter that is missing in the data, and an insertion is an extra letter in the data.

`--explain` uses the same names:

```
$ echo 'ACGT[1,2,3]' > p.pat
$ scan_for_matches --explain p.pat | sed -n 3,4p
unit  line  text         meaning
1     1     ACGT[1,2,3]  the letters ACGT, with up to 1 mismatch, 2 deletions, 3 insertions
```

More watch-outs:

* Errors make the search slower.
* The program is greedy: it takes the first alignment it finds and does not try all alignments. So it can miss a match when it uses deletions and insertions together (`--lint` gives a note).
* To avoid doubt about direction, use only mismatches `[m,0,0]` when you can.
* The same `[m,d,i]` works after a name (`p1[1,0,0]`) and after a reverse complement (`~p1[1,0,0]`).

## Ranges

A range `a...b` matches any `a` to `b` letters. Write three dots. The letters can be anything (a gap).

```
$ echo 'ATG 3...6 TAA' > p.pat
$ printf '>s\nATGcccTAAgggTAA\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,9]
ATG ccc TAA 
```

Here the gap is 3 letters. The range tries 3 first, so the first `TAA` is used.

The range can be empty. `0...3` matches 0, 1, 2 or 3 letters.

```
$ echo 'GGG 0...5 ATG' > p.pat
$ printf '>s\nGGGccATGaaaaGGGATG\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,8]
GGG cc ATG 
>s:[13,18]
GGG  ATG 
```

In the second hit the gap is empty, so there are two spaces in the output.

A range with a name (`p1=4...7`) remembers the letters it matched. See the next section.

Watch out: write the small number first. `10...3` is accepted, but it matches only length 10. Write `3...10` (see [Common mistakes](#common-mistakes)). A fixed length is written with the same number twice, for example `8...8`.

## Names and repeats

Put a name before a unit: `p0=` to `p50=` (the letter `p` and a number from 0 to 50; here the number is not the base `N`). Later you can use the name to match the same letters again.

```
p1=3...3 5...5 p1
```

This pattern means: 3 letters (call them `p1`), 5 any letters, then `p1` again. This is a direct repeat. A target site duplication (TSD) is a direct repeat:

```
$ echo 'p1=3...3 5...5 p1' > p.pat
$ printf '>s\nGATnnnnnGAT\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,11]
GAT nnnnn GAT 
```

Rules for names:

* Define a name before you use it. `p1=...` must come first, then `p1`.
* A name can be defined only once in a pattern.
* A bare `p1` must match the same letters again, exactly.
* You can allow errors: `p1[1,0,0]` is a repeat with one mismatch.

An example with an error:

```
$ echo 'p1=5...5 3...3 p1[1,0,0]' > p.pat
$ printf '>s\nACGTAnnnACCTAcc\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,13]
ACGTA nnn ACCTA 
```

The repeat `ACCTA` differs from `ACGTA` in one letter.

Names can be used again and again. This pattern from the original text matches a repeat of a repeat:

```
$ echo 'p1=4...8 0...3 p2=6...8 p1 0...3 p2' > p.pat
$ printf '>s\nATCTGTCTTTATCTTGTCTTT\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,21]
ATCT  GTCTTT ATCT T GTCTTT 
```

Watch out: a name can catch `N` letters (assembly gaps). Then `p1` and `<p1` accept any letter there, unless you use `--strict-n` (see [mistake 12](#12-n-assembly-gaps-inside-a-name)). A bare `p1` repeats the letters in the same direction. For a reverse complement, use `~p1`. For a reverse, use `<p1`.

## Reverse complement

A DNA strand has two ends, called 5' and 3'. A sequence is written from the 5' end to the 3' end. This is the reading direction. The program reads each sequence from left to right in this direction. The two strands of DNA run in opposite directions. Their bases pair by Watson-Crick pairs: A-T and G-C (in RNA: A-U and G-C).

`~p1` matches the reverse complement of `p1`. This finds inverted repeats: hairpin stems and transposon ends (TIRs).

```
5'-  C G T A A C C  ...  G G T T A C G  -3'
     p1                  ~p1
```

Test:

```
$ echo 'p1=5...5 3...3 ~p1' > e.pat
$ printf '>b\nACGTAnnnTACGT\n' > b.fa
$ scan_for_matches e.pat < b.fa
>b:[1,13]
ACGTA nnn TACGT 
```

`ACGTA` reversed is `ATGCA`. Complemented it is `TACGT`.

You can allow errors after `~p1`:

```
$ echo 'p1=5...5 3...3 ~p1' > e.pat
$ printf '>b\nACGTAnnnTACCTcc\n' > b.fa
$ scan_for_matches e.pat < b.fa
$
$ echo 'p1=5...5 3...3 ~p1[1,0,0]' > e.pat
$ scan_for_matches e.pat < b.fa
>b:[1,13]
ACGTA nnn TACCT 
```

The first pattern has no hit, because `TACCT` is not exactly `TACGT`. The second pattern allows one mismatch. In a stem, one mismatch means one pair that does not pair: `~p1[1,0,0]` allows 1 bad pair, `~p1[2,0,0]` allows 2.

Deletions and insertions work the same way as for words:

```
$ echo 'p1=5...5 3...3 ~p1[0,1,0]' > e.pat
$ printf '>b\nACGTAnnnTAGTcc\n' > b.fa      # TACGT with the C missing
$ scan_for_matches e.pat < b.fa
>b:[1,12]
ACGTA nnn TAGT 
$ echo 'p1=5...5 3...3 ~p1[0,0,1]' > e.pat
$ printf '>b\nACGTAnnnTACCGTcc\n' > b.fa    # TACGT with one extra C
$ scan_for_matches e.pat < b.fa
>b:[1,14]
ACGTA nnn TACCGT 
```

### Rule sets for RNA

In RNA, G can also pair with U. This is the G-U wobble pair (it is not a Watson-Crick pair). Define your own pairing rule with `r0={...}` to `r50={...}`. Each pair has two letters, for example `gu`. The rule set itself matches nothing. Use for example `r1~p1` to match the reverse complement of `p1` with rule set `r1`.

```
$ echo 'r1={au,ua,gc,cg,gu,ug} p1=4...4 3...3 r1~p1' > r.pat
$ printf '>z\nGGCUnnnGGCC\n' > z.fa
$ scan_for_matches r.pat < z.fa
>z:[1,11]
GGCU nnn GGCC 
```

The plain `~p1` does not match here:

```
$ echo 'p1=4...4 3...3 ~p1' > r0.pat
$ scan_for_matches r0.pat < z.fa
$
```

Why: `p1` is `GGCU`. Read it backward: `U C G G`. Each letter pairs with the data letter at the same place:

| place | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| `p1` read backward | U | C | G | G |
| data `GGCC` | G | G | C | C |
| pair (name letter, data letter) | `ug` | `cg` | `gc` | `gc` |

All four pairs are in `r1`. Plain `~p1` would need `AGCC`.

**The order inside a pair matters.** `xy` means: letter `x` in the name (`p1`) and letter `y` in the data at the paired place. So `gu` and `ug` are different. With `r1={au,ua,gc,cg,ug}` the stem `GGG ... UCC` is not found, because G of `p1` against U in the data needs `gu`. With `gu` it is found:

```
$ echo 'r1={au,ua,gc,cg,ug} p1=3...3 3...3 r1~p1' > ra.pat
$ echo 'r1={au,ua,gc,cg,gu} p1=3...3 3...3 r1~p1' > rb.pat
$ printf '>g\nGGGnnnUCC\n' > g.fa
$ scan_for_matches ra.pat < g.fa
$ scan_for_matches rb.pat < g.fa
>g:[1,9]
GGG nnn UCC
```

To allow G-U both ways, write both `gu` and `ug`. Use `--explain` to see the rule set:

```
$ scan_for_matches --explain r.pat | sed -n 7p
rule set r1 (line 1): letter pairs a-t, t-a, g-c, c-g, g-t, t-g
```

Rules for rule sets:

* `r1={AT,UA,gc,cg}` is the standard rule. U and T are the same.
* Mismatches work too: `r1~p1[1,0,0]`.
* A rule set must be defined before it is used.

Watch out: you must write the other end of the hairpin yourself. The program does not guess it. If a motif is at the left end, you must write its reverse complement at the right end by hand. For `CACTA` the right end is `TAGTG` (see [Worked example CACTA transposon](#worked-example-cacta-transposon)).

## Reverse

`<p1` matches `p1` read backward, but **not complemented**. If `p1` matched `GCAT`, then `<p1` matches `TACG`. This is a palindrome of letters (the same letters backward), not a hairpin.

```
$ echo 'p1=6...6 <p1' > p.pat
$ printf '>s\nGCATTAATTACG\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,12]
GCATTA ATTACG 
```

Compare `~p1` and `<p1`:

| unit | `p1` is `GAT`, the unit matches |
|---|---|
| `p1` | `GAT` |
| `~p1` | `ATC` |
| `<p1` | `TAG` |

```
$ echo 'p1=3...3 5...5 ~p1' > p.pat
$ printf '>s\nGATnnnnnATC\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,11]
GAT nnnnn ATC 
$ echo 'p1=3...3 5...5 <p1' > p.pat
$ printf '>s\nGATnnnnnTAG\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,11]
GAT nnnnn TAG 
```

Watch out: for hairpins and TIRs you want `~p1`, not `<p1`. A biologist's palindrome in DNA (for example the site `GAATTC`) is a sequence that is equal to its own reverse complement. Find it with `~p1`, not with `<p1`:

```
$ echo 'p1=3...3 ~p1' > pal.pat
$ printf '>s\nggGAATTCcc\n' > pal.fa
$ scan_for_matches pal.pat < pal.fa
>s:[3,8]
GAA TTC
```

`<p1` finds a letter palindrome such as `GAAAAG`, not `GAATTC`. `--lint` gives a note when you use `<p1`.

## Alternatives

**Do not use the same name in both branches.** A name can be defined only once. The example in [README.original](README.original) uses `p1=` in both branches. In this program that is an error (see [Common mistakes](#common-mistakes)).

`( A | B )` matches A or B. Write a space before and after `|`.

```
$ echo '( GAGA | GCGCA )' > p.pat
$ printf '>s\naaGCGCAttGAGAcc\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[10,13]
GAGA
```

Look at this result. The data holds `GCGCA` at 3 and `GAGA` at 10. The program reports only `GAGA`. The program searches branch 1 over the whole record first. After a hit, the search goes on after its end. So a branch 2 site that lies before a branch 1 hit, or inside it, is not reported. If you need every site of both branches, put each branch in its own pattern file and run both. If the pattern is one whole-record match (`^ ... $`), this does not matter.

Each branch can be a list of units:

```
$ echo '( p1=3...3 3...8 ~p1 | p2=5...5 4...4 ~p2 GGG )' > p.pat
$ printf '>s\naaaGATnnnATCaa\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[4,12]
GAT nnn ATC 
```

An alternative can sit inside a longer pattern:

```
$ echo '( GAGA | CCCC ) TTT' > p.pat
$ printf '>s\nGAGATTTggCCCCTTT\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,7]
GAGA TTT 
>s:[10,16]
CCCC TTT 
```

Rules:

* There are exactly **two** branches in one pair of brackets. For three, put brackets inside: `( GAGA | ( GCGCA | TTCGA ) )`.
* Write spaces around `|`. Without them the pattern cannot be read.

Example with three branches inside `^ ... $`: a whole record that is a start codon, 30 to 60 letters, and one of the stop codons TAA, TAG or TGA. Record `orf3` has no stop codon, so it gives no hit.

```
$ printf '>orf1\nATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATAG\n>orf2\nATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATGA\n>orf3\nATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAATTT\n' > orf.fa
$ echo '^ ATG 30...60 ( TAA | ( TAG | TGA ) ) $' > orf.pat
$ scan_for_matches orf.pat < orf.fa
>orf1:[1,45]
 ATG AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAA TAG
>orf2:[1,45]
 ATG AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAA TGA
```

The IUPAC letter `R` (A or G) makes it shorter: `TAA` and `TAG` are `TAR`. The result is the same:

```
$ echo '^ ATG 30...60 ( TAR | TGA ) $' > orf2.pat
$ scan_for_matches orf2.pat < orf.fa
>orf1:[1,45]
 ATG AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAA TAG
>orf2:[1,45]
 ATG AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAA TGA
```

## Length limit

`length(p1+p2) < N` checks the sum of the lengths of named units. It matches no letters. It only says yes or no. You can add as many names as you want, with `+`. Only `+` and `<` are allowed.

```
$ echo 'p1=1...5 GCA p2=1...6 length(p1+p2) < 6' > p.pat
$ printf '>s\naaaGCAtttt\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,7]
aaa GCA t  
```

`p1` has 3 letters and `p2` has 1, so the sum is 4, which is less than 6. With `< 3` the program must find a smaller sum, so it starts later:

```
$ echo 'p1=1...5 GCA p2=1...6 length(p1+p2) < 3' > p.pat
$ scan_for_matches p.pat < s.fa
>s:[3,7]
a GCA t  
```

Watch out: the units must be defined before the `length` unit. It is "less than", not "less than or equal".

## Start and end of the sequence

`^` matches the start of each FASTA record. `$` matches the end of each record. (A record is one sequence with its `>` line. They do not mean the start or end of the whole file.) They match no letters.

```
$ echo '^ AAA' > p.pat
$ printf '>s\nAAAtAAA\n' > s.fa
$ scan_for_matches p.pat < s.fa
>s:[1,3]
 AAA 
$ echo 'AAA $' > p.pat
$ scan_for_matches p.pat < s.fa
>s:[5,7]
AAA  
```

A run with two records shows this:

```
$ printf '>r1\nAAAtttAAA\n>r2\nAAAgggAAA\n' > two.fa
$ echo '^ AAA' > st.pat
$ scan_for_matches st.pat < two.fa
>r1:[1,3]
 AAA
>r2:[1,3]
 AAA
$ echo 'AAA $' > en.pat
$ scan_for_matches en.pat < two.fa
>r1:[7,9]
AAA
>r2:[7,9]
AAA
```

Use `^` and `$` to test a pattern on a very short sequence. Then random hits from other places cannot disturb your test.

## Weight matrices

A weight matrix gives each position its own score for A, C, G and T. The program adds up the scores along the data. The match works when the sum is larger than a cutoff.

Syntax:

```
{(a,c,g,t),(a,c,g,t),...} > T
```

Each `(a,c,g,t)` is one position. The numbers are integers (you may use negative numbers). The order is A, C, G, T. The sum must be **more than** `T`.

You can add an upper limit in front: `U > {...} > T`. Then the sum must be less than `U` and more than `T`. (`U` is a number. Do not write the word `MAX`.)

The example comes from the original text. It has 8 positions. Data `AACACCGC` scores 558. Data `GACACCGC` scores 626.

```
$ cat w.pat
{(16,0,84,0),(57,10,29,4),(0,80,0,20),(95,0,0,5),(0,100,0,0),(18,60,20,2),(0,0,100,0),(0,50,50,0)} > 450
$ printf '>a\nttAACACCGCtt\n>b\nttGACACCGCtt\n>c\nttTTTTTTTTtt\n' > w.fa
$ scan_for_matches w.pat < w.fa
>a:[3,10]
AACACCGC 
>b:[3,10]
GACACCGC 
```

Now with the upper limit 600:

```
$ cat w.pat
600 > {(16,0,84,0),(57,10,29,4),(0,80,0,20),(95,0,0,5),(0,100,0,0),(18,60,20,2),(0,0,100,0),(0,50,50,0)} > 450
$ scan_for_matches w.pat < w.fa
>a:[3,10]
AACACCGC 
```

Sequence `b` scores 626, which is more than 600, so it is not reported.

Watch out:

* The program does not print the score. It reports the first window above the cutoff.
* An `N` or IUPAC letter in the data gets an average score. Words never match an `N` in the data (but see [Common mistakes](#common-mistakes), 12, for names).
* The line must fit in the pattern file (about 31 999 characters).
* You rarely type a matrix by hand. The script `tools/jaspar2sfm.py` turns a JASPAR matrix into a weight unit. (JASPAR is a public database of transcription-factor binding profiles.) See [README.md](README.md), section "JASPAR matrix to weight unit".

## Protein patterns

Use `-p` when the FASTA file holds proteins. You cannot use `-c` or the reverse complement units (`~p1`, `r1~p1`) with `-p`. `-c` with `-p` is an error.

In protein mode the letters are the one-letter amino acid codes. For example `N` is asparagine, not "any base". The codes:

| | | | |
|---|---|---|---|
| A alanine | C cysteine | D aspartate | E glutamate |
| F phenylalanine | G glycine | H histidine | I isoleucine |
| K lysine | L leucine | M methionine | N asparagine |
| P proline | Q glutamine | R arginine | S serine |
| T threonine | V valine | W tryptophan | Y tyrosine |

`X` means any amino acid.

In a protein pattern:

* A word is a list of amino acid letters, for example `MKV`.
* `any(HQD)` matches one letter: H, Q or D.
* `notany(HK)` matches one letter that is not H and not K.
* `X` matches any one amino acid.
* Names, repeats, ranges, `^`, `$`, alternatives, `length()` and `[m,d,i]` work as in DNA.
* Upper and lower case are the same.

```
$ echo 'p1=0...4 any(HQD) 1...3 notany(HK) p1' > p.pat
$ printf '>s\nAAYWVDAACYWVKK\n' > s.fa
$ scan_for_matches -p p.pat < s.fa
>s:[3,12]
YWV D AA C YWV 
```

```
$ echo 'M any(ST) X X notany(P) K' > p.pat
$ printf '>s\nGMSAAGKLL\n' > s.fa
$ scan_for_matches -p p.pat < s.fa
>s:[2,7]
M S A A G K 
```

Watch out: `-p` and `-c` together stop the program with an error.

## Comments and pattern files

* A `%` starts a comment. The comment goes to the end of the line.
* A pattern file holds **one pattern**.
* You may split the pattern over several lines.
  The program joins all lines into one line with spaces. So break a line only where a space is allowed (between units).
* Blank lines are fine.

```
$ cat m.pat
% a hairpin
p1=4...7   % stem
3...8
~p1
$ scan_for_matches --explain m.pat | sed -n 2p
joined pattern: p1=4...7 3...8 ~p1
```

Watch out: two lines are one pattern. If you write two patterns in one file, the program does not search for two things. It searches for one thing that contains both.

```
$ printf 'AAA\nGGG\n' > j2.pat
$ printf '>s\nAAAtttGGG\nTTAAAGGGTT\n' > j2.fa
$ scan_for_matches j2.pat < j2.fa
>s:[12,17]
AAA GGG 
```

The file holds `AAA` and `GGG` on two lines, so the pattern is `AAA GGG`. It matches `AAAGGG` with no letters between. It does not match `AAAtttGGG`. Put each pattern in its own file. You can give several pattern files to one run with `--format` (see [README.md](README.md)).

## Diagrams

### A hairpin

The pattern `p1=7...7 8...8 ~p1` matches this hairpin. The stem has two halves. The second half is written backward and below the first half, so each letter sits under its partner:

```
5'-  C G T A A C C  -.
     | | | | | | |     G G T T A A C C   (loop = the gap, 8 letters)
3'-  G C A T T G G  -'
```

The top row is `p1`. The bottom row is `~p1` read from right to left (`GGTTACG` backward is `GCATTGG`). One pair per column: C-G, G-C, T-A, A-T, A-T, C-G, C-G. Read on one line, the stem, the loop and the stem look like this:

```
CGTAACC GGTTAACC GGTTACG
p1      loop     ~p1
```

### A TIR transposon

A TIR transposon has two ends that are reverse complements of each other (TIRs). It sits between two copies of a short host sequence (TSDs).

```
5'  TSD   TIR                                   TIR    TSD  3'
    GAT   CACTA ACGTTGC  ... inside ...  GCAACGT TAGTG  GAT
    p1    CACTA p2                          ~p2   TAGTG   p1
          |------ p2 pairs with ~p2 ------|
          |--- CACTA pairs with TAGTG ----|
```

* `p1` is a **direct repeat**: the same letters on both sides (`p1` ... `p1`).
* `p2` and `~p2` are an **inverted repeat**: the right side is the reverse complement of the left side.
* The pattern is: `p1=3...3 CACTA p2=7...7 500...15000 ~p2 TAGTG p1`.
* The part `500...15000` is only the inside. It is not the whole element.

### Both strands and --dedup

The TIR pattern reads the same on both strands. The reverse complement of the element is again a valid element.

```
forward  5'  ... GAT CACTA ... TAGTG GAT ...  3'   hit [201,830]
reverse  5'  ... ATC CACTA ... TAGTG ATC ...  3'   hit [830,201]
                 (the same piece of DNA, read on the other strand)
```

With `-c` the program finds the same element twice. With `--format` and `--dedup`, the second hit is dropped. The strand of the element is then `.` (unknown), because both strands fit. You can see this in the next section.

## Worked example CACTA transposon

The CACTA pattern in `tir_scan_patterns/labelled/00_DTC_published_2025.pat` finds CACTA transposons. (DTC is the name for the CACTA superfamily.) These elements have:

* a 3 letter TSD,
* a TIR that starts with `CACTA`,
* a TIR at the other end that ends with `TAGTG`.

We build the pattern in steps. The test sequence is the file `examples/te.fa` in this repository (980 letters). It has random letters, with one planted element at 201 to 830. Run the commands in the folder `examples/`, or give the path `examples/te.fa`. The element is `GAT CACTA ACGTTGC` + 600 random letters + `GCAACGT TAGTG GAT`.

**Step 1.** The TSD (3 letters) and the first TIR letters.

```
$ echo 'p1=3...3 CACTA' > s.pat
$ scan_for_matches s.pat < te.fa
>chr1:[153,160]
CTA CACTA 
>chr1:[201,208]
GAT CACTA 
```

There are two hits. The one at 153 is a random hit. Short patterns find many random hits.

**Step 2.** Add 7 more letters of the left TIR, and name them `p2`.

```
$ echo 'p1=3...3 CACTA p2=7...7' > s.pat
$ scan_for_matches s.pat < te.fa
>chr1:[153,167]
CTA CACTA ACTTGAA 
>chr1:[201,215]
GAT CACTA ACGTTGC 
```

**Step 3.** Add the inside (a gap of 500 to 15000 letters).

```
$ echo 'p1=3...3 CACTA p2=7...7 500...15000' > s.pat
$ scan_for_matches s.pat < te.fa | cut -c1-70
>chr1:[153,667]
CTA CACTA ACTTGAA CGCCTAGTGGTCAAAGAGTACTGGTAATCGTCGGATCACTAACGTTGCGCTA
```

A range always tries the shortest length first, so the gap takes 500 letters here, and this random hit is 515 letters long. It is not our element yet.

**Step 4.** Add the reverse complement of `p2`, then the end `TAGTG`.

```
$ echo 'p1=3...3 CACTA p2=7...7 500...15000 ~p2 TAGTG' > s.pat
$ scan_for_matches s.pat < te.fa | cut -c1-70
>chr1:[201,827]
GAT CACTA ACGTTGC GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGT
```

The hit is now our element, without the right TSD.

**Step 5.** Add the TSD again (`p1`).

```
$ echo 'p1=3...3 CACTA p2=7...7 500...15000 ~p2 TAGTG p1' > cacta1.pat
$ scan_for_matches cacta1.pat < te.fa | cut -c1-30
>chr1:[201,830]
GAT CACTA ACGTTGC GCTAAAGACAAT
```

The element is at 201 to 830, as planted. The pattern checks itself with `--lint`:

```
$ scan_for_matches --lint cacta1.pat
cacta1.pat: note: the pattern reads the same on both strands: with -c, every element is found twice; use --dedup
cacta1.pat:1: note: `p1`: an N that p1 caught from the data matches any letter here (false hits at assembly gaps); add --strict-n for genomes with gaps
0 errors, 0 warnings, 2 notes
```

The second note is about assembly gaps: see [mistake 12](#12-n-assembly-gaps-inside-a-name).

**Step 6.** Add labels to get GFF3.

A label is a comment that starts with `%@`. Old versions of the program ignore it. `%@element` names the whole element. `%@ TYPE` at the end of a line names the parts that start on that line.

```
$ cat cacta.pat
%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC
p1=3...3               %@ target_site_duplication
CACTA p2=7...7         %@ five_prime_terminal_inverted_repeat
500...15000
~p2 TAGTG              %@ three_prime_terminal_inverted_repeat
p1                     %@ target_site_duplication
```

Check the labels with `--explain`:

```
$ scan_for_matches --explain cacta.pat | sed -n 3,10p
unit  line  text         meaning
1     2     p1=3...3     any 3 letters; call them p1   [label: target_site_duplication]
2     3     CACTA        the letters CACTA, exactly   [label: five_prime_terminal_inverted_repeat]
3     3     p2=7...7     any 7 letters; call them p2   [label: five_prime_terminal_inverted_repeat]
4     4     500...15000  a gap of 500 to 15000 letters (the shortest gap that works is used)
5     5     ~p2          the reverse complement of p2, exactly   [label: three_prime_terminal_inverted_repeat]
6     5     TAGTG        the letters TAGTG, exactly   [label: three_prime_terminal_inverted_repeat]
7     6     p1           the same letters as p1, exactly   [label: target_site_duplication]
```

Run with `-c --dedup --strict-n --format gff3`. (Only columns 1 and 3 to 9 are shown, and the last column is shortened.)

```
$ scan_for_matches -c --dedup --strict-n --format gff3 cacta.pat < te.fa | cut -f1,3-9 | cut -c1-135
##gff-version 3
# scan_for_matches 0.4.1 pattern=cacta.pat
##sequence-region chr1 1 980
chr1	repeat_region	201	830	.	.	.	ID=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;...
chr1	target_site_duplication	201	203	.	.	.	ID=DTC1.lTSD;Parent=DTC1;Name=DTC1;...
chr1	CACTA_TIR_transposon	204	827	.	.	.	ID=DTC1.te;Parent=DTC1;Name=DTC1;...;TSD=GAT_GAT;TIR=CACTAACGTTGC_GCAACGTTAGTG
chr1	five_prime_terminal_inverted_repeat	204	215	.	.	.	ID=DTC1.lTIR;Parent=DTC1.te;...
chr1	three_prime_terminal_inverted_repeat	816	827	.	.	.	ID=DTC1.rTIR;Parent=DTC1.te;...
chr1	target_site_duplication	828	830	.	.	.	ID=DTC1.rTSD;Parent=DTC1;Name=DTC1;...
###
```

This is the real output with the first column and the source column shortened. The element (`.te`) is the hit without its two TSDs: 204 to 827. `repeat_region` covers TSD to TSD: 201 to 830. Column 7 (strand) is `.` because of `--dedup`.

For all output formats (GFF3, BED, JSON lines), see [README.md](README.md). For the 22 ready-made TIR patterns and how they were built, see [tir_scan_patterns/README.md](tir_scan_patterns/README.md).

## Common mistakes

Each mistake below lists the wrong pattern, what happens, and the right pattern.

### 1. A reversed range

Wrong: `10...3 AAA`. The program accepts it, but it matches only length 10. It does not try 3 to 10.

```
$ echo '10...3 AAA' > p.pat
$ scan_for_matches --lint p.pat
p.pat:1: warning: `10...3`: the first number is larger, so only length 10 is matched; write 3...10
0 errors, 1 warnings, 0 notes
```

Right: `3...10 AAA`.

### 2. A name used before it is defined, or never defined

Wrong: `p1=3...3 5...5 ~p2`. There is no `p2`. The program stops with "failed to parse pattern".

```
$ echo 'p1=3...3 5...5 ~p2' > p.pat
$ scan_for_matches --lint p.pat
p.pat:1: error: p2 is used but not defined (write p2=... before it)
1 errors, 0 warnings, 0 notes
```

Wrong: `~p1 p1=3...3`. The name is used before it is defined. The program does not stop. It uses the match from an earlier try, which is rarely what you want.

```
$ echo '~p1 p1=3...3' > p.pat
$ scan_for_matches --lint p.pat
p.pat: note: the pattern cannot be split for threads (a unit uses a match from an earlier attempt): with -t, each long sequence uses one thread
p.pat:1: warning: `~p1` uses p1 before it is defined; it then uses a match from an earlier attempt or sequence
0 errors, 1 warnings, 1 notes
```

Right: define first, use later: `p1=3...3 5...5 ~p1`.

### 3. `p1`, `~p1` and `<p1` mixed up

The three units give three different results. For a TIR or a hairpin you need `~`. In this test the planted element has an inverted repeat. These patterns find nothing:

```
$ echo 'p1=3...3 CACTA p2=7...7 500...15000 p2 TAGTG p1' > w1.pat
$ scan_for_matches w1.pat < te.fa | wc -l
0
$ echo 'p1=3...3 CACTA p2=7...7 500...15000 <p2 TAGTG p1' > w2.pat
$ scan_for_matches w2.pat < te.fa | wc -l
0
$ echo 'p1=3...3 CACTA p2=7...7 500...15000 ~p1 TAGTG p1' > w3.pat
$ scan_for_matches w3.pat < te.fa | wc -l
0
```

* `p2` asks for the same letters again (a direct repeat).
* `<p2` asks for the letters reversed, but not complemented.
* `~p1` is the reverse complement of the TSD, not of the TIR.

Right: `~p2`.

### 4. Forgetting to reverse-complement the motif by hand

The program does not know that the right end of a CACTA element is `TAGTG`. You must write it. If you write `CACTA` at the right end, nothing is found:

```
$ echo 'p1=3...3 CACTA p2=7...7 500...15000 ~p2 CACTA p1' > w4.pat
$ scan_for_matches w4.pat < te.fa | wc -l
0
```

Right: the reverse complement of `CACTA` is `TAGTG`. Reverse `CACTA` to get `ATCAC`, then complement to get `TAGTG`.

Tip: instead of typing the motif twice, use a name: `p3=CACTA ... ~p3`.

### 5. The gap is not the element length

In `p1=3...3 CACTA p2=7...7 500...15000 ~p2 TAGTG p1`, the range `500...15000` is only the inside. The whole hit is longer. Here the whole hit is 630 letters, because: 3 + 5 + 7 + 600 + 7 + 5 + 3.

```
$ echo 'p1=3...3 CACTA p2=7...7 0...400 ~p2 TAGTG p1' > w5.pat
$ scan_for_matches w5.pat < te.fa | wc -l
0
$ echo 'p1=3...3 CACTA p2=7...7 590...610 ~p2 TAGTG p1' > w6.pat
$ scan_for_matches w6.pat < te.fa | cut -c1-30
>chr1:[201,830]
GAT CACTA ACGTTGC GCTAAAGACAAT
```

The gap `0...400` is too small for an inside of 600 letters. Add the TSDs and the TIRs when you count the length of the whole element. `--explain` prints "hit length" for the whole pattern.

### 6. `-o -c`

`-o` takes a value. If you write `-o -c`, the program reads `-c` as the value of `-o`. It does not search the reverse strand. It does not print an error.

```
$ scan_for_matches --lint -o -c e.pat
e.pat: warning: -o takes a value: `-c` is read as the value of -o, not as an option, so the reverse strand is not searched; write -o 1
e.pat: note: the pattern reads the same on both strands: with -c, every element is found twice; use --dedup
0 errors, 1 warnings, 1 notes
```

Right: `-o 1 -c`. Watch out: with `-o`, the program shows overlapping hits on the forward strand only. On the reverse strand, hits never overlap.

### 7. A symmetric pattern with `-c`

A pattern that reads the same on both strands (hairpins, TIR patterns) is found twice with `-c`.

```
$ scan_for_matches --lint -c hairpin.pat
hairpin.pat: warning: the pattern reads the same on both strands: with -c, every element is found twice; use --dedup with --format
hairpin.pat:1: note: `p1=4...7`: the shortest length that works is used, not the longest
0 errors, 1 warnings, 1 notes
```

Right: use `--format gff3` (or `bed6`, `jsonl`) with `--dedup`.

```
$ scan_for_matches -c --dedup --format gff3 hairpin.pat < first.fa
...
```

Or skip `-c`: a symmetric pattern needs only one strand.

### 8. A wide gap before an inexact unit

A wide range, like `500...15000`, followed by an exact word is fast. The program skips gap lengths that cannot work. A wide range followed by a unit with errors can be slow.

```
$ echo 'p1=3...3 CACTA p2=7...7 500...15000 TAGTG[1,0,0] ~p2 p1' > w9.pat
$ scan_for_matches --lint --strict-n w9.pat
w9.pat:1: note: slow: each gap length of `500...15000` is tried (`TAGTG[1,0,0]` after it cannot be used to skip gap lengths)
0 errors, 0 warnings, 1 notes
```

Right: put the exact words first, and put errors only on the `~p2[1,0,0]` unit. It is fast when `p2` has at least 10 letters (5 more per extra mismatch: 15 for `[2,0,0]`), or when 5 exact letters follow it, as `TAGTG` here. Or make the gap narrower. `--lint` tells you when it is slow.

### 9. Two lines joined into one pattern

Two patterns in one file are one pattern. See [Comments and pattern files](#comments-and-pattern-files). Put each pattern in its own file.

### 10. The same name in both branches of an alternative

Wrong: `( p1=3...3 3...8 ~p1 | p1=5...5 4...4 ~p1 GGG )`. This is the example in the original text, but a name can only be defined once.

```
$ echo '( p1=3...3 3...8 ~p1 | p1=5...5 4...4 ~p1 GGG )' > p.pat
$ printf '>s\naaaGATnnnATCaa\n' > s.fa
$ scan_for_matches p.pat < s.fa
failed to parse pattern: ( p1=3...3 3...8 ~p1 | p1=5...5 4...4 ~p1 GGG ) 
$ scan_for_matches --lint p.pat
p.pat:1: error: p1 is defined a second time; a name can be defined only once
1 errors, 0 warnings, 0 notes
```

Right: use `p1` in one branch and `p2` in the other.

### 11. Nested or overlapping hits are lost

After a hit, the search goes on after the end of the hit. Overlapping hits are lost unless you use `-o 1`.

```
$ echo 'p1=3...3 p1' > n.pat
$ printf '>n\nGATGATGATGAT\n' > n.fa
$ scan_for_matches n.pat < n.fa
>n:[1,6]
GAT GAT 
>n:[7,12]
GAT GAT 
$ scan_for_matches -o 1 n.pat < n.fa
>n:[1,6]
GAT GAT 
>n:[2,7]
ATG ATG 
>n:[3,8]
TGA TGA 
>n:[4,9]
GAT GAT 
...
```

Watch out: with `-o` you may get a very large output for a loose pattern.

### 12. `N` (assembly gaps) inside a name

An `N` in the data never matches a word. But a name can catch `N`
letters, because a range takes any letter.  When the name is used
again as `p1` or `<p1`, each caught `N` works like the pattern letter `N`:
it accepts any letter.  So a TSD made of `N` "matches" anything, and you
get false hits next to assembly gaps.  The same is true for the other
IUPAC codes (`R`, `Y`, ...) in the data.  (`~p1` is safe: it never
matches a name that caught `N`.)

```
$ echo 'p1=8...8 20...50 p1' > n.pat
$ printf '>gap\nNNNNNNNNNNACGTACGTACGTACGTACGTGGCATTGACC\n' > n.fa
$ scan_for_matches n.pat < n.fa
>gap:[1,36]
NNNNNNNN NNACGTACGTACGTACGTAC GTGGCATT 
$ scan_for_matches --strict-n n.pat < n.fa | wc -l
0
```

Right: use `--strict-n` for genomes with gaps.  The rule has one
sentence: **with `--strict-n`, each letter that a name caught and that is
not A, C, G or T counts as 1 mismatch when `p1` or `<p1` uses the name
again.**  The mismatches the pattern allows (`[m,0,0]`) decide if the
hit stays.  There is no other limit.

| pattern | the name caught | the data there | without `--strict-n` | with `--strict-n` |
|---|---|---|---|---|
| `p1=3...3 7...7 p1` (exact TSD) | `GAT` | `GAT` | hit | hit |
| `p1=3...3 7...7 p1` (exact TSD) | `GNT` (1 N) | `GAT` | hit | no hit |
| `p1=4...4 6...6 p1[1,0,0]` | `GNTA` (1 N) | `GCTA` | hit | hit (the N uses the 1 mismatch) |
| `p1=4...4 6...6 p1[1,0,0]` | `GNNA` (2 N) | `GCTA` | hit | no hit |
| `p1=4...4 6...6 p1[1,0,0]` | `GNTA` (1 N) | `GCTT` (1 more mismatch) | hit | no hit |
| `p1=4...4 6...6 p1[1,0,0]` | `GRTA` (1 R) | `GCTA` | hit | hit (R is like N) |
| `p1=4...4 6...6 <p1` | `GNTA` (1 N) | `ATCG` | hit | no hit |
| `p1=4...4 6...6 ~p1[1,0,0]` | `GNTA` (1 N) | `TAGC` | no hit | no hit (as before) |
| `p1=4...4 6...6 p1` | `GCTA` | `GCTA`, `N` in the gap | hit | hit (a gap may contain N) |

`--strict-n` works inside the search: a false hit does not hide a real
element next to it.  Without `--strict-n`, the output is the same as in
the original C program.  `--lint` gives a note for each `p1` / `<p1` that
uses a range name, until you add `--strict-n`.

## Options used in this guide

| option | meaning |
|---|---|
| `-c` | also search the other strand (the reverse complement) |
| `-o 1` | also report overlapping and nested hits (the number is not used, but it must be there) |
| `-m N` | stop after N hits |
| `-p` | the sequences and the pattern are protein |
| `-t N` | use N threads (CPU cores); the output does not change |
| `--format gff3` | write GFF3 (also `bed6`, `bed12`, `jsonl`) instead of the original output |
| `--dedup` | with `-c` and `--format`: report an element found on both strands once |
| `--strict-n` | an `N` (or other non-ACGT letter) that a name caught is 1 mismatch when `p1` / `<p1` uses the name again; use it for genomes with gaps ([mistake 12](#12-n-assembly-gaps-inside-a-name)) |
| `--explain`, `--lint` | describe or check a pattern; no sequence is read |

All options: `scan_for_matches --help` and [README.md](README.md).

## Checking a pattern

Check a pattern before a long run. No FASTA file is read, so these commands are fast.

**`--explain`** describes each unit in plain words. It shows the source line, the meaning, the hit length and the labels. Read it and ask: "Is this what I meant?" Look at the "hit length" line. It tells you the shortest and longest possible hit.

```
$ scan_for_matches --explain cacta1.pat
...
hit length: 530 to 15030 letters
```

**`--lint`** finds traps. It prints errors, warnings and notes, and a count. The exit status is 1 when there is an error. Give `--lint` the same options as your real run (`-p`, `-c`, `--dedup`, `--strict-n`, `-o`), because some checks depend on them.

```
$ scan_for_matches --lint -c --dedup --strict-n --format gff3 cacta.pat
```

An error means the pattern does not work. A warning means the pattern runs but probably does not do what you want. A note is information.

**Test on a small planted sequence first.**

1. Make a short FASTA file with the thing you look for, like `te.fa` above. Add some random letters on both sides.
2. Run your pattern on it. Does it find the piece, at the right place?
3. Change one letter to make the piece wrong. Run again. The hit must disappear.
4. Only then run it on a genome.

Use `^` and `$` to test on a very short sequence. Use `-m 1` to stop after the first hit. For a genome with assembly gaps (`N`), add `--strict-n`: else a name that caught `N` can give false hits (see [mistake 12](#12-n-assembly-gaps-inside-a-name)).

## Glossary

- **5' and 3' ends**: the two ends of a DNA strand. A sequence is written and read from 5' to 3'.
- **Alignment**: a way to write two sequences side by side so that equal letters are in the same column.
- **Amino acid**: a building block of a protein. Written with one letter, for example M, K, V.
- **Complement**: the partner letter in a base pair. A pairs with T (or U in RNA). C pairs with G.
- **Deletion**: a letter of the pattern that is missing in the data (the data is shorter). In `[m,d,i]`, the 2nd number allows it.
- **Direct repeat**: the same sequence twice, in the same direction. Example: `GAT ... GAT`.
- **FASTA**: a text file format for sequences. A line that starts with `>` gives the name. The lines after it give the sequence.
- **Gap**: a part of the sequence between two units. You match it with a range such as `3...8`.
- **GFF3**: a table of 9 columns: sequence name, source, type, start, end, score, strand, phase, attributes.
- **Greedy**: the program takes the first alignment it finds and does not try all alignments.
- **G-U wobble pair**: a G-U pair in RNA. It is not a Watson-Crick pair, but it is common.
- **Hairpin**: a piece of RNA or DNA that folds back on itself. A stem, then a loop, then the reverse complement of the stem.
- **Hit**: a piece of the sequence that fits the pattern.
- **Insertion**: an extra letter in the data (the data is longer). In `[m,d,i]`, the 3rd number allows it.
- **Inverted repeat**: the same sequence twice, but the second copy is the reverse complement of the first.
- **IUPAC code**: a one-letter code for one or more bases (IUPAC = International Union of Pure and Applied Chemistry). For example R means A or G, and N means any base.
- **JASPAR**: a public database of transcription-factor binding profiles.
- **Loop**: the unpaired part in the middle of a hairpin.
- **Mismatch** (also called substitution): a position where the data has a different letter than the pattern.
- **Palindrome**: in DNA, a biologist's palindrome (for example `GAATTC`) is a sequence equal to its own reverse complement; `~p1` finds it. A letter palindrome (the same letters backward, for example `GAAAAG`) is found by `<p1`.
- **Pattern**: the text that describes what to search for.
- **Reverse**: the sequence read from the last letter to the first. `GCAT` becomes `TACG`.
- **Reverse complement**: reverse the sequence, then complement each letter. `GCAT` becomes `ATGC`. This is the other strand of the DNA, read 5' to 3'.
- **Stem**: the paired part of a hairpin.
- **Strand**: one of the two chains of DNA. The second strand is the reverse complement of the first.
- **Thread**: one CPU core at work. `-t N` uses N cores.
- **TIR (terminal inverted repeat)**: the two inverted-repeat sequences at the ends of a transposon (each is the reverse complement of the other).
- **Transposon**: a piece of DNA that can move in the genome.
- **TSD (target site duplication)**: a short piece of host DNA that is copied on both sides of a transposon when it inserts. The two copies are a direct repeat.
- **Unit**: one part of a pattern, such as a word, a range or a name.
- **Watson-Crick pair**: the normal base pairs A-T and G-C (A-U in RNA).
- **Weight matrix**: a table that gives a score for each letter at each position of a motif.
