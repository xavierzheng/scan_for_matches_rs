# TIR candidate patterns for Scan for Matches (after Wicker et al. 2007)

Traditional Chinese version: [README_zh-TW.txt](README_zh-TW.txt).
This file is an English translation of it (compiled 2026-09-30). The
last section, "Use with this Rust port", is new and is not in the
Chinese version.

## Important: what these patterns are

This folder has 22 `.pat` files; each one runs on its own.
`00_DTC_published_2025.pat` was supplied by the user and can be checked
against Zheng et al. (2025). All other files are **candidate search
patterns designed from structural features described in the
literature**. They are not complete patterns published by the original
authors, and they have not been benchmarked for sensitivity or precision
on real genomes.

The supplied `Wicker2007_TEclassify.pdf` was read: Figure 1 (p. 974),
Figure 2 (p. 977), Figure 3 (p. 979), Figure 4 (p. 980), and the text of
pp. 976, 977 and 980. The scope is the 2007 classification of that
paper; it is not claimed to cover the current classification.

## 1. Classification in the figure and search structure

One DNA strand, written 5' → 3':

```
TSD | left terminal seed | intervening sequence | reverse-complement seed | TSD
```

A TIR is the pair of reverse-complementary sequences at the two inner
ends of the transposon; a TSD is the host sequence repeated in the same
direction on both outer sides. `p1` stands for the same sequence read
earlier; `~p2` is the reverse complement of `p2` (not only reversed).
The TSD must touch the outer edge of the transposon: no gap may be added
between a TSD and the terminal seed.

The nine superfamilies of Figure 1:

- **DTT** Tc1-Mariner: TA TSD.
- **DTA** hAT: 8 bp TSD; TIRs of 5–27 bp, no general diagnostic motif.
- **DTM** Mutator: 9–11 bp TSD; ends often G...C; TIRs can be very long,
  very short, or not detectable.
- **DTE** Merlin: 8–9 bp TSD; TIRs up to tens or hundreds of bp.
- **DTR** Transib: 5 bp TSD.
- **DTP** P: 8 bp TSD.
- **DTB** piggyBac: TTAA TSD.
- **DTH** PIF-Harbinger: 3 bp TSD; the text mentions a preference for TAA.
- **DTC** CACTA: in plants 3 bp TSD, ends CACTA or CACTG; the animal /
  fungal type placed here by the paper can be CCC...GGG with a 2 bp TSD.

Maverick (DMM) has TIRs in the figure, but it is in subclass 2, order
Maverick, not a tenth superfamily of the TIR order. Helitrons have no
paired terminal TIRs, and the hairpin near their 3' end cannot be used as
the two ends of a TIR transposon. DIRS also shows inverted terminal
structures in the figure but belongs to class I. **An inverted-repeat
structure does not mean that a TIR-order TE has been identified.**
Figure 2 shows that one CACTA family can include derivatives with
different internal deletions: a complete structure does not mean the
element can transpose autonomously. Figure 3 is a reminder to combine
nucleotide homology, protein homology and structure.

## 2. How to read these patterns

Each `.pat` file holds one pattern only, with no title or comments.
**Do not join the patterns of the folder into one input file**: Scan for
Matches reads all consecutive text as one pattern.

- `p1=3...3`: read and remember 3 bp; the final `p1` requires exactly the
  same sequence in the same direction.
- `p2=12...12`: read a 12 bp terminal seed; `~p2` requires its reverse
  complement at the other end.
- `CACTA[0,0,0]`: CACTA must match exactly.
- `~p2[1,0,0]`: the other end may have at most 1 substitution relative to
  the reverse complement of `p2`; no indels.
- A repeated reference without error settings must match exactly.

Here the **seed** is the terminal piece checked by the search, not the
definition of the full TIR length of a superfamily. For example, if the
real TIR is 200 bp long, the outermost 20 bp can be the seed and the rest
of the TIR falls in the middle part. After a candidate is found, extend
the alignment to measure the true TIR length.

Except for the original DTC pattern and the DMM supplement, `50...30000`
is the middle-part search window chosen for this set; it is not the full
element size range claimed by the literature. Seed lengths of 8, 12, 15
and 20 bp are also search settings. This window holds some small
deletion derivatives and larger elements, but it can be slow, and it does
not cover all very short elements, very long elements, or elements with
large nested insertions.

With seed length *k*, middle part length *g* and no indels:

- transposon length (without the two TSDs) = 2*k* + *g*
- full hit length = 2 × TSD length + 2*k* + *g*

Original DTC: `p2` is only 7 bp, but with the 5 bp of CACTA each end
checks 12 bp. `500...15000` is the middle part, so the transposon length
is 524–15024 bp and the hit length with both TSDs is 530–15030 bp. It
misses TIRs shorter than 12 bp, mutated ends, the CACTG type, the CCC
type and elements outside the window. The literature describes CACTA
TIRs as short as 10 bp, so seed10 versions are also given. To limit the
**transposon body** strictly to 500–15000 bp with this 12 bp seed, change
the middle part to `476...14976`.

## 3. Choosing files, and limits of classification

- `01_DTT_seed12`: TA-TSD candidates without a Stowaway motif.
- `01_DTT_Stowaway_seed12`: the CTCCTCCC...GGGAGGAG subtype; it does not
  represent all DTT.
- `02_DTA_seed8`: 8 bp TSD + 8 bp seed; also matches other elements with
  an 8 bp TSD.
- `02_DTA_short5to7`: recovers very short TIRs; random and low-complexity
  background hits increase greatly. To narrow the window after Wicker's
  description of hAT bodies under 4 kb, the seed8 version can change
  `50...30000` to `50...3983`; this is an extra filter and can miss
  derivatives with insertions.
- `03_DTM_G_seed20`: G/C-terminal type, 9–11 bp TSD, 20 bp exact seed.
- `03_DTM_unanchored_seed20`: supplementary search without the G/C ends.
  Mutators without TIRs cannot be found by these IR patterns; they need
  homology and other evidence.
- `04_DTE_seed20`: 8–9 bp TSD, 20 bp seed; a longer seed cannot exclude
  hAT / P.
- `05_DTR_CAC_seed15`: CAC...GTG + 5 bp TSD. CAC comes from the Transib
  end alignment of Kapitonov & Jurka (2005), not from a sequence listed
  directly in Wicker Figure 1.
- `05_DTR_unanchored_seed15`: supplementary pattern without CAC.
- `06_DTP_seed12`: 8 bp TSD, 12 bp seed; only a P-compatible candidate,
  not a P-specific pattern.
- `07_DTB_seed12`: TTAA is the outer TSD; do not use it as the inner TIR.
- `08_DTH_seed12`: any identical 3 bp TSD; can also hit CACTA and other
  candidates.
- `08_DTH_TWA_seed12`: only the preferred TAA or TTA sites. `p1=TWA`
  remembers the triplet actually matched, so both ends must be TAA or both
  TTA. Do not change it to independent TWA on each side, or a left TAA
  with a right TTA could be accepted by mistake.
- `08_DTH_Tourist_seed8`: after the Tourist candidate criteria of RiTE:
  G/C ends, seed of at least 8 bp, 3 bp TSD. `p3=S` remembers the actual
  G or C and `~p3` requires its complement; it is not an independent S
  on each end. This pattern is our transcription; the RiTE paper does not
  give this complete pattern.
- `09_DTC_CACTA_seed12`, `09_DTC_CACTG_seed12`: the two plant end types,
  searched separately.
- `09_DTC_CACTA_seed10`, `09_DTC_CACTG_seed10`: also find shorter exact
  TIRs.
- `09_DTC_CCC_seed12`: the animal / fungal CCC type described by Wicker
  2007, with a 2 bp TSD.
- `09_DTC_CACTA_relaxed1`: allows one substitution in the `p2` part only;
  TSD and CACTA stay exact.
- `10_DMM_Maverick_seed30`: extra structural search for order Maverick;
  6 bp TSD, 30 bp seed; `9940...19940` gives a body length of
  10000–20000 bp. Long TIRs and homology to POLB, integrase etc. must
  still be checked; a hit alone does not confirm a Maverick.

DTA, DTP and DTE candidates overlap; seed length is a sensitivity choice,
not a reliable boundary between classes. DTH broad is not specific
either. Without evidence that separates them, keep candidates as
unclassified.

## 4. Running and post-processing

(This section describes the original C program from the SEED site; see
the last section for this Rust port.)

Run one `.pat` at a time where scan_for_matches is installed:

```sh
scan_for_matches -o 1 09_DTC_CACTA_seed12.pat < genome.fa > DTC.hits.fa
```

The version tested was the C source from the SEED site; `-o 1` keeps
overlapping and alternative matches. Options can differ between
distributions, so check your own version. By default overlaps are
skipped, which can miss nested or overlapping candidates.

The complete-element patterns of this set have the same structure after
reverse complement (DTH TWA already covers TAA/TTA), so scanning the
input strand finds elements in both orientations, and `-c` is not needed
for them. This does not tell the transcription direction of the
transposase. Only if you change a pattern into an asymmetric,
family-specific one must you handle both strands. In tests, `-c` of this
old version did not keep all overlapping matches on the reverse strand;
for a complete list, scan the forward and reverse-complement FASTA
separately, then convert coordinates and remove duplicates.

Post-processing should include at least:

1. Remove doubtful hits caused by N or low complexity in the TSD or
   terminal seeds; check for assembly gaps.
2. The candidate transposon body is the output without the two TSDs; for
   variable-length TSDs, read the actual matched length.
3. Extend inward from both outer edges and align the full TIRs; do not
   report the seed length as the biological TIR length.
4. Remove duplicates with the same coordinates that come from different
   seed splits; keep truly different nested candidates.
5. Combine multi-copy boundaries, empty insertion sites (if any), known TE
   libraries, transposase homology or protein domains. A low-complexity
   IR or a chance TSD alone does not prove a TE.
6. For large genomes, test a small region first, then scan by chromosome
   or overlapping windows; the window overlap must cover at least the
   longest allowed full hit, and coordinates must be converted back and
   duplicates removed. The longest full hit of the general patterns here
   is about 30.1 kb, so a 31 kb overlap can be used. Nested insertions
   longer than the window are still missed.

That source version fixes the size of one input sequence at 250,000,000
nt; split longer chromosomes into overlapping pieces or use a fixed
version. Splitting does not remove the cost of the wide backtracking
search.

## 5. Literature check: what is original and what is not

1. Wicker et al. 2007. A unified classification system for eukaryotic
   transposable elements. *Nature Reviews Genetics* 8:973–982. The
   supplied PDF is the main basis for the classification here; its
   figures and text were read. <https://doi.org/10.1038/nrg2165>
2. Dsouza, Larsen & Overbeek. 1997. Searching for patterns in genomic
   data. *Trends in Genetics* 13:497–498. The original paper of the tool;
   the citation can be checked, PubMed has no abstract. This citation is
   not evidence that patterns for the nine TE superfamilies were
   published. <https://pubmed.ncbi.nlm.nih.gov/9433140/>
3. Ashok Aiyar, 1999-01-20, bionet.software thread "How to find
   (imperfect) repeats in DNA sequences?". Recommends Scan for Matches to
   find inverted repeats; gives no diagnostic patterns for the nine TE
   superfamilies. <https://groups.google.com/g/bionet.software/c/YeMC40OL6_I>
4. The SEED Team, 2010-07-16. Scan For Matches official description.
   Gives the syntax for `p1`, `~p1`, errors, variable gaps and the
   backtracking search.
   <https://blog.theseed.org/servers/2010/07/scan-for-matches.html>
   Source package used for the tests:
   <https://www.theseed.org/servers/downloads/scan_for_matches.tgz>
5. Wicker et al. 2003. CACTA Transposons in Triticeae. A Diverse Family of
   High-Copy Repetitive Elements. *Plant Physiology* 132:52–63. Uses BLAST
   and dot plots, not a verified Scan for Matches method. Supports plant
   CACTA TIRs of 10–28 bp and deletion derivatives of all sizes; has a
   274 bp small element and discusses the 23 kb Candystripe1, so
   500–15000 is not the full range.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC166951/>
6. Nicolas et al. 2005. Suffix-tree analyser (STAN): looking for
   nucleotidic and peptidic patterns in chromosomes. *Bioinformatics*
   21:4408–4410. Compares STAN, PatScan and GenLang on AtREP3. AtREP3 is
   a Helitron-type TE, so its near-3' hairpin pattern must not be
   described as a paired-end pattern of a TIR-order element. Also notes
   the PatScan limits of that time for overlapping alternative solutions.
   <https://academic.oup.com/bioinformatics/article/21/24/4408/179826>
7. Biostars: Searching Repeats And Palindromic Sequences In Dna
   Sequences. The discussion recommends scan_for_matches but gives no
   complete per-superfamily patterns. <https://www.biostars.org/p/79567/>
8. Copetti et al. 2015. RiTE database: a resource database for genus-wide
   rice genomics and evolutionary biology. *BMC Genomics* 16:538. The
   Stowaway / Tourist section supports specific TIR / TSD criteria but
   does not print complete SFM patterns for them. The SFM pattern it does
   name and print is for the Helitron 3' hairpin:
   `p1=7...10 2...4 ~p1[1,0,0] 6...10 CTRRT`. It is only a method
   reference and is not one of the TIR `.pat` files here.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC4508813/>
9. Zheng et al. 2025. Transposable elements drive evolution and perturb
   gene expression in *Brassica rapa* and *B. oleracea*. *The Plant
   Journal* 123:e70452. The Methods confirm the supplied DTC pattern; file
   00 keeps it unchanged. <https://doi.org/10.1111/tpj.70452>,
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC12401559/>
10. Kapitonov & Jurka. 2005. RAG1 Core and V(D)J Recombination Signal
    Sequences Were Derived from Transib Transposons. *PLoS Biology*
    3:e181. Supports the conserved CAC at Transib ends and the 5 bp TSD;
    not a source of the SFM method.
    <https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030181>
11. Jurka & Kapitonov. 2001. PIFs meet Tourists and Harbingers: A
    superfamily reunion. *PNAS* 98:12315–12316. Supports the 3 bp TSD and
    the TTA preference of PIF / Tourist.
    <https://pmc.ncbi.nlm.nih.gov/articles/PMC60043/>

No verifiable early paper or discussion was found that printed complete
SFM patterns for each of Wicker's nine superfamilies. So the designs here
are not presented as published patterns, and no claim is made that all
early studies were found.

## 6. Scope of validation

The 22 patterns were run with the program compiled from the official
source; each passed: a synthetic positive case, a negative case with a
wrong TSD, a negative case with more TIR substitutions than allowed, and
a check that the unanchored pattern finds the expected coordinates.
Negative cases use `^` and `$` to isolate the structure under test, so
chance matches elsewhere do not interfere. Also passed: 14 tests of TWA,
G/C, 9/10/11 bp TSD, 8/9 bp TSD, a single mismatch, and the lower and
upper bounds of the DTC middle part, plus a check that `-o 1` lists three
overlapping AA on the forward strand. These results are **not** a recall,
specificity or speed benchmark on real data.

Extra test of the indel meaning in the old version:

- `ACGT[0,1,0]` matches AGT but not ACAGT.
- `ACGT[0,0,1]` matches ACAGT but not AGT.

So, described for "the sequence being matched relative to the pattern",
the actual order in this source is [substitutions, deletions,
insertions]. The variable names inside the source code can make this
confusing; all patterns here set the last two values to 0, to avoid
ambiguity in indel direction and end coordinates.

`catalog.json` stores the purpose of each pattern; `validation.json`
stores the test results and the SHA-256 of the source package. The
download did not include program binaries, and no user genome was run.

## Use with this Rust port (new; not in the Chinese version)

- The `.pat` files above are unchanged.
- `labelled/` holds copies of the 22 patterns with `%@` labels for
  `--format gff3|bed6|bed12|jsonl`; the labels name the element, the TSDs
  and the TIRs (see the main [README](../README.md), "Output formats").
  Without `--format` they give the same output as the plain files.
- The 250,000,000 nt limit of the C program is removed in this port
  (0.1.0).
- With `-c`, a symmetric pattern (all patterns here) reports each element
  on both strands; with `--format`, add `--dedup` to keep one copy (strand
  `.`).
- Example (one pattern, bgzip or plain FASTA, 8 threads):

  ```sh
  scan_for_matches -t 8 -c --dedup --format gff3 --name-prefix DTC \
      --input genome.fna.gz --output DTC.gff3 \
      tir_scan_patterns/labelled/00_DTC_published_2025.pat
  ```

  Give each pattern its own `--name-prefix` (or `--name-start`), so Names
  do not repeat between runs.
- Speed of all 22 patterns on the 1.0 Gb *B. napus* genome: see
  `CHANGELOG.md`.
