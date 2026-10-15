English | [繁體中文](USER_GUIDE_zh.md)

# User guide: scan_for_matches for beginners

This guide is for a biologist who has never used `scan_for_matches`.
You learn how to start, how to run a real genome, and how to read the result.
Every command here was run, and the output shown is the real output
(long lines are cut with `...`).

## Table of contents

0. [Words used in this guide](#words-used-in-this-guide)
1. [What the tool does](#1-what-the-tool-does)
2. [Get the program](#2-get-the-program)
3. [A first run in 1 minute](#3-a-first-run-in-1-minute)
4. [A real genome walk-through (B. oleracea)](#4-a-real-genome-walk-through-b-oleracea)
5. [Look at the results](#5-look-at-the-results)
6. [Run many patterns at once](#6-run-many-patterns-at-once)
7. [Write your own pattern](#7-write-your-own-pattern)
8. [All options](#8-all-options)
9. [Speed and memory](#9-speed-and-memory)
10. [Common problems (FAQ)](#10-common-problems-faq)

## Words used in this guide

Hard words are explained here. Click a word in the text to come back to this table.

| word | meaning |
|---|---|
| transposon | a piece of DNA that can copy or move itself in a genome |
| TIR transposon | a transposon with a TIR at each end (a DNA transposon) |
| TIR | terminal inverted repeat. The two ends of the element have the same letters, one end reverse-complemented |
| TSD | target site duplication. A short piece of host DNA, copied once on each side of the element when it jumped in |
| CACTA, DTC | CACTA is a family of plant TIR transposons. Its ends start with `CACTA`. DTC is its short code. Other codes: DTA, DTH, DTM, DTT ... are other families |
| assembly gap | a run of `N` in the genome: the assembly has no real sequence there |
| FASTA | a text file of sequences. A line starting with `>` is the name; the lines after it are the letters |
| GFF3 | a text table of features with 9 tab-separated columns (section 4.5) |
| BED | a simpler table: one line per region, start counted from 0 (section 5.2) |
| JSONL | JSON lines: one JSON object per line, one line per hit (section 5.3) |
| PATH | the list of folders where the shell looks for programs. `export PATH=/my/folder:$PATH` adds a folder |
| stdin, stdout | the "in" and "out" streams of a program. `< file` sends a file to stdin. `> file` sends stdout to a file |
| exit status | a number a program gives when it ends. `0` means no problem. See it with `echo $?` |
| thread, core | a core is one CPU unit. A thread is one worker that runs on a core. Use `-t N` threads on N cores |
| Slurm, sbatch | the job manager of a cluster. `sbatch script.sh` puts a script in the queue |
| login node | the computer where you log in. Do not run big jobs on it |
| `#SBATCH` lines | lines at the top of a script. They tell Slurm what you need (cores, memory, name) |
| bgzip | gzip made of many small blocks. The program can read it, and read it with many threads |
| IGV | a genome browser program. It shows BED and GFF3 files on the genome |
| Sequence Ontology (SO) | a list of standard names and numbers for sequence features, for example `SO:0000657` is `repeat_region` |
| greedy alignment | the matcher takes the first alignment that fits and does not try others. `--align` tries the others |
| gap skipping | for a wide gap such as `500...15000`, the program tries only the gap lengths where the next unit can match. This is much faster |
| Mb, Gb | million and billion bases (sequence length); not MB, GB (file or memory size) |
| 5-mer index | a table of the places of every 5-letter word in a record. The program uses it for gap skipping |

## 1. What the tool does

`scan_for_matches` reads sequences in [FASTA](#words-used-in-this-guide) format (DNA or protein).
It looks for places that fit a pattern that you write.
A pattern can describe a short motif, a hairpin, a direct repeat, or a whole
[TIR transposon](#words-used-in-this-guide) (two end repeats, the [TIRs](#words-used-in-this-guide), and two
target site duplications, the [TSDs](#words-used-in-this-guide)).
The result is a list of hits. You can get it as [GFF3, BED, JSONL](#words-used-in-this-guide)
(JSON lines: one JSON object per hit), or as the plain text of the original program.

It is a Rust port of an old C program by Ross Overbeek.
The pattern language is the same as in the old program.

## 2. Get the program

* macOS (Apple silicon): download `scan_for_matches-<version>-macos-arm64.tar.gz`
  from the GitHub Releases page:
  <https://github.com/xavierzheng/scan_for_matches_rs/releases>
* Linux and other systems: build it. You need Rust (<https://rustup.rs>).
  See [README.md](README.md#build).

  ```sh
  git clone https://github.com/xavierzheng/scan_for_matches_rs.git
  cd scan_for_matches_rs
  cargo build --release      # the program: target/release/scan_for_matches
  ```

Put the program in a folder that is in your `PATH`, or add its folder:

```sh
export PATH=$PWD/target/release:$PATH
scan_for_matches --help
```

The first line of the help text is the version. This guide was run with
`scan_for_matches 0.4.3`.

`git clone` makes the folder `scan_for_matches_rs/`. It holds the ready-made
patterns (`tir_scan_patterns/`) and this guide.
**The commands in sections 4 to 6 are run from one work folder.** That folder
holds your genome file and the folder `scan_for_matches_rs/` side by side:

```
$ ls
genome.fa.gz
scan_for_matches_rs
```

(Here both are links to the real files. A real copy works the same.) Output
files are written to this work folder. Run `sbatch` from this folder too: Slurm
starts the job in the folder where you ran `sbatch`.

## 3. A first run in 1 minute

Make a small FASTA file with `printf`. The file has two sequences.
`seq1` has the word `GATTACA` once.

```
$ printf '>seq1\nCCCCGATTACAAAAAAGTAATCGGGGG\n>seq2\nCCCCCCCCCCCCCCCCCCCCCC\n' > tiny.fa
```

Make a pattern file with a one-line pattern:

```
$ echo 'GATTACA' > motif.pat
```

Run it. The FASTA file goes in with `<`. The hits go to the screen.

```
$ scan_for_matches motif.pat < tiny.fa
>seq1:[5,11]
GATTACA 
```

How to read it:

* `>seq1:[5,11]` is the hit. It is in `seq1`, from letter 5 to letter 11.
  The first letter is number 1.
* The next line is the matched letters.
* `seq2` has no hit, so nothing is printed for it.

Now a pattern with more meaning: a hairpin. A hairpin is a stem, a loop, and
the reverse complement of the stem. Here the stem has 4 to 7 letters and the
loop has 3 to 8 letters.

```
$ echo 'p1=4...7 3...8 ~p1' > hairpin.pat
$ scan_for_matches hairpin.pat < tiny.fa
>seq1:[3,24]
CCGATTA CAAAAAAG TAATCGG 
```

`p1=4...7` takes 4 to 7 letters and calls them `p1`. `~p1` is the reverse
complement of `p1`. (Reverse the letters, then change A to T and C to G.)
`CCGATTA` becomes `TAATCGG`.

Ask the program to say in plain words what a pattern does. No FASTA file is
read.

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

Ask for GFF3 instead of the original output:

```
$ scan_for_matches --format gff3 hairpin.pat < tiny.fa
##gff-version 3
# scan_for_matches 0.4.3 pattern=hairpin.pat
##sequence-region seq1 1 27
seq1	scan_for_matches	sequence_motif	3	24	.	+	.	ID=sfm1;Name=sfm1;Method=structural;Sequence_ontology=SO:0001683
###
```

## 4. A real genome walk-through (B. oleracea)

The example genome is *Brassica oleracea* (file `genome.fa.gz`).
It has 289 records: the chromosomes `II_C01` to `II_C09`, and small contigs.
The file is [bgzip](#words-used-in-this-guide). The program reads it directly: you do not need to unpack it.
Run the commands of this section from the work folder (section 2).

We search for CACTA transposons (DTC, a TIR transposon family).
We use the ready-made pattern `00_DTC_published_2025.pat`.
All 22 ready-made patterns are in
[`tir_scan_patterns/labelled/`](tir_scan_patterns/README.md).

### 4.1 Look at the pattern

```
$ cat scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
% 00_DTC_published_2025.pat: labelled copy for --format (scan_for_matches 0.3.0)
%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC
p1=3...3               %@ target_site_duplication
CACTA[0,0,0] p2=7...7  %@ five_prime_terminal_inverted_repeat
500...15000
~p2 TAGTG[0,0,0]       %@ three_prime_terminal_inverted_repeat
p1                     %@ target_site_duplication
```

* A line that starts with `%` is a comment.
* `%@element ...` names the whole element, and gives the Name prefix (`DTC`).
* `%@ type` at the end of a line gives a name to the units on that line.
  These names become the features in the GFF3.
* `p1=3...3` takes any 3 letters (the TSD, target site duplication). The last
  `p1` needs the same 3 letters again.
* `CACTA[0,0,0] p2=7...7` is the 5' end of the element: `CACTA` with no
  mismatch, then 7 letters called `p2`.
* `500...15000` is the middle part, 500 to 15000 letters of anything.
* `~p2 TAGTG[0,0,0]` is the 3' end: the reverse complement of `p2`, then `TAGTG`.

### 4.2 Check it: `--lint` and `--explain`

`--lint` looks for traps. Give it the same options as the real run.

```
$ scan_for_matches --lint -c scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat: warning: the pattern reads the same on both strands: with -c, every element is found twice; use --dedup with --format
scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat:7: note: `p1`: an N that p1 caught from the data matches any letter here (false hits at assembly gaps); add --strict-n for genomes with gaps
0 errors, 1 warnings, 1 notes
```

The program tells us two things. With `-c`, add `--dedup`. For a genome with
assembly gaps (runs of `N`), add `--strict-n`. Check again with both:

```
$ scan_for_matches --lint -c --dedup --strict-n --format gff3 scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
0 errors, 0 warnings, 0 notes
```

`--explain` shows each unit (here the first lines):

```
$ scan_for_matches --explain scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
pattern file: scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat (DNA)
joined pattern: p1=3...3 CACTA[0,0,0] p2=7...7 500...15000 ~p2 TAGTG[0,0,0] p1
unit  line  text          meaning
1     3     p1=3...3      any 3 letters; call them p1   [label: target_site_duplication]
...
hit length: 530 to 15030 letters
```

The line `hit length` tells you the shortest and the longest element that the
pattern can find.

### 4.3 Run on ONE chromosome first

Always test on a small piece first. Cut chromosome `II_C09` out of the genome:

```
$ zcat genome.fa.gz | awk '/^>/{p=($1==">II_C09")} p' > II_C09.fa
```

Run the pattern on it:

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format gff3 --name-prefix DTC \
      --input II_C09.fa --output dtc_C09.gff3 \
      scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
$ head -4 dtc_C09.gff3 | cut -c1-200
##gff-version 3
# scan_for_matches 0.4.3 pattern=scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
##sequence-region II_C09 1 65908049
II_C09	scan_for_matches	repeat_region	854082	860374	.	.	.	ID=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657
$ grep -c repeat_region dtc_C09.gff3
396
```

This chromosome (65.9 Mb) has 396 elements. The run took a few seconds.
Looks good? Then do the whole genome.

### 4.4 The whole genome

What each option does:

| option | why |
|---|---|
| `-t 20` | use 20 [threads](#words-used-in-this-guide) (one thread per CPU core). Ask Slurm for the same number of cores |
| `-c` | search both strands |
| `--dedup` | the same element found on both strands is written once (strand `.`) |
| `--strict-n` | no false hits at [assembly gaps](#words-used-in-this-guide) (runs of `N`) |
| `--format gff3` | write GFF3 |
| `--name-prefix DTC` | the elements are named `DTC1`, `DTC2`, ... |
| `--progress` | write a progress line to stderr every minute |
| `--input`, `--output` | the FASTA file in, the hit file out |

Do not run a big job on the [login node](#words-used-in-this-guide) of a cluster. Send it with
[Slurm](#words-used-in-this-guide). The `#SBATCH` lines tell Slurm what you need.
Put this in a file `run_dtc_genome.sh` in the work folder:

```bash
#!/bin/bash
#SBATCH --job-name=sfm_dtc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem=20G
#SBATCH --output=sfm_dtc-%j.out

set -euo pipefail
export PATH=/path/to/folder_with_binary:${PATH}

GENOME=genome.fa.gz
PAT=scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat

time scan_for_matches -t 20 -c --dedup --strict-n --format gff3 \
    --name-prefix DTC --progress \
    --input ${GENOME} --output dtc_II.gff3 ${PAT}
```

Change `/path/to/folder_with_binary` to the folder of your program (after a
build: `scan_for_matches_rs/target/release`; if unsure: `which scan_for_matches`).
Some clusters also need `#SBATCH -p` (partition) or `#SBATCH -A` (account); ask
your cluster admin. Send it from
the work folder (the paths in the script are relative to it):

```
$ sbatch run_dtc_genome.sh
```

The cores in `-c 20` and the number after `-t` must be the same.
`--mem=20G` is a guess with room to spare. (README.md reports 2.5 GB at
maximum for this pattern on the 2.3 Gb maize genome.)
Slurm writes the progress lines to `sfm_dtc-<jobid>.out`:

```
scan_for_matches: [0:00:00] 0 records (0.0 Mb) done, 0 hits; started
scan_for_matches: [0:00:06] 289 records (563.4 Mb) done, 3272 hits; finished

real	0m6.650s
user	0m41.698s
sys	0m4.075s
```

The genome has 563.4 Mb. The job took 7 seconds with 20 threads. It found 3 272 elements.

**How much memory did it really use?** Put `/usr/bin/time -v` before the program
in the script (`/usr/bin/time -v scan_for_matches ...`). At the end of the
`.out` file look for "Maximum resident set size". For this job (run again with
`/usr/bin/time -v`):

```
$ grep 'Maximum resident' sfm_mem-580692.out
	Maximum resident set size (kbytes): 1026796
```

That is about 1.0 GB. Do not trust `sacct -j JOBID -o MaxRSS` for a plain
script: it shows only the small batch step. Real output of the same job:

```
$ sacct -j 580692 -o JobID,Elapsed,MaxRSS,State
JobID           Elapsed     MaxRSS      State 
------------ ---------- ---------- ---------- 
580692         00:00:40             COMPLETED 
580692.batch   00:00:40     12340K  COMPLETED 
580692.exte+   00:00:40      1048K  COMPLETED 
580692.0       00:00:32    885964K  COMPLETED 
```

Line `.batch` (12 MB) is wrong for the program. Line `.0` (886 MB) is right: it
exists because this script also started the program with `srun`. So `sacct`
and `seff` are right only if you start the program with `srun scan_for_matches ...`.

### 4.5 Read the GFF3

GFF3 is a table with 9 columns, separated by tabs. These are the first
lines of `dtc_II.gff3` (column 9 is cut):

```
##gff-version 3
# scan_for_matches 0.4.3 pattern=scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
##sequence-region II_C01 1 52350442
II_C01	scan_for_matches	repeat_region	2874799	2878133	.	.	.	ID=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657
II_C01	scan_for_matches	target_site_duplication	2874799	2874801	.	.	.	ID=DTC1.lTSD;Parent=DTC1;Name=DTC1;...
II_C01	scan_for_matches	CACTA_TIR_transposon	2874802	2878130	.	.	.	ID=DTC1.te;Parent=DTC1;Name=DTC1;...;Sequence_ontology=SO:0002285;TSD=TAA_TAA;TIR=CACTACAAGAAA_TTTCTTGTAGTG
II_C01	scan_for_matches	five_prime_terminal_inverted_repeat	2874802	2874813	.	.	.	ID=DTC1.lTIR;Parent=DTC1.te;Name=DTC1;...
II_C01	scan_for_matches	three_prime_terminal_inverted_repeat	2878119	2878130	.	.	.	ID=DTC1.rTIR;Parent=DTC1.te;Name=DTC1;...
II_C01	scan_for_matches	target_site_duplication	2878131	2878133	.	.	.	ID=DTC1.rTSD;Parent=DTC1;Name=DTC1;...
###
II_C01	scan_for_matches	repeat_region	3398540	3401290	.	.	.	ID=DTC2;Name=DTC2;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657
...
```

The columns:

| column | name | in this file |
|---|---|---|
| 1 | seqid | the FASTA record: `II_C01` |
| 2 | source | the program: `scan_for_matches` |
| 3 | type | the kind of feature (table below) |
| 4, 5 | start, end | position on the forward strand. The first letter is 1. Start is not more than end, also for hits on the other strand |
| 6 | score | `.` (the program gives no score) |
| 7 | strand | `+` or `-`. `.` means the element was found on both strands (see `--dedup`) |
| 8 | phase | `.` (not used) |
| 9 | attributes | `key=value` pairs, separated by `;` |

With the labelled TIR patterns, one element is a block of 6 lines. The block ends with `###`.
The features:

| type (column 3) | what it is | position in the example |
|---|---|---|
| `repeat_region` | the whole hit, from the first TSD to the last TSD. It has the bare name as ID (`DTC1`) | 2874799-2878133 |
| `target_site_duplication` (TSD) | the host DNA that is copied on both sides of the element. `lTSD` is the left one, `rTSD` the right one | 2874799-2874801 and 2878131-2878133 (`TAA`) |
| `CACTA_TIR_transposon` | the element without its TSDs. This name comes from `%@element` in the pattern | 2874802-2878130 |
| `five_prime_terminal_inverted_repeat` (TIR) | the repeat at the left end (`lTIR`). TIR = terminal inverted repeat | 2874802-2874813 |
| `three_prime_terminal_inverted_repeat` (TIR) | the repeat at the right end (`rTIR`). It is the reverse complement of the left one | 2878119-2878130 |

In column 9:

* `ID` is the name of the feature. `Parent` says which feature it belongs to.
  The TIRs belong to the element (`DTC1.te`). The TSDs belong to `DTC1`.
* `Name` is the same for all lines of one element. Names are unique in one run.
* `Classification`, `Method`, `Sequence_ontology` come from the pattern label
  (`Sequence_ontology` is the SO number of the type).
* On the element line (`.te`): `TSD=TAA_TAA` gives the two TSDs, and
  `TIR=CACTACAAGAAA_TTTCTTGTAGTG` gives the two TIRs (left_right, forward-strand letters).

## 5. Look at the results

### 5.1 Count and measure

Every element has one `repeat_region` line. Count them:

```
$ awk '$3=="repeat_region"' dtc_II.gff3 | wc -l
3272
```

Count the lines of each type:

```
$ grep -v '^#' dtc_II.gff3 | cut -f3 | sort | uniq -c
   3272 CACTA_TIR_transposon
   3272 five_prime_terminal_inverted_repeat
   3272 repeat_region
   6544 target_site_duplication
   3272 three_prime_terminal_inverted_repeat
```

Elements on each chromosome (the first 4 lines):

```
$ awk '$3=="repeat_region"{n[$1]++} END{for(c in n) print c, n[c]}' dtc_II.gff3 | sort | head -4
II_C01 345
II_C02 392
II_C03 448
II_C04 426
```

Number, shortest, longest and mean length of the element (without TSDs):

```
$ awk '$3=="CACTA_TIR_transposon"{l=$5-$4+1; s+=l; if(min==""||l<min)min=l; if(l>max)max=l; n++} END{print n, min, max, s/n}' dtc_II.gff3
3272 538 14991 6935.45
```

Strand (column 7) of the `repeat_region` lines:

```
$ awk '$3=="repeat_region"' dtc_II.gff3 | cut -f7 | sort | uniq -c
     28 +
     28 -
   3216 .
```

`.` means found on both strands. `+` or `-` means only one strand had this span.

### 5.2 BED for IGV or another genome browser

BED is a simpler format: one line for each element. Let the program write it:

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format bed6 --name-prefix DTC \
      --input II_C09.fa --output dtc_C09.bed \
      scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
$ head -3 dtc_C09.bed
II_C09	854081	860374	DTC1	0	.
II_C09	2427868	2430086	DTC2	0	.
II_C09	2466892	2469357	DTC3	0	.
```

BED6 columns: record, start (counted from 0), end, Name, score (0), strand.
Note: BED starts at 0 and GFF3 starts at 1. So `854081` in BED is `854082` in GFF3.
Most genome browsers (IGV for example) can open BED and GFF3 files.

`--format bed12` also writes the parts (TSD, TIR) as blocks:

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format bed12 --name-prefix DTC \
      --input II_C09.fa scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat | head -2
II_C09	854081	860374	DTC1	0	.	854084	860371	0	4	3,12,12,3,	0,3,6278,6290,
II_C09	2427868	2430086	DTC2	0	.	2427871	2430083	0	4	3,12,12,3,	0,3,2203,2215,
```

You can also make a BED from an existing GFF3 with awk (here from the whole-genome
file of section 4.4):

```
$ awk -F'\t' '$3=="repeat_region"{split($9,a,";"); sub("ID=","",a[1]); print $1"\t"$4-1"\t"$5"\t"a[1]"\t0\t"$7}' dtc_II.gff3 | head -3
II_C01	2874798	2878133	DTC1	0	.
II_C01	3398539	3401290	DTC2	0	.
II_C01	4397540	4400232	DTC3	0	.
```

### 5.3 JSON lines for Python

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format jsonl --name-prefix DTC \
      --input II_C09.fa --output dtc_C09.jsonl \
      scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
$ head -1 dtc_C09.jsonl | cut -c1-200
{"seq":"II_C09","strand":".","start":854082,"end":860374,"name":"DTC1","type":"CACTA_TIR_transposon","units":[{"line":3,"start":854082,"end":854084,"text":"TAT","label":"target_site_duplication"},...
```

Each line is one hit. `units` are the parts of the hit, in pattern order.
Python with only the standard library (`read_jsonl.py`):

```python
import json, sys
n = 0
total = 0
for line in open(sys.argv[1]):
    hit = json.loads(line)
    n += 1
    total += hit["end"] - hit["start"] + 1
print(n, "hits; mean length", round(total / n))
tsd = [u["text"] for u in hit["units"] if u.get("label") == "target_site_duplication"]
print("last hit", hit["seq"], hit["start"], hit["end"], hit["name"], "TSD", tsd)
```

```
$ python read_jsonl.py dtc_C09.jsonl
396 hits; mean length 6935
last hit II_C09 7202370 7211624 DTC396 TSD ['TAA', 'TAA']
```

(A unit without a label, such as the middle part, has no `label` key. That is why the
code uses `u.get("label")`.) If you have pandas: `pd.read_json(f, lines=True)`.

### 5.4 Join the results of several runs: `--merge`

Run two patterns in two separate runs. Both runs start the Names at 1:

```
$ L=scan_for_matches_rs/tir_scan_patterns/labelled
$ scan_for_matches -t 4 -c --dedup --strict-n --format gff3 --input II_C09.fa --output a.gff3 $L/00_DTC_published_2025.pat
$ scan_for_matches -t 4 -c --dedup --strict-n --format gff3 --input II_C09.fa --output b.gff3 $L/09_DTC_CACTA_seed10.pat
$ grep -c -P '\trepeat_region\t' a.gff3 b.gff3
a.gff3:396
b.gff3:452
```

Both files have an element `DTC1`. They are not the same element. `--merge` joins the
files, removes the elements with the same span (same record, start and end), and gives
new unique Names. The files are in priority order: the element of the first file is kept.

```
$ scan_for_matches --merge --output ab.gff3 a.gff3 b.gff3
merge: 848 elements read, 389 duplicates removed, 459 written
$ grep -v '^#' ab.gff3 | awk '$3=="repeat_region"' | cut -f1,3,4,5,9 | head -3
II_C09	repeat_region	212663	236410	ID=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657
II_C09	repeat_region	854082	860374	ID=DTC2;Name=DTC2;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657;Merged=b.gff3:DTC2
II_C09	repeat_region	2427869	2430086	ID=DTC3;Name=DTC3;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657;Merged=b.gff3:DTC3
```

Why is the element at 854082 now `DTC2`? In `a.gff3` it was `DTC1` (the first
element of `a.gff3` by position). `--merge` sorts all elements by position
and numbers them again, so **old Names change**. Here another element, at 212663 (only in
`b.gff3`), comes first and becomes `DTC1`.
`Merged=b.gff3:DTC2` lists the removed element by its OLD name: the
element `DTC2` in `b.gff3` had the same span as the kept one.

Useful merge options:

* `--overlap 0.9` also counts two elements as the same when each overlaps the
  other by at least 90 %:

  ```
  $ scan_for_matches --merge --overlap 0.9 --output ab90.gff3 a.gff3 b.gff3
  merge: 848 elements read, 391 duplicates removed, 457 written
  ```

  (391 duplicates removed instead of 389.)
* `--name-prefix P` gives all elements the prefix `P` (for example `--name-prefix II_CACTA_` gives `II_CACTA_1`).
* Give the files one by one. Do not join them with `cat`: the IDs would repeat.
  The program stops with an error:

```
$ cat a.gff3 b.gff3 | scan_for_matches --merge - > /dev/null
merge: -:2779: ID 'DTC1' is used twice (line 4); IDs must be unique in one file
```

### 5.5 Compare with an old annotation (optional)

`II.CC.TE.gff3` is just an example of an old annotation of this genome that you
may have. It has 3 244 `intact_CACTA` lines. We compare them with our 3 272 elements.
`bedtools` is a program for regions. On this cluster it is a module (check with
`module avail bedtools`). `bedtools2/2.30.0` is the module name here and can be
different on your system. Nothing is installed: we only load it.

```
$ module load bedtools2/2.30.0
$ awk -F'\t' '$3=="repeat_region"{split($9,a,";"); sub("ID=","",a[1]); print $1"\t"$4-1"\t"$5"\t"a[1]}' dtc_II.gff3 | sort -k1,1 -k2,2n > new.bed
$ awk -F'\t' '$3=="intact_CACTA"{print $1"\t"$4-1"\t"$5"\t"$3}' II.CC.TE.gff3 | sort -k1,1 -k2,2n > old.bed
$ wc -l new.bed old.bed
  3272 new.bed
  3244 old.bed
  6516 total
$ bedtools intersect -a new.bed -b old.bed -f 0.9 -r -u | wc -l     # same element: 90 % overlap of each
3255
$ bedtools intersect -a new.bed -b old.bed -f 0.9 -r -v | wc -l     # only in the new run
17
```

The flags of `bedtools intersect`:

| flag | meaning |
|---|---|
| `-a`, `-b` | the new regions (a) and the old regions (b) |
| `-f 0.9 -r` | keep a pair only when the overlap is at least 90 % of the length of a AND of b (`-r` = reciprocal) |
| `-u` | write each a-region once if it has a partner in b |
| `-v` | write only the a-regions that have NO partner in b |

So 3 255 of our 3 272 elements have a partner in the old file, and 17 do not.

## 6. Run many patterns at once

Give several pattern files, and say `--format`. The FASTA file is read once.
This runs all 22 ready-made TIR patterns on the whole genome. It is a longer job,
so use Slurm (`run_all22.sh`):

```bash
#!/bin/bash
#SBATCH --job-name=sfm_all22
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem=24G
#SBATCH --output=sfm_all22-%j.out

set -euo pipefail
export PATH=/path/to/folder_with_binary:${PATH}

time scan_for_matches -t 20 -c --dedup --strict-n --format gff3 --progress \
    --input genome.fa.gz --output tir_all.gff3 \
    scan_for_matches_rs/tir_scan_patterns/labelled/*.pat
```

The end of the Slurm output:

```
scan_for_matches: [0:01:00] 1 records (52.4 Mb) done, 41092 hits; now II_C02 (63.4 Mb), pattern 7/22, strand -, at 21.0 Mb
...
scan_for_matches: [0:05:58] 289 records (563.4 Mb) done, 440329 hits; finished

real	5m58.453s
```

(The shell expands `*.pat` to the 22 files.) The output file is large:
440 329 elements, 518 MB.

About the Names: each pattern file has its own prefix in its `%@element` line
(`DTC`, `DTH`, `DTT`, ...). Patterns with the same prefix share one counter.
So all Names are unique in this file. The elements per prefix:

```
$ awk -F'\t' '!/^#/ && $9 !~ /Parent=/{split($9,a,";"); sub("ID=","",a[1]); sub(/[0-9]+$/,"",a[1]); n[a[1]]++} END{for(k in n) print k, n[k]}' tir_all.gff3 | sort
DMM 1782
DTA 135967
DTB 13439
DTC 19645
DTE 7366
DTH 164419
DTM 6711
DTP 20226
DTR 8553
DTT 62221
```

The same element can be found by two patterns. In one run, each record is searched by
the patterns in the order of the command line, and the hits are written in that order.
So for the same span, the pattern given first on the command line comes first in the
GFF3. `--merge` keeps the first one it reads (README.md, `--merge`: "the element of the
earlier file, then the earlier line, is kept"). Merge the one file:

```
$ scan_for_matches --merge --output tir_merged.gff3 tir_all.gff3
merge: 440329 elements read, 29679 duplicates removed, 410650 written
$ grep -m2 'Merged=' tir_merged.gff3 | cut -f1,3,4,5,9 | cut -c1-300
II_C01	repeat_region	37068	38350	ID=DTT6;Name=DTT6;Classification=TIR/DTT;Method=structural;Sequence_ontology=SO:0000657;Merged=tir_all.gff3:DTC1539
II_C01	repeat_region	58913	59427	ID=DTE1;Name=DTE1;Classification=TIR/DTE;Method=structural;Sequence_ontology=SO:0000657;Merged=tir_all.gff3:DTP2
```

Check the first line. In `tir_all.gff3` the same span 37068-38350 on `II_C01` is
there twice:

```
$ awk -F'\t' '$1=="II_C01" && $4==37068 && $3=="repeat_region"{split($9,a,";"); print NR, $4,$5,a[1]}' tir_all.gff3
44097 37068 38350 ID=DTT5955
284099 37068 38350 ID=DTC1539
```

`DTT5955` comes from `01_DTT_seed12.pat` and is written first (line 44097).
`DTC1539` comes from `09_DTC_...` and is written later (line 284099). The shell
sorts `*.pat` in name order, so `01_` is before `09_`. `--merge` keeps the DTT
element (new name `DTT6`) and `Merged=` has the OLD name of the removed one
(`DTC1539`). The merge took 6 seconds.

Three patterns on one chromosome (the `# scan_for_matches` line of the GFF3 lists
each pattern):

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format gff3 --input II_C09.fa --output C09_three.gff3 \
      $L/00_DTC_published_2025.pat $L/09_DTC_CACTA_seed10.pat $L/09_DTC_CACTG_seed10.pat
```

Name rules, in short:

* No option: the prefix is the `Name=` of the `%@element` line, else `sfm`.
* `--name-prefix P` sets the prefix for all patterns of the run.
* `--name-start N` sets the first number.
* Separate runs repeat the Names. Give each run its own `--name-prefix`, or join the runs with `--merge`.
* `-m N` is one limit for the whole run: at most N hits of all patterns together.

## 7. Write your own pattern

Put the pattern in a text file. One pattern per file.
Always check it with `--explain` and `--lint`, and test it on a tiny FASTA with
the target planted in it. The full language is in [PATTERNS.md](PATTERNS.md).

### 7.1 A motif with ambiguity codes

`W` means A or T. This pattern finds TATA-box-like words:

```
$ echo 'TATAWAW' > tata.pat
$ printf '>demo1\nGGGCTATAAAAGGCCCTATATATGGCCGGGATTACAGGACCTTAAGCTAGCTAGCCGATTACAGGG\n>demo2\nAAAAAAAAAAAAAAAAAAAA\n' > demo.fa
$ scan_for_matches tata.pat < demo.fa
>demo1:[5,11]
TATAAAA 
>demo1:[17,23]
TATATAT 
```

### 7.2 A direct repeat

`p1=6...8` takes 6 to 8 letters. `10...30` is a gap. The last `p1` is the same letters again.

| pattern | input (`demo1`, part) | result |
|---|---|---|
| `p1=6...8 10...30 p1` | `GATTAC AGGACCTTAAGCTAGCTAGCC GATTAC` | hit at 30-62 |

```
$ echo 'p1=6...8 10...30 p1' > direct.pat
$ scan_for_matches direct.pat < demo.fa
>demo1:[30,62]
GATTAC AGGACCTTAAGCTAGCTAGCC GATTAC 
```

### 7.3 A hairpin where one letter may differ

`[1,0,0]` means: 1 mismatch, 0 deletions, 0 insertions. This stem has 1 mismatch
(`CTGTAAGC` is not exactly the reverse complement of `GATTACAG`, which is `CTGTAATC`).

| pattern | result on `stem.fa` |
|---|---|
| `p1=8...8 4...10 ~p1` (exact) | no hit |
| `p1=8...8 4...10 ~p1[1,0,0]` | hit at 5-24 |

```
$ printf '>stem\nAAAAGATTACAGTTTTCTGTAAGCAAAA\n' > stem.fa
$ echo 'p1=8...8 4...10 ~p1' > hp0.pat
$ echo 'p1=8...8 4...10 ~p1[1,0,0]' > hp1.pat
$ scan_for_matches hp0.pat < stem.fa
$ scan_for_matches hp1.pat < stem.fa
>stem:[5,24]
GATTACAG TTTT CTGTAAGC 
```

**Trap: the shortest stem wins.** `p1=5...8` uses the shortest length that works
(`--explain` says so). With a mismatch allowed, a short false stem is often found
first, and the search goes on after it, so the real hairpin is lost:

| pattern | result on `AAAAACCGATTACTTTTGTAAACGGAAAA` |
|---|---|
| `p1=5...8 3...6 ~p1[1,0,0]` | `[2,17]`: `AAAAC CGATTA CTTTT` (false: a 5-letter stem with 1 mismatch) |
| `p1=8...8 3...6 ~p1[1,0,0]` | `[6,25]`: `CCGATTAC TTTT GTAAACGG` (the real hairpin) |

Use a fixed stem length, or run one pattern for each length.

### 7.4 A small TIR transposon with labels

Labels (`%@`) give the parts a name in GFF3. This is a small CACTA-like element with
a gap of 20 to 100 letters. Put it in `mini.pat`:

```
%@element mini_TIR_transposon Name=mini
p1=3...3               %@ target_site_duplication
CACTA p2=7...7         %@ five_prime_terminal_inverted_repeat
20...100
~p2 TAGTG              %@ three_prime_terminal_inverted_repeat
p1                     %@ target_site_duplication
```

```
$ printf '>mini\nGGCATCGGATGCACTACAAGAAAATTGGCCTTAACCGGATCGCTAGCTAGGGATTTCTTGTAGTGATGGCGGTCATGGC\n' > mini.fa
$ scan_for_matches -c --dedup --strict-n --format gff3 mini.pat < mini.fa | cut -c1-260
##gff-version 3
# scan_for_matches 0.4.3 pattern=mini.pat
##sequence-region mini 1 79
mini	scan_for_matches	repeat_region	9	68	.	.	.	ID=mini1;Name=mini1;Method=structural;Sequence_ontology=SO:0000657
mini	scan_for_matches	target_site_duplication	9	11	.	.	.	ID=mini1.lTSD;Parent=mini1;Name=mini1;Method=structural;Sequence_ontology=SO:0000434
mini	scan_for_matches	mini_TIR_transposon	12	65	.	.	.	ID=mini1.te;Parent=mini1;Name=mini1;Method=structural;TSD=ATG_ATG;TIR=CACTACAAGAAA_TTTCTTGTAGTG
mini	scan_for_matches	five_prime_terminal_inverted_repeat	12	23	.	.	.	ID=mini1.lTIR;Parent=mini1.te;Name=mini1;Method=structural;Sequence_ontology=SO:0000420
mini	scan_for_matches	three_prime_terminal_inverted_repeat	54	65	.	.	.	ID=mini1.rTIR;Parent=mini1.te;Name=mini1;Method=structural;Sequence_ontology=SO:0000421
mini	scan_for_matches	target_site_duplication	66	68	.	.	.	ID=mini1.rTSD;Parent=mini1;Name=mini1;Method=structural;Sequence_ontology=SO:0000434
###
```

Break the TSD (the last `ATG` becomes `TTG`) and the hit goes away:

```
$ printf '>mini_broken\nGGCATCGGATGCACTACAAGAAAATTGGCCTTAACCGGATCGCTAGCTAGGGATTTCTTGTAGTGTTGGCGGTCATGGC\n' | scan_for_matches --format gff3 mini.pat | cut -c1-100
##gff-version 3
# scan_for_matches 0.4.3 pattern=mini.pat
```

### 7.5 A protein pattern

With `-p` the letters are amino acids. `any(ST)` is S or T. `X` is any amino acid.
`notany(P)` is any letter except P. You cannot use `-c` with `-p`.

```
$ echo 'M any(ST) X X notany(P) K' > p.pat
$ printf '>s\nGMSAAGKLL\n' > s.fa
$ scan_for_matches -p p.pat < s.fa
>s:[2,7]
M S A A G K 
```

### 7.6 More help

* The full language (ranges, names, alternatives, rule sets for RNA, weight
  matrices, length limits, common mistakes): [PATTERNS.md](PATTERNS.md).
* A JASPAR matrix into a weight unit: `python3 -I tools/jaspar2sfm.py MATRIX.jaspar --pvalue 1e-4`
  ([README.md](README.md#jaspar-matrix-to-weight-unit-toolsjaspar2sfmpy)).
* An AI assistant can write patterns for you: see
  [`skills/sfm-patterns/SKILL.md`](skills/sfm-patterns/SKILL.md).
  Check what it writes with `--explain` and `--lint`.

## 8. All options

All options of `scan_for_matches --help`:

| option | meaning | when to use it |
|---|---|---|
| `-c` | also search the opposite strand | DNA, when the target can be on either strand. Not with `-p` |
| `-p` | the sequences are protein | protein FASTA |
| `-n N` | stop (exit status 1) after N sequences with no hit in total (not only in a row) | rarely; stop a search that finds nothing |
| `-m N` | write at most N hits in total (all patterns together) | a quick test: `-m 1` |
| `-o N` | a leftover from the old program: show overlapping hits. Leave it out if you do not need it. If you use it, write `-o 1`: a number must follow, but it does nothing | when hits can overlap or be nested |
| `-i FILE` | file with the ids of records to skip | skip the organelle or contig records |
| `-t N` | N threads. Default 1. The output is the same for any N | every big run. Ask Slurm for N cores too |
| `--input FILE` | read FASTA from FILE (`-` is stdin) | instead of `< file` |
| `--output FILE` | write hits to FILE (`-` is stdout) | instead of `> file` |
| `--progress` | every minute, write a progress line to stderr | long runs; under Slurm it is in the `.out` file |
| `-h`, `--help` | show the help | |
| `--format F` | `gff3`, `bed6`, `bed12` or `jsonl` instead of the original output | always, for genome work. Needed for several patterns |
| `--name-prefix P` | name of each hit = P + number | to give each run its own names |
| `--name-start N` | first number (default 1) | to continue the numbers |
| `--type T` | column 3 of the element (default: the `%@element` type, else `sequence_motif`) | to give a pattern without labels a better type |
| `--strict-n` | an `N` (or other non-ACGT letter) caught by a name counts as 1 mismatch when the name is used again | every genome with assembly gaps. DNA only |
| `--align` | `~p1[m,d,i]` with deletions or insertions tries every alignment length | the `indel/` patterns, or any `~pN[m,d,i]` with `d` or `i` more than 0 |
| `--dedup` | with `-c` (and `--format`): the same span on both strands is written once, strand `.` | TIR patterns and hairpins with `-c` |
| `--lint` | check the pattern for traps. No input is read | before every long run |
| `--explain` | describe the pattern in plain words. No input is read | when you do not trust a pattern |
| `--merge` | join GFF3 files (its own options below) | after many runs |

Options of `--merge` (`scan_for_matches --merge --help`):

| option | meaning | when to use it |
|---|---|---|
| `--overlap F` | also the same element when the overlap is at least F (0 to 1) of the length of each of the two | patterns whose spans differ a little, for example `--overlap 0.9` |
| `--name-prefix P` | prefix for all new Names | one clean set of names |
| `--name-start N` | first number | |
| `--output FILE` | write to FILE instead of stdout | |

Notes:

* Long options can be written `--format gff3` or `--format=gff3`.
* `--dedup`, `--name-prefix`, `--name-start` and `--type` need `--format`.
* `-i` and `-o` are options of the original program (ids to skip, overlapping
  hits). They are **not** input and output files.

## 9. Speed and memory

The numbers below come from [README.md](README.md) (sections Quick start,
Several patterns in one run, Speed) and
[`tir_scan_patterns/README.md`](tir_scan_patterns/README.md).

* **Threads (`-t N`).** More threads, less time. One pattern
  (02_DTA_short5to7, 81 514 hits) on a 1.0 Gb genome took 349 s with `-t 1`,
  184 s with `-t 2`, 125 s with `-t 4` and 110 s with `-t 8`; the memory
  grew from 0.6 GB to 4.2 GB (README.md, Speed; Apple M4). That machine is
  limited by memory speed, so more than 4 threads gained little there.
  Ask Slurm for the same number of cores (`-c N`).
* **bgzip input.** With `-t N`, the blocks of a bgzip file are decompressed by
  N threads. Plain gzip uses one thread (README.md, Usage). You can run on
  `genome.fa.gz` directly. To make a bgzip file: `bgzip genome.fa`
  (part of htslib/samtools).
* **Many patterns, one run.** 22 TIR patterns on a 1.0 Gb *B. napus* genome,
  `-t 20`: 648 s in one run, 762 s as 22 runs (README.md,
  "Several patterns in one run"). About 11 min is the figure in the Quick start.
* **Gap skipping.** A wide range followed by an exact word, for example
  `50...30000 ~p1`, skips the gap lengths that cannot work: up to 1000 times
  faster on TIR patterns (README.md, Quick start). Keep the unit right after
  a wide gap exact (see the next point).
* **The `slow` note.** `--lint` says `slow` when the unit right after a wide gap
  (`500...15000`) cannot be used to skip gap lengths. Then the program tries every gap
  length one by one. Fix it in this order:
  1. Put an exact word or an exact `~pN` right after the gap.
  2. If the unit after the gap needs mismatches, write them on `~pN[m,0,0]` and
     make `pN` at least 5 x (m + 1) letters long (10 letters for 1 mismatch),
     or put an exact word of 5 letters after it.
  3. With deletions or insertions (`~pN[m,d,i]`), `pN` needs at least
     4 x (m + d + i + 1) letters (20 for `[2,1,1]`). Exact letters after it do not help.
  4. Or make the gap narrower.

  (Numbers: README.md, 0.4.2 and 0.4.3 lines; PATTERNS.md, mistake 8.)
  Real output, a slow pattern and a fast one:

  ```
  $ echo 'p1=3...3 CACTA p2=7...7 500...15000 TAGTG[1,0,0] ~p2 p1' > w9.pat
  $ scan_for_matches --lint --strict-n w9.pat
  w9.pat:1: note: slow: each gap length of `500...15000` is tried (`TAGTG[1,0,0]` after it cannot be used to skip gap lengths)
  0 errors, 0 warnings, 1 notes
  $ echo 'p1=3...3 CACTA p2=7...7 500...15000 ~p2[1,0,0] TAGTG p1' > w12.pat
  $ scan_for_matches --lint --strict-n w12.pat
  w12.pat: note: the pattern reads the same on both strands: with -c, every element is found twice; use --dedup with --format
  0 errors, 0 warnings, 1 notes
  ```

  The second pattern has the mismatch on `~p2[1,0,0]` and 5 exact letters (`TAGTG`)
  after it, so there is no `slow` note. (Without `TAGTG`, `p2` of 7 letters is too short:
  the note comes back.)
* **`mm2/` and `indel/` patterns.** `tir_scan_patterns/mm2/` allow 2 mismatches
  in the right TIR. On about 100 Mb of *B. napus* (`-t 40`) they took 4 s,
  27 s, 23 s and 34 s (from 3 s, 14 s, 10 s and 4 s for the exact versions).
  `tir_scan_patterns/indel/` (use `--align`) took 6 min, 2.7 h, 1.8 h and
  5 min on the same data (tir_scan_patterns/README.md). Run indel patterns with
  `sbatch`, and check the hits: the random share is higher.
* **`--progress`** shows how far a long run is (records and Mb done, hits, the
  record searched now).
* **Memory.** For long gaps the program builds a [5-mer index](#words-used-in-this-guide) of the record: about
  4 bytes per base (a 74 Mb chromosome: about 300 MB). Maize (2.3 Gb), `-t 20`,
  pattern 00_DTC_published_2025: 2.5 GB at maximum; 02_DTA_short5to7: 5.1 GB
  (README.md, Speed). A pattern that gives millions of hits makes a big output
  file: the 22-pattern run in section 6 wrote 518 MB. To measure your own run,
  see section 4.4.
* **Start small.** Test on one chromosome or a cut before the full run.

## 10. Common problems (FAQ)

**I used `-o out.txt` and there is no output file.**
`-o` is NOT the output file (see the `-o 1` row in section 8). The program reads
`out.txt` as the value of `-o`, writes the hits to the screen, and makes no file.
Use `--output out.txt`.

```
$ scan_for_matches -o out.txt hp1.pat < stem.fa
>stem:[5,24]
GATTACAG TTTT CTGTAAGC 
$ ls out.txt
ls: cannot access 'out.txt': No such file or directory
```

**`-o -c` does not search the other strand.**
`-o` takes the next word as its value, so `-c` is eaten. Write `-o 1 -c`.
`--lint` warns you:

```
hp1.pat: warning: -o takes a value: `-c` is read as the value of -o, not as an option, so the reverse strand is not searched; write -o 1
```

**Every element is found twice with `-c`.**
A pattern that reads the same on both strands (hairpins, TIR patterns) gives each
element twice. Add `--dedup` (with `--format`). Strand becomes `.`.

| command | result on `stem.fa` (BED6) |
|---|---|
| `scan_for_matches -c --format bed6 hp1.pat` | `stem 4 24 sfm1 0 +` and `stem 4 24 sfm2 0 -` |
| `scan_for_matches -c --dedup --format bed6 hp1.pat` | `stem 4 24 sfm1 0 .` |

`--dedup` without `--format` is an error:
`scan_for_matches: --dedup, --name-prefix, --name-start and --type need --format`.

**False hits where the assembly has gaps (`N`).**
Add `--strict-n`. Without it, an `N` that a name caught matches any letter when the
name is used again, so a TSD made of `N` gives a false element. With it, each such
letter is 1 mismatch. See [README.md](README.md#assembly-gaps---strict-n). `--lint`
reminds you with a note until you add it.

**Indel patterns (`~p2[2,1,1]`) miss elements.**
Add `--align`. Without it the matcher takes one alignment length and loses some
hits. `--lint` says so, and says nothing when `--align` is there:

```
$ scan_for_matches --lint -c --dedup --strict-n --format gff3 scan_for_matches_rs/tir_scan_patterns/indel/labelled/04_DTE_seed20_indel.pat
...:3: note: `p1=8...9`: the shortest length that works is used, not the longest
...:6: note: `~p2[2,1,1]` has deletions and insertions: the matcher takes one greedy alignment and can miss some matches (add --align to try every length)
0 errors, 0 warnings, 2 notes
$ scan_for_matches --lint -c --dedup --strict-n --align --format gff3 scan_for_matches_rs/tir_scan_patterns/indel/labelled/04_DTE_seed20_indel.pat
...:3: note: `p1=8...9`: the shortest length that works is used, not the longest
0 errors, 0 warnings, 1 notes
```

**Names collide between runs (two files have `DTC1`).**
Each run starts at 1. Give each run its own `--name-prefix`, or join the files
with `--merge` (section 5.4). Do not join GFF3 files with `cat`.

**"several pattern files need --format".**
More than one pattern file needs `--format gff3` (or another format). The old output
cannot say which pattern made a hit. Exit status is 2.

**"-c cannot be used with -p".**
Protein has no reverse strand. Remove `-c`. Exit status is 2.

**A pattern file with two patterns.**
All text in one file is ONE pattern. Use one file per pattern.

**No hits and no error.**
Check the pattern with `--explain` (is this what you meant?) and test it on a tiny
FASTA with the target planted in it. Check that the FASTA ids are not in an `-i` file.
Remember: a pattern needs the letters in the order of the forward strand.

**`-n 1` stopped my run.**
`-n N` stops (exit status 1) after N sequences with no hit. Leave `-n` out.

**`--lint` says `slow`.**
See the `slow` note in section 9 (four fixes, in order).

**How do I watch a running job?**
Two commands. `squeue` shows the queue: `PD` = waiting, `R` = running. When the job
is not in the list any more, it has ended. `tail` shows the last lines of the Slurm
output file; with `--progress` there is a new line every minute.

```
$ squeue -u $USER -n sfm_mem
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
            580692       CPU  sfm_mem    pxzhe PD       0:00      1 (None)
$ grep 'scan_for_matches:' sfm_all22-580689.out | tail -n 2 | cut -c1-170
scan_for_matches: [0:05:00] 8 records (467.2 Mb) done, 367678 hits; now II_C09 (65.9 Mb), pattern 4/22, strand -, at 65.9 Mb
scan_for_matches: [0:05:58] 289 records (563.4 Mb) done, 440329 hits; finished
```

(Use `tail -n 2 sfm_dtc-<jobid>.out` the same way; `<jobid>` is the number that
`sbatch` printed. To stop a job: `scancel <jobid>`.)

**Where to get help.**
`scan_for_matches --help`, [PATTERNS.md](PATTERNS.md), [README.md](README.md),
[CHANGELOG.md](CHANGELOG.md), and the 22 ready-made patterns in
[`tir_scan_patterns/`](tir_scan_patterns/README.md).
