[English](USER_GUIDE.md) | 繁體中文

# 使用者指南：給初學者的 scan_for_matches

這份指南寫給從沒用過 `scan_for_matches` 的生物學家。
你會學到怎麼開始、怎麼跑一個真實的基因體，以及怎麼看懂結果。
這裡每個指令都實際執行過，顯示的輸出也是真實的輸出
（太長的行用 `...` 截斷）。

## 目錄

0. [本指南用到的詞](#本指南用到的詞)
1. [這個工具做什麼](#1-這個工具做什麼)
2. [取得程式](#2-取得程式)
3. [1 分鐘完成第一次執行](#3-1-分鐘完成第一次執行)
4. [真實基因體操作範例（B. oleracea）](#4-真實基因體操作範例b-oleracea)
5. [查看結果](#5-查看結果)
6. [一次執行多個 pattern](#6-一次執行多個-pattern)
7. [自己寫 pattern](#7-自己寫-pattern)
8. [所有選項](#8-所有選項)
9. [速度與記憶體](#9-速度與記憶體)
10. [常見問題（FAQ）](#10-常見問題faq)

## 本指南用到的詞

難懂的詞在這裡說明。點內文中的詞，就會回到這張表。

| 詞 | 意思 |
|---|---|
| 轉位子（transposon） | 一段能在基因體中複製或移動自己的 DNA |
| TIR 轉位子（TIR transposon） | 兩端各有一個 TIR 的轉位子（一種 DNA 轉位子） |
| TIR | 末端反向重複（terminal inverted repeat）。元件的兩端有相同的字母，其中一端是反向互補 |
| TSD | 目標位點重複（target site duplication）。元件跳入時，一小段宿主 DNA 被複製，元件兩側各有一份 |
| CACTA、DTC | CACTA 是植物 TIR 轉位子的一個家族。它的兩端以 `CACTA` 開頭。DTC 是它的簡碼。其他簡碼：DTA、DTH、DTM、DTT ... 是其他家族 |
| 組裝缺口（assembly gap） | 基因體中一串連續的 `N`：組裝在那裡沒有真實序列 |
| FASTA | 存放序列的文字檔案。以 `>` 開頭的那一行是名稱；後面的行是字母 |
| GFF3 | 一種特徵（feature）的文字表格，有 9 個以 tab 分隔的欄位（第 4.5 節） |
| BED | 較簡單的表格：每個區域一行，起點從 0 開始算（第 5.2 節） |
| JSONL | JSON lines：每行一個 JSON 物件，每個 hit 一行（第 5.3 節） |
| PATH | shell 尋找程式的資料夾清單。`export PATH=/my/folder:$PATH` 會加入一個資料夾 |
| stdin、stdout | 程式的「輸入」與「輸出」資料流。`< file` 把檔案送進 stdin。`> file` 把 stdout 送到檔案 |
| 結束狀態（exit status） | 程式結束時給出的數字。`0` 代表沒有問題。用 `echo $?` 查看 |
| 執行緒、核心（thread、core） | 核心是一個 CPU 單位。執行緒是在核心上工作的一個工作者。在 N 個核心上用 `-t N` 個執行緒 |
| Slurm、sbatch | 叢集的工作管理程式。`sbatch script.sh` 把腳本放進佇列 |
| 登入節點（login node） | 你登入的那台電腦。不要在上面跑大型工作 |
| `#SBATCH` 行 | 腳本最上方的幾行。它們告訴 Slurm 你需要什麼（核心、記憶體、名稱） |
| bgzip | 由許多小區塊組成的 gzip。程式可以讀它，也可以用多個執行緒讀它 |
| IGV | 一個基因體瀏覽器程式。它把 BED 和 GFF3 檔案顯示在基因體上 |
| Sequence Ontology（SO） | 序列特徵的標準名稱與編號清單，例如 `SO:0000657` 是 `repeat_region` |
| 貪婪比對（greedy alignment） | 比對器採用第一個符合的比對，不再嘗試其他比對。`--align` 會嘗試其他比對 |
| 間隔跳躍（gap skipping） | 對於像 `500...15000` 這樣的寬間隔，程式只嘗試下一個單元可能符合的間隔長度。這樣快很多 |
| Mb、Gb | 百萬鹼基、十億鹼基（序列長度），不是 MB、GB（檔案或記憶體大小） |
| 5-mer 索引（5-mer index） | 一張表，記錄一筆紀錄中每個 5 字母字詞的位置。程式用它來做間隔跳躍 |

## 1. 這個工具做什麼

`scan_for_matches` 讀取 [FASTA](#本指南用到的詞) 格式的序列（DNA 或蛋白質）。
它尋找符合你所寫 pattern 的位置。
一個 pattern 可以描述一個短的模體（motif）、一個髮夾結構、一個正向重複，或一整個
[TIR 轉位子](#本指南用到的詞)（兩個末端重複，也就是 [TIR](#本指南用到的詞)，以及兩個
目標位點重複，也就是 [TSD](#本指南用到的詞)）。
結果是一份 hit 清單。你可以得到 [GFF3、BED、JSONL](#本指南用到的詞)
（JSON lines：每個 hit 一個 JSON 物件），或原始程式的純文字輸出。

它是 Ross Overbeek 一個舊 C 程式的 Rust 移植版本。
pattern 語言和舊程式相同。

## 2. 取得程式

* macOS（Apple silicon）：從 GitHub Releases 頁面下載 `scan_for_matches-<version>-macos-arm64.tar.gz`：
  <https://github.com/xavierzheng/scan_for_matches_rs/releases>
* Linux 和其他系統：自己編譯。你需要 Rust（<https://rustup.rs>）。
  請看 [README.md](README.md#build)。

  ```sh
  git clone https://github.com/xavierzheng/scan_for_matches_rs.git
  cd scan_for_matches_rs
  cargo build --release      # the program: target/release/scan_for_matches
  ```

把程式放到一個在 `PATH` 裡的資料夾，或把它的資料夾加進去：

```sh
export PATH=$PWD/target/release:$PATH
scan_for_matches --help
```

說明文字的第一行是版本。這份指南是用
`scan_for_matches 0.4.3` 執行的。

`git clone` 會建立資料夾 `scan_for_matches_rs/`。裡面有現成的
pattern（`tir_scan_patterns/`）和這份指南。
**第 4 到 6 節的指令都在同一個工作資料夾中執行。** 這個資料夾
裡並排放著你的基因體檔案和資料夾 `scan_for_matches_rs/`：

```
$ ls
genome.fa.gz
scan_for_matches_rs
```

（這裡兩者都是指向真實檔案的連結。用真正的複本也一樣可以。）輸出
檔案會寫到這個工作資料夾。`sbatch` 也要在這個資料夾執行：Slurm
會在你執行 `sbatch` 的資料夾中啟動工作。

## 3. 1 分鐘完成第一次執行

用 `printf` 做一個小的 FASTA 檔案。這個檔案有兩條序列。
`seq1` 含有一次 `GATTACA` 這個字詞。

```
$ printf '>seq1\nCCCCGATTACAAAAAAGTAATCGGGGG\n>seq2\nCCCCCCCCCCCCCCCCCCCCCC\n' > tiny.fa
```

做一個 pattern 檔案，內容是一行 pattern：

```
$ echo 'GATTACA' > motif.pat
```

執行它。FASTA 檔案用 `<` 送進去。hit 顯示在螢幕上。

```
$ scan_for_matches motif.pat < tiny.fa
>seq1:[5,11]
GATTACA 
```

怎麼看：

* `>seq1:[5,11]` 是 hit。它在 `seq1` 中，從第 5 個字母到第 11 個字母。
  第一個字母是第 1 號。
* 下一行是符合的字母。
* `seq2` 沒有 hit，所以不會為它輸出任何東西。

現在換一個更有意義的 pattern：髮夾結構（hairpin）。髮夾結構是一個莖（stem）、一個環（loop），以及
莖的反向互補。這裡莖有 4 到 7 個字母，
環有 3 到 8 個字母。

```
$ echo 'p1=4...7 3...8 ~p1' > hairpin.pat
$ scan_for_matches hairpin.pat < tiny.fa
>seq1:[3,24]
CCGATTA CAAAAAAG TAATCGG 
```

`p1=4...7` 取 4 到 7 個字母，並把它們叫做 `p1`。`~p1` 是 `p1` 的反向
互補。（把字母順序倒過來，再把 A 換成 T、C 換成 G。）
`CCGATTA` 變成 `TAATCGG`。

請程式用簡單的話說明一個 pattern 做什麼。這不會讀取
FASTA 檔案。

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

要求輸出 GFF3，而不是原始輸出：

```
$ scan_for_matches --format gff3 hairpin.pat < tiny.fa
##gff-version 3
# scan_for_matches 0.4.3 pattern=hairpin.pat
##sequence-region seq1 1 27
seq1	scan_for_matches	sequence_motif	3	24	.	+	.	ID=sfm1;Name=sfm1;Method=structural;Sequence_ontology=SO:0001683
###
```

## 4. 真實基因體操作範例（B. oleracea）

範例基因體是 *Brassica oleracea*（檔案 `genome.fa.gz`）。
它有 289 筆紀錄：染色體 `II_C01` 到 `II_C09`，以及一些小的重疊群（contig）。
這個檔案是 [bgzip](#本指南用到的詞) 格式。程式可以直接讀它：你不需要先解壓縮。
請在工作資料夾（第 2 節）執行本節的指令。

我們要找 CACTA 轉位子（DTC，一個 TIR 轉位子家族）。
我們使用現成的 pattern `00_DTC_published_2025.pat`。
全部 22 個現成 pattern 都在
[`tir_scan_patterns/labelled/`](tir_scan_patterns/README.md)（中文說明：[README_zh-TW.md](tir_scan_patterns/README_zh-TW.md)）。

### 4.1 看看這個 pattern

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

* 以 `%` 開頭的行是註解。
* `%@element ...` 為整個元件命名，並給出 Name 的前綴（`DTC`）。
* 行尾的 `%@ type` 為這一行上的單元命名。
  這些名稱會成為 GFF3 中的特徵。
* `p1=3...3` 取任意 3 個字母（TSD，目標位點重複）。最後的
  `p1` 需要再次出現相同的 3 個字母。
* `CACTA[0,0,0] p2=7...7` 是元件的 5' 端：`CACTA`，不允許
  錯配，接著是 7 個叫做 `p2` 的字母。
* `500...15000` 是中間部分，500 到 15000 個任意字母。
* `~p2 TAGTG[0,0,0]` 是 3' 端：`p2` 的反向互補，接著是 `TAGTG`。

### 4.2 檢查它：`--lint` 和 `--explain`

`--lint` 會找出陷阱。給它和正式執行時相同的選項。

```
$ scan_for_matches --lint -c scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat: warning: the pattern reads the same on both strands: with -c, every element is found twice; use --dedup with --format
scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat:7: note: `p1`: an N that p1 caught from the data matches any letter here (false hits at assembly gaps); add --strict-n for genomes with gaps
0 errors, 1 warnings, 1 notes
```

程式告訴我們兩件事。使用 `-c` 時，要加 `--dedup`。基因體有
組裝缺口（一串 `N`）時，要加 `--strict-n`。兩個都加上，再檢查一次：

```
$ scan_for_matches --lint -c --dedup --strict-n --format gff3 scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
0 errors, 0 warnings, 0 notes
```

`--explain` 會列出每個單元（這裡只顯示前幾行）：

```
$ scan_for_matches --explain scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
pattern file: scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat (DNA)
joined pattern: p1=3...3 CACTA[0,0,0] p2=7...7 500...15000 ~p2 TAGTG[0,0,0] p1
unit  line  text          meaning
1     3     p1=3...3      any 3 letters; call them p1   [label: target_site_duplication]
...
hit length: 530 to 15030 letters
```

`hit length` 這一行告訴你，這個 pattern 能找到的最短與最長元件
是多長。

### 4.3 先在「一條」染色體上執行

一定要先在小片段上測試。從基因體中切出染色體 `II_C09`：

```
$ zcat genome.fa.gz | awk '/^>/{p=($1==">II_C09")} p' > II_C09.fa
```

在它上面執行這個 pattern：

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

這條染色體（65.9 Mb）有 396 個元件。執行只花了幾秒。
看起來沒問題？那就跑整個基因體。

### 4.4 整個基因體

每個選項的作用：

| 選項 | 原因 |
|---|---|
| `-t 20` | 使用 20 個[執行緒](#本指南用到的詞)（每個 CPU 核心一個執行緒）。向 Slurm 要求相同數量的核心 |
| `-c` | 搜尋兩條鏈 |
| `--dedup` | 在兩條鏈都找到的同一個元件只寫一次（鏈為 `.`） |
| `--strict-n` | 在[組裝缺口](#本指南用到的詞)（一串 `N`）不產生假 hit |
| `--format gff3` | 輸出 GFF3 |
| `--name-prefix DTC` | 元件命名為 `DTC1`、`DTC2`、... |
| `--progress` | 每分鐘輸出一行進度到 stderr |
| `--input`、`--output` | 輸入的 FASTA 檔案、輸出的 hit 檔案 |

不要在叢集的[登入節點](#本指南用到的詞)上跑大型工作。請用
[Slurm](#本指南用到的詞) 送出。`#SBATCH` 行告訴 Slurm 你需要什麼。
把下面內容放進工作資料夾中的檔案 `run_dtc_genome.sh`：

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

把 `/path/to/folder_with_binary` 改成你程式所在的資料夾（自己編譯的話是
`scan_for_matches_rs/target/release`；不確定時用 `which scan_for_matches` 查）。
有些叢集還要加 `#SBATCH -p`（partition）或 `#SBATCH -A`（帳號），請問叢集管理員。從
工作資料夾送出（腳本中的路徑是相對於這個資料夾）：

```
$ sbatch run_dtc_genome.sh
```

`-c 20` 中的核心數和 `-t` 後面的數字必須相同。
`--mem=20G` 是留有餘裕的估計值。（README.md 記錄這個 pattern 在
2.3 Gb 的玉米基因體上最多用 2.5 GB。）
Slurm 把進度行寫到 `sfm_dtc-<jobid>.out`：

```
scan_for_matches: [0:00:00] 0 records (0.0 Mb) done, 0 hits; started
scan_for_matches: [0:00:06] 289 records (563.4 Mb) done, 3272 hits; finished

real	0m6.650s
user	0m41.698s
sys	0m4.075s
```

這個基因體有 563.4 Mb。這個工作用 20 個執行緒花了 7 秒。它找到 3,272 個元件。

**它實際用了多少記憶體？** 在腳本中把 `/usr/bin/time -v` 放在程式
前面（`/usr/bin/time -v scan_for_matches ...`）。在
`.out` 檔案的結尾找 "Maximum resident set size"。以這個工作為例（用
`/usr/bin/time -v` 再跑一次）：

```
$ grep 'Maximum resident' sfm_mem-580692.out
	Maximum resident set size (kbytes): 1026796
```

這大約是 1.0 GB。對於一般腳本，不要相信 `sacct -j JOBID -o MaxRSS`：
它只顯示小小的 batch 步驟。同一個工作的真實輸出：

```
$ sacct -j 580692 -o JobID,Elapsed,MaxRSS,State
JobID           Elapsed     MaxRSS      State 
------------ ---------- ---------- ---------- 
580692         00:00:40             COMPLETED 
580692.batch   00:00:40     12340K  COMPLETED 
580692.exte+   00:00:40      1048K  COMPLETED 
580692.0       00:00:32    885964K  COMPLETED 
```

`.batch` 這一行（12 MB）對這個程式來說是錯的。`.0` 這一行（886 MB）才是對的：它
之所以存在，是因為這個腳本也用 `srun` 啟動了程式。所以只有當你用
`srun scan_for_matches ...` 啟動程式時，`sacct` 和 `seff` 才是對的。

### 4.5 看懂 GFF3

GFF3 是一個有 9 個欄位、以 tab 分隔的表格。以下是
`dtc_II.gff3` 的前幾行（第 9 欄被截斷）：

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

各欄位：

| 欄 | 名稱 | 在這個檔案中 |
|---|---|---|
| 1 | seqid | FASTA 紀錄：`II_C01` |
| 2 | source | 程式：`scan_for_matches` |
| 3 | type | 特徵的種類（見下表） |
| 4, 5 | start, end | 在正鏈上的位置。第一個字母是 1。start 不會大於 end，另一條鏈上的 hit 也一樣 |
| 6 | score | `.`（程式不給分數） |
| 7 | strand | `+` 或 `-`。`.` 代表元件在兩條鏈都有找到（見 `--dedup`） |
| 8 | phase | `.`（不使用） |
| 9 | attributes | `key=value` 配對，以 `;` 分隔 |

使用有標籤的 TIR pattern 時，一個元件是一個 6 行的區塊。區塊以 `###` 結束。
各個特徵：

| type（第 3 欄） | 它是什麼 | 在範例中的位置 |
|---|---|---|
| `repeat_region` | 整個 hit，從第一個 TSD 到最後一個 TSD。它的 ID 是不加後綴的名稱（`DTC1`） | 2874799-2878133 |
| `target_site_duplication`（TSD） | 在元件兩側被複製的宿主 DNA。`lTSD` 是左邊那個，`rTSD` 是右邊那個 | 2874799-2874801 和 2878131-2878133（`TAA`） |
| `CACTA_TIR_transposon` | 不含 TSD 的元件。這個名稱來自 pattern 中的 `%@element` | 2874802-2878130 |
| `five_prime_terminal_inverted_repeat`（TIR） | 左端的重複（`lTIR`）。TIR = 末端反向重複 | 2874802-2874813 |
| `three_prime_terminal_inverted_repeat`（TIR） | 右端的重複（`rTIR`）。它是左端那個的反向互補 | 2878119-2878130 |

在第 9 欄中：

* `ID` 是特徵的名稱。`Parent` 說明它屬於哪個特徵。
  TIR 屬於元件（`DTC1.te`）。TSD 屬於 `DTC1`。
* 同一個元件的所有行，`Name` 都相同。在一次執行中，Name 不會重複。
* `Classification`、`Method`、`Sequence_ontology` 來自 pattern 的標籤
  （`Sequence_ontology` 是這個 type 的 SO 編號）。
* 在元件那一行（`.te`）：`TSD=TAA_TAA` 給出兩個 TSD，
  `TIR=CACTACAAGAAA_TTTCTTGTAGTG` 給出兩個 TIR（左_右，正鏈字母）。

## 5. 查看結果

### 5.1 計數與測量

每個元件有一行 `repeat_region`。計算它們的數量：

```
$ awk '$3=="repeat_region"' dtc_II.gff3 | wc -l
3272
```

計算每種 type 的行數：

```
$ grep -v '^#' dtc_II.gff3 | cut -f3 | sort | uniq -c
   3272 CACTA_TIR_transposon
   3272 five_prime_terminal_inverted_repeat
   3272 repeat_region
   6544 target_site_duplication
   3272 three_prime_terminal_inverted_repeat
```

每條染色體上的元件數（前 4 行）：

```
$ awk '$3=="repeat_region"{n[$1]++} END{for(c in n) print c, n[c]}' dtc_II.gff3 | sort | head -4
II_C01 345
II_C02 392
II_C03 448
II_C04 426
```

元件（不含 TSD）的數量、最短、最長與平均長度：

```
$ awk '$3=="CACTA_TIR_transposon"{l=$5-$4+1; s+=l; if(min==""||l<min)min=l; if(l>max)max=l; n++} END{print n, min, max, s/n}' dtc_II.gff3
3272 538 14991 6935.45
```

`repeat_region` 行的鏈（第 7 欄）：

```
$ awk '$3=="repeat_region"' dtc_II.gff3 | cut -f7 | sort | uniq -c
     28 +
     28 -
   3216 .
```

`.` 代表在兩條鏈都有找到。`+` 或 `-` 代表只有一條鏈有這個範圍。

### 5.2 給 IGV 或其他基因體瀏覽器用的 BED

BED 是較簡單的格式：每個元件一行。讓程式來寫：

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format bed6 --name-prefix DTC \
      --input II_C09.fa --output dtc_C09.bed \
      scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
$ head -3 dtc_C09.bed
II_C09	854081	860374	DTC1	0	.
II_C09	2427868	2430086	DTC2	0	.
II_C09	2466892	2469357	DTC3	0	.
```

BED6 欄位：紀錄、start（從 0 開始算）、end、Name、score（0）、strand。
注意：BED 從 0 開始，GFF3 從 1 開始。所以 BED 中的 `854081` 就是 GFF3 中的 `854082`。
大多數基因體瀏覽器（例如 IGV）都能開啟 BED 和 GFF3 檔案。

`--format bed12` 也會把各部分（TSD、TIR）寫成區塊（block）：

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format bed12 --name-prefix DTC \
      --input II_C09.fa scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat | head -2
II_C09	854081	860374	DTC1	0	.	854084	860371	0	4	3,12,12,3,	0,3,6278,6290,
II_C09	2427868	2430086	DTC2	0	.	2427871	2430083	0	4	3,12,12,3,	0,3,2203,2215,
```

你也可以用 awk 從現有的 GFF3 做出 BED（這裡用第 4.4 節的整個基因體
檔案）：

```
$ awk -F'\t' '$3=="repeat_region"{split($9,a,";"); sub("ID=","",a[1]); print $1"\t"$4-1"\t"$5"\t"a[1]"\t0\t"$7}' dtc_II.gff3 | head -3
II_C01	2874798	2878133	DTC1	0	.
II_C01	3398539	3401290	DTC2	0	.
II_C01	4397540	4400232	DTC3	0	.
```

### 5.3 給 Python 用的 JSON lines

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format jsonl --name-prefix DTC \
      --input II_C09.fa --output dtc_C09.jsonl \
      scan_for_matches_rs/tir_scan_patterns/labelled/00_DTC_published_2025.pat
$ head -1 dtc_C09.jsonl | cut -c1-200
{"seq":"II_C09","strand":".","start":854082,"end":860374,"name":"DTC1","type":"CACTA_TIR_transposon","units":[{"line":3,"start":854082,"end":854084,"text":"TAT","label":"target_site_duplication"},...
```

每一行是一個 hit。`units` 是 hit 的各個部分，依 pattern 中的順序排列。
只用標準函式庫的 Python（`read_jsonl.py`）：

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

（沒有標籤的單元，例如中間部分，沒有 `label` 這個鍵。所以
程式碼使用 `u.get("label")`。）如果你有 pandas：`pd.read_json(f, lines=True)`。

### 5.4 合併多次執行的結果：`--merge`

用兩次分開的執行跑兩個 pattern。兩次執行的 Name 都從 1 開始：

```
$ L=scan_for_matches_rs/tir_scan_patterns/labelled
$ scan_for_matches -t 4 -c --dedup --strict-n --format gff3 --input II_C09.fa --output a.gff3 $L/00_DTC_published_2025.pat
$ scan_for_matches -t 4 -c --dedup --strict-n --format gff3 --input II_C09.fa --output b.gff3 $L/09_DTC_CACTA_seed10.pat
$ grep -c -P '\trepeat_region\t' a.gff3 b.gff3
a.gff3:396
b.gff3:452
```

兩個檔案都有一個元件 `DTC1`。它們不是同一個元件。`--merge` 會合併
這些檔案，移除範圍相同（相同紀錄、start 和 end）的元件，並給出
新的、不重複的 Name。檔案依優先順序排列：保留第一個檔案中的元件。

```
$ scan_for_matches --merge --output ab.gff3 a.gff3 b.gff3
merge: 848 elements read, 389 duplicates removed, 459 written
$ grep -v '^#' ab.gff3 | awk '$3=="repeat_region"' | cut -f1,3,4,5,9 | head -3
II_C09	repeat_region	212663	236410	ID=DTC1;Name=DTC1;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657
II_C09	repeat_region	854082	860374	ID=DTC2;Name=DTC2;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657;Merged=b.gff3:DTC2
II_C09	repeat_region	2427869	2430086	ID=DTC3;Name=DTC3;Classification=TIR/DTC;Method=structural;Sequence_ontology=SO:0000657;Merged=b.gff3:DTC3
```

為什麼位於 854082 的元件現在叫 `DTC2`？它在 `a.gff3` 裡原本叫 `DTC1`。
`--merge` 會把所有元件依位置排序，再重新編號，所以**舊的 Name 會改變**。
位於 212663 的元件只在 `b.gff3` 裡，它的位置最前面，所以變成新的 `DTC1`。

`Merged=b.gff3:DTC2` 的意思是：`b.gff3` 裡有一個元件（舊名 `DTC2`）跟這個元件範圍完全相同，
所以它被移除了。

好用的合併選項：

* `--overlap 0.9` 在兩個元件彼此重疊至少 90 % 時，也把它們算作
  同一個元件：

  ```
  $ scan_for_matches --merge --overlap 0.9 --output ab90.gff3 a.gff3 b.gff3
  merge: 848 elements read, 391 duplicates removed, 457 written
  ```

  （移除了 391 個重複，而不是 389 個。）
* `--name-prefix P` 讓所有元件使用前綴 `P`（例如 `--name-prefix II_CACTA_` 會得到 `II_CACTA_1`）。
* 把檔案一個一個給它。不要用 `cat` 把它們接在一起：ID 會重複。
  程式會以錯誤停止：

```
$ cat a.gff3 b.gff3 | scan_for_matches --merge - > /dev/null
merge: -:2779: ID 'DTC1' is used twice (line 4); IDs must be unique in one file
```

### 5.5 和舊的註解比較（選用）

`II.CC.TE.gff3` 只是一個範例，代表你可能已有的這個基因體的舊註解。
它有 3,244 行 `intact_CACTA`。我們把它們和我們的 3,272 個元件比較。
`bedtools` 是一個處理區域的程式。在這個叢集上它是一個模組（module）（用
`module avail bedtools` 檢查）。`bedtools2/2.30.0` 是這裡的模組名稱，在你的系統上
可能不同。這裡沒有安裝任何東西：我們只是載入它。

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

`bedtools intersect` 的旗標（flag）：

| 旗標 | 意思 |
|---|---|
| `-a`、`-b` | 新的區域（a）和舊的區域（b） |
| `-f 0.9 -r` | 只有當重疊長度至少是 a 的長度「且」b 的長度的 90 % 時，才保留這一對（`-r` = 互相，reciprocal） |
| `-u` | 如果 a 區域在 b 中有對應，就把它寫出一次 |
| `-v` | 只寫出在 b 中「沒有」對應的 a 區域 |

所以我們的 3,272 個元件中，有 3,255 個在舊檔案中有對應，17 個沒有。

## 6. 一次執行多個 pattern

給它多個 pattern 檔案，並指定 `--format`。FASTA 檔案只會讀取一次。
下面的例子在整個基因體上跑全部 22 個現成的 TIR pattern。這是比較長的工作，
所以要用 Slurm（`run_all22.sh`）：

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

Slurm 輸出的結尾：

```
scan_for_matches: [0:01:00] 1 records (52.4 Mb) done, 41092 hits; now II_C02 (63.4 Mb), pattern 7/22, strand -, at 21.0 Mb
...
scan_for_matches: [0:05:58] 289 records (563.4 Mb) done, 440329 hits; finished

real	5m58.453s
```

（shell 會把 `*.pat` 展開成這 22 個檔案。）輸出檔案很大：
440,329 個元件，518 MB。

關於 Name：每個 pattern 檔案在它的 `%@element` 行中有自己的前綴
（`DTC`、`DTH`、`DTT`、...）。前綴相同的 pattern 共用一個計數器。
所以在這個檔案中，所有 Name 都不重複。每個前綴的元件數：

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

同一個元件可能被兩個 pattern 找到。在一次執行中，每筆紀錄會依照
命令列上的順序被各個 pattern 搜尋，hit 也依這個順序寫出。
所以對於相同的範圍，命令列上先給的 pattern 會在
GFF3 中排在前面。`--merge` 保留它先讀到的那一個（README.md，`--merge`：「the element of the
earlier file, then the earlier line, is kept」）。合併這一個檔案：

```
$ scan_for_matches --merge --output tir_merged.gff3 tir_all.gff3
merge: 440329 elements read, 29679 duplicates removed, 410650 written
$ grep -m2 'Merged=' tir_merged.gff3 | cut -f1,3,4,5,9 | cut -c1-300
II_C01	repeat_region	37068	38350	ID=DTT6;Name=DTT6;Classification=TIR/DTT;Method=structural;Sequence_ontology=SO:0000657;Merged=tir_all.gff3:DTC1539
II_C01	repeat_region	58913	59427	ID=DTE1;Name=DTE1;Classification=TIR/DTE;Method=structural;Sequence_ontology=SO:0000657;Merged=tir_all.gff3:DTP2
```

檢查第一行。在 `tir_all.gff3` 中，`II_C01` 上相同的範圍 37068-38350
出現了兩次：

```
$ awk -F'\t' '$1=="II_C01" && $4==37068 && $3=="repeat_region"{split($9,a,";"); print NR, $4,$5,a[1]}' tir_all.gff3
44097 37068 38350 ID=DTT5955
284099 37068 38350 ID=DTC1539
```

`DTT5955` 來自 `01_DTT_seed12.pat`，先被寫出（第 44097 行）。
`DTC1539` 來自 `09_DTC_...`，較晚被寫出（第 284099 行）。shell
會依名稱順序排列 `*.pat`，所以 `01_` 在 `09_` 前面。`--merge` 保留 DTT
元件（新名稱 `DTT6`），而 `Merged=` 記錄被移除那個元件的「舊」名稱
（`DTC1539`）。這次合併花了 6 秒。

在一條染色體上跑三個 pattern（GFF3 的 `# scan_for_matches` 行會列出
每個 pattern）：

```
$ scan_for_matches -t 4 -c --dedup --strict-n --format gff3 --input II_C09.fa --output C09_three.gff3 \
      $L/00_DTC_published_2025.pat $L/09_DTC_CACTA_seed10.pat $L/09_DTC_CACTG_seed10.pat
```

命名規則，簡單來說：

* 不加選項：前綴是 `%@element` 行的 `Name=`，沒有的話就是 `sfm`。
* `--name-prefix P` 為這次執行的所有 pattern 設定前綴。
* `--name-start N` 設定第一個編號。
* 分開的執行會產生重複的 Name。給每次執行各自的 `--name-prefix`，或用 `--merge` 合併這些執行的結果。
* `-m N` 是整次執行共用的一個上限：所有 pattern 加起來最多 N 個 hit。

## 7. 自己寫 pattern

把 pattern 放在一個文字檔案中。每個檔案一個 pattern。
一定要用 `--explain` 和 `--lint` 檢查它，並在一個放入目標序列的小 FASTA 上
測試它。完整的語言說明在 [PATTERNS.md](PATTERNS.md)。

### 7.1 使用模糊碼的模體

`W` 代表 A 或 T。這個 pattern 會找出類似 TATA box 的字詞：

```
$ echo 'TATAWAW' > tata.pat
$ printf '>demo1\nGGGCTATAAAAGGCCCTATATATGGCCGGGATTACAGGACCTTAAGCTAGCTAGCCGATTACAGGG\n>demo2\nAAAAAAAAAAAAAAAAAAAA\n' > demo.fa
$ scan_for_matches tata.pat < demo.fa
>demo1:[5,11]
TATAAAA 
>demo1:[17,23]
TATATAT 
```

### 7.2 正向重複

`p1=6...8` 取 6 到 8 個字母。`10...30` 是間隔。最後的 `p1` 是再次出現的相同字母。

| pattern | 輸入（`demo1` 的一部分） | 結果 |
|---|---|---|
| `p1=6...8 10...30 p1` | `GATTAC AGGACCTTAAGCTAGCTAGCC GATTAC` | hit 位置 30-62 |

```
$ echo 'p1=6...8 10...30 p1' > direct.pat
$ scan_for_matches direct.pat < demo.fa
>demo1:[30,62]
GATTAC AGGACCTTAAGCTAGCTAGCC GATTAC 
```

### 7.3 允許一個字母不同的髮夾結構

`[1,0,0]` 的意思是：1 個錯配、0 個缺失、0 個插入。這個莖有 1 個錯配
（`CTGTAAGC` 不完全是 `GATTACAG` 的反向互補，反向互補是 `CTGTAATC`）。

| pattern | 在 `stem.fa` 上的結果 |
|---|---|
| `p1=8...8 4...10 ~p1`（完全相符） | 沒有 hit |
| `p1=8...8 4...10 ~p1[1,0,0]` | hit 位置 5-24 |

```
$ printf '>stem\nAAAAGATTACAGTTTTCTGTAAGCAAAA\n' > stem.fa
$ echo 'p1=8...8 4...10 ~p1' > hp0.pat
$ echo 'p1=8...8 4...10 ~p1[1,0,0]' > hp1.pat
$ scan_for_matches hp0.pat < stem.fa
$ scan_for_matches hp1.pat < stem.fa
>stem:[5,24]
GATTACAG TTTT CTGTAAGC 
```

**陷阱：最短的莖會勝出。** `p1=5...8` 會使用可行的最短長度
（`--explain` 會這樣說明）。允許錯配時，常會先找到一個短的假莖，
而搜尋會從它後面繼續，於是真正的髮夾結構就被漏掉了：

| pattern | 在 `AAAAACCGATTACTTTTGTAAACGGAAAA` 上的結果 |
|---|---|
| `p1=5...8 3...6 ~p1[1,0,0]` | `[2,17]`：`AAAAC CGATTA CTTTT`（假的：一個有 1 個錯配的 5 字母莖） |
| `p1=8...8 3...6 ~p1[1,0,0]` | `[6,25]`：`CCGATTAC TTTT GTAAACGG`（真正的髮夾結構） |

請使用固定的莖長度，或為每個長度各跑一個 pattern。

### 7.4 帶標籤的小型 TIR 轉位子

標籤（`%@`）在 GFF3 中為各部分命名。這是一個類似 CACTA 的小元件，
間隔為 20 到 100 個字母。把它放進 `mini.pat`：

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

破壞 TSD（最後的 `ATG` 變成 `TTG`），hit 就消失了：

```
$ printf '>mini_broken\nGGCATCGGATGCACTACAAGAAAATTGGCCTTAACCGGATCGCTAGCTAGGGATTTCTTGTAGTGTTGGCGGTCATGGC\n' | scan_for_matches --format gff3 mini.pat | cut -c1-100
##gff-version 3
# scan_for_matches 0.4.3 pattern=mini.pat
```

### 7.5 蛋白質 pattern

使用 `-p` 時，字母是胺基酸。`any(ST)` 是 S 或 T。`X` 是任意胺基酸。
`notany(P)` 是 P 以外的任意字母。`-c` 不能和 `-p` 一起使用。

```
$ echo 'M any(ST) X X notany(P) K' > p.pat
$ printf '>s\nGMSAAGKLL\n' > s.fa
$ scan_for_matches -p p.pat < s.fa
>s:[2,7]
M S A A G K 
```

### 7.6 更多說明

* 完整的語言說明（範圍、名稱、替代選項、RNA 規則集、權重
  矩陣、長度限制、常見錯誤）：[PATTERNS.md](PATTERNS.md)。
* 把 JASPAR 矩陣轉成權重單元：`python3 -I tools/jaspar2sfm.py MATRIX.jaspar --pvalue 1e-4`
  （[README.md](README.md#jaspar-matrix-to-weight-unit-toolsjaspar2sfmpy)）。
* AI 助理可以幫你寫 pattern：請看
  [`skills/sfm-patterns/SKILL.md`](skills/sfm-patterns/SKILL.md)。
  用 `--explain` 和 `--lint` 檢查它寫的內容。

## 8. 所有選項

`scan_for_matches --help` 的所有選項：

| 選項 | 意思 | 什麼時候用 |
|---|---|---|
| `-c` | 也搜尋相反的那一條鏈 | DNA，目標可能在任一條鏈上時。不能和 `-p` 一起用 |
| `-p` | 序列是蛋白質 | 蛋白質 FASTA |
| `-n N` | 累計 N 條序列沒有 hit 後停止（不必連續；結束狀態 1） | 很少用；停止一個找不到任何東西的搜尋 |
| `-m N` | 總共最多寫出 N 個 hit（所有 pattern 加起來） | 快速測試：`-m 1` |
| `-o N` | 舊程式留下來的選項：顯示重疊的 hit。不需要時就不要寫。要用時一律寫 `-o 1`：後面一定要接一個數字，但數字沒有作用 | hit 可能重疊或巢狀時 |
| `-i FILE` | 列出要跳過的紀錄 id 的檔案 | 跳過胞器或 contig 紀錄 |
| `-t N` | N 個執行緒。預設 1。不論 N 是多少，輸出都相同 | 每次大型執行。也要向 Slurm 要求 N 個核心 |
| `--input FILE` | 從 FILE 讀取 FASTA（`-` 是 stdin） | 取代 `< file` |
| `--output FILE` | 把 hit 寫到 FILE（`-` 是 stdout） | 取代 `> file` |
| `--progress` | 每分鐘輸出一行進度到 stderr | 長時間的執行；在 Slurm 下它在 `.out` 檔案中 |
| `-h`、`--help` | 顯示說明 | |
| `--format F` | 用 `gff3`、`bed6`、`bed12` 或 `jsonl` 取代原始輸出 | 處理基因體時一律使用。多個 pattern 時必須使用 |
| `--name-prefix P` | 每個 hit 的名稱 = P + 編號 | 讓每次執行有自己的名稱 |
| `--name-start N` | 第一個編號（預設 1） | 接續編號 |
| `--type T` | 元件的第 3 欄（預設：`%@element` 的 type，沒有的話是 `sequence_motif`） | 為沒有標籤的 pattern 給一個更好的 type |
| `--strict-n` | 被名稱抓到的 `N`（或其他非 ACGT 字母），在這個名稱再次使用時算作 1 個錯配 | 每個有組裝缺口的基因體。只限 DNA |
| `--align` | 有缺失或插入的 `~p1[m,d,i]` 會嘗試每一種比對長度 | `indel/` pattern，或任何 `d` 或 `i` 大於 0 的 `~pN[m,d,i]` |
| `--dedup` | 搭配 `-c`（和 `--format`）：兩條鏈上相同的範圍只寫一次，鏈為 `.` | 搭配 `-c` 的 TIR pattern 和髮夾結構 |
| `--lint` | 檢查 pattern 中的陷阱。不讀取輸入 | 每次長時間執行之前 |
| `--explain` | 用簡單的話描述 pattern。不讀取輸入 | 當你不確定一個 pattern 是否正確時 |
| `--merge` | 合併 GFF3 檔案（它自己的選項見下方） | 多次執行之後 |

`--merge` 的選項（`scan_for_matches --merge --help`）：

| 選項 | 意思 | 什麼時候用 |
|---|---|---|
| `--overlap F` | 兩個元件重疊的長度，同時達到兩者各自長度的 F（0 到 1）時，也算作同一個元件 | 範圍稍有不同的 pattern，例如 `--overlap 0.9` |
| `--name-prefix P` | 所有新 Name 的前綴 | 一套乾淨的名稱 |
| `--name-start N` | 第一個編號 | |
| `--output FILE` | 寫到 FILE，而不是 stdout | |

注意事項：

* 長選項可以寫成 `--format gff3` 或 `--format=gff3`。
* `--dedup`、`--name-prefix`、`--name-start` 和 `--type` 需要 `--format`。
* `-i` 和 `-o` 是原始程式的選項（要跳過的 id、重疊的
  hit）。它們**不是**輸入和輸出檔案。

## 9. 速度與記憶體

下面的數字來自 [README.md](README.md)（Quick start、
Several patterns in one run、Speed 這幾節）和
[`tir_scan_patterns/README.md`](tir_scan_patterns/README.md)。

* **執行緒（`-t N`）。** 執行緒越多，時間越短。一個 pattern
  （02_DTA_short5to7，81,514 個 hit）在 1.0 Gb 的基因體上，`-t 1` 花 349 s，
  `-t 2` 花 184 s，`-t 4` 花 125 s，`-t 8` 花 110 s；記憶體
  從 0.6 GB 增加到 4.2 GB（README.md，Speed；Apple M4）。那台機器
  受記憶體速度限制，所以在那裡超過 4 個執行緒幾乎沒有幫助。
  向 Slurm 要求相同數量的核心（`-c N`）。
* **bgzip 輸入。** 使用 `-t N` 時，bgzip 檔案的區塊由
  N 個執行緒解壓縮。一般的 gzip 只用一個執行緒（README.md，Usage）。你可以直接在
  `genome.fa.gz` 上執行。製作 bgzip 檔案：`bgzip genome.fa`
  （屬於 htslib/samtools）。
* **多個 pattern，一次執行。** 在 1.0 Gb 的 *B. napus* 基因體上跑 22 個 TIR pattern，
  `-t 20`：一次執行花 648 s，分成 22 次執行花 762 s（README.md，
  「Several patterns in one run」）。Quick start 中的數字約為 11 分鐘。
* **間隔跳躍。** 一個寬範圍後面接一個完全相符的字詞，例如
  `50...30000 ~p1`，會跳過不可能成立的間隔長度：在 TIR pattern 上最多
  快 1000 倍（README.md，Quick start）。讓緊接在寬間隔後面的單元
  保持完全相符（見下一點）。
* **`slow` 提示。** 當緊接在寬間隔（`500...15000`）後面的單元
  無法用來跳過間隔長度時，`--lint` 會說 `slow`。這時程式會逐一嘗試
  每個間隔長度。依下列順序修正：
  1. 在間隔後面緊接一個完全相符的字詞，或一個完全相符的 `~pN`。
  2. 如果間隔後面的單元需要錯配，把錯配寫在 `~pN[m,0,0]` 上，並且
     讓 `pN` 至少有 5 x (m + 1) 個字母長（1 個錯配就是 10 個字母），
     或在它後面放一個 5 個字母、完全相符的字詞。
  3. 有缺失或插入時（`~pN[m,d,i]`），`pN` 至少需要
     4 x (m + d + i + 1) 個字母（`[2,1,1]` 需要 20 個）。後面接完全相符的字母沒有幫助。
  4. 或把間隔縮窄。

  （數字來源：README.md，0.4.2 與 0.4.3 的條目；PATTERNS.md，mistake 8。）
  真實輸出，一個慢的 pattern 和一個快的 pattern：

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

  第二個 pattern 把錯配放在 `~p2[1,0,0]` 上，後面接 5 個完全相符的字母（`TAGTG`），
  所以沒有 `slow` 提示。（如果沒有 `TAGTG`，7 個字母的 `p2` 就太短：
  提示會再出現。）
* **`mm2/` 和 `indel/` pattern。** `tir_scan_patterns/mm2/` 允許右側 TIR 有 2 個錯配。
  在大約 100 Mb 的 *B. napus* 上（`-t 40`），它們花了 4 s、
  27 s、23 s 和 34 s（完全相符的版本是 3 s、14 s、10 s 和 4 s）。
  `tir_scan_patterns/indel/`（要用 `--align`）在相同資料上花了 6 min、2.7 h、1.8 h 和
  5 min（tir_scan_patterns/README.md）。用 `sbatch` 執行 indel pattern，
  並檢查 hit：隨機 hit 的比例比較高。
* **`--progress`** 顯示長時間執行進行到哪裡（已完成的紀錄數和 Mb、hit 數、
  目前正在搜尋的紀錄）。
* **記憶體。** 遇到長間隔時，程式會為紀錄建立一個 [5-mer 索引](#本指南用到的詞)：大約
  每個鹼基 4 bytes（一條 74 Mb 的染色體：大約 300 MB）。玉米（2.3 Gb），`-t 20`，
  pattern 00_DTC_published_2025：最多 2.5 GB；02_DTA_short5to7：5.1 GB
  （README.md，Speed）。產生數百萬個 hit 的 pattern 會產生很大的輸出
  檔案：第 6 節的 22 個 pattern 執行寫出了 518 MB。要測量你自己的執行，
  請看第 4.4 節。
* **從小開始。** 在完整執行之前，先在一條染色體或一小段序列上測試。

## 10. 常見問題（FAQ）

**我用了 `-o out.txt`，卻沒有輸出檔案。**
`-o` 不是輸出檔案（見第 8 節的 `-o 1` 那一列）。程式把
`out.txt` 讀成 `-o` 的值，把 hit 寫到螢幕上，不會建立檔案。
請用 `--output out.txt`。

```
$ scan_for_matches -o out.txt hp1.pat < stem.fa
>stem:[5,24]
GATTACAG TTTT CTGTAAGC 
$ ls out.txt
ls: cannot access 'out.txt': No such file or directory
```

**`-o -c` 沒有搜尋另一條鏈。**
`-o` 會把下一個字當成它的值，所以 `-c` 被吃掉了。請寫 `-o 1 -c`。
`--lint` 會警告你：

```
hp1.pat: warning: -o takes a value: `-c` is read as the value of -o, not as an option, so the reverse strand is not searched; write -o 1
```

**使用 `-c` 時，每個元件都被找到兩次。**
在兩條鏈上讀起來都一樣的 pattern（髮夾結構、TIR pattern）會讓每個
元件出現兩次。加上 `--dedup`（搭配 `--format`）。鏈會變成 `.`。

| 指令 | 在 `stem.fa` 上的結果（BED6） |
|---|---|
| `scan_for_matches -c --format bed6 hp1.pat` | `stem 4 24 sfm1 0 +` 和 `stem 4 24 sfm2 0 -` |
| `scan_for_matches -c --dedup --format bed6 hp1.pat` | `stem 4 24 sfm1 0 .` |

`--dedup` 沒有搭配 `--format` 是錯誤：
`scan_for_matches: --dedup, --name-prefix, --name-start and --type need --format`。

**組裝有缺口（`N`）的地方出現假 hit。**
加上 `--strict-n`。沒有它時，被名稱抓到的 `N` 在這個名稱再次使用時
會符合任意字母，所以由 `N` 組成的 TSD 會產生假元件。加上它之後，每個這樣的
字母算作 1 個錯配。請看 [README.md](README.md#assembly-gaps---strict-n)。在你加上它之前，
`--lint` 會用一則提示提醒你。

**indel pattern（`~p2[2,1,1]`）漏掉元件。**
加上 `--align`。沒有它時，比對器只採用一種比對長度，會漏掉一些
hit。`--lint` 會這樣說，而有 `--align` 時就不會說：

```
$ scan_for_matches --lint -c --dedup --strict-n --format gff3 scan_for_matches_rs/tir_scan_patterns/indel/labelled/04_DTE_seed20_indel.pat
...:3: note: `p1=8...9`: the shortest length that works is used, not the longest
...:6: note: `~p2[2,1,1]` has deletions and insertions: the matcher takes one greedy alignment and can miss some matches (add --align to try every length)
0 errors, 0 warnings, 2 notes
$ scan_for_matches --lint -c --dedup --strict-n --align --format gff3 scan_for_matches_rs/tir_scan_patterns/indel/labelled/04_DTE_seed20_indel.pat
...:3: note: `p1=8...9`: the shortest length that works is used, not the longest
0 errors, 0 warnings, 1 notes
```

**不同執行之間的名稱衝突（兩個檔案都有 `DTC1`）。**
每次執行都從 1 開始。給每次執行各自的 `--name-prefix`，或用
`--merge` 合併這些檔案（第 5.4 節）。不要用 `cat` 接合 GFF3 檔案。

**"several pattern files need --format"。**
多於一個 pattern 檔案時，需要 `--format gff3`（或其他格式）。舊的輸出
無法說明是哪個 pattern 產生了某個 hit。結束狀態是 2。

**"-c cannot be used with -p"。**
蛋白質沒有反鏈。請移除 `-c`。結束狀態是 2。

**一個 pattern 檔案裡有兩個 pattern。**
一個檔案中的所有文字是「一個」pattern。每個 pattern 用一個檔案。

**沒有 hit，也沒有錯誤。**
用 `--explain` 檢查 pattern（這是你想要的嗎？），並在一個放入目標序列的小
FASTA 上測試它。檢查 FASTA 的 id 沒有出現在 `-i` 檔案中。
記住：pattern 中的字母必須依照正鏈的順序。

**`-n 1` 讓我的執行停止了。**
`-n N` 會在累計 N 條序列沒有 hit 後停止（不必連續；結束狀態 1）。不要使用 `-n`。

**`--lint` 說 `slow`。**
請看第 9 節的 `slow` 提示（四種修正方式，依順序進行）。

**怎麼查看正在執行的工作？**
兩個指令。`squeue` 顯示佇列：`PD` = 等待中，`R` = 執行中。當工作
不再出現在清單中，就表示它已經結束。`tail` 顯示 Slurm
輸出檔案的最後幾行；使用 `--progress` 時，每分鐘會有新的一行。

```
$ squeue -u $USER -n sfm_mem
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
            580692       CPU  sfm_mem    pxzhe PD       0:00      1 (None)
$ grep 'scan_for_matches:' sfm_all22-580689.out | tail -n 2 | cut -c1-170
scan_for_matches: [0:05:00] 8 records (467.2 Mb) done, 367678 hits; now II_C09 (65.9 Mb), pattern 4/22, strand -, at 65.9 Mb
scan_for_matches: [0:05:58] 289 records (563.4 Mb) done, 440329 hits; finished
```

（用同樣方式使用 `tail -n 2 sfm_dtc-<jobid>.out`；`<jobid>` 是
`sbatch` 輸出的那個數字。要停止工作：`scancel <jobid>`。）

**哪裡可以取得幫助。**
`scan_for_matches --help`、[PATTERNS.md](PATTERNS.md)、[README.md](README.md)、
[CHANGELOG.md](CHANGELOG.md)，以及
[`tir_scan_patterns/`](tir_scan_patterns/README.md)（中文說明：[README_zh-TW.md](tir_scan_patterns/README_zh-TW.md)）中的 22 個現成 pattern。
