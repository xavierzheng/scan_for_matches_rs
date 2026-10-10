Scan for Matches：依 Wicker et al. (2007) 整理的 TIR 候選搜尋式
整理日期：2026-09-30

重要定位
本資料夾提供 22 個可各自執行的 .pat 檔。除 00_DTC_published_2025.pat
是使用者提供、且可在 Zheng et al. (2025) 核對的原式外，其餘皆為依文獻
結構特徵設計的候選搜尋式，不是原作者發表過的完整 pattern，也未經真實
基因組的 sensitivity/precision benchmark。

已實際檢視使用者提供的 Wicker2007_TEclassify.pdf：Figure 1 (p.974)、
Figure 2 (p.977)、Figure 3 (p.979)、Figure 4 (p.980)，並核對 pp.976、
977、980 的文字。範圍以該文 2007 年分類為準，不主張是現行分類的全集。

一、圖中分類與搜尋結構

同一條 DNA 以 5' -> 3' 表示：
TSD | left terminal seed | intervening sequence | reverse-complement seed | TSD

TIR 是轉座子內兩端的反向互補序列；TSD 是轉座子外兩側同向重複的宿主
序列。p1 代表前面實際讀到的同一序列；~p2 代表 p2 的反向互補，不是單純
倒序。TSD 必須接在轉座子外緣，不能在 TSD 與末端間任意增加空隙。

Figure 1 的九類依序為：
DTT Tc1-Mariner：TA TSD。
DTA hAT：8 bp TSD；本文描述 TIR 為 5-27 bp，且缺乏通用診斷 motif。
DTM Mutator：9-11 bp TSD；常見末端 G...C，TIR 可很長、很短或無法辨識。
DTE Merlin：8-9 bp TSD；TIR 可達數十至數百 bp。
DTR Transib：5 bp TSD。
DTP P：8 bp TSD。
DTB piggyBac：TTAA TSD。
DTH PIF-Harbinger：3 bp TSD；正文提到 TAA 偏好。
DTC CACTA：植物 3 bp TSD，末端 CACTA 或 CACTG；本文歸入此類的
動物／真菌型可為 CCC...GGG，搭配 2 bp TSD。

圖中 Maverick (DMM) 雖有 TIR，卻位於 subclass 2、order Maverick，
不是上述 TIR order 的第十個 superfamily。Helitron 無成對末端 TIR，
其 3' 附近髮夾也不能當成 TIR transposon 的兩端。DIRS 圖中也有反向
末端結構，但屬 class I。看到 IR 結構不等於已辨識出 TIR-order TE。
Figure 2 顯示同一 CACTA family 可含不同內部缺失衍生物；結構完整
不等於具有自主轉座能力。Figure 3 提醒應整合核酸、蛋白同源性與結構。

二、如何讀取本套 pattern

所有 .pat 檔只有一個 pattern，不含標題或註解。不要把整個資料夾的
pattern 接成一個輸入檔：Scan for Matches 會將連續文字當成同一 pattern。

p1=3...3：讀取並記住 3 bp；末尾 p1 要求完全相同的同向序列。
p2=12...12：讀取 12 bp 的末端 seed；~p2 要求另一端反向互補。
CACTA[0,0,0]：要求 CACTA 完全匹配。
~p2[1,0,0]：允許另一端相對於 p2 的反向互補有至多 1 個替換；不允許 indel。
所有未加誤差修飾的重複參照，預設為完全匹配。

此處 seed 是搜尋時檢查的末端片段，不是該 superfamily 的完整 TIR 長度
定義。例如 TIR 真正長 200 bp，仍可用最外端 20 bp 作 seed，剩下的 TIR
序列落入中間區段。找到候選後必須延伸比對，以量測真正的 TIR 長度。

除原始 DTC 式與 DMM 補充式外，50...30000 都是本次選擇的中間區段
搜尋窗口，不是文獻宣稱的全部元素大小；8、12、15、20 bp 等 seed 長度
也屬搜尋設定。這個窗口同時容納部分小型缺失衍生物和較大元素，但可能
很耗時，也不涵蓋所有極短元素、超長元素或含大型巢狀插入的元素。

若 seed 長度為 k，中間區段長度為 g，且無 indel：
轉座子長度（排除兩側 TSD）= 2*k + g。
完整命中長度 = 2*TSD長度 + 2*k + g。

原始 DTC：p2 只有 7 bp，但加上 CACTA 的 5 bp，各端檢查共 12 bp。
500...15000 是中間區段，所以轉座子長度是 524-15024 bp；含兩側 TSD
的命中長度是 530-15030 bp。它會漏掉不足 12 bp 的 TIR、末端突變、
CACTG 型、CCC 型以及窗口外的元素。CACTA 文獻描述的 TIR 可短至 10 bp，
因此另提供 seed10 版本。若須嚴格限制「轉座子本體」500-15000 bp，
此 12 bp seed 寫法的中間區段應改為 476...14976。

三、檔案選擇與分類限制

01_DTT_seed12：不限定 Stowaway motif 的 TA-TSD 候選。
01_DTT_Stowaway_seed12：CTCCTCCC...GGGAGGAG 子型；不能代表所有 DTT。

02_DTA_seed8：8-bp TSD + 8-bp seed；也會匹配其他 8-bp-TSD 元素。
02_DTA_short5to7：補找極短 TIR；隨機及低複雜度背景會大幅增加。
若希望依 Wicker 文中 hAT 本體小於 4 kb 的描述縮窗，seed8 版本可將
50...30000 改為 50...3983；這是額外篩選，可能漏掉有插入的衍生物。

03_DTM_G_seed20：末端 G/C 型、9-11-bp TSD；有 20-bp exact seed。
03_DTM_unanchored_seed20：不要求 G/C 末端的補充搜尋。
無 TIR 的 Mutator 不能靠這些 IR pattern 找出，須依同源性與其他證據。

04_DTE_seed20：8-9-bp TSD、20-bp seed；不能用較長 seed 就排除 hAT/P。
05_DTR_CAC_seed15：CAC...GTG + 5-bp TSD；CAC 來自 Kapitonov & Jurka
(2005) 的 Transib 末端比對，而非 Wicker Figure 1 直接列出的序列。
05_DTR_unanchored_seed15：不要求 CAC 的補充式。
06_DTP_seed12：8-bp TSD、12-bp seed，只是與 P 相容的候選，不是 P 專屬式。
07_DTB_seed12：TTAA 是外側 TSD，不能誤用為內側 TIR。

08_DTH_seed12：任意相同的 3-bp TSD；也可能命中 CACTA 等其他候選。
08_DTH_TWA_seed12：只找 TAA 或 TTA 偏好位點；p1=TWA 記住實際匹配的
triplet，所以兩端必須同為 TAA 或同為 TTA。不可改成左右獨立的 TWA，
否則可能錯誤允許左 TAA、右 TTA。
08_DTH_Tourist_seed8：依 RiTE 的 Tourist 候選條件，末端 G/C、seed 至少
8 bp、3-bp TSD。p3=S 會記住實際 G 或 C，~p3 要求正確互補；不是兩端
各自任選 S。本式是本次轉寫，RiTE 原文沒有提供這一行完整 pattern。

09_DTC_CACTA_seed12、09_DTC_CACTG_seed12：分開找植物兩種末端。
09_DTC_CACTA_seed10、09_DTC_CACTG_seed10：補找較短的 exact TIR。
09_DTC_CCC_seed12：Wicker 2007 所述動物／真菌 CCC 型，TSD 為 2 bp。
09_DTC_CACTA_relaxed1：僅放寬 p2 部分的一個替換；TSD 與 CACTA 保持 exact。

10_DMM_Maverick_seed30：order Maverick 的額外結構搜尋；6-bp TSD、
30-bp seed，9940...19940 使本體長度為 10000-20000 bp。仍須驗證長 TIR
與 POLB、integrase 等同源性；單獨命中不等於確診 Maverick。

DTA、DTP、DTE 的候選有交集；seed 長度是敏感度選擇，不是可靠的分類
分界。DTH broad 也不具專一性。缺乏可區分證據時應保留未定類候選。

四、執行與後處理

在已安裝 scan_for_matches 的環境中，一次執行一個 .pat：
scan_for_matches -o 1 09_DTC_CACTA_seed12.pat < genome.fa > DTC.hits.fa

本次測試的是 SEED 官網提供的 C 原始碼版本，-o 1 可保留重疊及替代
匹配。各發行版選項可能不同，請以自己的版本核對；預設略過重疊可能
漏掉巢狀或重疊候選。

本套完整元素 pattern 在反向互補後具有相同結構（DTH 的 TWA 已同時
涵蓋 TAA/TTA），因此掃描輸入鏈即可捕捉兩種元素朝向，不必為它們再加
-c。這不會判定 transposase 的轉錄方向。若自行改成非對稱、家族特定的
pattern，才須另處理兩條鏈。實測此舊版 -c 的反向鏈不保留全部重疊
匹配；需完整列舉時可分別掃描正向和反向互補 FASTA 後轉換座標、去重。

後處理至少包括：
1. 排除 TSD/末端 seed 含 N 或低複雜度造成的可疑命中，檢查 assembly gap。
2. 從輸出移除兩側 TSD 才是候選轉座子本體；變長 TSD 請讀實際匹配長度。
3. 從兩個外緣向內延伸、比對完整 TIR；不要把 seed 長度寫成生物學 TIR 長度。
4. 去除同座標、不同 seed 拆法造成的重複；保留真正不同的巢狀候選。
5. 整合多拷貝邊界、空插入位點（若有）、已知 TE library、transposase
   同源性或蛋白 domain。低複雜度 IR 或偶然 TSD 不能單獨證明是 TE。
6. 大型基因組先小區域測試，再按染色體或重疊窗口掃描；窗口 overlap
   應至少覆蓋允許的最長完整命中，並將座標還原和去重。本套一般式最長
   完整命中為約 30.1 kb，可用 31 kb overlap。巢狀插入超過窗口仍會漏。

此來源版本把單條輸入序列的大小固定為 250,000,000 nt；較長染色體應
先切成重疊區段或使用已修正的版本。分段不會消除廣泛回溯搜尋的成本。

五、文獻查核：哪些是原文、哪些不是

[1] Wicker et al. 2007. A unified classification system for eukaryotic
transposable elements. Nature Reviews Genetics 8:973-982.
使用者提供的 PDF 是本次主要分類依據；實際閱讀圖與正文。
https://doi.org/10.1038/nrg2165

[2] Dsouza, Larsen & Overbeek. 1997. Searching for patterns in genomic data.
Trends in Genetics 13:497-498. 工具原始文獻；可核對書目，PubMed 無摘要。
不能將此書目本身當成九種 TE pattern 都已發表的證據。
https://pubmed.ncbi.nlm.nih.gov/9433140/

[3] Ashok Aiyar, 1999-01-20, bionet.software 討論串：
How to find (imperfect) repeats in DNA sequences?
推薦 Scan for Matches 尋找 inverted repeats；未提供九類 TE 診斷式。
https://groups.google.com/g/bionet.software/c/YeMC40OL6_I

[4] The SEED Team, 2010-07-16. Scan For Matches 官方說明。
提供 p1、~p1、誤差、變長間隔與回溯搜尋語法。
https://blog.theseed.org/servers/2010/07/scan-for-matches.html
本次測試用的來源套件：
https://www.theseed.org/servers/downloads/scan_for_matches.tgz

[5] Wicker et al. 2003. CACTA Transposons in Triticeae. A Diverse Family
of High-Copy Repetitive Elements. Plant Physiology 132:52-63.
使用 BLAST 和 dot plot，不是已核實的 Scan for Matches 方法。
支持植物 CACTA 的 10-28 bp TIR 與各種大小的缺失衍生物；有 274 bp
的小型元素，也討論 23 kb 的 Candystripe1，故 500-15000 不是全集界線。
https://pmc.ncbi.nlm.nih.gov/articles/PMC166951/

[6] Nicolas et al. 2005. Suffix-tree analyser (STAN): looking for
nucleotidic and peptidic patterns in chromosomes. Bioinformatics 21:4408-4410.
以 AtREP3 比較 STAN、PatScan、GenLang；AtREP3 是 Helitron 型 TE，
不可將其近 3' 髮夾 pattern 說成 TIR-order 元素的成對末端 pattern。
也指出當時 PatScan 對重疊替代解的限制。
https://academic.oup.com/bioinformatics/article/21/24/4408/179826

[7] Biostars: Searching Repeats And Palindromic Sequences In Dna Sequences.
討論有推薦 scan_for_matches，但沒有逐 superfamily 的完整 pattern。
https://www.biostars.org/p/79567/

[8] Copetti et al. 2015. RiTE database: a resource database for genus-wide
rice genomics and evolutionary biology. BMC Genomics 16:538.
Stowaway/Tourist 小節支持特定 TIR/TSD 條件，未明確刊出其完整 SFM 式。
明確點名 SFM 並刊出的 pattern 則用於 Helitron 的 3' 髮夾：
p1=7...10 2...4 ~p1[1,0,0] 6...10 CTRRT
此式僅為方法來源對照，不列入本套 TIR 搜尋 .pat 檔。
https://pmc.ncbi.nlm.nih.gov/articles/PMC4508813/

[9] Zheng et al. 2025. Transposable elements drive evolution and perturb
gene expression in Brassica rapa and B. oleracea. The Plant Journal 123:e70452.
Methods 可核對使用者提供的 DTC pattern。本套 00 檔原樣保留。
https://doi.org/10.1111/tpj.70452
https://pmc.ncbi.nlm.nih.gov/articles/PMC12401559/

[10] Kapitonov & Jurka. 2005. RAG1 Core and V(D)J Recombination Signal
Sequences Were Derived from Transib Transposons. PLoS Biology 3:e181.
支持 Transib 末端 CAC 的保守性及 5-bp TSD；本文非 SFM 方法來源。
https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030181

[11] Jurka & Kapitonov. 2001. PIFs meet Tourists and Harbingers:
A superfamily reunion. PNAS 98:12315-12316.
支持 PIF/Tourist 的 3-bp TSD 與 TTA 偏好。
https://pmc.ncbi.nlm.nih.gov/articles/PMC60043/

本次沒有找到可核實、早期就逐一刊出 Wicker 九類完整 SFM pattern 的
論文或討論，因此不將本次設計冒稱文獻原式，也不宣稱窮盡所有早期研究。

六、驗證範圍

以官網原始碼編譯後實際執行 22 個交付 pattern，每個均通過：
合成正例、錯誤 TSD 反例、超過容許數目的 TIR 替換反例，以及未加
anchor 的交付 pattern 是否能找到預期座標。反例使用 ^ 和 $ 隔離待測
結構，避免其他位置的偶然匹配干擾測試。
另通過 14 個 TWA、G/C、9/10/11 bp TSD、8/9 bp TSD、單一 mismatch、
DTC 中間區段上下界的測試，並確認 -o 1 可列出正向的三個重疊 AA。
這些結果不構成真實資料上的召回率、專一性或執行速度 benchmark。

額外實測舊版 indel 語意：
ACGT[0,1,0] 可匹配 AGT，但不能匹配 ACAGT。
ACGT[0,0,1] 可匹配 ACAGT，但不能匹配 AGT。
所以以「相對於 pattern 的待匹配序列」描述，此來源實際順序為
[替換, 刪除, 插入]。原始碼內部變數命名可能使解讀混淆；本套 pattern
全部把後兩項設為 0，避免 indel 方向與末端座標歧義。

catalog.json 保存每個 pattern 的用途；validation.json 保存測試結果
與來源套件 SHA-256。下載包不包含軟體可執行檔，也未執行使用者基因組。

七、搭配本 Rust 版使用（新增）

上面的 .pat 檔都沒有改。下面是本 Rust 版（scan_for_matches_rs）才有的
功能，與「四、執行與後處理」中 C 原版的限制不同。

labelled/：22 個 pattern 加上 %@ 標籤的複本，給 --format gff3、bed6、
bed12、jsonl 用；標籤寫出元素、TSD 和 TIR 的名稱（見主 README 的
"Output formats"）。不加 --format 時，結果和沒有標籤的檔案一樣。

C 原版每條序列 250,000,000 nt 的上限已拿掉（0.1.0），不必再切段。

-c：本套 pattern 兩條鏈的結構相同，加 -c 會讓每個元素在兩條鏈各找到
一次；配合 --format 時加 --dedup，只留一份（strand 為 "."）。

--strict-n（0.4.1）：沒有它時，TSD（p1）若從 assembly gap 抓到 N，
另一端的 p1 會接受任何字母，所以 gap 旁會出現假的元素。加上它，抓到
的每個非 A、C、G、T 字母算 1 個 mismatch；本套 TSD 不容許 mismatch，
所以這種 hit 會被丟掉。有 gap 的基因組都應該加。

--progress（0.4.3）：每分鐘在 stderr 寫一行進度（做完幾條序列、幾 Mb、
幾個 hit、現在在找哪一條；用 -t N 時還有哪一股、找到第幾 Mb）。結果
（stdout 或 --output）不變。在 Slurm 下這些行會出現在 slurm-*.out。

範例（一個 pattern，bgzip 或一般 FASTA，8 個執行緒）：
scan_for_matches -t 8 -c --dedup --strict-n --format gff3 --name-prefix DTC \
    --input genome.fna.gz --output DTC.gff3 \
    tir_scan_patterns/labelled/00_DTC_published_2025.pat
每個 pattern 用不同的 --name-prefix（或 --name-start），Name 才不會重複。

一次跑全部 pattern（0.4.0；基因組只讀一次），再合成一個檔案：
scan_for_matches -t 20 -c --dedup --strict-n --format gff3 \
    --input genome.fna.gz --output tir_all.gff3 tir_scan_patterns/labelled/*.pat
scan_for_matches --merge --output tir_merged.gff3 tir_all.gff3
--merge 在兩個元素範圍相同時留第一個（在前面的 pattern）；
--overlap 0.9 也會合併彼此重疊 90 % 的元素。
修改 pattern 後，可用 scan_for_matches --lint FILE.pat 檢查、
--explain FILE.pat 說明。

22 個 pattern 在 1.0 Gb 甘藍型油菜（B. napus）基因組上的速度：見
CHANGELOG.md。

八、容許 mismatch 的版本：mm2/（新增）

老的元素有突變，兩端 TIR 已不是完全相同的複本。TIR 較長（19-30 bp）
的 4 個 pattern 在 mm2/（無標籤）和 mm2/labelled/ 各有一份複本，右端
TIR 容許最多 2 個替換：~p2 改成 ~p2[2,0,0]。TSD 和 G/C 末端仍要完全
相同。原本 22 個檔案沒有改，labelled/*.pat 也不包含這些複本。

mm2/ 的檔案                    真實 hit   真實 hit   打亂 hit   時間
                               exact      mm2        mm2        exact → mm2
03_DTM_G_seed20_mm2              885      1 285        0        3 s → 4 s
03_DTM_unanchored_seed20_mm2   1 523      2 565        0        14 s → 27 s
04_DTE_seed20_mm2              1 859      3 297        1        10 s → 23 s
10_DMM_Maverick_seed30_mm2       290        568        0        4 s → 34 s

測試資料：B. napus 染色體 A1-A3（約 100 Mb），-t 40 -c --dedup
--strict-n --format gff3。「打亂」：同樣的染色體，每 10 kb 內把字母
打亂（真的元素已經不在），所以在那裡找到的 hit 是隨機的；exact 版本
是 0 個。有標籤和沒標籤的檔案結果相同。

TIR 較短（5-15 bp）的 pattern 不要加 mismatch：同一個測試中，1 個
mismatch 讓隨機 hit 的比例升到 8-33 %（12-15 bp），或讓 hit 大多是隨機
的（5-8 bp）。mm2 檔可以和 exact 檔一起跑，或取代它；一起跑時，用
--merge 去掉被找到兩次的元素。
