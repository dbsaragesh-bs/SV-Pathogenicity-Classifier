# Study guide: the whole project from zero (with viva answers)

Written for someone with a computer-science / ML background and no biology. Every part builds on the one before. Numbers marked (data) come from your own files; design reasons marked (LOG Dn) come from `LOG.md`'s decision register.

---

# PART 0. The project in one paragraph

DNA can be treated as a very long text. Some edits to that text (large deletions or duplications) cause disease; others are harmless. Doctors want software that looks at such an edit and says "harmful" or "harmless". Edits can be detected from two kinds of sequencing machines: **short-read** (Illumina) and **long-read** (PacBio HiFi). Long-read machines are much more expensive, so the question is: **does feeding long-read information to the classifier make it measurably better?** We answered this with a controlled computer simulation and found: **no, not once the edit has been detected**. Long reads mainly help *find* edits that short reads miss.

---

# PART 1. Biology from zero

## 1.1 DNA, genome, chromosome
- **DNA** is a molecule that stores information as a sequence of four letters: **A, C, G, T**. Think of it as a text file with a 4-symbol alphabet.
- A **genome** is the full DNA of an organism. The human genome is about **3.1 billion letters**.
- It is split into **23 pairs of chromosomes** (chr1 ... chr22, plus X and Y). Think of them as 24 big files (chr1 is about 248 million letters).
- You have **two copies** of each autosome (chr1-22): one from each parent. These two copies are called **haplotypes**. This matters: a variant can be on one copy only (**heterozygous, "het"**) or on both (**homozygous, "hom"**). All variants in this project are het.

## 1.2 Genes, exons
- Only about 1-2 percent of the text is **genes** (instructions to build proteins). A gene is made of **exons** (the parts that code for the protein) separated by introns (which are cut out).
- If an edit damages an exon, the protein can break. This is the main way a deletion causes disease.
- **CDS** (coding sequence) means the exon letters that actually get translated into protein.

## 1.3 Variants
A **variant** is a place where a person's DNA differs from the standard reference.
- **SNV**: one letter changed. Tiny.
- **Indel**: a few letters inserted or deleted.
- **Structural variant (SV)**: a change of **at least 50 letters** (50 is the conventional cut-off). Types:
  - **DEL** (deletion): a chunk is missing.
  - **DUP** (duplication): a chunk is repeated.
  - INS (insertion) and INV (inversion) exist but are **not** used here (LOG D4).

## 1.4 Pathogenic vs benign
- **Pathogenic**: causes or strongly contributes to disease.
- **Benign**: harmless.
- Clinical labs classify variants into five tiers (pathogenic, likely pathogenic, uncertain, likely benign, benign). We use only the extremes and merge "likely" into the main class: P + LP = pathogenic (label 1), B + LB = benign (label 0). "Uncertain" is discarded because we have no reliable label.

## 1.5 Why SVs matter
SVs affect far more letters than SNVs, so one SV can delete a whole gene. They are also much harder to detect than SNVs, which is why the technology choice matters.

---

# PART 2. How DNA is read: sequencing

## 2.1 The core problem
No machine can read a whole 3-billion-letter genome end to end. All technologies work like this:
1. Break the DNA into fragments.
2. Read each fragment, producing a **read** (a string of letters plus a quality score per letter).
3. Repeat until each position in the genome has been read many times. The average number of times each position is covered is the **coverage** (or depth), written like **30x**.
4. A computer puts the reads back in place (alignment, Part 4).

Analogy: shred 30 copies of the same book and read every scrap, then reassemble.

## 2.2 Illumina (short reads)
- **How it works (sequencing by synthesis):** DNA fragments are attached to a glass slide and copied into dense clusters. Fluorescent letters are added one cycle at a time, and a camera records the colour for each cluster each cycle. After 150 cycles you have 150 letters from each cluster.
- **Paired-end:** each fragment (insert) is read from **both ends**, giving two reads that face each other. The middle of the fragment (roughly 100 letters here) is not read. The distance between the two reads is the **insert size** (about 400 here, set as `frag_mean: 400`, `frag_sd: 50`). That distance is a clue for SV detection.
- **Read length:** about 150 letters (here 2 x 150).
- **Accuracy:** very high, about 0.1 percent errors, mostly substitutions.
- **Cost:** cheapest per letter. This is why it is the clinical standard.
- **Weakness:** a 150-letter read cannot cross a big SV, and it cannot be placed uniquely in repeated regions.
- **Why we chose it:** it is the most widely used clinical technology; it is the baseline that long reads must beat.

## 2.3 PacBio HiFi (long reads)
- **How it works:** DNA fragments are turned into circles. A polymerase goes around the same circle many times, reading the same molecule over and over. The multiple passes are combined into one **consensus** read (called **circular consensus sequencing, CCS**). Random errors cancel out.
- **Read length:** 10-25 thousand letters (here median about 14.8 kb, range 3-40 kb).
- **Accuracy:** above 99.9 percent (here 0.4 percent error is simulated, which is conservative). "Hi-Fi" means high fidelity.
- **Cost:** several times more expensive per letter than Illumina.
- **Strength:** a read can cross a whole SV and its flanking unique sequence, so it shows the SV directly and also maps uniquely in repeats.
- **Why HiFi, not ONT (Oxford Nanopore)?** ONT reads are longer but have higher error (about 5-6 percent in this config). HiFi is accurate and a popular clinical long-read choice, and SV callers are best validated on it. The simulator has an ONT mode configured, but ONT was **not run** (LOG: "HiFi only").

## 2.4 Side by side

| | Illumina short | PacBio HiFi long |
|---|---|---|
| Read length | 150 | about 15,000 |
| Error rate | about 0.1 percent | about 0.1-0.4 percent |
| Cost per base | low | high |
| Sees whole SV in one read | no | yes (if SV < read length) |
| Unique mapping in repeats | poor | good |
| Typical use | routine clinical sequencing | hard regions, SVs |

---

# PART 3. Reference genome and public databases

## 3.1 The reference genome: GRCh38
- A **reference genome** is one agreed-upon, assembled human sequence that everyone compares against. Coordinates like "chr1:12,007,233" mean "position 12,007,233 on chromosome 1 of the reference".
- **GRCh38** = "Genome Reference Consortium Human build 38" (released 2013). "hg38" is the same sequence with the UCSC name.
- The older build is **GRCh37 / hg19**. The two give different coordinates for the same location. **Never mix builds.**
- **Why GRCh38 (LOG D5):** every source we use (ClinVar, gnomAD v4.1, GENCODE v50) provides GRCh38 coordinates.
- **What we used as the file (LOG D7):** the Broad Institute's `Homo_sapiens_assembly38.fasta`, restricted to chr1-22, X, Y, M. We removed the alt/decoy contigs, because those extra sequences create ambiguous alignments (low mapping quality artefacts).
- **FASTA** is the file format for such a reference (see Part 8).

## 3.2 Databases used (what each one is for)

| Source | What it is | Used for | Why this one |
|---|---|---|---|
| **ClinVar** (NCBI) | Public archive where labs submit "this variant is pathogenic/benign" with evidence and a **review status** (0-4 stars) | pathogenic label and some benign | The main open source of expert-labelled variants. We require at least 1 star (criteria provided). |
| **dbVar** (NCBI), file `nstd102` | The structural-variant counterpart of ClinVar | checking how precise the ClinVar coordinates are (confidence interval widths), used only for filtering | ClinVar SV coordinates often come from low-resolution arrays, so wide intervals are dropped (`max_ci_frac: 0.10`). |
| **gnomAD-SV v4.1** | Catalogue of SVs found in tens of thousands of presumably healthy genomes, with **allele frequency (AF)** | benign label (PASS, AF at least 1 percent) | A variant that is common in the healthy population is very unlikely to cause severe disease. |
| **GENCODE v50** | Gene annotation: where genes, exons, CDS are | gene-context features | The standard human gene catalogue. |
| **UCSC tables** (rmsk, genomicSuperDups, simpleRepeat, gap, centromere) | Locations of repeats, segmental duplications, tandem repeats, gaps, centromeres | repeat and "difficult region" features | Repeats and segmental duplications are exactly where short reads fail. |
| **gnomAD constraint** (pLI, LOEUF) | How intolerant each gene is to losing function | gene sensitivity features | A deletion in a loss-intolerant gene is more likely harmful. |

**Important detail:** allele frequency (AF) *defines* the benign label for gnomAD variants, so it must **never** be a feature, otherwise the model would just read the answer (LOG D2). This is called **leakage**.

## 3.3 Why the sources are somewhat biased (a viva favourite)
- Pathogenic examples all come from ClinVar; 84 percent of benign come from gnomAD. A model could in principle learn "is it a ClinVar-style variant" instead of "is it harmful".
- ClinVar pathogenic SVs were submitted because someone suspected disease, so they sit in gene-rich regions.
- Mitigations: size and type matching; a "source swap" test (train on gnomAD-benign, test on ClinVar-benign: AUROC 0.89 / 0.96); a ClinVar-only labels control.

---

# PART 4. Alignment and variant calling

## 4.1 Alignment
Take a read, find where on the reference it fits best, and record that. The output is a **BAM** file.
- **bwa-mem** (short reads): the standard aligner for Illumina. It finds exact seeds then extends them.
- **minimap2** (long reads, preset `map-hifi`): the standard long-read aligner, built for long reads and tuned for HiFi.
- **MAPQ** (mapping quality, 0-60): the aligner's confidence that the read is in the right place. Low in repeats.
- **CIGAR string**: a compact description of how the read matches the reference. `150M` = 150 matches. `687M1I599M1D153M...` = matches with a 1-letter insertion (`I`) and 1-letter deletion (`D`). A big deletion looks like `500M1187D700M`.
- **Soft clip (S)**: part of the read that doesn't match here (typically because it belongs on the other side of a breakpoint).

## 4.2 What clues does an SV leave in the alignments?
| Clue | Name | Technology |
|---|---|---|
| Pair lands further apart than normal insert size (about 400) | **discordant pair** | short |
| One read is split: first half aligns here, second half far away | **split read / clipped read** | short and long |
| Coverage halves inside a deletion | **read depth** change | both |
| One long read contains a gap compared to the reference | **in-read deletion** | long |

## 4.3 Variant calling
A **caller** scans the alignments for clues and reports suspected SVs as a **VCF** file.

| Caller | For | Idea |
|---|---|---|
| **Manta** | short | Illumina's caller; combines discordant pairs and split reads, then assembles the breakpoint. |
| **Delly** | short | Classic caller using discordant pairs and split reads. |
| **Sniffles2** | long | Reads signals directly from long-read alignments (split reads and in-read gaps). |
| **cuteSV** | long | Collects signatures from long reads and clusters them. |

**Why two callers per technology:** they use different algorithms, so combining them improves detection and lets us create a feature "how many callers found it". It also avoids the conclusion depending on one tool's quirks.

Real example from your data, SV001822 (truth: 12,007,233-12,008,420, 1,187 bp):
- Delly (short reads) reported `DEL END=12008421 SVLEN=-1188 PE=15 SR=3 PRECISE`: found, within 1 base of the truth.
- Sniffles2 (long reads) reported `Sniffles2.DEL` at position 12,007,233: found.

## 4.4 Matching calls to truth
A caller's output must be matched back to the planted SV. We use **reciprocal overlap (RO)**: two intervals match if each covers at least X percent of the other. This is also used to deduplicate and group the truth set. We wrote our own matching code (plan used Truvari; LOG deviations). Detection was stable at RO 0.3 / 0.5 / 0.8.

---

# PART 5. Why simulation (the central design decision)

## 5.1 The problem
We want the **same SV** observed by both technologies with a **known true answer**.
- Pathogenic SVs in ClinVar are just coordinates: **no reads exist**.
- Real benign reads: short and long reads come from partly different samples and labs.
- So a fair real-data comparison is impossible (LOG D1).

## 5.2 The solution
Build a fake genome piece with the SV planted, then generate both read types from it. Technology is then the only thing that differs.

## 5.3 Step by step

**(a) Window.** Take the reference from 50 kb before to 50 kb after the SV (`flank: 50000`, LOG D8). About 101,000 letters for SV001822. Why a window, not the whole genome? Compute: 4,030 whole genomes is impossible; 100 kb is enough context for aligners and for computing depth normalisation.

**(b) Two haplotypes.** Make `h0` (reference plus tiny random typos) and `h1` (the same, plus the deletion). Background typos (about 1 per 1,000 letters, 90 percent SNV, 10 percent indel) are added identically regardless of label so the model cannot use them.

**(c) Short reads.** ART (`art_illumina`) is an established Illumina simulator with an empirical error profile (`HS25` = HiSeq 2500, 150 bp). It cuts 400-letter fragments and reads 150 from each end. Coverage about 30x; we simulate 20 percent extra then discard reads from abnormal-GC fragments to mimic real GC bias (`overshoot: 1.2`, `gc_mu: 0.45`).

**(d) Long reads.** An in-house simulator (LOG D6): random start, length from a lognormal distribution (median about 14.8 kb), random strand, uniform errors (0.4 percent: 25 sub / 35 ins / 40 del), constant quality. Alternative PBSIM3 was rejected for speed, control and seeding.

**(e) Independence.** Short and long reads use different seeds (`stable_seed(sv_id, rep, tech)`), so they are independent samples of the same haplotypes. **Short reads are not derived from long reads.**

**(f) Batches (pseudo-samples).** Aligners and callers expect a whole sample, e.g. Manta estimates the insert-size distribution from all reads. So up to 150 windows (at least 1 Mb apart) are combined into one pseudo-sample (batch) of 27-33x coverage. 29 batches cover all 4,030 SVs.

**(g) Label-blindness.** All simulation parameters are identical across classes. Test (KS test on parameters, all p > 0.27) confirmed this.

## 5.4 What simulation costs us
Conclusions only hold for the simulated world (Part 11).

---

# PART 6. Building the dataset (truth set)

Goal: a balanced, clean, leak-free list of 4,030 labelled SVs.

1. Start: 247,113 ClinVar DEL/DUP rows on GRCh38 autosomes.
2. Keep P/LP vs B/LB only; at least 1 review star.
3. Add gnomAD-SV: PASS, AF at least 1 percent, DEL/DUP (benign).
4. Drop wide-CI coordinates (imprecise).
5. Keep size 50 bp to 1 Mb (`min_svlen`, `max_svlen`; 1 Mb is a compute cap).
6. Window must lie inside the chromosome and have under 1 percent N (unknown bases).
7. Remove **conflicts**: a pathogenic and a benign SV overlapping at least 50 percent (they cannot both be right).
8. **Deduplicate** near-identical SVs (RO at least 0.8).
9. **Cap local density** (at most 12 SVs per overlapping region, LOG D10) because ClinVar hotspots such as 22q11 would dominate.
10. **Match** classes on type x size bins (bins are 0.25 decades wide, i.e. log10 size) and keep only the overlapping range. Why: otherwise "bigger = pathogenic" is a free shortcut. After matching, median sizes are 977 vs 955 bp, KS D = 0.028 (p = 0.42).
11. Subsample to about 2,000 per class (compute, LOG D9): final **2,015 + 2,015 = 4,030**, 3,766 DEL and 264 DUP.

**Group id:** SVs overlapping at least 50 percent share a `group_id`, and a group is never split between train and test.

---

# PART 7. Features (what the model sees)

One row per SV, 245 columns, in families (registry: `src/features/registry.py`):

| Prefix | Meaning | Arm |
|---|---|---|
| `a_`, `s_` | short-read caller fields (quality, PE/SR counts), short alignment evidence (supporting reads, clipped, discordant, depth ratio, MAPQ), plus `a_ann_*` gene/repeat annotation recomputed from the **short-read called** coordinates | A |
| `b_`, `l_` | long-read caller fields, long alignment evidence (support, breakpoint dispersion, allele fraction), microhomology at long-read breakpoints, plus `b_ann_*` annotation from the **long-read called** coordinates | B adds these to A |
| `o_` | the same annotation from the **true planted** coordinates | Oracle |

Annotation features include: exons hit, CDS overlap, fraction of CDS affected, number of genes disrupted, gene constraint (LOEUF, pLI), repeat fractions (Alu, L1, simple repeats, segmental duplications), GC content, distance to gaps and TSS.

**Excluded on purpose** (FORBIDDEN list): allele frequency, source, review stars, cluster/batch/replicate ids, absolute coverage, size-derived IDs, and truth-derived QC fields. Reason: each either encodes the label or a simulation artefact.

**Why annotation is recomputed per technology (LOG D3):** if both arms used true coordinates, a long-read benefit from better breakpoints could never show up.

---

# PART 8. File formats (also likely in your viva; Task 1 of your course)

| Format | What it holds | Structure | Where in this project |
|---|---|---|---|
| **FASTA** | sequences without quality | `>name` line, then letters | the reference genome |
| **FASTQ** | reads **with** quality | 4 lines per read: `@name`, sequence, `+`, quality string | simulator output (deleted after calling to save disk, `--rm-fastq`) |
| **SAM/BAM** | aligned reads | tab-separated: name, flag, chrom, pos, MAPQ, CIGAR, mate info, sequence, qualities, tags. BAM is the compressed binary form | `short.bam`, `hifi.bam` |
| **VCF** | variants | header `##...`, then one line per variant: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample | caller outputs |
| **Parquet** | the feature table | columnar compressed table | `features_all.parquet` |

Real short-read BAM line (your data):
`SV001822_h1_9674  97  chr1  12006923  60  150M  =  12007...  tlen=1651` (flag 97 = paired, mate reverse, first in pair; MAPQ 60; the pair spans 1,651 letters instead of about 400, a discordant pair).

Real VCF line (Delly): `chr1 12007233 DEL00000000 A <DEL> 1046 PASS PRECISE;SVTYPE=DEL;END=12008421;SVLEN=-1188;PE=15;SR=3 ...`.

**Storage:** FASTQ is the biggest (text), BAM is compressed (about 2-5x smaller, random-access with an index `.bai`), VCF is tiny (only differences). In your data, one batch BAM pair was about 480 MB.

---

# PART 9. Machine learning design

## 9.1 Task
Binary classification: feature row -> pathogenic (1) or benign (0).

## 9.2 The experiment: ablation
An **ablation** removes a component to see its value.
- **Arm A**: short-read features only.
- **Arm B**: A + long-read features.
- **Oracle**: true coordinates (an upper bound for "perfect breakpoints").
- Key quantity: **dAUPRC = score(B) - score(A)**, with confidence intervals.

## 9.3 Models and why several
| Model | Why included |
|---|---|
| **XGBoost** (gradient-boosted trees, on GPU) | Best general choice for tabular data; main model. |
| **Logistic regression** (elastic net) | Simple linear baseline; shows result is not due to complexity. |
| **Random forest** | A second, different tree ensemble. |
| **MLP** (neural net, PyTorch, on GPU) | Checks whether a neural model changes the conclusion. |

All four agree: no model shows a positive gain with CI excluding zero. Agreement across model families is what makes the null result credible.

## 9.4 Fair comparison rules
- **Same data splits** for every arm and model.
- **Same hyperparameter search budget** per arm (Optuna, inner 3-fold CV, objective AUPRC). Otherwise B could win just by being tuned more.
- **Preprocessing inside training folds only** (median imputation plus "missing" indicators plus quantile transform), to avoid test leakage.
- **Missing values:** a short caller may not detect an SV, leaving NaN. The model gets a missing indicator, since missingness itself is information.

## 9.5 Splits
- **Outer:** 5 folds of **held-out chromosome groups**. Train on four groups of chromosomes, test on the unseen fifth.
- Why by chromosome and by group: SVs near each other or overlapping are near duplicates; random splits would let the model memorise.
- **Inner:** grouped, stratified 3-fold CV for tuning.

## 9.6 Metrics
- **AUPRC** (area under the precision-recall curve): primary. Focuses on how well the positive class is ranked. 1.0 is perfect.
- **AUROC**: secondary; the 50/50 balance makes it meaningful too (0.5 = chance).
- Also MCC, F1, recall at precision 0.9, Brier score (calibration).
- **Cluster bootstrap** (4,000 resamples, on GPU): resample whole overlap-groups, giving confidence intervals that respect the dependence between near-duplicates. **Paired** differences compare A and B on the same samples.

## 9.7 Interpretation
- **SHAP (TreeSHAP):** attributes each prediction to features. Top features: exons hit, CDS overlap, fraction of CDS affected, constraint, repeat content.

## 9.8 Why a GPU
RTX 4060: XGBoost (`device=cuda`), MLP training, bootstrap statistics and TreeSHAP are all massively parallel.

## 9.9 Infrastructure choices
- **WSL2 (Linux on Windows):** Manta, bwa, minimap2, Delly and others are Linux-only.
- **micromamba envs:** reproducible tool versions (Manta needs Python 2.7, Delly clashes on boost, so separate envs).
- **Snakemake:** workflow manager, reruns only missing steps.
- **Seeds from SV ids:** the same SV always gets the same random draws, so runs are reproducible.

---

# PART 10. Results and why

**Main (n = 3,510, detected by both):** AUPRC A 0.984, B 0.985, Oracle 0.986. dAUPRC(B-A) = +0.0004 [-0.0021, +0.0038], one-sided p = 0.42.

**Why no gain:**
1. Pathogenicity here is mostly "does it cut an exon": 97.3 percent of pathogenic hit an exon vs 13 percent of benign (data).
2. Exon hit depends only on the start and end coordinates.
3. Both technologies locate them within about 1 bp in the simulation.
4. Exon status changed for **0 percent** of SVs between short and long coordinates; Oracle only +0.002.
5. Arm A is already at 0.984, leaving little room.

**Where long reads do help: detection.** HiFi detects 94.5 percent vs short 91.6 percent of planted SVs; the gap is mostly benign SVs at segmental duplications (short reads detect only 52 percent there vs 93 percent elsewhere). Keeping undetected SVs in gives dAUPRC +0.0036 [+0.0021, +0.0054] (the only robustly positive effect).

**Depth sweep:** HiFi detection 97.8 / 97.4 / 94.8 / 76.4 percent at 30 / 20 / 10 / 5x; classifier AUPRC stays flat. Doubling short reads to 60x raises short detection only 91.5 -> 92.9 percent.

**Leakage audit:** size-only, type-only, chromosome-only, position-only, coverage/batch-only, label-shuffled and pipeline-column models all score AUROC 0.48-0.55 (chance).

---

# PART 11. Limitations (state these before your professor does)

1. **Simulation is cleaner than reality.** On real HG002 data, a classifier separates real from simulated evidence at AUROC 0.97-0.98. Real long-read breakpoints wobble about 5-10 bp; simulated ones do not. Real data might give long reads a larger advantage.
2. **Easy benchmark.** Gene overlap nearly decides the label, so AUPRC about 0.98 leaves no headroom.
3. **Label biases:** pathogenic only from ClinVar; benign mostly gnomAD.
4. **Scope:** DEL and DUP only (264 DUP), heterozygous only, one replicate, HiFi only (no ONT), no phasing.
5. **Large benign CNVs are rare**, so large pathogenic CNVs are under-represented.
6. **No real pathogenic reads exist** to test on.
7. **Small real-data check:** 104 deletions, only 11 short-read detections.
8. **Not clinical.**

---

# PART 12. Likely viva questions with answers

**Q1. Why simulate instead of real data?**
Real pathogenic SVs have no reads, and real benign short/long reads come from different samples; a matched comparison is impossible. Simulation gives the same SV seen by both technologies with a known truth.

**Q2. How do you know the simulation is realistic? What did you lose?**
We tested it on 104 real HG002 deletions: real vs simulated evidence is easily separable (AUROC 0.97), mainly because real data has noisier breakpoints and lower MAPQ. So it is optimistic about breakpoint precision. We lose realism of errors and artefacts, so conclusions hold for the simulated regime.

**Q3. Why GRCh38 and not GRCh37 or T2T?**
All annotation sources (ClinVar, gnomAD v4.1, GENCODE v50) are available in GRCh38 coordinates. Mixing builds would put SVs in the wrong places. T2T is newer but these databases are not natively on it.

**Q4. Why Illumina and HiFi? Why not ONT?**
Illumina is the clinical short-read standard. HiFi has high accuracy and well-tested callers. ONT has higher errors and different caller settings; its simulator profile is in the config but it was out of scope.

**Q5. Why two callers per technology?**
Different algorithms catch different SVs; it increases recall and gives a "number of callers" feature; it avoids dependence on one tool.

**Q6. Why not use allele frequency as a feature?**
Because the benign label from gnomAD is defined by AF at least 1 percent. The feature would reveal the label (leakage).

**Q7. Why split by chromosome and cluster?**
Overlapping or neighbouring SVs are near-duplicates; a random split would leak them across train and test.

**Q8. Why AUPRC and not accuracy?**
AUPRC measures ranking quality for the positive class and is robust under class imbalance. Accuracy depends on the threshold. (Here classes are balanced by design, but the metric is standard for rare-positive problems and we also report AUROC.)

**Q9. Why match on size?**
So size cannot be used as a shortcut. Size-only models on held-out chromosomes score AUROC 0.48-0.53.

**Q10. Why is long-read gain null? Is it because long reads are useless?**
No. Once an SV is detected, both technologies give about the same coordinates, and pathogenicity here is determined by gene overlap that comes from coordinates. Long reads help detection (94.5 vs 91.6 percent) and would likely help more on real data, where short-read breakpoints are fuzzier.

**Q11. Why is the Oracle arm included?**
It gives the ceiling for any breakpoint-precision benefit. Even perfect coordinates give only +0.002, showing no headroom.

**Q12. Why did you use several models?**
To show the conclusion does not depend on the learner. XGBoost, LR, RF and MLP all agree.

**Q13. Why run a depth sweep?**
To answer the practical question: is it worth paying for more long-read coverage? Detection saturates around 10-20x; classifier AUPRC is flat.

**Q14. Why batches of windows rather than one window per BAM?**
Aligners and callers estimate sample-wide statistics (e.g. insert-size distribution), so they need realistic multi-SV samples.

**Q15. What would you do next?**
Re-test with real long-read call sets; add INS/INV and ONT; make the simulator more realistic (chimeras, homopolymer errors, quality); try a harder benchmark not dominated by gene overlap; add replicates.

---

# One-page memory aid

```
ClinVar + gnomAD-SV  ->  4,030 labelled DEL/DUP (balanced, size-matched, cleaned)
        |
        v   plant each SV into a 100 kb GRCh38 window (two haplotypes, one has SV)
   h0, h1 strings
        |-- ART (Illumina 2x150, 30x)   -> bwa-mem  -> Manta + Delly
        |-- own HiFi simulator (15kb,30x)-> minimap2 -> Sniffles2 + cuteSV
        v
 245-column feature table (short / long / gene context / oracle)
        v
 Arm A (short)  vs  Arm B (short+long)  vs  Oracle
 XGBoost / LR / RF / MLP, chromosome-held-out CV, cluster-bootstrap CIs
        v
 dAUPRC(B-A) = +0.0004 (null). Long reads help detection (94.5 vs 91.6 %), not classification.
```
