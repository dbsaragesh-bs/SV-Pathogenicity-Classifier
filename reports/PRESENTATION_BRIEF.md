# Presentation brief: the whole project, and why long reads did not win

Everything below comes from files in this repo. Numbers marked **[new]** come from `scripts/diagnose_why_no_gain.py` (outputs in `reports/diagnostics/`), run on 2026-10-14 to answer "why is the long-read gain zero?". Everything else is from `reports/REPORT.md`, `LOG.md`, `reports/data_audit/`.

---

## 0. The 30-second story

1. A **structural variant (SV)** is a large DNA edit (>= 50 letters): a deletion (DEL) or a duplication (DUP) here.
2. Goal: a classifier that says *pathogenic* or *benign* for an SV. We train it two ways: **Arm A** sees only what **short-read** (Illumina) sequencing gives; **Arm B** sees A **plus** what **long-read** (PacBio HiFi) sequencing gives.
3. Real data cannot do this fairly (pathogenic SVs have no reads), so each SV is **planted into a simulated 100 kb piece of the human genome**, and both read types are simulated from it. Technology is the only thing that differs.
4. Result: B is **not** better than A (dAUPRC = +0.0004, CI -0.002 to +0.004). Long reads help *find* more SVs (94.5 % vs 91.6 %), not *classify* them better.
5. Reason: in this dataset the label is almost decided by **"does the SV cut an exon?"**, a fact that depends only on the two breakpoint coordinates, and in the simulation both technologies give those within about 1 bp. When we make short-read breakpoints realistically sloppy, the long-read gain appears (section 6).

---

## 1. The data: where the 4,030 rows come from

| Piece | Source | Role |
|---|---|---|
| Pathogenic SVs | ClinVar (Pathogenic / Likely pathogenic, >= 1 review star) | label = 1 |
| Benign SVs | gnomAD-SV v4.1 (PASS, allele freq >= 1 %) + ClinVar Benign/Likely benign | label = 0 |
| Reference genome | GRCh38 (chr1-22, X, Y, M) | what reads are simulated from and aligned to |
| Genes / exons / CDS | GENCODE v50 | gene-context features |
| Repeats, segdups, gaps | UCSC tables | repeat features |
| Gene constraint | gnomAD v4.1 (pLI, LOEUF) | "how bad is losing this gene" features |

**Filtering waterfall** (`data/processed/truthset_manifest.json`): 247,113 ClinVar DEL/DUP rows → 4,841 pathogenic + 32,496 benign of usable size (50 bp to 1 Mb) → remove benign/pathogenic conflicts (overlap >= 50 %) and duplicates (overlap >= 80 %) → cap hotspots (<= 12 overlapping windows) → **match size and type** → **2,015 pathogenic + 2,015 benign = 4,030** (3,766 DEL, 264 DUP; 132 DUP per class).

What the final table looks like (`data/processed/truthset_matched.parquet`, 4,030 x 23): `sv_id, chrom, start, end, svtype, svlen, label, source, stars, group_id, zygosity (all het) ...`.

Composition you should know (data):
- Benign: 1,698 gnomAD + 317 ClinVar. Pathogenic: 100 % ClinVar.
- Size is matched on purpose: median 977 bp (pathogenic) vs 955 bp (benign), KS D = 0.028, p = 0.42.
- 50/50 balanced by design (so AUPRC baseline = 0.5 here).

---

## 2. The pipeline, step by step (what happens to one SV)

```
ClinVar/gnomAD row (chr10:119,676,453-119,677,500 DEL, label 1)
   |  plant into window = SV +/- 50 kb of GRCh38; two haplotypes: h0 = reference + 0.1 % random SNV/indel, h1 = same + the deletion
   v
   |-- short reads : ART Illumina 2x150, fragment 400+-50, ~30x (GC-biased)          --> bwa-mem   --> Manta + Delly
   |-- long reads  : in-house HiFi simulator, median ~15 kb, 0.4 % error, ~30x       --> minimap2  --> Sniffles2 + cuteSV
   |   (different random seeds, so short reads are NOT derived from long reads)
   v
 up to 150 windows pooled into one "pseudo-sample" (batch) so callers see a realistic sample (29 batches)
   v
 match each caller call back to the planted SV (reciprocal overlap) -> detected_short / detected_hifi flags
   v
 features (245 columns): short evidence, long evidence, gene/repeat annotation computed from each tech's CALLED coordinates, oracle annotation from TRUE coordinates
   v
 models: XGBoost / LogReg / RandomForest / MLP, 5 outer folds = held-out chromosome groups
```

Key design points to say out loud:
- **Label-blind simulation**: every simulation parameter is drawn from the same distribution for both classes (KS tests, all p > 0.27). The simulator never knows the label.
- **Annotation recomputed per technology (decision D3)**: Arm A computes "which genes/exons does this SV hit" from the *short-read called* start/end; Arm B additionally from the *long-read called* start/end. This is the only way a better breakpoint can turn into a benefit.
- **Oracle arm**: same annotation from the true planted coordinates; the ceiling for "perfect breakpoints".
- **Not used as features** (forbidden list): allele frequency, source, review stars, batch, replicate, absolute coverage, truth-derived QC columns.

---

## 3. What the model is actually given (training data)

One row per SV; label in the last column; 245 columns in the file, of which the model arms use:

| Arm | # columns used | What is in them |
|---|---|---|
| A (short) | 87 | caller fields (Manta/Delly quality, PE/SR counts, CI width), alignment evidence (discordant pairs, split reads, clipped reads, MAPQ, depth ratio inside vs flank), **40 annotation columns from short-called coordinates** |
| B (short+long) | 167 | all of A + long-read caller fields (Sniffles/cuteSV support, allele fraction), long-read evidence (supporting reads, breakpoint dispersion, strand balance, MAPQ, depth ratio), microhomology at junction, **40 annotation columns from long-called coordinates** |
| Oracle | 125 | evidence from both + 40 annotation columns from the true coordinates |

Real rows (same size ~1.1 kb deletion, both detected by both technologies):

| feature | pathogenic SV001488 | benign SV004747 |
|---|---|---|
| exons hit (short-called coords) | **1** | **0** |
| CDS bases hit (log) | 6.7 | 0 |
| max fraction of a gene's CDS removed | 0.47 | 0 |
| min LOEUF (gene constraint, low = intolerant) | 0.54 | 1.0 |
| Manta/Delly callers agreeing | 2 | 2 |
| Delly paired-end / split-read support | 8 / 5 | 8 / 0 |
| short-read mean MAPQ | 60 | 52.4 |
| Sniffles supporting reads | 18 | 12 |
| long-read allele fraction | 0.48 | 0.49 |
| short / HiFi breakpoint error (bp) | -1 / -1 | 0 / -1 |

Look at the last rows: long-read and short-read evidence look the same for both rows. Only the gene-context rows differ. That is the whole story in one table.

**Missing values**: if a caller misses an SV, its columns are NaN. Trees use NaN directly; other models get a median fill + a "was missing" indicator, fitted **inside the training fold only**. Heavy-tailed columns are quantile-transformed inside the fold. Large SVs are never dropped as outliers.

**Training protocol**: 5 outer folds, each = a group of whole chromosomes held out (never random rows; overlapping SVs share a `group_id` and stay together). Hyper-parameters tuned by Optuna in an inner grouped CV with the same budget for every arm. Metric: AUPRC (+ AUROC, MCC, calibration). Confidence intervals: cluster-bootstrap (resample groups), paired A vs B.

---

## 4. Results in five numbers

| | A short | B short+long | Oracle |
|---|---|---|---|
| AUPRC (XGBoost, n = 3,510 both-detected) | 0.984 | 0.985 | 0.986 |
| dAUPRC B - A | +0.0004 [-0.0021, +0.0038] | | |
| Detection of planted SVs | 91.6 % | 94.5 % (97.8 % at 30x in the depth sweep subset) | |
| Short 60x vs long 30x detection (8 batches) | 92.9 % | 97.8 % | |
| Real HG002 realism check (104 deletions) | short 10.6 % found | HiFi 64.4 % found | simulation is clearly cleaner than reality |

All four model families (XGB, LR, RF, MLP) agree: no positive gain with CI above zero.

---

## 5. Why long reads did not win: the actual reasons (verified)

Your instinct is right in principle: long reads capture SVs better. The experiment just measures something else. Seven reasons, in order of importance.

### Reason 1. The label is almost "does it hit an exon" **[new numbers]**
- Share of SVs overlapping a coding exon: **97.3 % of pathogenic vs 13.0 % of benign**.
- One binary feature, "hits an exon", alone gives AUROC 0.935 / AUPRC 0.905.
- The 40 annotation columns alone give **AUPRC 0.984**; adding all 127 other columns gives 0.985.
- Why the dataset is like this: ClinVar pathogenic SVs were labelled pathogenic *because* they remove genes/exons; common gnomAD SVs are common precisely because they avoid genes. The label is partly *defined* by gene overlap, so the task is easy and near-circular. That is a property of the real databases (it is also how StrVCTVRE-type tools get high scores), not a bug of the pipeline, but it leaves no headroom.

### Reason 2. Exon-hit depends only on two numbers: start and end
Evidence features (read counts, MAPQ, depth) do not enter the exon calculation. The only way a long read can help the classifier is by giving better start/end coordinates.

### Reason 3. In the simulation both technologies already get the breakpoints right
Breakpoint error among detected SVs **[new]**: short median 1 bp, p90 2 bp, p99 47-52 bp; HiFi median 1 bp, p99 8-10 bp. The exon-hit status computed from short-called coordinates differs from truth in **4 of 3,693** SVs (0.1 %). Long: 4 of 3,808.
Why so clean: the ClinVar coordinates are planted exactly, a simulated read has an exact breakpoint, and Delly/Manta then find split reads at exactly that base. Real data is noisier (section 7).

### Reason 4. Perfect breakpoints (Oracle) add almost nothing
Oracle - A = +0.002 (CI includes 0). If even true coordinates cannot help, nothing a long read provides about coordinates can.

### Reason 5. Evidence features are about "was this SV detected well", not "is it harmful" **[new, corrects an over-simple statement in REPORT.md]**
- Long-read evidence alone (no annotation): AUPRC **0.60** (near chance for a 50/50 set).
- Short-read evidence alone: AUPRC **0.79**. This is not pathogenicity. It is *genomic context seen through mappability*: short-read mean MAPQ alone has AUROC 0.72, because benign gnomAD SVs sit more in repeats/segmental duplications (repeat fraction AUROC 0.23 for pathogenic, i.e. benign have more) where short reads map badly, and pathogenic SVs sit in unique gene-rich regions where MAPQ is 60.
- Irony worth telling: **long reads map uniquely in repeats, so they erase that context shortcut** (0.60 vs 0.79). Long-read evidence is *less* label-informative than short-read evidence. The extra 80 long-read columns therefore add mostly noise on top of an already saturated annotation signal.
- Careful wording for the professor: the leakage audit tested coverage/batch/size/etc. and they are at chance; it did not test mappability features alone. This 0.79 is a real, explainable correlation of genome context with label, not leakage of the label, but it is a context shortcut inside Arm A, and it should be stated.

### Reason 6. Everything is saturated and the test set is small
AUPRC ~0.984 leaves at most 0.016 of headroom; at n = 3,510 the smallest detectable difference is about +-0.003 overall and about 0.02 inside a 250-SV stratum. Even a true +0.005 gain would be invisible.

### Reason 7. The one place long reads help is *detection*, and we mostly removed it
The headline cohort is "detected by both" (n = 3,510), which throws away exactly the SVs where long reads shine: **288 benign and 49 pathogenic SVs are missed by short reads** (benign at segmental-duplication breakpoints: 52 % detected vs 93 % elsewhere). Keep them (cohort C_all with detection flags, n = 4,030) and dAUPRC becomes **+0.0036 [+0.0021, +0.0054]**, the only robustly positive effect.

(Side note that may come up: short-read detection by label is 85.7 % benign vs 97.6 % pathogenic. This is the segdup effect, not label leakage, since the simulation is label-blind, but it is why both-detected is a slightly selected cohort.)

---

## 6. The experiment that proves the mechanism **[new]** (`reports/diagnostics/noise_stress.png`)

Question: how sloppy would short-read breakpoints have to be before long reads help? Method: keep everything the same, but add Gaussian noise of SD *s* bp to the short arm's coordinates before recomputing its gene annotation; Arm B keeps its ~1 bp long-read coordinates. Same fixed XGBoost for every point, same held-out-chromosome folds, paired cluster bootstrap.

| short-read breakpoint noise SD | exon status wrong in short arm | AUPRC A | AUPRC B | B - A [95 % CI] |
|---|---|---|---|---|
| 0 bp (our simulation) | 0.1 % | 0.985 | 0.985 | -0.0004 [-0.0020, +0.0012] |
| 25 bp | 2.1 % | 0.980 | 0.984 | +0.0041 [+0.0005, +0.0079] |
| 100 bp | 5.4 % | 0.974 | 0.985 | +0.0106 [+0.0055, +0.0158] |
| 250 bp | 8.3 % | 0.962 | 0.985 | +0.0234 [+0.0160, +0.0313] |
| 1,000 bp | 13.5 % | 0.946 | 0.987 | +0.0401 [+0.0308, +0.0502] |
| 5,000 bp | 17.4 % | 0.932 | 0.986 | +0.0540 [+0.0436, +0.0642] |

Reading: the long-read benefit is not zero in principle. It is a **monotone function of short-read breakpoint imprecision**, and it becomes significant as soon as ~2 % of SVs get the wrong exon status (SD about 25 bp, which is within what imprecise short-read-only calls such as read-depth or discordant-pair-only calls produce in real data). Our simulation sits at the left edge because it is too clean.

**Label this honestly: it is a synthetic what-if, not a measurement on real data.** It does not say real data would sit at 100 bp; it says where the effect switches on.

---

## 7. Is the "null" a design flaw or a finding? (be ready for this)

It is a **finding conditional on the simulated regime**, and the report says so. Arguments for the professor:

- **What the experiment is built to measure**: value of long-read *information* for *classification*, holding the SV set fixed. It is not a test of "are long reads better at sequencing", which is already known (and which the detection numbers confirm: 94.5 vs 91.6 %, 97.8 % vs 92.9 % even when doubling short depth).
- **Why the gain is structurally capped here**: classification uses gene context; long reads improve *where the SV is* and *whether it is seen*, not *what genes are there*.
- **Where the design is weak (say it before they do)**:
  1. Simulation is cleaner than reality: on 104 real HG002 deletions an adversary separates real from simulated evidence at AUROC 0.97-0.98; real HiFi breakpoint dispersion is ~5-10 bp vs 0 in simulation; Manta crashed on the real BAM so only Delly was used.
  2. Easy benchmark: label ~ exon overlap (reason 1). A harder benchmark (benign SVs that also hit exons, e.g. from ClinVar only, or non-coding regulatory pathogenic SVs) would give more room. ClinVar-only control: +0.001 [-0.001, +0.004], still null.
  3. Pathogenic and benign come from different databases (84 % gnomAD benign vs 100 % ClinVar pathogenic). Mitigated by size/type matching, source-swap test (AUROC 0.89 / 0.96), ClinVar-only control. Not fully removable.
  4. ClinVar coordinates are array-resolution and were planted as exact; truth coordinates are therefore "nicer" than reality.
  5. No insertions/inversions (long reads' biggest advantage is resolved insertion sequence), no phasing, no ONT, DUP only 264, one replicate, het only.
  6. Realistic short callers also have **breakpoint homology and repeat ambiguity**, which our single-substitution error model underestimates.
- **What *is* solid**: leakage audit (nuisance, size, type, chromosome, position, label-shuffle all at chance), identical folds for all arms, tuning only inside training folds, four model families agree, cluster-bootstrap CIs, and the mechanism is demonstrated, not guessed (sections 5-6).

---

## 8. Suggested slide order (12 slides)

1. Title + one-sentence question.
2. What is an SV; why it matters (deletion cutting a gene).
3. Short vs long reads (150 bp vs 15 kb; cost; what each can see).
4. Why real data cannot answer the question → simulation-first design.
5. Dataset: sources, waterfall, 2,015 + 2,015, size-matched.
6. Pipeline diagram (plant → simulate → align → call → match → features).
7. What the model sees: Arm A / B / Oracle, example row table (section 3).
8. Training protocol: held-out chromosomes, grouped CV, AUPRC, no leakage.
9. Main result (PR curves, dAUPRC ~ 0): `reports/figures/pr_curves.png`.
10. Why null: exon-hit 97 % vs 13 %; breakpoint error 1 bp; oracle +0.002.
11. Where long reads do help: detection (94.5 / 91.6 / 97.8 vs 92.9) + noise-stress curve: `reports/diagnostics/noise_stress.png`.
12. Limitations + next steps (real long-read callsets, INS/INV, harder benchmark).

An existing deck is at `reports/SV_ablation_slides.pptx`; it does not contain the section 5-6 material yet.

---

## 9. Likely questions, short answers

- **"Long reads see SVs better, so why no gain?"** They do see them better (detection +3 points, 5x depth detects 76 %, 10x already 95 %). But the classifier's input is gene context computed from two coordinates, and in the simulation both technologies give those within 1 bp. The effect switches on once short-read breakpoints are imprecise (table in section 6).
- **"Is something wrong with your training?"** Checked: same folds, same tuning budget, preprocessing inside folds, label-shuffle at chance, four model families agree, single-feature and family-only analyses reproduce the same ceiling. Nothing indicates a training defect; the ceiling is the data.
- **"Did you give the long-read arm a chance?"** Yes: B includes annotation from long coordinates, long evidence, microhomology; Oracle gives perfect coordinates; ladder and leave-one-family-out show no single long-read family helps.
- **"Why is short-read evidence alone 0.79 if the sim is label-blind?"** Label-blind simulation does not mean label-independent genome. Pathogenic SVs sit in unique gene-rich DNA, benign in repeats; mappability is a property of the DNA. Long reads remove it (0.60).
- **"So what is the contribution?"** A leakage-audited, matched framework that quantifies long-read value separately as *detection* vs *description*; a demonstration that the classification gain is bounded by gene-context saturation and short-read breakpoint precision; a depth answer (10-20x long saturates detection; doubling short depth barely helps).
