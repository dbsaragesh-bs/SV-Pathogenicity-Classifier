# plan.md — Step-by-Step Build Plan

**Project:** SV Pathogenicity Classifier — Short-Read vs Short+Long-Read Ablation (simulation-first)
**Read first:** `context.md`. **Log everything in:** `LOG.md`. **Learning guide:** `LEARN.md`.
**Priority (set by owner):** dataset quality above all. Every phase ends with a quality gate. A phase is not done until its gate passes and `LOG.md` is updated.

---

## 0. How to read this plan

- Each phase has: **Goal, Steps, Outputs, Quality gate, Log/learn actions**.
- Gate failures are not skipped. Fix, or record a decision in `LOG.md` and ask the owner if scope changes.
- Numbers marked **(verify)** are expectations or assumptions to confirm in Phase 0, not facts.
- Defaults live in `configs/*.yaml`; no thresholds hard-coded in code.

### Phase map

| Phase | Name | Rough effort (part-time) |
|---|---|---|
| P0 | Feasibility, environment, pilot | 4–5 days |
| P1 | Truth-set curation (labels) | 4–5 days |
| P2 | Simulation engine | 5–6 days |
| P3 | Alignment and SV calling | 4–5 days |
| P4 | Truth matching and dataset QC | 3–4 days |
| P5 | Feature engineering | 5–6 days |
| P6 | Dataset freeze and leakage audit | 3 days |
| P7 | Modelling | 4–5 days |
| P8 | Ablation, depth sweep, interpretation | 4–5 days |
| P9 | Real-data realism check | 3 days |
| P10 | Report, slides, learning guide | 4–5 days |

Total about 6–8 weeks part-time. **Minimum viable path** (if time is short): DEL only, one caller per technology, two held-out chromosomes, no depth sweep. Order of cuts if time runs short: (1) ONT and robustness settings, (2) depth sweep grid reduced to long {10, 30} × short {30}, (3) second caller per technology, (4) INS/INV, (5) phasing, (6) P9 reduced to a feature-distribution comparison on a few regions. Never cut: label-blind simulation, leakage audit (P6), grouped splits, Arm A/B/Oracle comparison.

---

## 1. Principles (senior-ML-engineer rules for this project)

1. **Data first.** Model choice moves metrics by points; a leak or label error moves them by tens of points.
2. **Label-blind pipeline.** Nothing upstream of the features may know the label.
3. **One source of truth per fact.** Configs, manifests, and seeds define everything. No manual edits to data.
4. **Everything deterministic.** Seeds derived from `hash(sv_id, replicate, tech)`; re-running gives identical bytes.
5. **Test the pipeline, not only the model.** Unit tests for SV application, coordinates, matching; integration test on a tiny dataset; sanity checks on every output.
6. **Assume leakage until disproved.** Run adversarial and shortcut tests before trusting any metric.
7. **Report uncertainty.** Cluster-bootstrap CIs, paired differences, effect sizes. No single-number claims.
8. **Log as you go.** If it is not in `LOG.md`, it did not happen.

---

## 2. Open decisions (defaults used unless the owner changes them)

| # | Question | Default | Needed by |
|---|---|---|---|
| Q1 | Compute: cores, RAM, free disk on the Windows machine? Cloud/HPC allowed? | Assume 8 cores, 16 GB RAM, 300 GB disk, no cloud | P0 |
| Q2 | Deadline (date)? Sets scope | None given; plan assumes ~8 weeks | P0 |
| Q3 | Target dataset size | ~3,000 pathogenic + ~3,000 benign unique DEL/DUP after cleaning, 2 replicates each | P1 |
| Q4 | Include INS/INV in main analysis? | No (exploratory only) | P1 |
| Q5 | Long-read technology | HiFi primary, ONT as secondary sweep | P2 |
| Q6 | Is the report graded on a real-data result as well? | Real data only as realism check (P9) | P9 |
| Q7 | Team size / who owns what? | Treat as one person plus AI help | P0 |

---

## 3. Repository layout

```
D:\Projects\DNA\
├── context.md  plan.md  LOG.md  LEARN.md  dna.md
├── configs/
│   ├── data_sources.yaml        # verified URLs, versions, builds, checksums
│   ├── truthset.yaml            # label rules, filters, size cap
│   ├── simulation.yaml          # window, coverage, read models, SNP rate
│   ├── calling.yaml             # tool versions + params
│   ├── features.yaml            # feature definitions and arm membership
│   └── experiment.yaml          # splits, models, seeds, metrics
├── envs/environment.yml  envs/environment.lock.txt
├── workflow/Snakefile  workflow/rules/*.smk
├── src/
│   ├── truthset/    # download, parse, harmonise, liftover, clean, cluster
│   ├── sim/         # apply_sv.py, make_window.py, simulate_reads.py
│   ├── calling/     # wrappers, merging
│   ├── matching/    # truvari wrappers, label attachment
│   ├── features/    # annotation.py, short_evidence.py, long_evidence.py
│   ├── validation/  # schema, QC, leakage tests
│   ├── modelling/   # splits, pipelines, tuning, evaluation
│   └── utils/
├── tests/           # unit + integration
├── data/
│   ├── raw/         # untouched downloads + checksums
│   ├── interim/     # truth sets, windows, reads, BAMs, VCFs (large, not committed)
│   └── processed/   # frozen feature tables (parquet), manifests
├── reports/         # data_audit/, ablation/, figures/, model_cards/
└── notebooks/       # exploration only; no logic lives only here
```

---

## P0 — Feasibility, environment, pilot

**Goal:** prove the approach works on a tiny scale and confirm every external dependency before building anything big.

**Steps**
1. **Environment.** Install WSL2 Ubuntu. Create the conda/mamba env from bioconda: `samtools bcftools bedtools minimap2 bwa-mem2 art pbsim3 manta delly sniffles cutesv truvari survivor whatshap pyfaidx pysam cyvcf2 pybedtools snakemake`; Python libs: `pandas numpy scikit-learn xgboost lightgbm shap optuna pandera pyarrow matplotlib`. Record all versions in `envs/environment.lock.txt`. Check which tools actually install (verify; Manta needs Python 2.7-era dependencies in some builds).
2. **Resource check (Q1).** CPU/RAM/disk numbers into `LOG.md`; this fixes dataset size.
3. **Source verification.** For each source below record URL, version, genome build, size, date accessed, SHA-256:
   - Reference: GRCh38 no-alt analysis set + indexes (bwa-mem2, minimap2, faidx).
   - Pathogenic: ClinVar SV records / dbVar ClinVar SV study (accession to verify, not assumed).
   - Benign: gnomAD-SV (GRCh38 release to verify), optionally ClinVar benign SV, optionally a long-read-resolved population callset (HGSVC / 1000G long-read; verify build).
   - Annotation: GENCODE (match build), RepeatMasker, segmental duplications, simple repeats, mappability, gaps/centromeres/blacklist, gene constraint (LOEUF/pLI), haploinsufficiency/triplosensitivity predictions.
   - Realism-check data: GIAB HG002 short-read and long-read data and SV benchmark (verify availability, size; may be too big, then use a region subset).
4. **Counts feasibility (no reads yet).** Parse sources, apply draft filters, report counts by label, type, size band, chromosome. Thresholds: proceed only if DEL+DUP has at least about 2,000 pathogenic and 2,000 benign after cleaning (verify; otherwise ask owner: lower target, add sources, or widen types).
5. **Pilot simulation (20 SVs).** Run the full chain on 10 pathogenic-like and 10 benign DEL/DUP: plant → simulate both read types → align → call → match. This exposes tool failures early.
6. **Windows vs batches vs whole chromosome.** Pilot compares: (a) one window per job, (b) 100 windows pooled into one pseudo-sample, (c) a whole small chromosome (chr21). Metric: recall of planted SVs, runtime, caller stability. Decision recorded in `LOG.md` (D-entry). Reason: callers such as Manta estimate global statistics (insert size, depth) and may behave badly on sparse input.
7. **Calibrate simulators against real data.** Pull read-length, insert-size, error, GC-bias profiles from HG002 (or any real run you can access) and set ART/PBSIM3 parameters to match (verify). Record in `simulation.yaml`.

**Outputs:** `configs/data_sources.yaml` (verified), `envs/environment.lock.txt`, pilot report `reports/pilot.md`, decisions D6+ in `LOG.md`.

**Gate P0**
- All tools run end to end on the pilot.
- Pilot recall of planted DEL > 90% for long reads at 30x and > 70% for short reads (verify, expected from literature; if far lower, the simulation is broken, not the idea).
- Counts feasibility meets thresholds, or owner approves a revised target.
- Build compatibility table filled: every source on one build, or a liftover plan with expected loss.

**Log/learn:** write the resource table, the verified-source table, the pilot decision, and a short "what each tool does" note into `LEARN.md` §2/§6.

---

## P1 — Truth-set curation (the labels)

**Goal:** a small, clean, unambiguous, deduplicated list of SVs with labels. This is the highest-leverage step for dataset quality.

**Steps**
1. **Harmonise** every source into one schema (`sv_id, chrom, start, end, svtype, svlen, label, source, source_id, evidence_flags, ci_start, ci_end, build`). Normalise chromosome names, SV types, coordinate base (VCF 1-based vs BED 0-based; unit-tested).
2. **Build lift** if any source differs from GRCh38. Record fails and multi-maps; drop them; log counts before/after.
3. **Pathogenic filters.**
   - Significance = Pathogenic or Likely pathogenic. Exclude VUS, conflicting, benign classes.
   - Review status at least one star with assertion criteria (verify exact field names). Optionally a stricter "gold" tier (expert panel / multiple submitters).
   - Types: DEL and DUP (copy-number gain/loss) primary.
   - Coordinate quality: drop records with large start/stop uncertainty (array-resolution). Threshold set from the data distribution, logged. Reason: planted coordinates must be a believable truth.
   - Size: 50 bp to cap (default 5 Mb; decide in P0).
   - Remove records overlapping assembly gaps, centromeres, or acrocentric short arms.
4. **Benign filters.**
   - Primary benign: gnomAD-SV (or equivalent) with AF above threshold (default 1%) AND PASS AND not flagged low quality (verify fields). Add ClinVar Benign / Likely benign SVs at one star or better. Record `benign_basis` (AF-based vs clinical) for analysis only; never a feature.
   - Same type and size and gap filters as pathogenic.
5. **Conflict and noise removal.**
   - Remove any benign SV with at least 50% reciprocal overlap with a pathogenic SV (and vice versa); log both removed sets.
   - Remove SVs that are multiallelic or flagged complex.
6. **Deduplicate and cluster.**
   - Merge exact and near-duplicate records (same type, reciprocal overlap at least 80%, breakpoint distance threshold) with Truvari collapse or equivalent.
   - Recurrent CNVs (e.g., 22q11.2) appear many times in ClinVar; keep one representative plus a `cluster_id` and `cluster_size`. Cluster size is NOT a feature.
   - Final `cluster_id` is the unit for later splitting.
7. **Distribution matching (anti-shortcut).** Pathogenic SVs are larger and are mostly DEL/DUP; benign are smaller. Build a **benign sampling plan** stratified by (type × size band × chromosome-group) so that size/type alone cannot separate classes. Keep two versions:
   - `matched`: benign downsampled to match the pathogenic type×size distribution (main modelling set).
   - `natural`: unmatched, for evaluation under realistic prevalence only.
8. **Genotype assignment (label-independent).** Assign zygosity (het or hom-alt) to each SV by a random draw from the same distribution for both classes (default: 100% het for the main run; hom as a sensitivity run), using a fixed seed. Reason: ClinVar SVs are mostly het, population SVs are mixed; real differences must not leak.
9. **Spacing for batches.** Place SVs on pseudo-samples so SVs in a batch are far apart (default ≥ 2 Mb, no shared window) and chromosome/label composition is balanced per batch.
10. **Insertion/inversion handling (exploratory only).** Insertions: inserted sequence for BOTH classes drawn from the same sequence pool (e.g., real mobile-element consensus plus real benign long-read insertion sequences matched by length) so sequence composition cannot reveal the label. Inversions: only if counts allow.

**Outputs:** `data/processed/truthset_v1.parquet`, `data/processed/truthset_manifest.json` (counts at every filter, thresholds, tool versions, hashes), `reports/data_audit/truthset.md` with plots (size by label, type by label, chromosome by label, cluster-size histogram).

**Gate P1**
- Pandera schema passes (no nulls in keys, valid types, start < end, svlen consistent with coordinates, chromosome in allowed set).
- Zero duplicates under the clustering rule across classes.
- Zero benign–pathogenic conflicts remaining.
- After matching: pathogenic and benign type and size distributions are statistically similar (e.g., KS test p > 0.05 or standardized mean difference < 0.1 per size band; logged).
- Size-only, type-only, chromosome-only baseline models on the truth set give AUPRC near class prevalence (test with cross-validation by cluster). If a trivial baseline is strong, the dataset still has a shortcut; fix the sampling before continuing.
- Counts meet the Q3 target, or owner is informed.

**Log/learn:** record the full filter waterfall (n after each step) in `LOG.md`; write the "why each label rule exists" table into `LEARN.md`.

---

## P2 — Simulation engine

**Goal:** generate reads that are identical in process for every SV, so that technology is the only varying factor.

**Steps**
1. **SV application module** (`src/sim/apply_sv.py`). Own, small, deterministic code (pysam/pyfaidx) to construct the mutated haplotype for DEL, DUP (tandem), INV, INS. Unit tests with hand-computed expected sequences, edge cases (SV at window edge, adjacent to N, overlapping repeat), VCF-coordinate off-by-one tests. Cross-check a random 50 SVs against an independent simulator (VISOR or Sim-it) and require identical resulting sequences.
2. **Window construction.** Window = SV ± flank (default 100 kb) from GRCh38; haplotype 1 = reference, haplotype 2 = reference + SV (het) or both mutated (hom). Add **background variation** identically for all windows: random heterozygous SNVs and small indels at realistic density (default 0.1%) with fixed seeds. Required for realistic alignment and for the phasing feature. Assumption: DUP is simulated as tandem duplication (orientation unknown in ClinVar; state in limitations).
3. **Read simulation.**
   - Short: ART Illumina, 2×150, fragment 400 ± 50 (calibrated), ~30x diploid-total; GC bias and error profile from the calibrated model.
   - Long: PBSIM3 HiFi/CCS model (mean about 15 kb, accuracy about Q30+) at ~30x; ONT model as a secondary sweep.
   - Coverage is drawn per pseudo-sample from a label-independent distribution (e.g., 25–35x) to avoid a fixed-depth artifact.
   - Seeds = `hash(sv_id, replicate, tech)`.
   - Replicates: R=2 by default (same SV, different read noise); all replicates share `cluster_id` and group.
4. **Sweep arms for later.** Generate extra long-read depths (5, 10, 15, 20x) and short-read depths (15, 30, 60x) by simulating at 60x and sub-sampling reads deterministically (so lower depths are true subsets of higher depths).
5. **Negative controls.** Generate reference-only windows (no SV) matching the same count distribution (about 10% of windows). Used to measure false-positive calls and to train/evaluate the callers' error behavior.
6. **Sanity checks per batch** (automated, saved as `reports/data_audit/sim_qc.md`):
   - Achieved coverage within ±10% of target (counted from FASTQ).
   - Insert-size and read-length distributions match config.
   - Fraction of reads spanning the planted breakpoint matches expectation (het 50%).
   - Per-class distributions of every simulation parameter identical (KS tests) — label-blindness proof.

**Outputs:** `data/interim/windows/*.fa`, `reads/{short,long}/*.fastq.gz`, `sim_manifest.parquet` (seed, parameters, coverage, haplotype, hash).

**Gate P2**
- Unit and cross-simulator tests pass.
- Simulation parameter distributions do not differ between pathogenic and benign (all p > 0.05 or SMD < 0.05).
- Re-running a sample of 20 windows reproduces identical FASTQ hashes.
- Realism: simulated read-length, insert-size, quality, GC distributions within tolerance of the calibration targets.

**Log/learn:** save the parameter table and the "what the simulator models and does not model" list in `LEARN.md` (limits: no sample-prep artifacts, no real-chimeric reads beyond model, uniform GC bias approximations).

---

## P3 — Alignment and SV calling

**Goal:** produce call sets for each technology with the same procedure for every SV.

**Steps**
1. **Align to the FULL GRCh38** (not just the window) so multi-mapping in repeats and segmental duplications behaves realistically. Short: bwa-mem2 + samtools sort/markdup; long: minimap2 (`map-hifi`/`map-ont`). Index once.
2. **QC BAMs:** mapped fraction, MAPQ distribution, coverage inside the window, insert-size distribution. Fail the batch if outside tolerance.
3. **Call SVs** with two callers per technology (short: Manta + DELLY; long: Sniffles2 + cuteSV), run per batch with fixed parameters from `calling.yaml`. Keep all calls (including low quality) and their evidence fields.
4. **Merge callers** per technology (Truvari/SURVIVOR) into a consensus plus per-caller columns. Keep "supported by k callers" as a quality field (a legitimate feature of the technology arm).
5. **Haplotagging (long):** call small variants in the window (e.g., from the simulated SNV truth or a lightweight caller), phase with WhatsHap, tag reads. Gives phased flag, haplotype-block id, and allele-specific support.
6. **Restrict analysis to window:** discard calls outside the SV window (batch artifacts).

**Outputs:** BAMs, per-caller VCFs, merged VCFs, `calls_{short,long}.parquet`.

**Gate P3**
- Negative-control windows: false-positive call rate documented (calls per Mb), and no systematic dependency on the label (labels do not exist there; this is checked via later nuisance tests).
- Recall of planted SVs by technology and type within pilot expectations; if not, investigate simulator/caller parameters before continuing.
- Same tool versions and parameters for all batches (manifest check).

**Log/learn:** table of caller recall/precision by type and size (also an interesting result for the report), and a one-page "how each caller decides it saw an SV" note in `LEARN.md`.

---

## P4 — Truth matching and dataset QC

**Goal:** attach exact ground truth to each call and quantify detection, without contaminating features with the truth.

**Steps**
1. **Match calls to planted SVs** with Truvari bench (type-aware; thresholds on size similarity, sequence/breakpoint distance, reciprocal overlap; sensitivity at 3 stringency levels).
2. **Per-SV outcome table:** `detected_short`, `detected_long`, `n_calls_matched`, `best_call_id`, `breakpoint_error_short/long` (for QC and analysis only, never a feature).
3. **Cohort definition:**
   - `C_both`: detected by both technologies (headline comparison).
   - `C_all`: all planted SVs (detection component).
4. **Detection-vs-label audit.** Detection rate by label within type and size bands. If detection differs by label after conditioning on type and size, investigate (e.g., pathogenic SVs sit in more repetitive regions), report, and stratify analysis.
5. **Unmatched calls (false positives):** keep in a separate table; excluded from the classifier, but counted and reported (honest scope: classifier assumes a true SV call).
6. **Truth-label integrity check:** manually inspect (IGV or pileup plot) a random sample of about 30 SVs per class/type/technology, including a few failures.

**Outputs:** `data/interim/matching/*.parquet`, `reports/data_audit/matching.md`.

**Gate P4**
- Matching stable across the 3 stringency settings (changes in matched counts small and documented).
- Detection by label difference explained or removed.
- Manual inspection finds no systematic planting/matching error.

---

## P5 — Feature engineering

**Goal:** build clean, leakage-free, well-documented features for Arm A, Arm B, and Oracle.

**Principles:** features come from alignments, calls, reference, and public annotations only. No truth, no AF, no source, no clinical text, no cluster size, no replicate id, no batch id, no coverage knobs from the simulator.

**Feature families**

| Family | Examples | Arm A | Arm B | Oracle |
|---|---|---|---|---|
| Call geometry | svtype, log(called length), position relative to chromosome ends | short call | long call | truth |
| Breakpoint quality | CI width (short), spread of supporting-read breakpoints (long) | CI | spread | 0 |
| Short evidence | discordant pairs, split reads, depth ratio, soft-clip consensus, MAPQ of supporting reads, GC bias | yes | yes | n/a |
| Long evidence | supporting reads, allele fraction, mean MAPQ, strand balance, consensus junction sequence quality | no | yes | n/a |
| Caller agreement | number of callers, caller-specific flags | yes | yes | n/a |
| Local sequence context | GC flank, repeat fractions (Alu/L1/LTR/simple/tandem), segdup overlap, mappability, distance to gap | from short coords | from long coords | truth |
| Junction mechanism | microhomology length, homology/repeat at breakpoints | no (not resolvable) | yes | truth |
| Gene annotation | n genes, n exons, CDS overlap fraction, TSS distance, gene-dosage/constraint scores (LOEUF/pLI, haploinsufficiency/triplosensitivity predictions) | short coords | long coords | truth |
| Insertion sequence (exploratory) | inserted length, GC, repeat class | limited | yes | n/a |
| Phasing | phased flag, haplotype-block size | no | yes (needs SNVs) | n/a |

**Steps**
1. **Annotation recomputation:** for each arm, compute gene/exon/CDS/regulatory overlap from that arm's coordinates (D3). For Oracle, from the truth. Use pybedtools/bedtools with explicit coordinate-convention tests.
2. **Gene-list safety:** do not use ClinGen/ClinVar-derived pathogenic gene lists as features (circular). Allowed: constraint metrics derived from population data. Run a sensitivity analysis with and without gene-dosage predictors.
3. **Evidence features from BAMs** with pysam (independent of caller where possible): counts of discordant/split reads, depth ratio between SV interior and flanks, supporting-read breakpoint dispersion.
4. **Missingness policy:** two kinds are distinguished and encoded separately.
   - *Structural* (feature not applicable, e.g., long-read evidence for an SV the long reads missed): keep NaN + `detected_X` indicator; trees handle NaN; for logistic regression use fixed constant + indicator, fitted inside the training fold.
   - *Technical* (tool failed): audit and fix, do not impute silently.
5. **Transformations (inside pipelines only):** log1p for counts/length, rank-gauss or standardization for LR, winsorize at training quantiles. Never drop large SVs as "outliers" (they are the biology).
6. **Feature registry:** `configs/features.yaml` lists every feature with definition, source, arm membership, type, allowed range, missingness rule. A test enforces that the code outputs exactly these columns.
7. **Collinearity and redundancy:** correlation clusters, VIF for LR; keep feature groups intact for ablation (do not drop within a group by correlation).

**Outputs:** `data/processed/features_armA.parquet`, `features_armB.parquet`, `features_oracle.parquet`, plus `feature_dictionary.md`.

**Gate P5**
- Schema/range tests pass; no inf; missingness map matches the policy.
- **Forbidden-feature test:** automated check that no column is in the forbidden list and no column has |correlation| > 0.95 with the label or with truth-derived QC fields.
- Feature distributions between labels compared and visualised; nothing implausible (e.g., a feature that separates classes perfectly).

---

## P6 — Dataset freeze and leakage audit

**Goal:** freeze a dataset only after trying hard to break it.

**Tests** (each saved in `reports/data_audit/leakage.md`)
1. **Nuisance-only model:** predict label from coverage, depth knob, batch id, window id, replicate, read count, chromosome, flank GC, caller versions. Expect AUPRC ≈ prevalence; AUROC ≈ 0.5.
2. **Shortcut baselines:** size-only, type-only, chromosome-only, position-only, "number of callers" only. Must be near chance.
3. **Label-shuffle test:** permute labels within cluster-balanced strata; full pipeline AUROC ≈ 0.5.
4. **Adversarial validation:** classifier that distinguishes pathogenic rows from benign rows using only simulation/pipeline-derived columns (no biology columns). Must fail. Also classifier simulated-vs-real (against HG002 calls, in P9).
5. **Duplicate/leak scan:** check no cluster spans folds; check no near-identical rows (cosine/Hamming on feature vectors) across folds.
6. **Temporal/source split sanity:** e.g., train on one benign basis (AF-based) and test on the other (ClinVar-benign) to see if performance collapses (source shortcut).
7. **Feature-importance sanity:** quick tree model; top features must be explainable biologically or mechanically; anything surprising is investigated before proceeding.

**Steps after tests pass**
- Freeze: write `data/processed/dataset_v1/` with parquet files, `manifest.json` (all hashes, config hashes, seeds, tool versions, counts), `DATA_CARD.md` (provenance, intended use, known limitations, label policy, collection process, known biases).
- Version with Git (small files) + DVC or checksums (big files). Tag `dataset_v1`.
- **Open the final test partition only once**; before that it is represented only by hashes.

**Gate P6**
- All six leakage tests pass or have documented, accepted explanations.
- Frozen dataset reproducible from configs (spot rebuild of 20 SVs yields identical rows).
- `DATA_CARD.md` complete. Owner signs off in `LOG.md`.

---

## P7 — Modelling

**Goal:** a trustworthy, comparable training harness, then models.

**Design**
- **Splits.** Outer test = held-out chromosomes chosen to balance counts by label and type (e.g., pick 3–4 chromosomes; list in `experiment.yaml`; sex chromosomes excluded). Inner = StratifiedGroupKFold (groups = `cluster_id`, stratified by label × type) for tuning. Repeat the outer split over several different chromosome sets (leave-chromosome-group-out) for variance estimates. **Identical folds for A, B, Oracle.**
- **Imbalance handling (priority order).**
  1. Matched benign sampling (done in P1) to remove size/type shortcuts.
  2. Class weights (`class_weight`, `scale_pos_weight`).
  3. Replicate-level augmentation already in the dataset (legitimate, grouped).
  4. Threshold chosen on validation folds for a target (e.g., precision 0.9 or max F-beta); never on test.
  5. Evaluate under both matched and natural prevalence; show prevalence-adjusted PPV.
  6. No SMOTE by default (synthetic points in this feature space have no biological meaning); only as a labelled sensitivity run inside training folds.
- **Models.** Logistic Regression (L2/elastic net), Random Forest, XGBoost and/or LightGBM. Same hyperparameter budget for every arm (Optuna, nested in inner CV, fixed trials, fixed seed).
- **Pipelines.** scikit-learn Pipelines: imputation/indicators, transforms, model — fitted only on training folds. No leakage through preprocessing.
- **Calibration.** Isotonic/Platt fit on a held-out calibration slice of the training data; report reliability curves and Brier score.
- **Metrics.** Primary AUPRC. Also AUROC, MCC, F1/F2 at the chosen threshold, recall at fixed precision, precision at fixed recall, confusion matrix, calibration. Per type, per size band, per repeat context (segdup/non-segdup).
- **Uncertainty.** Cluster-bootstrap CIs (resample clusters, not rows); across-fold mean ± std.
- **Experiment tracking.** One JSON per run: git hash, dataset hash, config hash, seed, fold ids, metrics, artifacts. (MLflow optional.)

**Steps**
1. Reproduce shortcut baselines (should be near chance) on the final split.
2. Baseline models on Arm A. Sanity: does it beat baselines by a believable margin?
3. Tune within inner CV, retrain on full training set, evaluate on outer folds.
4. Repeat for Arm B and Oracle with the same code path (only the feature-set config changes).
5. Save models, predictions per row (with fold ids), and a model card.

**Gate P7**
- Pipeline passes unit tests; permutation of labels gives chance performance.
- Arm A beats shortcut baselines. If performance is suspiciously high (e.g., AUROC > 0.98 with weak features), stop and re-run leakage audit.
- Variance across outer splits reported; no tuning on test data (verified by code review/log).

---

## P8 — Ablation, depth sweep, interpretation (the scientific result)

**Goal:** answer "how much, and why".

**Experiments**
1. **Main comparison on C_both:** Arm A vs Arm B vs Oracle. ΔAUPRC and ΔAUROC with cluster-bootstrap CIs and paired bootstrap p-values; effect sizes; the same for MCC and recall@precision. Correct for multiple comparisons across models/metrics (Holm or just state primary metric = AUPRC, model = pre-specified best-by-inner-CV).
2. **Ablation ladder:** A → A + long-read coordinates only (annotation recomputed) → + long-read evidence → + junction/microhomology → + phasing → full B → Oracle. Shows which ingredient carries the gain.
3. **Leave-one-feature-group-out** from B.
4. **Detection component:** run on C_all with detection flags; report how much of long-read value is "finding more SVs" vs "describing found SVs better".
5. **Depth sweep:** long depth {5, 10, 15, 20, 30x} × short depth {15, 30, 60x}. Plot ΔAUPRC vs added long-read coverage; compare to a short-read-only 60x arm to answer "is long-read worth it versus simply more short reads". Include a rough relative-cost note (state assumptions; do not present as precise pricing).
6. **Stratified analysis:** by SV type, size band, repeat/segdup context, zygosity, gene-overlap status, and by whether the breakpoint shift changed exon overlap.
7. **Technology robustness:** ONT instead of HiFi; second short-read caller set; hom-alt vs het. Conclusions should not hinge on one setting.
8. **Why analysis (flip table).** For test SVs where B is right and A is wrong (and vice versa), tabulate the cause: exon status changed by breakpoint shift, gene count changed, repeat context, evidence-based features, other. SHAP (TreeExplainer) for global and local explanations on best Arm B model; grouped importance by feature family. Distinguish "model uses it" (SHAP) from "removing it hurts" (ablation).
9. **Confound controls:** repeat with size-removed models (and size-adjusted metrics), without gene-dosage predictors, with the stricter pathogenic "gold" tier, and with the clinical-benign basis only.
10. **Negative result handling:** if ΔAUPRC CI includes zero, report the null with power discussion (minimum detectable effect given the sample size).

**Outputs:** `reports/ablation/*.md`, figures: PR curves A/B/O; ablation ladder bars; depth-sweep curve; stratified Δ plots; SHAP summary; flip-table; calibration.

**Gate P8**
- Every claim in the report maps to a saved table/figure and a run id.
- Results stable across at least 2 robustness settings (or instability is reported).
- No test-set peeking after the freeze (log audit).

---

## P9 — Real-data realism check (HG002)

**Goal:** confirm that the simulator is not teaching the model something unreal.

**Steps**
1. Take real HG002 short-read and HiFi/ONT data over the same genomic windows used in simulation (or a subset; size permitting). Run the **same** alignment, calling, and feature code.
2. **Feature-distribution comparison:** simulated vs real for matched benign SVs (GIAB benchmark truth): KS/Wasserstein per feature; adversarial validation AUC (simulated vs real) on technology-evidence features. Report which features differ and by how much.
3. **Detection comparison:** recall/precision of callers on real vs simulated at the same depth.
4. **Apply the trained models** to real HG002 benign SVs: report predicted-pathogenic fraction (a false-positive proxy); expect low if the model is not simulator-specific. Real pathogenic SVs at scale are not available; say so.
5. If large gaps exist, tune simulator parameters (not the model) and re-run P2–P6, logging a new dataset version.

**Gate P9:** domain gap documented; any feature with severe shift is either fixed in the simulator, removed with justification, or flagged as a limitation.

---

## P10 — Report, slides, learning guide

1. **Report** (structure from `dna.md` §47, adapted): problem, related work (StrVCTVRE, simulation benchmarks, SV callers), data and simulation design, label policy, leakage controls, models, results, why-analysis, limits, conclusion that leads with the measured gain.
2. **Slides:** 12–14; most time on the main result and the "why".
3. **DATA_CARD.md, MODEL_CARD.md** finalized.
4. **LEARN.md** completed: results, re-run commands with timings, personal notes, final viva Q&A.
5. **Reproducibility test:** fresh clone → environment → `snakemake` on a 50-SV mini config reproduces the mini results; document runtime and disk.
6. **Final LOG.md entry** summarizing the project, decisions, deviations from this plan.

**Gate P10:** another person can reproduce the mini run from the README alone; report numbers match saved artifacts.

---

## 4. Quality gates summary (what "clean and accurate dataset" means here)

| Dimension | How we ensure it |
|---|---|
| Label correctness | Strict ClinVar significance + review status; conflict removal; coordinate-quality filter; manual spot checks |
| No duplicates / near-duplicates | Clustering at 80% reciprocal overlap; cluster-level splits |
| No shortcut | Matched type×size sampling; trivial-baseline tests |
| No leakage | AF/source/clinical text/gene-list exclusion; nuisance, shuffle, adversarial, source-swap tests |
| Label-blind generation | Same parameters for both classes; KS tests on every sim parameter |
| Realism | Simulator calibration to HG002; P9 distribution comparison |
| Correct ground truth | Unit-tested SV application; cross-simulator agreement; manual inspection |
| Traceability | Seeds, manifests, hashes, config hashes, DATA_CARD |
| Balance without distortion | Matched sampling, class weights, replicate augmentation; no default SMOTE |
| Missing data | Typed missingness, indicators, fold-internal imputation |
| Reproducibility | Snakemake, lockfile, deterministic seeds, mini-run test |

## 5. Risk register

| Risk | Impact | Mitigation |
|---|---|---|
| Too few clean pathogenic DEL/DUP | Weak statistics | Relax review status tier with sensitivity analysis; add sources; accept DEL-only |
| Sparse-window input breaks callers | Bad recall | Pilot (P0) decides batch vs whole-chromosome design |
| Simulated reads unrealistic | Result not credible | Calibration + P9 check; robustness to ONT/HiFi and callers |
| Compute/disk too large | Delays | Windows not whole genome; delete FASTQ after BAM (keep seeds to regenerate); store CRAM; downscale N |
| Callers fail on large CNVs (> few Mb) | Biased detection | Size cap; separate depth-based CNV caller if needed (e.g., CNVnator/Control-FREEC) logged as extension |
| Detection differs by label | Confounded | P4 audit and stratification |
| Pathogenic coordinates imprecise (arrays) | Label/coordinate noise | CI filter; gold tier; limitation statement |
| Tool version drift | Irreproducible | Lockfile, container option (Docker/Singularity) |
| Gene-level annotation circularity | Inflated performance | Exclude ClinGen-like lists; sensitivity run |
| Null result | Anticlimax | Pre-registered as acceptable; power analysis; value in "why/depth" analysis |

## 6. Logging protocol (applies to every phase)

1. After each task: append a `LOG.md` entry (template inside) with exact numbers and output paths.
2. After each decision: add a row to the Decision register.
3. After each phase: update the Status board, record the gate result, add concepts learned to `LEARN.md`.
4. Every figure and table in the report cites a run id from the log.

## 7. Definition of done

- [ ] Dataset v1 frozen with manifest, DATA_CARD, and all P6 leakage tests passed
- [ ] Arm A, Arm B, Oracle models trained under identical protocol
- [ ] ΔAUPRC (and other metrics) with CIs, ablation ladder, depth sweep, stratified results
- [ ] Why-analysis (flip table + SHAP + ablation)
- [ ] Real-data realism check documented
- [ ] Report, slides, MODEL_CARD
- [ ] LOG.md complete; LEARN.md complete
- [ ] Fresh-clone mini reproduction passes

## 8. Immediate next actions

1. Owner answers §2 Q1, Q2 (compute, deadline); others use defaults.
2. Start P0: WSL2 + environment; source verification; pilot of 20 SVs.
3. Update `LOG.md` after each step.
