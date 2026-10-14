# context.md — Authoritative Project Context (v2, simulation-first)

**Project:** SV Pathogenicity Classifier — Short-Read vs Short+Long-Read Ablation on a Matched Simulated Dataset
**Course:** DNA Sequencing Technology (semester project)
**Status:** Planning complete, implementation not started (see `LOG.md` for live status)
**Supersedes:** `dna.md` (v1 context pack). Keep `dna.md` for history; where the two disagree, THIS file wins.
**Last updated:** 2026-10-03

Any AI agent or teammate: read this file, then `plan.md`, then the latest entries in `LOG.md`, before doing anything.

---

## 1. One-paragraph summary

We train two classifiers that label a structural variant (SV) as pathogenic or benign. Model A sees only features derivable from short-read sequencing. Model B sees those features plus features that only long reads provide. Both models are trained on the same SVs, with reads that are simulated from the same synthetic genomes, so the only thing that differs is the sequencing technology. The headline result is the measured gain of B over A (ΔAUPRC with confidence interval), why it occurs, and how it changes with sequencing depth.

## 2. Why simulation (the decision that shaped v2)

Real data cannot support the design:
- Pathogenic SVs in ClinVar/dbVar are catalog records (coordinates plus a label). No reads exist for them.
- Real benign samples have short and long reads, but only ~20% of the short reads are derived from the same sample material as the long reads. That makes the comparison confounded.

Fix: plant every SV, benign and pathogenic, into a synthetic genome window. Simulate long reads and short reads from that same window with the same label-blind procedure. Ground truth is then exact and the two technologies are matched by construction.

Cost of this decision (state it honestly in the report): the result is a **controlled-experiment** result. It depends on the simulators' error models and on the SV callers. Real HG002 data is used only as a realism check, never for training.

## 3. Research question and hypotheses

**RQ:** How much does adding long-read-derived information improve SV pathogenicity classification over short-read-derived information alone, and why?

- **H1:** AUPRC(B) > AUPRC(A) on the matched cohort.
- **H2:** Most of the gain comes from more precise breakpoints, which change gene/exon overlap calls (not from the classifier seeing more columns).
- **H3:** The gain is larger for SVs in repeat/segmental-duplication regions and for small SVs (50 bp to 1 kb).
- **H4:** Gain from extra long-read depth plateaus (roughly 10–15x is expected from the literature; we measure it ourselves).
- **Null result is acceptable.** Pathogenicity is driven mostly by gene context, which both arms share. A small or zero gain is a valid, reportable finding.

## 4. Corrections to v1 (`dna.md`) — read these

| # | v1 said | v2 says | Why |
|---|---|---|---|
| 1 | Match real short/long callsets by overlap | Simulate matched reads from planted SVs | Pathogenic reads do not exist; benign mismatch |
| 2 | `population_af` is a feature | `population_af` is **NOT** a feature | Benign label is defined by AF > threshold, so AF leaks the label |
| 3 | Genome build GRCh37 | **GRCh38** (verify in Phase 0) | gnomAD-SV v4, current ClinVar/dbVar, GENCODE are GRCh38 |
| 4 | Annotation features use the VCF coordinates | Annotation features are computed from each technology's **called** coordinates | This is how breakpoint precision turns into a measurable benefit |
| 5 | Phasing as a feature | Phasing only if background heterozygous SNPs are simulated | Without SNPs there is nothing to phase |
| 6 | Source-database features excluded | Also exclude clinical-significance-derived gene lists (e.g., ClinGen HI genes) as features | Circular: those lists are built from pathogenic CNV evidence |
| 7 | Primary types DEL/DUP/INS/INV | Primary: **DEL and DUP**. INS/INV exploratory | ClinVar pathogenic is almost all DEL/DUP; ClinVar insertions have no sequence |

## 5. Definitions

- **SV:** genomic variant ≥ 50 bp (DEL, DUP, INS, INV; TRA out of scope).
- **Unit of analysis:** one planted SV (a row = one SV in one replicate). Replicates of the same SV share a group id and are never split across train/test.
- **Label:** 1 = pathogenic (ClinVar/dbVar Pathogenic or Likely pathogenic, review status passes the quality filter); 0 = benign (common population SV and/or ClinVar Benign/Likely benign, see `plan.md` §Phase 1).
- **Arm A:** features from short-read alignments and short-read SV calls only (plus shared annotation computed from short-read coordinates).
- **Arm B:** Arm A + long-read evidence + annotation recomputed from long-read coordinates.
- **Arm O (oracle):** annotation computed from the true planted coordinates. Upper bound on what better breakpoints could give.
- **Cohort C_both:** SVs detected by both technologies. Used for the headline fair comparison.
- **Cohort C_all:** all planted SVs, with detection flags. Used to quantify the detection component of long-read value.

## 6. Non-negotiable rules

1. Label-blind pipeline: no simulation, alignment, calling, or feature parameter may depend on the label.
2. Nuisance factors (coverage, batch, window id, replicate id, caller version) must be unable to predict the label (tested, not assumed).
3. No feature may be derived from the label, the source database, AF, or clinical text.
4. Split before any fitting (imputation, scaling, feature selection, threshold, calibration).
5. Replicates and near-duplicate SVs stay in the same fold.
6. Never tune on the final test set. The test set is opened once.
7. Never report accuracy as a headline metric. Primary: AUPRC.
8. Never present simulated results as real-data results.
9. Never claim clinical utility.
10. Every filtering step logs counts in `LOG.md` and `reports/data_audit/`.
11. Verify every URL, accession, and version before using it; record access date and checksum.
12. Nothing is fabricated: no invented data, metrics, accessions.

## 7. Key technical choices (defaults, adjustable in `configs/`)

| Item | Default | Notes |
|---|---|---|
| Reference | GRCh38 no-alt analysis set | Align to the full genome even though reads come from windows, to keep realistic multi-mapping |
| Environment | WSL2 Ubuntu + conda/mamba | Manta etc. are Linux-only |
| Workflow | Snakemake | Deterministic, resumable |
| Window | SV ± 100 kb flank, SV length cap 5 Mb | Cap decided in Phase 0 |
| Short reads | ART Illumina, 2×150, ~30x | Calibrated to real HG002 short-read stats |
| Long reads | PacBio HiFi via PBSIM3 (~15 kb mean, ~30x) | ONT as secondary |
| Short callers | Manta + DELLY (merged via Truvari/SURVIVOR) | |
| Long callers | Sniffles2 + cuteSV (merged) | |
| Truth matching | Truvari bench | |
| Models | Logistic Regression, Random Forest, XGBoost/LightGBM | No deep learning |
| Split | Outer: held-out chromosomes. Inner: StratifiedGroupKFold by overlap cluster | Same folds for all arms |
| Primary metric | AUPRC with cluster-bootstrap CI | |

## 8. File map

```
context.md   this file (what and why)
plan.md      step-by-step build plan with quality gates (how)
LOG.md       running log: what was done, decisions, numbers (update after every task)
LEARN.md     learning guide: concepts, reading order, how to re-run, viva Q&A
dna.md       v1 context pack (historical)
```

Project code layout is defined in `plan.md` §3.
