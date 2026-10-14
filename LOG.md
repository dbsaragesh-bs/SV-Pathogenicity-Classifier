# LOG.md — Project Log

**How to use:** append a new entry at the TOP of "Entries" after every task or decision. Never delete old entries; correct them with a new entry that references the old one.

**Entry template**
```
### YYYY-MM-DD — Phase X.Y — short title
- Did: what was done (commands/scripts, versions)
- Why: reason / decision
- Numbers: counts, metrics, file sizes (exact)
- Outputs: file paths + checksum (if data)
- Gate: pass / fail / n/a (and which gate)
- Next: what follows
- Learned: one line worth remembering for LEARN.md
```

## Status board

| Phase | Name | Status | Gate result |
|---|---|---|---|
| P0 | Feasibility and environment | done | gate passed / documented |
| P1 | Truth-set curation (labels) | done | gate passed / documented |
| P2 | Simulation engine | done | gate passed / documented |
| P3 | Alignment and SV calling | done | gate passed / documented |
| P4 | Truth matching and dataset QC | done | gate passed / documented |
| P5 | Feature engineering (Arm A / B / O) | done | gate passed / documented |
| P6 | Dataset freeze and leakage audit | done | gate passed / documented |
| P7 | Modelling | done | gate passed / documented |
| P8 | Ablation, depth sweep, interpretation | done | gate passed / documented |
| P9 | Real-data realism check (HG002) | done (domain gap documented) | gate passed / documented |
| P10 | Report, slides, learning guide | done | gate passed / documented |

## Decision register

| ID | Date | Decision | Alternatives rejected | Reason |
|---|---|---|---|---|
| D1 | 2026-10-03 | Simulate matched short+long reads for both classes | Match real callsets; chop real long reads for benign only | Pathogenic has no reads; real benign short reads are not derived from the long-read sample |
| D2 | 2026-10-03 | Drop `population_af` as a feature | Keep with confound ablation | Benign label defined by AF, so the feature leaks the label |
| D3 | 2026-10-03 | Annotation features computed from called (not true) coordinates | Use true coordinates for both arms | Otherwise breakpoint-precision benefit cannot appear |
| D4 | 2026-10-03 | Primary scope DEL + DUP | Include INS/INV in primary | ClinVar pathogenic is almost all DEL/DUP; ClinVar INS lacks sequence |
| D5 | 2026-10-03 | GRCh38 | GRCh37 | Source availability (verified in P0) |
| D6 | 2026-10-03 | In-house long-read simulator (HiFi/ONT error models) instead of PBSIM3 | PBSIM3 | Controllable, fast (20 windows in 18 s), seeded; limits documented in LEARN.md |
| D7 | 2026-10-03 | Reference = Broad hg38 primary chromosomes only (no alt/decoy) | NCBI/UCSC analysis set (throttled downloads) | Same primary sequence; avoids alt-contig MAPQ artefacts |
| D8 | 2026-10-03 | SV length cap 1 Mb, flank 50 kb | 5 Mb cap, 100 kb flank | Compute: window volume scales with SV size; large benign CNVs are too rare to match anyway |
| D9 | 2026-10-03 | Matched set capped at ~2,000 per class | All 2,359 per class | Compute budget (~29 batches) |
| D10 | 2026-10-03 | Cap local density to <=12 overlapping windows | No cap | ClinVar hotspots (e.g. 22q11, NRXN1) dominate; also bounds the number of pseudo-sample batches needed |

## Entries

### 2026-10-03 — P8.5 follow-up — 60x short-read arm
- Did: simulated an independent second 30x short-read draw (replicate 1, same haplotypes) for batches B000-B007, merged with the 30x BAM, re-called (Manta+Delly), re-featurised.
- Numbers: short detection 91.5 % (30x) -> 92.9 % (60x) vs HiFi 30x 97.8 %; AUPRC A60 0.980, A30 0.979, B 0.981 (n=811).
- Gate: n/a. Learned: extra short depth does not recover segdup-breakpoint SVs.

### 2026-10-03 — P2-P10 — pipeline run end to end, results, real-data check
- Did: simulated 29 batches (4,030 windows; sim ~75 s/batch, minimap2 ~2 min, bwa ~2 min, callers < 1 min); stage-2 features ~4 min/batch. Leakage audit, main experiment (4 arms x 4 models x 5 chromosome-group folds, GPU XGBoost/MLP, 80 fits), ladder, leave-one-family-out, 6 controls, depth sweep (8 batches, L5/L10/L20/S15), SHAP, HG002 realism check (104 deletions, Delly-only for short reads because Manta crashed on the real BAM).
- Numbers: detection short 0.916 / HiFi 0.945; median breakpoint error 1 bp both. Main XGBoost AUPRC A 0.984, B 0.985, O 0.986; B-A +0.0004 [-0.0021, +0.0038]. C_all with detection flags B-A +0.0036 [+0.0021, +0.0054]. Depth: long detection 97.8/97.4/94.8/76.4 % at 30/20/10/5x, classifier AUPRC flat. HG002: adversarial AUROC 0.97 (short) / 0.98 (long); real detection short 0.106, HiFi 0.644.
- Gates: P1 pass, P2 pass, P3/P4 pass (detection stable across stringency; label-detection gap explained by segdup breakpoints), P6 pass (all shortcut/nuisance/shuffle tests at chance), P7 pass (identical folds, no test tuning), P8 done, P9 domain gap documented (simulation too clean), P10 report/cards written.
- Incidents: (1) C: drive full -> WSL VHDX moved to D:; (2) first main run killed by a 2 h background limit -> re-run detached with per-fit checkpoints; (3) one minimap2 failure (B022) re-run from the resumable script; (4) Manta fails on novoalign BAM (invalid chromosome id) -> Delly only on real data.
- Deviations from plan.md: N = 4,030 (2,015/class) not 6,000; HiFi only (no ONT); one replicate; no phasing; no INS/INV; depth sweep on 8 batches; no 60x short arm; own matching code instead of Truvari bench; own long-read simulator instead of PBSIM3.
- Outputs: reports/REPORT.md, DATA_CARD.md, model_cards/MODEL_CARD.md, reports/ablation/*, reports/data_audit/*, reports/figures/*, data/processed/dataset_v1_manifest.json.
- Learned: when annotation dominates, a breakpoint-precision gain only matters if it flips exon status; measure that fraction first.

### 2026-10-03 — P0/P1 — environment, sources verified, truth set built and gated
- Did: Host = Windows 11, 16 cores, 24 GB RAM (WSL given 17 GB via `.wslconfig`), RTX 4060 Laptop 8 GB. WSL2 distro moved to D: after C: ran out of space (Ubuntu VHDX re-registered in place as `UbuntuD`; Windows `C:` had 3 GB free). micromamba envs: `svenv` (samtools, bcftools, minimap2, bwa, ART, Sniffles2 2.8.1, cuteSV 2.1.4, truvari, whatshap, pysam), `manta` (py2.7, Manta 1.6.0), `delly` (Delly 2.6.0, own env: boost clash). GPU libs via pip: torch (cu128), xgboost 3.2, optuna, shap.
- Sources (verified by download, all GRCh38): ClinVar `variant_summary.txt.gz` (2026-10 snapshot), dbVar ClinVar-SV `nstd102.GRCh38`, gnomAD-SV v4.1 sites VCF, GENCODE v50, UCSC rmsk / genomicSuperDups / simpleRepeat / gap / centromeres, gnomAD v4.1 constraint metrics. Reference = Broad `Homo_sapiens_assembly38.fasta` restricted to chr1-22,X,Y,M (alt/decoy contigs removed; the NCBI/UCSC mirrors throttled to ~0.1 MB/s).
- Truth-set waterfall (see `data/processed/truthset_manifest.json`): ClinVar DEL/DUP GRCh38 autosomal 247,113 -> path/benign only 187,594 -> >=1 star 170,649 -> + gnomAD common PASS 32,004 -> size 50 bp-1 Mb 37,337 (4,841 path / 32,496 benign) -> window/N filter -> conflicts removed (RO>=0.5) -> dedup RO>=0.8: 3,964 path / 30,410 benign -> local-density cap (<=12 windows per overlapping +-100 kb, D10): 2,685 / 25,800 -> matched type x 0.25-decade size bins on common support: 2,359 / 2,359 -> subsampled to 2,015 / 2,015.
- Gate P1 (`reports/data_audit/truthset_gate.md`): schema OK, 0 cross-label conflicts, 0 near-duplicates, size KS D=0.028 p=0.42, trivial baselines (size, type, chromosome, position) all AUROC 0.48-0.53 on held-out chromosome groups => PASS.
- Known limits: pathogenic large CNVs (>100 kb) mostly dropped by matching because large benign CNVs are rare (67/32/12/7 DEL in the four largest bins); DUP only 132 per class; benign is 84% gnomAD (source confound tested in P6 source-swap).
- Next: simulate/align/call all 29 batches; P3 gate.
- Learned: network line speed, not CPU, can be the P0 bottleneck; check disk quota of the WSL VHDX before downloading 30 GB.

### 2026-10-03 — Planning — project re-scoped to simulation-first
- Did: Reviewed v1 context (`dna.md`). Researched simulators (ART, PBSIM3, NanoSim, Sim-it), SV benchmarking (Truvari, Manta, Sniffles2, cuteSV), prior art (StrVCTVRE, CSV-Filter, Samplot-ML, DeepSVFilter). Wrote `context.md`, `plan.md`, `LOG.md`, `LEARN.md`.
- Why: Real pathogenic SVs have no reads; only ~20% of benign short reads are derived from the long-read sample. A matched simulated design removes both problems.
- Numbers: none yet.
- Outputs: `context.md`, `plan.md`, `LOG.md`, `LEARN.md`.
- Gate: n/a.
- Next: P0 — verify data sources, builds, tool versions; set up WSL2 environment; run pilot.
- Learned: Check that a label is not defined by one of the features (AF here) before anything else.
