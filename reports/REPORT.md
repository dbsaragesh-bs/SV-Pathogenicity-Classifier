# Does adding long reads improve SV pathogenicity classification? A matched, simulated ablation

*Course project, DNA Sequencing Technology. All numbers come from saved artefacts under `reports/`; run ids are the `tag` of each experiment (`main`, `ladder`, `logo`, `controls/*`, `depth/*`). Results are from a controlled simulation, not clinical data, and make no clinical claim.*

## 1. Question and design
We train two classifiers that label a structural variant (SV) pathogenic or benign. **Arm A** sees only information obtainable from short reads: short-read call geometry, caller fields, alignment evidence, and gene/repeat annotation recomputed from the *short-read called coordinates*. **Arm B** adds long-read call fields, long-read alignment evidence, microhomology at the long-read breakpoints, and annotation recomputed from the *long-read coordinates*. **Oracle** replaces coordinate-derived columns with the true planted coordinates (an upper bound for better breakpoints).

Real data cannot support a matched comparison (pathogenic SVs have no reads; benign short and long reads come from partly different material), so every SV is **planted into a simulated GRCh38 window**, and reads of both technologies are simulated from the same haplotypes by the same label-blind procedure. Technology is then the only varying factor. Cost: conclusions depend on the simulators and callers (Sections 7-8).

## 2. Data
* Labels: ClinVar (variant_summary, GRCh38) Pathogenic/Likely pathogenic vs Benign/Likely benign DEL/DUP with at least one review star, plus common (AF >= 1 %) PASS gnomAD-SV v4.1 DEL/DUP as benign. 50 bp to 1 Mb, autosomes.
* Cleaning waterfall (`data/processed/truthset_manifest.json`): 247,113 ClinVar DEL/DUP rows -> 4,841 pathogenic / 32,496 benign in range -> conflicts (RO >= 0.5) and near-duplicates (RO >= 0.8) removed -> local-density cap (<= 12 overlapping windows) -> matched on type x 0.25-decade size bins -> **4,030 SVs, 2,015 per class** (3,766 DEL, 264 DUP).
* Gate P1: 0 cross-label conflicts, 0 near-duplicates, size KS D = 0.028 (p = 0.42); size-only, type-only, chromosome-only and position-only models score AUROC 0.48-0.53 on held-out chromosomes (`reports/data_audit/truthset_gate.md`).
* Simulation: 29 pseudo-sample batches of up to 150 windows (SV +/- 50 kb, heterozygous SV, background SNV/indel 0.1 %, 27-33x per batch). Short reads: ART HS25 2x150, GC-thinned. Long reads: in-house HiFi simulator (median about 15 kb, 0.4 % error). Aligned to full GRCh38 (bwa-mem, minimap2); called with Manta + Delly and Sniffles2 + cuteSV. Gate P2: simulation parameters identical across classes (all KS p > 0.27, `sim_qc.md`).

## 3. Detection and call quality (P4, `matching.md`)
Detection is stable across matching stringency (short 0.917/0.916/0.912, HiFi 0.946/0.945/0.944 at RO 0.3/0.5/0.8). Median breakpoint error is 1 bp for both technologies (p99: short 50 bp, HiFi 9 bp). Short reads miss many benign SVs at segmental-duplication breakpoints (52 % detected vs 93 % at unique breakpoints). This is real detection biology, but it makes the both-detected cohort slightly selected, so all-SV and non-segdup analyses are also reported.

## 4. Leakage audit (P6, `leakage.md`)
Nuisance-only (coverage, batch, read counts), size/type/chromosome/position/number-of-callers-only, label-shuffled Arm A, and adversarial pipeline-column models all score AUROC 0.49-0.55 (chance). Source swap: a model trained with gnomAD-benign scores AUROC 0.89 on ClinVar-benign (0.96 the other way), so the signal is not source-specific. No duplicated feature rows. The pipeline is label-blind.

## 5. Main result
Cohort = detected by both technologies (n = 3,510); 5 outer folds of held-out chromosome groups, identical for all arms; hyper-parameters tuned in inner grouped CV only.

| Arm (XGBoost, GPU) | AUPRC [95 % cluster-bootstrap CI] | AUROC |
|---|---|---|
| A short only | 0.984 [0.976, 0.990] | 0.987 |
| B short + long | 0.985 [0.977, 0.990] | 0.986 |
| Oracle | 0.986 [0.980, 0.991] | 0.986 |

Paired differences (`paired_deltas.csv`): **dAUPRC(B-A) = +0.0004 [-0.0021, +0.0038]**, one-sided p = 0.42; Oracle-A = +0.0017 [-0.0005, +0.0048]. Logistic regression, random forest and a GPU MLP agree (dAUPRC -0.003 to +0.0004; no model has a positive CI excluding 0). **H1 is not supported: the gain is null, with an upper bound of about +0.004 AUPRC.**

Why: Arm A already reaches about 0.98 because pathogenicity here is dominated by *gene context*. SHAP (GPU TreeSHAP) ranks annotation first (exons hit, CDS overlap, fraction of CDS affected). Both technologies localise breakpoints to about 1 bp in the simulation, so annotation recomputed from long-read coordinates equals that from short-read coordinates: exon-overlap status changed in **0 %** of the 3,510 SVs between A and B coordinates, and the Oracle (perfect coordinates) is itself only +0.002 above A. There is no headroom for breakpoint precision to exploit (H2 not supported). Flip analysis: 19 SVs are right in B and wrong in A, 31 the reverse, with no exon-status change in either set.

## 6. Ablations and controls
* **Ladder** (`ablation_ladder.png`): adding long coordinates and annotation +0.001, long evidence +0.003, junction/microhomology (= B) +0.003, Oracle +0.001 over A; every CI includes 0.
* **Leave-one-family-out from B** (`logo.png`): removing any single family changes AUPRC by at most 0.005; the largest drop is short-read depth features (-0.005). Information is redundant across families.
* **Controls** (`controls_summary.csv`, dAUPRC B-A): without gene-dosage predictors +0.0009 [-0.002, +0.004]; size removed -0.0003; ClinVar-only labels +0.001 [-0.001, +0.004]; non-segdup breakpoints -0.003 [-0.007, +0.0004]; exon-hit SVs only (n = 1,988) **+0.0019 [+0.0004, +0.0035]**; **all planted SVs with detection flags (n = 4,030) +0.0036 [+0.0021, +0.0054]**. The only robustly positive effects appear when SVs missed by short reads stay in the data. The long-read benefit is mostly *detecting* SVs that short reads miss (94.5 % vs 91.6 % detection), not *describing* detected ones better. Treat the exon-hit result as exploratory (one of six controls, no multiplicity correction).
* **Stratified dAUPRC** (`stratified_deltas.csv`): every stratum CI includes 0, including small SVs (50-300 bp: +0.002 [-0.002, +0.010]) and segdup context (-0.017 [-0.040, +0.005], n = 263). H3 is not supported at this sample size; the minimum detectable effect for a stratum of about 250 SVs is roughly 0.02 AUPRC.
* **Depth sweep** (`depth/`; 811 SVs detected at every depth, 8 batches): long-read detection 97.8 % (30x), 97.4 % (20x), 94.8 % (10x), 76.4 % (5x), while classifier AUPRC stays flat (A 0.979-0.983, B 0.974-0.982). The detection benefit plateaus at about 10-20x (H4 holds for detection). Short-read depth 15x lowers short detection to 86.8 %. **More short reads vs long reads:** doubling short-read depth to 60x (independent second 30x draw merged) raises short detection only from 91.5 % to 92.9 %, whereas adding 30x HiFi reaches 97.8 % on the same 8 batches (1,200 SVs); classifier AUPRC is unchanged either way (A at 60x 0.980 vs A at 30x 0.979, B 0.981; n = 811). Per added coverage unit, long reads recover far more of the SVs short reads miss (segdup breakpoints), and extra short depth does not.

## 7. Realism check on real HG002 data (P9, `p9_realdata.md`)
104 GIAB HG002 v5.0q GRCh38 deletions (52-6,640 bp) were streamed from real Illumina 2x250 and Revio HiFi alignments, subsampled to about 30x, and passed through the same calling and feature code (Manta crashed on the real novoalign BAM, so Delly only for short reads). Detection: short (Delly only) 0.106 real vs 0.850 simulated; HiFi 0.644 real vs 0.958 simulated. Adversarial validation separates real from simulated evidence features with AUROC 0.97 (short, n_real = 11) and 0.98 (long, n_real = 67). The most shifted features are the MAPQ scale, discordant-pair counts, long-read breakpoint dispersion (real about 5-10 bp vs 0 in simulation) and clip dispersion. Simulation-trained models call 0 of 11 real benign deletions pathogenic. **The simulation is clearly cleaner than real data.** Real breakpoint noise would plausibly give long reads a larger relative advantage than measured here, so the null result describes the simulated regime, and its bias favours finding a long-read gain rather than hiding one.

## 8. Limitations
Simulated errors only (uniform long-read errors, no homopolymer bias, no chimeras, constant quality); one tandem-duplication model; heterozygous SVs only; DEL-dominant (264 DUP); large benign CNVs are rare, so large pathogenic CNVs are under-represented after matching; the benign set is 84 % gnomAD while pathogenic is 100 % ClinVar (mitigated by size matching, tested by source-swap and ClinVar-only controls); ClinVar array-based coordinates are imprecise (planted as given); the gene-context shortcut makes the benchmark easy (AUPRC about 0.98), which lowers sensitivity to any gain; no real pathogenic reads exist, so sensitivity on real pathogenic SVs is untested; INS/INV, phasing and ONT were not run; single replicate; the depth and 60x analyses use only 8 batches (1,200 SVs); the P9 sample is small and short-read-limited (11 detections).

## 9. Conclusion
In a matched, leakage-audited simulation, long-read information does **not** measurably improve SV pathogenicity classification once an SV is detected (dAUPRC +0.0004, CI -0.002 to +0.004). The classifier is limited by gene context that both technologies resolve equally. The measurable long-read benefit is detection: HiFi recovers 94.5 % vs 91.6 % of planted SVs (mostly benign SVs at segmental duplications), and keeping undetected SVs in the analysis gives dAUPRC +0.004 [+0.002, +0.005]; detection saturates around 10-20x. The real-data check shows the simulation is optimistic about breakpoint precision, so re-testing with real long-read call sets is the key next step.

## 10. Reproduce
Scripts in `scripts/` (WSL), `workflow/Snakefile`, configs in `configs/`, seeds derived from SV ids. Unit tests: `pytest tests`. Compute: 16 CPU cores for simulation, alignment, calling and feature extraction; RTX 4060 GPU for XGBoost, MLP, cluster-bootstrap statistics and TreeSHAP.
