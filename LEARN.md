# LEARN.md — How to Understand This Project

Goal: after finishing the project you can explain every step, defend every choice, and re-run everything. This file is a skeleton now. Fill the "Your notes" and "Results" parts as phases finish (instruction repeated in `plan.md` phase gates).

## 1. Read in this order
1. `context.md` — what the project is and why it was designed this way (10 min).
2. `plan.md` — how it is built, phase by phase (30 min).
3. `LOG.md` — what actually happened, with numbers and decisions.
4. `reports/` — data audit, ablation results, figures (created during the project).
5. `src/` and `workflow/` — code, in pipeline order (simulate, align, call, match, features, model).

## 2. Concept primer (plain language)

**Structural variant (SV).** A change of 50 bp or more in the genome: deletion (piece missing), duplication (piece repeated), insertion (new piece), inversion (piece flipped). Some SVs are harmless, some cause disease by removing or doubling important genes.

**Short reads vs long reads.** A short read is about 150 letters of DNA. A long read is 10,000+ letters. A short read cannot span a big SV, so software infers it from clues: reads whose mates are too far apart (discordant pairs), reads cut in half at the breakpoint (split reads), and drops in coverage (depth). A long read can cross the whole breakpoint, so the exact position and the inserted letters are visible directly.

**Why simulate?** We need the same SV observed by both technologies, with a known correct answer. Real pathogenic SVs have no read data, so we build a synthetic genome with the SV inside and generate both kinds of reads from it. Because we know exactly what we planted, we can check every step.

**Pipeline in one line.** Plant SV into reference, simulate reads, align reads, call SVs, match calls to what was planted, compute features, train classifiers, compare.

**Alignment.** Placing each read at its best-matching location on the reference genome (minimap2 for long reads, bwa-mem2 for short).

**SV calling.** Software that looks at aligned reads and reports SVs (Manta, DELLY for short; Sniffles2, cuteSV for long).

**Truth matching.** Deciding which called SV corresponds to which planted SV (Truvari). Needed to attach the label to a call.

**Features.** Numbers describing an SV: type, size, how many genes/exons it hits, how repetitive the area is, how many reads support it, how precisely the breakpoints are known, etc.

**Arm A / Arm B / Oracle.** A = short-read information only. B = A + long-read information. Oracle = perfect breakpoints. Comparing them shows how much long reads add and the best that could be gained.

**AUPRC.** Area under the precision-recall curve. Better than accuracy when one class is rarer.

**Data leakage.** Information in the training data that would not be available at prediction time or that directly encodes the answer. Most "too good" results are leakage. Examples in this project: AF defines the benign label; the same SV appearing in train and test; simulation settings that differ between classes.

**Class imbalance.** Pathogenic SVs are rarer than benign. We handle it by matched sampling, class weights, replicate-level augmentation, and by judging with AUPRC.

## 3. Why each design choice exists (viva-ready)

| Choice | One-line reason |
|---|---|
| Simulation for both classes | Only way to get matched short and long reads for pathogenic SVs |
| Same pipeline, label-blind | Prevents the model learning "simulated vs real" or batch artifacts |
| AF not a feature | The benign label is defined by AF, so it would leak |
| Annotation from called coordinates | Shows the real effect of breakpoint precision |
| Oracle arm | Gives the ceiling of the long-read benefit |
| Depth sweep | Answers "is the extra sequencing worth it" |
| Cluster/chromosome-held-out splits | Stops near-duplicate SVs leaking from train to test |
| Replicates grouped | Same SV with different read noise must stay in one fold |
| Shortcut baselines (size-only etc.) | A good model must beat trivial rules |
| Label-shuffle and nuisance tests | Prove the pipeline has no hidden leak |
| Real HG002 check | Shows simulated feature distributions are realistic |

## 4. How to re-run (fill in when the pipeline exists)
```
# 1. create environment      -> envs/environment.yml
# 2. fetch + verify data     -> snakemake -s workflow/Snakefile fetch
# 3. build truth sets        -> snakemake truth
# 4. simulate + call         -> snakemake sim_call
# 5. build features          -> snakemake features
# 6. train + evaluate        -> snakemake model
# 7. produce report figures  -> snakemake report
```
Exact commands, runtimes, and disk use: to be filled at the end of P10.

## 5. Results (fill at the end)
- Final dataset size per class/type/size band:
- Arm A vs Arm B vs Oracle (AUPRC, CI):
- Depth sweep conclusion:
- Which long-read features mattered:
- Where long reads helped / did not help:
- Biggest limitations:

## 6. Your notes
(Write what confused you and how you resolved it, per phase.)

## 7. Likely viva questions
- Why simulate instead of using real data? What is lost?
- How do you know the simulated data is realistic?
- Why is AF not a feature?
- How did you prevent the model from learning the simulator?
- What exactly is "Arm B" and why is it a fair comparison?
- What does the Oracle arm show?
- How did you handle imbalance without SMOTE?
- Why held-out chromosomes plus clusters?
- What if long reads gave no gain? What would that mean?
- Is this clinically usable? (No. Why not?)
- What are the limits of treating "common = benign"?


---
## 8. Results summary and how to re-run (filled at project end)

**Headline:** B - A = +0.0004 AUPRC [-0.0021, +0.0038]; Oracle - A = +0.0017. Long reads help by *detecting* more SVs (94.5 % vs 91.6 %; all-SV analysis +0.0036 [+0.0021, +0.0054]), not by better description of detected SVs, because annotation (exons/CDS hit) already gives AUPRC ~0.98 and both technologies place breakpoints within ~1 bp in the simulation.

**Re-run (inside WSL, `source scripts/wsl_env.sh`):**
```bash
python -m src.truthset.build_truthset && python -m src.truthset.make_batches      # ~3 min
bash scripts/run_stage1.sh 0      # simulate + HiFi align/call, ~75 min for 29 batches
bash scripts/run_stage2.sh 0 base # short align/call + features, ~2 h
python -m src.features.assemble base features_all.parquet
python -m src.validation.sim_qc && python -m src.validation.leakage && python -m src.validation.matching_audit
bash scripts/run_main.sh          # GPU models, ~2.5 h
python -m src.modelling.analysis --pred reports/ablation/main/oof_main.parquet --metrics reports/ablation/main/metrics_main.csv --features data/processed/features_all.parquet --out reports/ablation/main --arms A,B,O,OA
bash scripts/run_abl.sh; bash scripts/run_sweep.sh; python -m src.modelling.depth_analysis; python -m src.modelling.report_figs
bash scripts/run_p9a.sh; bash scripts/run_p9b.sh; python -m src.validation.p9_compare
```

**Concepts worth remembering:** (1) a matched comparison needs the same SV on both technologies, so we simulate; (2) leakage tests are as important as models; (3) cluster-bootstrap resamples whole overlap groups, paired across arms; (4) AUPRC near 1 can mean an easy benchmark, not a good model; (5) a null result needs a confidence interval and a minimum detectable effect.

## 9. Viva answers (short)
* *Why simulate?* Pathogenic SVs have no reads; benign short/long reads come from different material.
* *Why is accuracy so high?* Gene overlap separates ClinVar pathogenic from gnomAD benign SVs; checked by leakage tests and ClinVar-only controls.
* *Why no gain from long reads?* Short-read breakpoints are already ~1 bp in the simulation; exon status never changed between arms (0/3,510).
* *What would change the answer?* Real breakpoint noise in repeats, harder (non-gene) pathogenic classes, INS/INV, lower short-read quality.
* *Where is the GPU used?* XGBoost (CUDA), MLP, bootstrap AUPRC/AUROC, TreeSHAP; alignment and calling are CPU-only tools.
