# MODEL_CARD — SV pathogenicity classifiers (Arm A, Arm B, Oracle)

**Purpose:** research instrument to measure the value of long-read information; **not a clinical classifier.** Trained and evaluated only on simulated reads of planted ClinVar/gnomAD SVs.

* **Models:** XGBoost (GPU, `device=cuda`), logistic regression (elastic net), random forest, MLP (PyTorch, GPU, 2-model ensemble). Same Optuna budget per arm (inner 3-fold grouped, stratified CV; AUPRC objective); preprocessing (median imputation + missing indicators, quantile transform) fitted inside training folds only.
* **Arms:** A = short-read features; B = A + long-read features; Oracle = truth-coordinate annotation/junction + shared evidence. Feature membership: `src/features/registry.py`.
* **Evaluation:** 5 outer folds of held-out chromosome groups (identical for every arm); groups = 50 %-overlap clusters; threshold (max F1) and isotonic calibration from inner OOF predictions; metrics AUPRC (primary), AUROC, MCC, F1, recall@precision 0.9, Brier; paired cluster bootstrap (4,000 resamples, GPU).
* **Performance (cohort detected by both, n = 3,510):** XGBoost AUPRC A 0.984, B 0.985, Oracle 0.986; B-A = +0.0004 [-0.0021, +0.0038].
* **Dominant signals (SHAP):** exons hit, CDS overlap, fraction of CDS affected, constraint (LOEUF/pLI), repeat content.
* **Limitations:** performance is high because gene overlap separates ClinVar pathogenic from gnomAD benign SVs; calibration reflects the 50/50 matched prevalence (not clinical prevalence); real HG002 benign deletions are all scored benign (0/11 pathogenic) but the simulation-to-real feature gap is large (adversarial AUROC ~0.97).
* **Ethics/use:** no patient data; public databases only; do not use for diagnosis.
