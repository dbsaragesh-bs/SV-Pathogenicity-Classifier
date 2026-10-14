"""P6: try hard to break the dataset.

python -m src.validation.leakage --features data/processed/features_all.parquet
Writes reports/data_audit/leakage.md and leakage.json
Every test uses the same outer folds (held-out chromosome groups) as the real experiment.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import average_precision_score, roc_auc_score

from src.features.registry import arm_columns
from src.modelling.splits import outer_folds
from src.utils.common import PROCESSED, REPORTS, ensure

XGB = dict(tree_method="hist", device="cuda", n_estimators=200, max_depth=3, learning_rate=0.05, subsample=0.8,
           colsample_bytree=0.8, verbosity=0, n_jobs=4)


def cv_score(df, cols, y=None, k=5, seed=0, model_kw=None):
    """Out-of-fold AUPRC/AUROC over held-out chromosome groups for a given column set."""
    y = df.label.values if y is None else y
    X = df[cols].to_numpy(float)
    oof = np.full(len(df), np.nan)
    for f, tr, te in outer_folds(df, k, seed):
        m = xgb.XGBClassifier(random_state=seed, **(model_kw or XGB)).fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    ok = ~np.isnan(oof)
    return {"auprc": float(average_precision_score(y[ok], oof[ok])), "auroc": float(roc_auc_score(y[ok], oof[ok])),
            "prevalence": float(y[ok].mean()), "n": int(ok.sum())}


def make_nuisance(df: pd.DataFrame) -> pd.DataFrame:
    nz = pd.read_parquet(PROCESSED / "nuisance.parquet")
    cols = ["sv_id", "replicate", "cov_short_target", "cov_long_target", "n_small", "short_cov_ach", "hifi_cov_ach",
            "short_pairs_kept", "hifi_reads", "win_len"]
    cols = [c for c in cols if c in nz.columns]
    out = df[["sv_id", "replicate", "batch", "chrom", "start"]].merge(nz[cols], on=["sv_id", "replicate"], how="left")
    out["batch_idx"] = out.batch.str[1:].astype(int)
    out["chrom_idx"] = out.chrom.str[3:].astype(int)
    return out


def run(df: pd.DataFrame) -> dict:
    res = {}
    df = df.reset_index(drop=True)
    nz = make_nuisance(df)
    nz_cols = [c for c in nz.columns if c not in ("sv_id", "replicate", "batch", "chrom")]
    d2 = pd.concat([df[["label", "chrom", "group_id", "svtype", "svlen"]], nz[nz_cols]], axis=1)
    d2["svlen_log"] = np.log10(d2.svlen)
    d2["svtype_dup"] = (d2.svtype == "DUP").astype(int)
    d2["start_pos"] = nz["start"]
    # 1. nuisance-only
    res["nuisance_only"] = cv_score(d2, [c for c in nz_cols if c not in ("start", "win_len")])
    # 2. shortcut baselines
    res["size_only"] = cv_score(d2, ["svlen_log"])
    res["type_only"] = cv_score(d2, ["svtype_dup"])
    res["chrom_only"] = cv_score(d2, ["chrom_idx"])
    res["position_only"] = cv_score(d2, ["start_pos"])
    ncall = [c for c in df.columns if c in ("a_n_callers", "b_n_callers")]
    if ncall:
        res["n_callers_only"] = cv_score(df.assign(), ncall)
    # 3. label shuffle within (type x size-bin) strata, full Arm A
    colsA = [c for c in arm_columns("A", df.columns) if df[c].nunique() > 1]
    rng = np.random.default_rng(1)
    ysh = df.label.values.copy()
    strata = df.svtype + "_" + pd.qcut(df.svlen, 10, duplicates="drop").astype(str)
    for s in strata.unique():
        idx = np.where(strata.values == s)[0]
        ysh[idx] = rng.permutation(ysh[idx])
    res["label_shuffle_armA"] = cv_score(df, colsA, y=ysh)
    res["real_armA_reference"] = cv_score(df, colsA)
    # 4. adversarial validation on pipeline-only (non-biological) columns
    pipe = [c for c in df.columns if c.endswith(("_reads", "reads_at_bp", "readlen_mean", "nm_rate"))
            and c.startswith(("s_", "l_"))]
    if pipe:
        res["adversarial_pipeline_cols"] = cv_score(df, pipe)
    # 6. source swap: benign basis AF vs ClinVar
    if "benign_basis" in df:
        pos = df.label == 1
        for tr_basis, te_basis in (("af", "clinvar"), ("clinvar", "af")):
            tr = pos & (df.chrom.map(lambda c: True)) | (df.benign_basis == tr_basis)
            te = pos | (df.benign_basis == te_basis)
            trm = (df.label == 1) | (df.benign_basis == tr_basis)
            tem = (df.label == 1) | (df.benign_basis == te_basis)
            if (df[trm].label == 0).sum() < 30 or (df[tem].label == 0).sum() < 30:
                continue
            X = df[colsA].to_numpy(float)
            # train on trm rows from chroms outside test fold, test on tem rows inside fold
            sc_te, y_te = [], []
            for f, tr_i, te_i in outer_folds(df, 5, 0):
                tr_i = np.intersect1d(tr_i, np.where(trm)[0])
                te_i = np.intersect1d(te_i, np.where(tem)[0])
                m = xgb.XGBClassifier(random_state=0, **XGB).fit(X[tr_i], df.label.values[tr_i])
                sc_te.append(m.predict_proba(X[te_i])[:, 1])
                y_te.append(df.label.values[te_i])
            s, y = np.concatenate(sc_te), np.concatenate(y_te)
            res[f"source_swap_train_{tr_basis}_test_{te_basis}"] = {
                "auprc": float(average_precision_score(y, s)), "auroc": float(roc_auc_score(y, s)),
                "prevalence": float(y.mean()), "n": int(len(y))}
    # 5. duplicates: identical feature vectors across rows with different group_id
    Xs = pd.util.hash_pandas_object(df[colsA].astype(float).round(6), index=False)
    dup = Xs.duplicated(keep=False)
    res["identical_feature_rows"] = int(dup.sum())
    res["identical_feature_rows_cross_group"] = int(df[dup].groupby(Xs[dup]).group_id.nunique().gt(1).sum())
    # 7. importance sanity
    X = df[colsA].to_numpy(float)
    m = xgb.XGBClassifier(random_state=0, **XGB).fit(X, df.label.values)
    imp = pd.Series(m.feature_importances_, index=colsA).sort_values(ascending=False)
    res["top_features_armA"] = {k: float(v) for k, v in imp.head(15).items()}
    return res


def write_report(res: dict, path):
    L = ["# Leakage and shortcut audit (P6)", "",
         "All scores are out-of-fold over held-out chromosome groups. Chance AUROC = 0.5; chance AUPRC = prevalence.", "",
         "| test | AUPRC | AUROC | prevalence | verdict |", "|---|---|---|---|---|"]
    for k, v in res.items():
        if isinstance(v, dict) and "auprc" in v:
            chance = abs(v["auroc"] - 0.5) < 0.08
            verdict = "PASS (≈chance)" if chance else ("reference" if "reference" in k else "REVIEW")
            L.append(f"| {k} | {v['auprc']:.3f} | {v['auroc']:.3f} | {v['prevalence']:.3f} | {verdict} |")
    L += ["", f"Identical feature rows: {res.get('identical_feature_rows')}; "
              f"across different groups: {res.get('identical_feature_rows_cross_group')}", "",
          "## Top Arm A features (sanity)", ""]
    for k, v in res["top_features_armA"].items():
        L.append(f"- {k}: {v:.3f}")
    path.write_text("\n".join(L))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default=str(PROCESSED / "features_all.parquet"))
    ap.add_argument("--cohort", default="both")
    a = ap.parse_args()
    df = pd.read_parquet(a.features)
    if a.cohort == "both":
        df = df[(df.detected_short == 1) & (df.detected_hifi == 1)]
    res = run(df)
    out = ensure(REPORTS / "data_audit")
    (out / "leakage.json").write_text(json.dumps(res, indent=1))
    write_report(res, out / "leakage.md")
    print(open(out / "leakage.md").read())
