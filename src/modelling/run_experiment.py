"""P7: train every arm with an identical protocol on identical outer folds.

python -m src.modelling.run_experiment --features data/processed/features_all.parquet --cohort both \
       --arms A,B,O --models xgb,lr,rf,mlp --trials 20 --out reports/ablation/main
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import (average_precision_score, brier_score_loss, f1_score, matthews_corrcoef,
                             precision_score, recall_score, roc_auc_score)

from src.features.registry import arm_columns
from src.modelling.models import fit_tuned
from src.modelling.splits import outer_folds
from src.utils.common import REPORTS, ensure


def cohort_mask(df, cohort):
    if cohort == "both":
        return (df.detected_short == 1) & (df.detected_hifi == 1)
    if cohort == "all":
        return np.ones(len(df), bool)
    if cohort == "short":
        return df.detected_short == 1
    raise ValueError(cohort)


def best_threshold(y, s):
    ts = np.unique(np.quantile(s, np.linspace(0.02, 0.98, 97)))
    f = [f1_score(y, s >= t) for t in ts]
    return float(ts[int(np.argmax(f))])


def recall_at_precision(y, s, p=0.9):
    o = np.argsort(-s)
    yy = y[o]
    tp = np.cumsum(yy)
    prec = tp / np.arange(1, len(yy) + 1)
    ok = np.where(prec >= p)[0]
    return float(tp[ok[-1]] / yy.sum()) if len(ok) else 0.0


def run(df: pd.DataFrame, arms, models, k_outer=5, trials=20, seed=0, out: Path | None = None, tag="", arm_cols: dict | None = None):
    y_all = df.label.to_numpy().astype(int)
    groups = df.group_id.to_numpy()
    recs, metrics = [], []
    part_dir = (out / f"parts_{tag}") if out is not None else None
    if part_dir is not None:
        ensure(part_dir)
    TR = {"xgb": trials, "lr": max(8, trials * 3 // 4), "rf": max(6, trials // 3), "mlp": max(5, trials // 4)}
    for fold, tr, te in outer_folds(df, k_outer, seed):
        for arm in arms:
            cols = arm_cols[arm] if arm_cols and arm in arm_cols else arm_columns(arm, df.columns)
            cols = [c for c in cols if df[c].nunique(dropna=True) > 1]
            X = df[cols].to_numpy(dtype=float)
            X = np.where(np.isfinite(X), X, np.nan)
            for mname in models:
                t0 = time.time()
                pf = part_dir / f"{fold}_{arm}_{mname}.pkl" if part_dir is not None else None
                if pf is not None and pf.exists():
                    import pickle
                    m, rdf = pickle.load(open(pf, "rb"))
                    metrics.append(m)
                    recs.append(rdf)
                    continue
                model, best, oof_in, cv_ap = fit_tuned(mname, X[tr], y_all[tr], groups[tr], TR.get(mname, trials), seed)
                thr = best_threshold(y_all[tr], oof_in)
                iso = IsotonicRegression(out_of_bounds="clip").fit(oof_in, y_all[tr])
                s = model.predict_proba(X[te])[:, 1]
                sc = iso.predict(s)
                yt = y_all[te]
                pred = (s >= thr).astype(int)
                m = dict(tag=tag, arm=arm, model=mname, fold=fold, n_train=len(tr), n_test=len(te), n_feat=len(cols),
                         pos_test=int(yt.sum()), cv_auprc=cv_ap,
                         auprc=average_precision_score(yt, s), auroc=roc_auc_score(yt, s) if 0 < yt.sum() < len(yt) else np.nan,
                         mcc=matthews_corrcoef(yt, pred), f1=f1_score(yt, pred), precision=precision_score(yt, pred, zero_division=0),
                         recall=recall_score(yt, pred), recall_at_p90=recall_at_precision(yt, s, 0.9),
                         brier=brier_score_loss(yt, np.clip(sc, 0, 1)), threshold=thr, seconds=time.time() - t0,
                         best_params=json.dumps(best))
                metrics.append(m)
                recs.append(pd.DataFrame({"sv_id": df.sv_id.values[te], "replicate": df.replicate.values[te], "label": yt,
                                          "group_id": groups[te], "chrom": df.chrom.values[te], "svtype": df.svtype.values[te],
                                          "arm": arm, "model": mname, "fold": fold, "score": s, "score_cal": sc,
                                          "thr": thr}))
                if pf is not None:
                    import pickle
                    pickle.dump((m, recs[-1]), open(pf, "wb"))
                print(f"[{tag}] fold {fold} arm {arm} {mname}: AUPRC {m['auprc']:.3f} AUROC {m['auroc']:.3f} "
                      f"(cv {cv_ap:.3f}, {m['seconds']:.0f}s)", flush=True)
    pred = pd.concat(recs, ignore_index=True)
    met = pd.DataFrame(metrics)
    if out is not None:
        ensure(out)
        pred.to_parquet(out / f"oof_{tag}.parquet", index=False)
        met.to_csv(out / f"metrics_{tag}.csv", index=False)
    return pred, met


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    ap.add_argument("--cohort", default="both")
    ap.add_argument("--arms", default="A,B,O")
    ap.add_argument("--models", default="xgb,lr,rf,mlp")
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=str(REPORTS / "ablation" / "main"))
    ap.add_argument("--tag", default="main")
    a = ap.parse_args()
    df = pd.read_parquet(a.features)
    df = df[cohort_mask(df, a.cohort)].reset_index(drop=True)
    print(f"cohort {a.cohort}: {len(df)} rows, {df.label.sum()} pathogenic")
    run(df, a.arms.split(","), a.models.split(","), a.k, a.trials, a.seed, Path(a.out), a.tag)
