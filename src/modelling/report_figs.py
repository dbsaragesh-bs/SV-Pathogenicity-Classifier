"""Figures/tables for the report: SHAP by family, ladder, LOGO, controls, depth sweep."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
import xgboost as xgb

from src.features.registry import arm_columns, feature_groups
from src.modelling.gpu_metrics import ci, cluster_bootstrap, point_metrics
from src.utils.common import PROCESSED, REPORTS, ensure

FIG = ensure(REPORTS / "figures")


def delta_vs(pred, ref_arm, arms, B=3000):
    p = pred.pivot_table(index=["sv_id", "replicate"], columns="arm", values="score")
    meta = pred.drop_duplicates(["sv_id", "replicate"]).set_index(["sv_id", "replicate"])[["label", "group_id"]]
    W = p.join(meta).dropna()
    bs = cluster_bootstrap(W.label.values, {a: W[a].values for a in arms}, W.group_id.values, B=B, seed=5)
    rows = []
    for a in arms:
        d = bs[a]["auprc"] - bs[ref_arm]["auprc"]
        rows.append(dict(arm=a, auprc=point_metrics(W.label.values, W[a].values)["auprc"],
                         delta=point_metrics(W.label.values, W[a].values)["auprc"] - point_metrics(W.label.values, W[ref_arm].values)["auprc"],
                         lo=ci(d)[0], hi=ci(d)[1]))
    return pd.DataFrame(rows)


def bar_delta(df, title, path, ref):
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(df) + 1.3))
    y = np.arange(len(df))
    ax.errorbar(df.delta, y, xerr=[df.delta - df.lo, df.hi - df.delta], fmt="o", capsize=2)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(df.arm, fontsize=8)
    ax.set_xlabel(f"Δ AUPRC vs {ref} (95% cluster-bootstrap CI)")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def shap_by_family():
    df = pd.read_parquet(PROCESSED / "features_all.parquet")
    df = df[(df.detected_short == 1) & (df.detected_hifi == 1)].reset_index(drop=True)
    cols = [c for c in arm_columns("B", df.columns) if df[c].nunique() > 1]
    X = df[cols].to_numpy(float)
    m = xgb.XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                          tree_method="hist", device="cuda", verbosity=0, random_state=0).fit(X, df.label.values)
    booster = m.get_booster()
    contrib = booster.predict(xgb.DMatrix(X, feature_names=cols), pred_contribs=True)[:, :-1]  # GPU TreeSHAP
    imp = pd.Series(np.abs(contrib).mean(0), index=cols).sort_values(ascending=False)
    fam = feature_groups(cols)
    inv = {c: f for f, cs in fam.items() for c in cs}
    fam_imp = imp.groupby(imp.index.map(lambda c: inv.get(c, "other"))).sum().sort_values(ascending=True)
    imp.head(25).to_csv(REPORTS / "ablation" / "main" / "shap_top_features.csv")
    fam_imp.to_csv(REPORTS / "ablation" / "main" / "shap_by_family.csv")
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    imp.head(15)[::-1].plot.barh(ax=ax[0])
    ax[0].set_title("Mean |SHAP|, top 15 features (Arm B, XGBoost, GPU TreeSHAP)", fontsize=8)
    fam_imp.plot.barh(ax=ax[1], color="C1")
    ax[1].set_title("Mean |SHAP| summed by feature family", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "shap_families.png", dpi=160)
    plt.close(fig)
    return imp, fam_imp


def main():
    out = {}
    lad = REPORTS / "ablation" / "ladder" / "oof_ladder.parquet"
    if lad.exists():
        pr = pd.read_parquet(lad)
        arms = list(pr.arm.unique())
        d = delta_vs(pr, "R0_A", arms)
        d.to_csv(REPORTS / "ablation" / "ladder" / "ladder_deltas.csv", index=False)
        bar_delta(d, "Ablation ladder (XGBoost)", FIG / "ablation_ladder.png", "R0_A")
        print(d.round(4).to_string())
    lg = REPORTS / "ablation" / "logo" / "oof_logo.parquet"
    if lg.exists():
        pr = pd.read_parquet(lg)
        arms = list(pr.arm.unique())
        d = delta_vs(pr, "B_full", arms)
        d.to_csv(REPORTS / "ablation" / "logo" / "logo_deltas.csv", index=False)
        bar_delta(d, "Leave-one-feature-family-out from Arm B", FIG / "logo.png", "B_full")
        print(d.round(4).to_string())
    ctl = REPORTS / "ablation" / "controls"
    rows = []
    for f in sorted(ctl.glob("oof_*.parquet")):
        tag = f.stem[4:]
        pr = pd.read_parquet(f)
        arms = [a for a in ("A", "B", "O") if a in set(pr.arm)]
        d = delta_vs(pr, "A", arms)
        d["control"] = tag
        d["n"] = pr.drop_duplicates(["sv_id", "replicate"]).shape[0]
        d["n_path"] = int(pr.drop_duplicates(["sv_id", "replicate"]).label.sum())
        rows.append(d)
    if rows:
        c = pd.concat(rows)
        c.to_csv(ctl / "controls_summary.csv", index=False)
        print(c.round(4).to_string())
    imp, fam = shap_by_family()
    print(fam)
    dp = REPORTS / "ablation" / "depth" / "depth_summary.csv"
    if dp.exists():
        t = pd.read_csv(dp)
        order = ["L5", "L10", "L20", "L30_S30"]
        t = t.set_index("config").loc[order].reset_index()
        x = [5, 10, 20, 30]
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
        ax[0].plot(x, t.A_auprc, "o-", label="Arm A (short only)")
        ax[0].plot(x, t.B_auprc, "s-", label="Arm B (short + long)")
        ax[0].set_xlabel("Long-read depth (x)")
        ax[0].set_ylabel("AUPRC (common cohort)")
        ax[0].legend(fontsize=8)
        ax[1].plot(x, t.det_long, "o-", label="all")
        ax[1].plot(x, t.det_long_path, "s--", label="pathogenic")
        ax[1].plot(x, t.det_long_benign, "^--", label="benign")
        ax[1].set_xlabel("Long-read depth (x)")
        ax[1].set_ylabel("Long-read detection rate")
        ax[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG / "depth_sweep.png", dpi=160)
        plt.close(fig)


if __name__ == "__main__":
    main()
