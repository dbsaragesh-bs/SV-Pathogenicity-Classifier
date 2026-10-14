"""P8: bootstrap tables, paired deltas, stratified deltas, flip table, figures.

python -m src.modelling.analysis --pred reports/ablation/main/oof_main.parquet \
        --metrics reports/ablation/main/metrics_main.csv --features data/processed/features_all.parquet \
        --out reports/ablation/main
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve

from src.modelling.gpu_metrics import ci, cluster_bootstrap, point_metrics
from src.utils.common import ensure

ARM_NAMES = {"A": "A: short-read only", "B": "B: short + long", "O": "Oracle: true coordinates",
             "OA": "Oracle (short evidence)"}


def pivot_scores(pred: pd.DataFrame, model: str, arms) -> pd.DataFrame:
    p = pred[pred.model == model]
    wide = p.pivot_table(index=["sv_id", "replicate"], columns="arm", values="score")
    meta = p.drop_duplicates(["sv_id", "replicate"]).set_index(["sv_id", "replicate"])[["label", "group_id", "chrom", "svtype", "fold"]]
    return wide[list(arms)].join(meta).dropna()


def best_model_per_arm(metrics: pd.DataFrame) -> dict:
    """Pre-specified selection rule: highest mean inner-CV AUPRC (never test AUPRC)."""
    m = metrics.groupby(["arm", "model"]).cv_auprc.mean().reset_index()
    return {a: g.sort_values("cv_auprc", ascending=False).model.iloc[0] for a, g in m.groupby("arm")}


def summary_table(pred, metrics, arms, B=2000):
    rows = []
    for model in sorted(pred.model.unique()):
        W = pivot_scores(pred, model, arms)
        bs = cluster_bootstrap(W.label.values, {a: W[a].values for a in arms}, W.group_id.values, B=B, seed=1)
        for a in arms:
            pm = point_metrics(W.label.values, W[a].values)
            rows.append(dict(model=model, arm=a, n=len(W), auprc=pm["auprc"], auprc_lo=ci(bs[a]["auprc"])[0],
                             auprc_hi=ci(bs[a]["auprc"])[1], auroc=pm["auroc"], auroc_lo=ci(bs[a]["auroc"])[0],
                             auroc_hi=ci(bs[a]["auroc"])[1],
                             cv_auprc=metrics[(metrics.arm == a) & (metrics.model == model)].cv_auprc.mean()))
    return pd.DataFrame(rows)


def paired_deltas(pred, arms_pairs, B=4000):
    rows = []
    for model in sorted(pred.model.unique()):
        arms = sorted({a for p in arms_pairs for a in p})
        W = pivot_scores(pred, model, arms)
        bs = cluster_bootstrap(W.label.values, {a: W[a].values for a in arms}, W.group_id.values, B=B, seed=2)
        for hi, lo in arms_pairs:
            for met in ("auprc", "auroc"):
                d = bs[hi][met] - bs[lo][met]
                pt = point_metrics(W.label.values, W[hi].values)[met] - point_metrics(W.label.values, W[lo].values)[met]
                rows.append(dict(model=model, contrast=f"{hi}-{lo}", metric=met, delta=pt, lo=ci(d)[0], hi=ci(d)[1],
                                 p_one_sided=float((d <= 0).mean()), p_two_sided=float(min(1, 2 * min((d <= 0).mean(), (d >= 0).mean())))))
    return pd.DataFrame(rows)


def stratified_deltas(pred, feats, model, hi="B", lo="A", strata: dict | None = None, B=1000, min_n=60):
    W = pivot_scores(pred, model, [hi, lo]).reset_index().merge(feats, on=["sv_id", "replicate"], how="left", suffixes=("", "_f"))
    rows = []
    for sname, series in (strata or {}).items():
        s = series(W)
        for level, g in W.groupby(s, observed=True):
            if len(g) < min_n or g.label.nunique() < 2 or g.label.sum() < 15 or (g.label == 0).sum() < 15:
                continue
            bs = cluster_bootstrap(g.label.values, {hi: g[hi].values, lo: g[lo].values}, g.group_id.values, B=B, seed=3)
            d = bs[hi]["auprc"] - bs[lo]["auprc"]
            pt = point_metrics(g.label.values, g[hi].values)["auprc"] - point_metrics(g.label.values, g[lo].values)["auprc"]
            rows.append(dict(stratum=sname, level=str(level), n=len(g), n_path=int(g.label.sum()), delta_auprc=pt,
                             lo=ci(d)[0], hi=ci(d)[1]))
    return pd.DataFrame(rows)


def default_strata():
    size_bins = [50, 300, 1000, 5000, 30000, 100000, 1e7]
    return {
        "svtype": lambda W: W.svtype,
        "size_band": lambda W: pd.cut(W.svlen, size_bins, labels=["50-300", "300-1k", "1k-5k", "5k-30k", "30k-100k", ">100k"]),
        "segdup_context": lambda W: np.where(W.a_ann_segdup_frac.fillna(0) > 0.1, "segdup", "non-segdup"),
        "repeat_at_breakpoint": lambda W: np.where(W.a_ann_bp_repeat_max.fillna(0) > 0.5, "repeat bp", "unique bp"),
        "gene_overlap": lambda W: np.where(W.o_ann_n_genes_pc > 0, "overlaps PC gene", "no PC gene"),
        "exon_overlap": lambda W: np.where(W.o_ann_n_exons_hit > 0, "hits exon", "no exon"),
    }


def flip_table(pred, feats, model, thr_col="thr"):
    """Rows where B and A disagree on the thresholded call; tabulate likely causes."""
    p = pred[pred.model == model]
    A = p[p.arm == "A"].set_index(["sv_id", "replicate"])
    Bm = p[p.arm == "B"].set_index(["sv_id", "replicate"])
    j = A[["label", "score", "thr"]].join(Bm[["score", "thr"]], lsuffix="_A", rsuffix="_B", how="inner").reset_index()
    j["pred_A"] = (j.score_A >= j.thr_A).astype(int)
    j["pred_B"] = (j.score_B >= j.thr_B).astype(int)
    j["B_right_A_wrong"] = (j.pred_B == j.label) & (j.pred_A != j.label)
    j["A_right_B_wrong"] = (j.pred_A == j.label) & (j.pred_B != j.label)
    f = j.merge(feats, on=["sv_id", "replicate"], how="left", suffixes=("", "_f"))
    # causes
    f["exon_status_changed"] = ((f.a_ann_n_exons_hit > 0) != (f.b_ann_n_exons_hit > 0))
    f["gene_count_changed"] = (f.a_ann_n_genes_pc != f.b_ann_n_genes_pc)
    f["breakpoint_in_exon_changed"] = (f.a_ann_bp_in_exon != f.b_ann_bp_in_exon)
    f["repeat_context"] = f.a_ann_segdup_frac.fillna(0) > 0.1
    f["lr_evidence_strong"] = f.l_support_total.fillna(0) >= 10
    rows = []
    for name, mask in (("B right, A wrong", f.B_right_A_wrong), ("A right, B wrong", f.A_right_B_wrong)):
        g = f[mask]
        rows.append(dict(group=name, n=len(g), exon_status_changed=float(g.exon_status_changed.mean()) if len(g) else np.nan,
                         gene_count_changed=float(g.gene_count_changed.mean()) if len(g) else np.nan,
                         bp_in_exon_changed=float(g.breakpoint_in_exon_changed.mean()) if len(g) else np.nan,
                         in_segdup=float(g.repeat_context.mean()) if len(g) else np.nan,
                         median_svlen=float(g.svlen.median()) if len(g) else np.nan))
    base = f
    rows.append(dict(group="all rows", n=len(base), exon_status_changed=float(base.exon_status_changed.mean()),
                     gene_count_changed=float(base.gene_count_changed.mean()),
                     bp_in_exon_changed=float(base.breakpoint_in_exon_changed.mean()),
                     in_segdup=float(base.repeat_context.mean()), median_svlen=float(base.svlen.median())))
    return pd.DataFrame(rows), f


def plot_pr(pred, best, out: Path, arms=("A", "B", "O")):
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for a in arms:
        p = pred[(pred.arm == a) & (pred.model == best[a])]
        if p.empty:
            continue
        pr, rc, _ = precision_recall_curve(p.label, p.score)
        ax.plot(rc, pr, label=f"{ARM_NAMES[a]} ({best[a]})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_ylim(0.4, 1.02)
    ax.legend(fontsize=8)
    ax.set_title("Precision-recall, held-out chromosomes")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_forest(df: pd.DataFrame, out: Path, title: str):
    fig, ax = plt.subplots(figsize=(6, 0.45 * len(df) + 1.2))
    y = np.arange(len(df))
    ax.errorbar(df.delta_auprc, y, xerr=[df.delta_auprc - df.lo, df.hi - df.delta_auprc], fmt="o", color="C0", capsize=2)
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.stratum}: {r.level} (n={r.n})" for r in df.itertuples()], fontsize=8)
    ax.set_xlabel("Δ AUPRC (B − A), 95% cluster-bootstrap CI")
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--arms", default="A,B,O")
    a = ap.parse_args()
    out = ensure(Path(a.out))
    arms = a.arms.split(",")
    pred = pd.read_parquet(a.pred)
    metrics = pd.read_csv(a.metrics)
    feats = pd.read_parquet(a.features)
    best = best_model_per_arm(metrics)
    print("pre-specified best model per arm:", best)
    summ = summary_table(pred, metrics, arms)
    summ.to_csv(out / "summary_table.csv", index=False)
    contrasts = [("B", "A")] + ([("O", "B"), ("O", "A")] if "O" in arms else [])
    dl = paired_deltas(pred, contrasts)
    dl.to_csv(out / "paired_deltas.csv", index=False)
    print(summ.round(3).to_string())
    print(dl.round(4).to_string())
    model = best.get("B", "xgb")
    st = stratified_deltas(pred, feats, model, strata=default_strata())
    st.to_csv(out / "stratified_deltas.csv", index=False)
    if len(st):
        plot_forest(st, out / "stratified_forest.png", f"Gain of long reads by stratum ({model})")
    ft, flips = flip_table(pred, feats, model)
    ft.to_csv(out / "flip_table.csv", index=False)
    plot_pr(pred, best, out / "pr_curves.png", [x for x in ("A", "B", "O") if x in arms])
    (out / "best_models.json").write_text(json.dumps(best))
    print(ft.round(3).to_string())


if __name__ == "__main__":
    main()
