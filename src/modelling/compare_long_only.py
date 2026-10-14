"""Compare A (short), B (short+long), L (long only), Oracle on identical rows and folds.

python -m src.modelling.compare_long_only
Reads reports/ablation/main/oof_main.parquet and reports/ablation/long_only/oof_*.parquet.
Writes comparison.csv, paired_deltas.csv, long_parts.csv, comparison.png under reports/ablation/long_only/.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, f1_score, matthews_corrcoef, precision_score, recall_score,
                             roc_auc_score)

from src.utils.common import REPORTS

OUT = REPORTS / "ablation" / "long_only"
NAMES = {"A": "A  short only", "B": "B  short+long", "L": "L  long only", "O": "Oracle (true coords)"}
BOOT = 1000


def load():
    m = pd.read_parquet(REPORTS / "ablation" / "main" / "oof_main.parquet")
    m = m[m.arm.isin(["A", "B", "O"])]
    l = pd.read_parquet(OUT / "oof_long_only.parquet")
    d = pd.concat([m, l], ignore_index=True)
    d["key"] = d.sv_id + "_" + d.replicate.astype(str)
    return d


def wacc(y, pred, W):
    return ((y == pred) * W).sum(1) / W.sum(1)


def main():
    d = load()
    keys = sorted(d.key.unique())
    ref = d[(d.arm == "A") & (d.model == "xgb")].set_index("key").loc[keys]
    y = ref.label.to_numpy().astype(int)
    grp = pd.factorize(ref.group_id)[0]
    G = grp.max() + 1
    rng = np.random.default_rng(7)
    cnt = rng.multinomial(G, np.ones(G) / G, size=BOOT)  # BOOT x G cluster multiplicities
    W = cnt[:, grp].astype(float)  # BOOT x n row weights
    print("rows", len(y), "pos", int(y.sum()), "groups", G, flush=True)

    rows, store = [], {}
    for (arm, model), g in d.groupby(["arm", "model"]):
        g = g.set_index("key").loc[keys]
        s = g.score.to_numpy()
        thr = g.thr.to_numpy()
        p_t = (s >= thr).astype(int)
        p_5 = (s >= 0.5).astype(int)
        ap_b = np.array([average_precision_score(y, s, sample_weight=w) for w in W])
        auc_b = np.array([roc_auc_score(y, s, sample_weight=w) for w in W])
        acc_b = wacc(y, p_t, W)
        store[(arm, model)] = dict(ap=ap_b, auc=auc_b, acc=acc_b)
        folds = g.fold.to_numpy()
        pf = [average_precision_score(y[folds == f], s[folds == f]) for f in np.unique(folds)]
        ci = lambda v: np.percentile(v, [2.5, 97.5])
        rows.append(dict(arm=arm, model=model, n=len(y),
                         auprc=average_precision_score(y, s), auprc_lo=ci(ap_b)[0], auprc_hi=ci(ap_b)[1],
                         auroc=roc_auc_score(y, s), auroc_lo=ci(auc_b)[0], auroc_hi=ci(auc_b)[1],
                         accuracy=(p_t == y).mean(), acc_lo=ci(acc_b)[0], acc_hi=ci(acc_b)[1],
                         accuracy_at_0p5=(p_5 == y).mean(),
                         precision=precision_score(y, p_t), recall=recall_score(y, p_t), f1=f1_score(y, p_t),
                         mcc=matthews_corrcoef(y, p_t), fold_auprc_mean=np.mean(pf), fold_auprc_sd=np.std(pf),
                         fold_auprc_min=np.min(pf)))
    res = pd.DataFrame(rows)
    res["arm"] = pd.Categorical(res.arm, ["A", "B", "L", "O"], ordered=True)
    res = res.sort_values(["model", "arm"]).reset_index(drop=True)
    res.to_csv(OUT / "comparison.csv", index=False)

    dl = []
    pairs = [("B", "A"), ("L", "A"), ("B", "L"), ("O", "L"), ("O", "A")]
    for model in sorted(d.model.unique()):
        for x, z in pairs:
            if (x, model) not in store or (z, model) not in store:
                continue
            for met in ("ap", "auc", "acc"):
                dv = store[(x, model)][met] - store[(z, model)][met]
                base = {"ap": "auprc", "auc": "auroc", "acc": "accuracy"}[met]
                pt = res[(res.arm == x) & (res.model == model)][base].iloc[0] - res[(res.arm == z) & (res.model == model)][base].iloc[0]
                lo, hi = np.percentile(dv, [2.5, 97.5])
                dl.append(dict(model=model, metric=base, contrast=f"{x} - {z}", delta=pt, lo=lo, hi=hi,
                               significant=bool(lo > 0 or hi < 0)))
    pd.DataFrame(dl).to_csv(OUT / "paired_deltas.csv", index=False)

    # annotation-vs-evidence split of the long-only arm (XGBoost)
    try:
        lp = pd.read_parquet(OUT / "oof_long_parts.parquet")
        out = []
        for arm, g in lp.groupby("arm"):
            g = g.assign(key=g.sv_id + "_" + g.replicate.astype(str)).set_index("key").loc[keys]
            s = g.score.to_numpy()
            out.append(dict(arm=arm, auprc=average_precision_score(y, s), auroc=roc_auc_score(y, s),
                            accuracy=((s >= g.thr.to_numpy()) == y).mean()))
        pd.DataFrame(out).to_csv(OUT / "long_parts.csv", index=False)
    except FileNotFoundError:
        pass

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import precision_recall_curve
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    models = ["xgb", "rf", "lr", "mlp"]
    w = 0.2
    cols = {"A": "#c0392b", "B": "#2471a3", "L": "#27ae60", "O": "#7f8c8d"}
    for i, a in enumerate(["A", "B", "L", "O"]):
        r = res[res.arm == a].set_index("model").reindex(models)
        ax[0].bar(np.arange(4) + (i - 1.5) * w, r.auprc, w, yerr=[r.auprc - r.auprc_lo, r.auprc_hi - r.auprc],
                  label=NAMES[a], color=cols[a], capsize=2)
    ax[0].set_xticks(range(4)); ax[0].set_xticklabels(models); ax[0].set_ylim(0.9, 1.0)
    ax[0].set_ylabel("AUPRC (held-out chromosomes, 95% cluster-bootstrap CI)"); ax[0].legend(fontsize=8)
    ax[0].set_title("AUPRC by arm and model")
    for a in ["A", "B", "L", "O"]:
        g = d[(d.arm == a) & (d.model == "xgb")].assign(key=lambda t: t.sv_id + "_" + t.replicate.astype(str)).set_index("key").loc[keys]
        p, r_, _ = precision_recall_curve(y, g.score.to_numpy())
        ax[1].plot(r_, p, label=f"{NAMES[a]}  AP={average_precision_score(y, g.score.to_numpy()):.3f}", color=cols[a])
    ax[1].set_xlabel("recall"); ax[1].set_ylabel("precision"); ax[1].set_ylim(0.8, 1.01); ax[1].legend(fontsize=8)
    ax[1].set_title("Precision-recall (XGBoost, out-of-fold)")
    plt.tight_layout(); plt.savefig(OUT / "comparison.png", dpi=160)
    print(res.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
