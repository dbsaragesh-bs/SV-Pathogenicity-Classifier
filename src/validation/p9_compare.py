"""P9 comparison: simulated vs real HG002 deletions (same code path)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.stats import ks_2samp, wasserstein_distance
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from src.features.registry import arm_columns
from src.utils.common import INTERIM, PROCESSED, REPORTS, ensure

XG = dict(tree_method="hist", device="cuda", n_estimators=200, max_depth=3, learning_rate=0.05, verbosity=0)


def main():
    real = pd.read_parquet(INTERIM / "feat" / "real" / "R000_r0.parquet")
    sim = pd.read_parquet(PROCESSED / "features_all.parquet")
    hi = real.svlen.max()
    simb = sim[(sim.svtype == "DEL") & (sim.label == 0) & (sim.svlen <= hi)].copy()
    simb["a_delly_only"] = (simb.a_delly_present == 1)
    L = ["# Real-data realism check (P9): simulated vs GIAB HG002 deletions", "",
         f"Real: {len(real)} HG002 v5.0q GRCh38 benchmark DEL, size {int(real.svlen.min())}-{int(hi)} bp; "
         f"Illumina 2x250 novoalign and Revio HiFi streamed from GIAB, subsampled to ~30x, run through the SAME Delly/Sniffles2/cuteSV "
         f"and feature code (Manta crashed on the real novoalign BAM -> Delly only for short reads).",
         f"Simulated comparison set: {len(simb)} benign DEL <= {int(hi)} bp.", "",
         "## Detection (recall of planted/benchmark SV)", "", "| tech | simulated | real |", "|---|---|---|",
         f"| short, Delly only | {simb.a_delly_only.mean():.3f} | {real.detected_short.mean():.3f} |",
         f"| short, all callers | {simb.detected_short.mean():.3f} | n/a |",
         f"| HiFi | {simb.detected_hifi.mean():.3f} | {real.detected_hifi.mean():.3f} |", ""]
    fp = pd.read_parquet(PROCESSED / "features_all_fp.parquet") if (PROCESSED / "features_all_fp.parquet").exists() else None
    res = {"detect": {"sim_short_delly": float(simb.a_delly_only.mean()), "real_short": float(real.detected_short.mean()),
                      "sim_hifi": float(simb.detected_hifi.mean()), "real_hifi": float(real.detected_hifi.mean())}}
    # feature shift on SVs detected in both
    def shift(prefixes, det):
        r = real[real[det] == 1]
        s = simb[simb[det] == 1]
        cols = [c for c in real.columns if c.startswith(prefixes) and c in s.columns and not c.endswith(("depth_in", "depth_flank"))
                and real[c].nunique() > 1 and s[c].nunique() > 1 and pd.api.types.is_numeric_dtype(real[c])
                and "manta" not in c and "_cutesv_" not in c or False]
        rows = []
        for c in cols:
            a, b = r[c].dropna().astype(float), s[c].dropna().astype(float)
            if len(a) < 10 or len(b) < 10:
                continue
            sd = pd.concat([a, b]).std() or 1.0
            rows.append(dict(feature=c, ks=ks_2samp(a, b).statistic, wass_std=wasserstein_distance(a, b) / sd,
                             real_median=a.median(), sim_median=b.median()))
        return pd.DataFrame(rows).sort_values("ks", ascending=False), cols, r, s
    out = ensure(REPORTS / "data_audit")
    for tag, pref, det in (("short", ("s_",), "detected_short"), ("long", ("l_",), "detected_hifi")):
        t, cols, r, s = shift(pref, det)
        t.to_csv(out / f"p9_shift_{tag}.csv", index=False)
        # adversarial validation sim vs real on those evidence features
        X = pd.concat([r[cols], s[cols]]).to_numpy(float)
        y = np.r_[np.ones(len(r)), np.zeros(len(s))]
        oof = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(X, y):
            oof[te] = xgb.XGBClassifier(random_state=0, **XG).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        auc = roc_auc_score(y, oof)
        res[f"adversarial_auc_{tag}"] = float(auc)
        L += [f"## {tag}-read evidence features: real vs simulated", "",
              f"Adversarial validation AUROC (real vs simulated, {len(cols)} features, n_real={len(r)}, n_sim={len(s)}): **{auc:.3f}** (0.5 = indistinguishable).", "",
              "Most shifted features (KS):", "", t.head(10).round(3).to_markdown(index=False), ""]
    # apply models trained on simulation to real SVs
    cohort = sim[(sim.detected_short == 1) & (sim.detected_hifi == 1)]
    both = real[(real.detected_short == 1) & (real.detected_hifi == 1)]
    L += ["## Simulation-trained models applied to real HG002 benign-proxy deletions", ""]
    for arm in ("A", "B"):
        cols = [c for c in arm_columns(arm, sim.columns) if cohort[c].nunique() > 1 and c in real.columns]
        m = xgb.XGBClassifier(random_state=0, n_estimators=300, max_depth=3, learning_rate=0.05, tree_method="hist", device="cuda", verbosity=0)
        m.fit(cohort[cols].to_numpy(float), cohort.label.values)
        rr = real[(real.detected_short == 1)] if arm == "A" else both
        if len(rr):
            p = m.predict_proba(rr[cols].to_numpy(float))[:, 1]
            res[f"real_pred_pathogenic_frac_{arm}"] = float((p >= 0.5).mean())
            res[f"n_real_scored_{arm}"] = int(len(rr))
            L.append(f"- Arm {arm}: {len(rr)} real SVs scored; predicted pathogenic (p>=0.5): {(p >= 0.5).mean():.3f}; mean p {p.mean():.3f}")
    L += ["", "Real HG002 SVs are germline variants of a healthy individual, so a low predicted-pathogenic fraction is expected. "
          "Real pathogenic SVs with reads are not available at scale; this check therefore tests false-positive behaviour and "
          "feature realism, not sensitivity."]
    (out / "p9_realdata.md").write_text("\n".join(L))
    (out / "p9_realdata.json").write_text(json.dumps(res, indent=1))
    print("\n".join(L))


if __name__ == "__main__":
    main()
