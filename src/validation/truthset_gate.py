"""P1 gate: schema, duplicates, conflicts, matched-distribution tests, trivial baselines."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score

from src.modelling.splits import outer_folds
from src.truthset.intervals import reciprocal_overlap
from src.utils.common import PROCESSED, REPORTS, ensure


def schema_checks(d: pd.DataFrame) -> dict:
    res = {}
    res["no_null_keys"] = bool(d[["sv_id", "chrom", "start", "end", "svtype", "label"]].notna().all().all())
    res["start_lt_end"] = bool((d.start < d.end).all())
    res["svlen_consistent"] = bool((d.svlen == d.end - d.start).all())
    res["chrom_allowed"] = bool(d.chrom.isin([f"chr{i}" for i in range(1, 23)]).all())
    res["label_values"] = sorted(d.label.unique().tolist()) == [0, 1]
    res["unique_sv_id"] = bool(d.sv_id.is_unique)
    # conflicts: any cross-label pair with RO >= 0.5 (same chrom/type)
    n_conf = 0
    for (c, t), g in d.groupby(["chrom", "svtype"]):
        p, n = g[g.label == 1], g[g.label == 0]
        if len(p) and len(n):
            ro = reciprocal_overlap(p.start.values[:, None], p.end.values[:, None], n.start.values[None], n.end.values[None])
            n_conf += int((ro >= 0.5).sum())
    res["cross_label_conflicts"] = n_conf
    # duplicates within label at RO >= 0.8
    n_dup = 0
    for (c, t, l), g in d.groupby(["chrom", "svtype", "label"]):
        if len(g) > 1:
            ro = reciprocal_overlap(g.start.values[:, None], g.end.values[:, None], g.start.values[None], g.end.values[None])
            np.fill_diagonal(ro, 0)
            n_dup += int((ro >= 0.8).sum() // 2)
    res["near_duplicates_ro80"] = n_dup
    return res


def trivial_baselines(d: pd.DataFrame) -> dict:
    d = d.reset_index(drop=True)
    d["svlen_log"] = np.log10(d.svlen)
    d["dup"] = (d.svtype == "DUP").astype(int)
    d["chrom_idx"] = d.chrom.str[3:].astype(int)
    d["pos"] = d.start
    out = {}
    for name, cols in {"size_only": ["svlen_log"], "type_only": ["dup"], "chrom_only": ["chrom_idx"],
                       "position_only": ["pos"], "size+type": ["svlen_log", "dup"]}.items():
        X = d[cols].to_numpy(float)
        oof = np.full(len(d), np.nan)
        for f, tr, te in outer_folds(d, 5, 0):
            m = GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=0).fit(X[tr], d.label.values[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]
        out[name] = {"auprc": float(average_precision_score(d.label, oof)), "auroc": float(roc_auc_score(d.label, oof)),
                     "prevalence": float(d.label.mean())}
    return out


def main():
    d = pd.read_parquet(PROCESSED / "truthset_matched.parquet")
    res = {"n": len(d), "n_path": int((d.label == 1).sum()), "schema": schema_checks(d)}
    a, b = d[d.label == 1].svlen, d[d.label == 0].svlen
    res["size_ks"] = {"stat": float(ks_2samp(np.log10(a), np.log10(b)).statistic), "p": float(ks_2samp(np.log10(a), np.log10(b)).pvalue)}
    res["type_counts"] = pd.crosstab(d.svtype, d.label).to_dict()
    res["source_counts"] = pd.crosstab(d.source, d.label).to_dict()
    res["baselines"] = trivial_baselines(d)
    ensure(REPORTS / "data_audit")
    (REPORTS / "data_audit" / "truthset_gate.json").write_text(json.dumps(res, indent=1, default=str))
    L = ["# Truth-set gate (P1)", "", f"Matched set: {res['n']} SVs ({res['n_path']} pathogenic).", "", "## Schema / integrity", ""]
    L += [f"- {k}: {v}" for k, v in res["schema"].items()]
    L += ["", f"Size distribution KS (log10 length): D={res['size_ks']['stat']:.3f}, p={res['size_ks']['p']:.3g}", "",
          "## Trivial baselines (held-out chromosome groups, expect ≈ chance)", "", "| baseline | AUPRC | AUROC |", "|---|---|---|"]
    L += [f"| {k} | {v['auprc']:.3f} | {v['auroc']:.3f} |" for k, v in res["baselines"].items()]
    L += ["", "## Type x label", "", pd.crosstab(d.svtype, d.label).to_markdown(), "", "## Source x label", "",
          pd.crosstab(d.source, d.label).to_markdown()]
    (REPORTS / "data_audit" / "truthset_gate.md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
