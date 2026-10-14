"""Depth sweep analysis: does extra long-read depth change the A-vs-B gap, and how does detection scale?"""
from __future__ import annotations

import glob
import re

import numpy as np
import pandas as pd

from src.modelling.run_experiment import run
from src.utils.common import INTERIM, PROCESSED, REPORTS, ensure


def main(trials=12):
    base = pd.read_parquet(PROCESSED / "features_all.parquet")
    batches = sorted({re.match(r"(B\d+)_r0", p.split("/")[-1]).group(1) for p in glob.glob(str(INTERIM / "feat" / "sweep" / "B*_r0_*.parquet")) if "_fp" not in p})
    sub = base[base.batch.isin(batches)].reset_index(drop=True)
    variants = {"L30_S30": sub}
    for kind, d in (("L", 20), ("L", 10), ("L", 5), ("S", 15), ("S", 60)):
        parts = [pd.read_parquet(INTERIM / "feat" / "sweep" / f"{b}_r0_{kind}{d}.parquet") for b in batches]
        variants[f"{kind}{d}"] = pd.concat(parts, ignore_index=True)
    key = ["sv_id", "replicate"]
    det = pd.DataFrame({n: v.set_index(key)[["detected_short", "detected_hifi"]].stack() for n, v in variants.items()})
    rows = []
    for n, v in variants.items():
        rows.append(dict(config=n, n=len(v), det_short=v.detected_short.mean(), det_long=v.detected_hifi.mean(),
                         det_long_path=v[v.label == 1].detected_hifi.mean(), det_long_benign=v[v.label == 0].detected_hifi.mean(),
                         med_bp_err_long=v.hifi_bp_err_start.abs().median(), med_bp_err_short=v.short_bp_err_start.abs().median()))
    det_tab = pd.DataFrame(rows)
    # common cohort: detected by both technologies in every configuration
    ok = None
    for n, v in variants.items():
        s = set(map(tuple, v[(v.detected_short == 1) & (v.detected_hifi == 1)][key].values))
        ok = s if ok is None else ok & s
    out = ensure(REPORTS / "ablation" / "depth")
    res = []
    for n, v in variants.items():
        vv = v[[tuple(x) in ok for x in v[key].values]].reset_index(drop=True)
        p, m = run(vv, ["A", "B"], ["xgb"], 5, trials, 0, out, f"depth_{n}")
        g = m.groupby("arm")[["auprc", "auroc"]].mean()
        res.append(dict(config=n, n=len(vv), A_auprc=g.loc["A", "auprc"], B_auprc=g.loc["B", "auprc"],
                        A_auroc=g.loc["A", "auroc"], B_auroc=g.loc["B", "auroc"]))
    tab = pd.DataFrame(res).merge(det_tab, on="config", suffixes=("", "_all"))
    tab.to_csv(out / "depth_summary.csv", index=False)
    det_tab.to_csv(out / "depth_detection.csv", index=False)
    print(tab.round(4).to_string())


if __name__ == "__main__":
    main()
