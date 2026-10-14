"""Long-read-only arm (L): same cohort, folds, models and tuning budget as the main experiment.

python -m src.modelling.long_only --trials 20
Writes reports/ablation/long_only/{oof,metrics}_long_only.* (+ long_parts for annotation-vs-evidence split).
Arm L = every column derived from long reads only: long caller fields, long alignment evidence, junction
microhomology, and gene/repeat annotation recomputed from the LONG-READ called coordinates. No short-read column.
"""
from __future__ import annotations

import argparse

import pandas as pd

from src.features.registry import arm_columns
from src.modelling.run_experiment import cohort_mask, run
from src.utils.common import PROCESSED, REPORTS


def long_columns(cols):
    L = [c for c in arm_columns("B", cols) if c.startswith(("b_", "l_"))]
    return {
        "L": L,
        "L_ann": [c for c in L if c.startswith("b_ann_")],
        "L_ev": [c for c in L if not c.startswith("b_ann_")],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default=str(PROCESSED / "features_all.parquet"))
    ap.add_argument("--trials", type=int, default=20)
    ap.add_argument("--cohort", default="both")
    a = ap.parse_args()
    df0 = pd.read_parquet(a.features)
    base = df0[cohort_mask(df0, a.cohort)].reset_index(drop=True)
    lc = long_columns(list(df0.columns))
    print({k: len(v) for k, v in lc.items()}, "rows", len(base), flush=True)
    out = REPORTS / "ablation" / "long_only"
    run(base, ["L"], ["xgb", "lr", "rf", "mlp"], 5, a.trials, 0, out, "long_only", arm_cols=lc)
    run(base, ["L_ann", "L_ev"], ["xgb"], 5, a.trials, 0, out, "long_parts", arm_cols=lc)
    print("LONG_ONLY_DONE", flush=True)


if __name__ == "__main__":
    main()
