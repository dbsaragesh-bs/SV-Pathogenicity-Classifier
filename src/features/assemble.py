"""Concatenate per-batch feature tables into data/processed/features_all.parquet (+ schema checks)."""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from src.utils.common import INTERIM, PROCESSED, ensure

def main(tag="base", out_name="features_all.parquet"):
    parts = sorted((INTERIM / "feat" / tag).glob("B*_r*.parquet"))
    parts = [p for p in parts if not p.name.endswith("_fp.parquet")]
    df = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    assert df.sv_id.notna().all()
    assert not df.duplicated(["sv_id", "replicate"]).any(), "duplicate sv/replicate rows"
    num = df.select_dtypes("number")
    assert not np.isinf(num.to_numpy(float)).any(), "inf values in features"
    ensure(PROCESSED)
    df.to_parquet(PROCESSED / out_name, index=False)
    fp = [pd.read_parquet(p) for p in sorted((INTERIM / "feat" / tag).glob("B*_fp.parquet"))]
    if fp:
        pd.concat(fp, ignore_index=True).to_parquet(PROCESSED / out_name.replace(".parquet", "_fp.parquet"), index=False)
    print(f"{len(df)} rows from {len(parts)} batch tables; columns {df.shape[1]}")
    for t in ("short", "hifi", "ont"):
        if f"detected_{t}" in df:
            print(f"  detected_{t}: {df[f'detected_{t}'].mean():.3f}")

if __name__ == "__main__":
    main(*(sys.argv[1:]))
