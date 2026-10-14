import pandas as pd, numpy as np
d = pd.read_parquet("/home/dbs/dna_data/interim/feat/base/B000_r0.parquet")
fp = pd.read_parquet("/home/dbs/dna_data/interim/feat/base/B000_r0_fp.parquet")
print(d.shape, "FP calls", len(fp), fp.groupby("tech").size().to_dict())
for t in ("short","hifi"):
    print(t, "det", d[f"detected_{t}"].mean(), {lv: d[f"det_{t}_{lv}"].mean().round(3) for lv in ("loose","std","strict")})
d["sz"]=pd.cut(d.svlen,[0,300,1000,5000,30000,1e7])
print(d.groupby(["sz","label"],observed=True)[["detected_short","detected_hifi"]].mean().round(2))
for t in ("short","hifi"):
    print(t, "bp err abs median", d[f"{t}_bp_err_start"].abs().median(), d[f"{t}_bp_err_end"].abs().median(), "p90", d[f"{t}_bp_err_start"].abs().quantile(.9))
