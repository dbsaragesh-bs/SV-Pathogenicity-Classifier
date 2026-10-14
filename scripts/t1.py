import pandas as pd, numpy as np
from src.truthset.parse_sources import parse_clinvar, parse_dbvar_ci
from src.utils.common import RAW, load_cfg
cfg = load_cfg("truthset")
cv = parse_clinvar(RAW/"variant_summary.txt.gz", cfg)
print(len(cv)); print(cv.label.value_counts())
ci = parse_dbvar_ci(RAW/"nstd102.GRCh38.variant_call.vcf.gz")
print(ci.head(), len(ci), ci.ci_width.describe())
cv = cv.merge(ci, on="variation_id", how="left")
print("with dbvar:", cv.imprecise.notna().mean())
cv = cv[(cv.label>=0)&cv.review.isin(cfg["clinvar_ok_review"])]
cv["lb"]=np.floor(np.log10(cv.svlen)*2)/2
print(cv.groupby(["lb","label"]).size().unstack().fillna(0))
print(cv.groupby(["svtype","label"]).size().unstack())
print(cv.groupby(["coord_exact","label"]).size().unstack())
print(cv.review.value_counts())
