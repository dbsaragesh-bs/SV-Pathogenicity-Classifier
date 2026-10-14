import numpy as np, pandas as pd
from src.truthset.build_truthset import cap_local_density
# approximate pre-cap set by re-reading natural (which is capped) is not possible; rebuild quickly
import src.truthset.build_truthset as bt
from src.truthset.parse_sources import parse_clinvar, parse_gnomad_sv, parse_dbvar_ci
from src.utils.common import RAW, load_cfg
cfg = load_cfg("truthset")
cv = parse_clinvar(RAW/"variant_summary.txt.gz", cfg)
cv = cv[(cv.label>=0)&cv.review.isin(cfg["clinvar_ok_review"])]
gn = parse_gnomad_sv(RAW/"gnomad.v4.1.sv.sites.vcf.gz", cfg, 0.01)
df = pd.concat([cv,gn],ignore_index=True)
df = df[(df.svlen>=50)&(df.svlen<=1_000_000)].drop_duplicates(["chrom","start","end","svtype","label"]).reset_index(drop=True)
print(len(df), df.label.value_counts().to_dict())
rng = np.random.default_rng(0)
for pad,k in [(100000,10),(100000,20),(250000,12),(250000,20),(50000,30)]:
    r = cap_local_density(df,pad,k,np.random.default_rng(0))
    print(pad,k,len(r),r.label.value_counts().to_dict())
