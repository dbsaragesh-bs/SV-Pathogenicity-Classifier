"""Tiny pilot batch (ClinVar only) to exercise the whole chain before the real truth set exists."""
import numpy as np, pandas as pd
from pyfaidx import Fasta
from src.truthset.parse_sources import parse_clinvar
from src.utils.common import RAW, PROCESSED, load_cfg
cfg = load_cfg("truthset")
cv = parse_clinvar(RAW / "variant_summary.txt.gz", cfg)
cv = cv[(cv.label >= 0) & cv.review.isin(cfg["clinvar_ok_review"]) & cv.coord_exact & cv.svlen.between(2000, 200000)]
fa = Fasta(str(RAW / "GRCh38_no_alt.fa"))
rng = np.random.default_rng(1)
picked = []
for lab, t, n in [(1, "DEL", 6), (1, "DUP", 4), (0, "DEL", 6), (0, "DUP", 4)]:
    g = cv[(cv.label == lab) & (cv.svtype == t)].sample(frac=1, random_state=int(rng.integers(1e6)))
    k = 0
    for _, r in g.iterrows():
        if any(r.chrom == p.chrom and abs(r.start - p.start) < 3_000_000 for p in picked):
            continue
        seq = str(fa[r.chrom][r.start - 100000:r.end + 100000]).upper()
        if len(seq) < r.svlen + 200000 or seq.count("N") / len(seq) > 0.001:
            continue
        picked.append(r); k += 1
        if k == n: break
df = pd.DataFrame(picked).reset_index(drop=True)
df["sv_id"] = [f"PV{i:03d}" for i in range(len(df))]
df["group_id"] = range(len(df)); df["zygosity"] = "het"; df["batch"] = "P000"
df.to_parquet(PROCESSED / "pilot_batches.parquet", index=False)
print(df[["sv_id", "chrom", "start", "end", "svtype", "svlen", "label"]])
