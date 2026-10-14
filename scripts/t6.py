import pandas as pd
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50); pd.set_option("display.max_rows", 200)
d = pd.read_parquet("/home/dbs/dna_data/interim/feat/pilot/P000_r0.parquet")
print(d.shape)
L = [c for c in d.columns if c.startswith(("l_","b_")) and not c.startswith("b_ann_")]
print(d[["sv_id","svtype","svlen","label"]+L].round(3).T.to_string())
print(d.isna().mean().sort_values(ascending=False).head(15))
