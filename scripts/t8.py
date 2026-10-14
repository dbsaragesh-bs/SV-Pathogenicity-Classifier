import pandas as pd
pd.set_option("display.width",250)
d = pd.read_parquet("/home/dbs/dna_data/interim/feat/base/B000_r0.parquet")
m = d[(d.svlen>5000)&(d.label==0)]
print(m[["sv_id","svtype","svlen","source","detected_short","detected_hifi","o_ann_segdup_frac","o_ann_rep_allrep_frac","o_ann_bp_l_segdup","o_ann_bp_r_segdup","o_ann_bp_repeat_max","o_ann_rep_satellite_frac","o_ann_gap_overlap_frac"]].round(2).to_string())
