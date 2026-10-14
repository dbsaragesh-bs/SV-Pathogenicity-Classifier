import pandas as pd
from src.utils.common import INTERIM, BATCH_FILE
from src.calling.parse_calls import load_batch_calls, merge_callers
from src.matching.match_truth import match_batch
d = INTERIM/"sim"/"P000_r0"
truth = pd.read_parquet(BATCH_FILE)
raw = load_batch_calls(d, "short")
print(raw.groupby(["caller","svtype"]).size())
m = merge_callers(raw)
print("merged", len(m))
for lv in ("loose","std","strict"):
    r = match_batch(truth, m, lv)
    print(lv, r.matched.mean(), r.groupby(truth.label.values).matched.mean().to_dict())
r = match_batch(truth, m, "std")
print(pd.concat([truth[["sv_id","chrom","start","end","svtype","svlen","label"]], r[["matched","ro","bp_err_start","bp_err_end","n_candidates"]]], axis=1).to_string())
print(m[["chrom","start","end","svtype","n_callers"]].head(40))
