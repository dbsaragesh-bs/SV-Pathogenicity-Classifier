import time, pandas as pd
from pyfaidx import Fasta
from src.features.annotation import Annotations
from src.utils.common import RAW
t=time.time()
ref = Fasta(str(RAW/"GRCh38_no_alt.fa"))
a = Annotations(ref)
print("load", time.time()-t, "genes", sum(len(g) for g in a.genes.values()))
t=time.time()
# 22q11.2 deletion region approx hg38 chr22:18,900,000-21,500,000 ; SHANK3 etc
for c,s,e in [("chr22",18_900_000,21_500_000),("chr1",1_000_000,1_000_500),("chr17",44_000_000,44_400_000)]:
    f = a.interval_features(c,s,e)
    print(c,s,e,{k:(round(v,3) if isinstance(v,float) else v) for k,v in f.items()})
print("per-interval", (time.time()-t)/3)
