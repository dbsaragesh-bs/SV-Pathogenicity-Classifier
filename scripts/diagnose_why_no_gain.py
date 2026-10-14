"""Diagnostics for the presentation (read-only wrt project data; writes to $DNA_TMP).
 (1) which feature FAMILY alone carries the label (fixed-param XGBoost, held-out chromosome groups)
 (2) stress test: add breakpoint noise of scale s bp to the SHORT-read arm's annotation and see when long reads start to help."""
import sys, time, json, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from pyfaidx import Fasta
from xgboost import XGBClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from src.features.annotation import Annotations
from src.features.registry import arm_columns
from src.modelling.splits import outer_folds
from src.utils.common import RAW
import os

f = pd.read_parquet('data/processed/features_all.parquet')
df = f[(f.detected_short==1)&(f.detected_hifi==1)].reset_index(drop=True)
print('cohort', len(df), flush=True)
y = df.label.values; groups = df.group_id.values

def oof(cols, X=None, seed=0):
    X = df[cols].to_numpy(float) if X is None else X
    s = np.zeros(len(df))
    for fold, tr, te in outer_folds(df, 5, 0):
        m = XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                          min_child_weight=2, reg_lambda=1.0, n_jobs=8, random_state=seed, tree_method='hist',
                          eval_metric='logloss')
        m.fit(X[tr], y[tr]); s[te] = m.predict_proba(X[te])[:, 1]
    return s

def boot_delta(s1, s0, n=1000, seed=1):
    rng = np.random.default_rng(seed); g = pd.factorize(groups)[0]; G = g.max()+1
    idx = [np.where(g==k)[0] for k in range(G)]
    d = []
    for _ in range(n):
        pick = rng.integers(0, G, G); ii = np.concatenate([idx[k] for k in pick])
        d.append(average_precision_score(y[ii], s1[ii]) - average_precision_score(y[ii], s0[ii]))
    return np.percentile(d, [2.5, 97.5])

allcols = [c for c in arm_columns('B', df.columns)]
famA_ann = [c for c in arm_columns('A', df.columns) if c.startswith('a_ann_')]
famA_geo_ev = [c for c in arm_columns('A', df.columns) if not c.startswith('a_ann_')]
famB_ann = [c for c in arm_columns('B', df.columns) if c.startswith('b_ann_')]
famB_ev = [c for c in arm_columns('B', df.columns) if c.startswith(('b_','l_')) and not c.startswith('b_ann_')]
o_ann = [c for c in df.columns if c.startswith('o_ann_')]
print('family sizes', {k: len(v) for k, v in dict(A_ann=famA_ann, A_geo_ev=famA_geo_ev, B_ann=famB_ann, B_ev=famB_ev, O_ann=o_ann).items()}, flush=True)

res = []
sets = {
 'short evidence only (no annotation)': famA_geo_ev,
 'long evidence only (no annotation)': famB_ev,
 'short+long evidence only (no annotation)': famA_geo_ev + famB_ev,
 'annotation only, short coords': famA_ann,
 'annotation only, long coords': famB_ann,
 'annotation only, TRUE coords': o_ann,
 'Arm A (all short)': arm_columns('A', df.columns),
 'Arm B (short+long)': arm_columns('B', df.columns),
}
S = {}
for k, cols in sets.items():
    cols = [c for c in cols if df[c].nunique(dropna=True) > 1]
    s = oof(cols); S[k] = s
    print(f'{k:45s} n_feat={len(cols):3d} AUPRC={average_precision_score(y, s):.3f} AUROC={roc_auc_score(y, s):.3f}', flush=True)
    res.append(dict(set=k, n_feat=len(cols), auprc=average_precision_score(y, s), auroc=roc_auc_score(y, s)))
ex = (df.o_ann_n_exons_hit.fillna(0) > 0).astype(float).values
print(f'single binary feature "hits an exon"         AUPRC={average_precision_score(y, ex):.3f} AUROC={roc_auc_score(y, ex):.3f}', flush=True)
res.append(dict(set='single feature: hits an exon', n_feat=1, auprc=average_precision_score(y, ex), auroc=roc_auc_score(y, ex)))
pd.DataFrame(res).to_csv(os.path.join(os.environ['DNA_TMP'], 'family_alone.csv'), index=False)

# ---- stress test ----
ref = Fasta(str(RAW / 'GRCh38_no_alt.fa'))
t0 = time.time(); ann = Annotations(ref); print('annotations loaded', round(time.time()-t0), 's', flush=True)
base_A = [c for c in arm_columns('A', df.columns) if c.startswith('a_ann_')]
colsA = [c for c in arm_columns('A', df.columns) if df[c].nunique(dropna=True) > 1]
colsB = [c for c in arm_columns('B', df.columns) if df[c].nunique(dropna=True) > 1]
out = []
for sd in [0, 25, 100, 250, 500, 1000, 2500, 5000]:
    rng = np.random.default_rng(42)
    dN = df.copy()
    st = df.start.values + rng.normal(0, sd, len(df)); en = df.end.values + rng.normal(0, sd, len(df))
    a = np.round(np.minimum(st, en)).astype(int); b = np.round(np.maximum(st, en)).astype(int)
    b = np.maximum(b, a + 1); a = np.maximum(a, 0)
    rows = [ann.interval_features(c, x, z, prefix='a_ann_') for c, x, z in zip(df.chrom.values, a, b)] if sd > 0 else None
    if rows is not None:
        nd = pd.DataFrame(rows)
        for c in nd.columns:
            if c in dN.columns:
                dN[c] = nd[c].values
    XA = dN[colsA].to_numpy(float); XB = dN[colsB].to_numpy(float)
    sA = oof(colsA, XA); sB = oof(colsB, XB)
    ap_a, ap_b = average_precision_score(y, sA), average_precision_score(y, sB)
    ci = boot_delta(sB, sA)
    exA = (dN['a_ann_n_exons_hit'].fillna(0) > 0).values; exT = (df.o_ann_n_exons_hit.fillna(0) > 0).values
    flips = float((exA != exT).mean())
    print(f'noise sd={sd:5d} bp | exon-status flips vs truth {flips*100:5.1f}% | AUPRC A {ap_a:.4f}  B {ap_b:.4f}  B-A {ap_b-ap_a:+.4f} [{ci[0]:+.4f},{ci[1]:+.4f}]', flush=True)
    out.append(dict(noise_sd_bp=sd, exon_flip_pct=flips*100, auprc_A=ap_a, auprc_B=ap_b, delta=ap_b-ap_a, lo=ci[0], hi=ci[1]))
pd.DataFrame(out).to_csv(os.path.join(os.environ['DNA_TMP'], 'noise_stress.csv'), index=False)
print('DONE', flush=True)
