import pandas as pd, numpy as np, xgboost as xgb
from sklearn.metrics import roc_auc_score
from src.features.registry import arm_columns
from src.modelling.splits import outer_folds
d = pd.read_parquet("data/processed/features_partial.parquet")
d = d[(d.detected_short==1)&(d.detected_hifi==1)].reset_index(drop=True)
cols = [c for c in arm_columns("A", d.columns) if d[c].nunique()>1]
X = d[cols].to_numpy(float); y = d.label.values
m = xgb.XGBClassifier(n_estimators=200,max_depth=3,device="cuda",verbosity=0).fit(X,y)
imp = pd.Series(m.feature_importances_,index=cols).sort_values(ascending=False)
print(imp.head(15))
print(pd.crosstab(d.source,d.label))
def cv(mask, label):
    sub = d[mask].reset_index(drop=True); Xs = sub[cols].to_numpy(float); ys=label(sub)
    oof=np.full(len(sub),np.nan)
    for f,tr,te in outer_folds(sub,3,0):
        if len(set(ys[tr]))<2 or len(set(ys[te]))<2: continue
        oof[te]=xgb.XGBClassifier(n_estimators=150,max_depth=3,device="cuda",verbosity=0).fit(Xs[tr],ys[tr]).predict_proba(Xs[te])[:,1]
    ok=~np.isnan(oof); return roc_auc_score(ys[ok],oof[ok]), ok.sum(), ys[ok].mean()
print("clinvar path vs clinvar benign", cv((d.source=="clinvar"), lambda s: s.label.values))
print("clinvar path vs gnomad benign", cv(((d.label==1)|(d.source=="gnomad_sv")), lambda s: s.label.values))
print("gnomad benign vs clinvar benign (source detect)", cv((d.label==0), lambda s:(s.source=="gnomad_sv").astype(int).values))
for bins,name in [((0,1000),"<1k"),((1000,30000),"1-30k"),((30000,2e6),">30k")]:
    print(name, cv(d.svlen.between(*bins), lambda s: s.label.values))
