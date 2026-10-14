"""Model zoo with an identical tuning budget per arm.  GPU: XGBoost (cuda) and a torch MLP."""
from __future__ import annotations

import warnings

import numpy as np
import optuna
import torch
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import QuantileTransformer, StandardScaler

from src.modelling.splits import inner_cv

optuna.logging.set_verbosity(optuna.logging.WARNING)
warnings.filterwarnings("ignore")
DEV = "cuda" if torch.cuda.is_available() else "cpu"


# ------------------------------------------------------------------ XGBoost (GPU)
def _xgb(params, seed, n_jobs=4):
    return xgb.XGBClassifier(
        tree_method="hist", device=DEV, n_jobs=n_jobs, random_state=seed, eval_metric="aucpr",
        max_depth=params["max_depth"], learning_rate=params["lr"], n_estimators=params["n_est"],
        subsample=params["subsample"], colsample_bytree=params["colsample"],
        min_child_weight=params["mcw"], reg_lambda=params["lambda"], verbosity=0)


def _xgb_space(t):
    return dict(max_depth=t.suggest_int("max_depth", 2, 6), lr=t.suggest_float("lr", 0.02, 0.3, log=True),
                n_est=t.suggest_int("n_est", 80, 500, step=20), subsample=t.suggest_float("subsample", 0.5, 1.0),
                colsample=t.suggest_float("colsample", 0.4, 1.0), mcw=t.suggest_int("mcw", 1, 10),
                **{"lambda": t.suggest_float("lambda", 1e-2, 10, log=True)})


# ------------------------------------------------------------------ Logistic regression
def _lr(params, seed):
    return Pipeline([("imp", SimpleImputer(strategy="median", add_indicator=True)),
                     ("qt", QuantileTransformer(n_quantiles=200, output_distribution="normal", random_state=seed)),
                     ("sc", StandardScaler()),
                     ("m", LogisticRegression(C=params["C"], penalty="elasticnet", l1_ratio=params["l1"],
                                              solver="saga", max_iter=3000, random_state=seed))])


def _lr_space(t):
    return dict(C=t.suggest_float("C", 1e-3, 10, log=True), l1=t.suggest_float("l1", 0.0, 0.8))


# ------------------------------------------------------------------ Random forest
def _rf(params, seed):
    return RandomForestClassifier(n_estimators=params["n_est"], min_samples_leaf=params["leaf"],
                                  max_features=params["mf"], n_jobs=8, random_state=seed)


def _rf_space(t):
    return dict(n_est=t.suggest_int("n_est", 200, 600, step=100), leaf=t.suggest_int("leaf", 1, 8),
                mf=t.suggest_float("mf", 0.1, 0.6))


# ------------------------------------------------------------------ MLP (GPU)
class MLPClf:
    def __init__(self, hidden=64, depth=2, drop=0.2, wd=1e-3, lr=2e-3, epochs=200, seed=0, n_ens=2):
        self.h, self.d, self.drop, self.wd, self.lr, self.ep, self.seed, self.n_ens = hidden, depth, drop, wd, lr, epochs, seed, n_ens
        self.pre = Pipeline([("imp", SimpleImputer(strategy="median", add_indicator=True)),
                             ("qt", QuantileTransformer(n_quantiles=200, output_distribution="normal", random_state=seed)),
                             ("sc", StandardScaler())])

    def _net(self, d_in):
        layers, d = [], d_in
        for _ in range(self.d):
            layers += [torch.nn.Linear(d, self.h), torch.nn.GELU(), torch.nn.Dropout(self.drop)]
            d = self.h
        layers.append(torch.nn.Linear(d, 1))
        return torch.nn.Sequential(*layers).to(DEV)

    def fit(self, X, y):
        Xp = torch.as_tensor(self.pre.fit_transform(X), dtype=torch.float32, device=DEV).clamp(-6, 6)
        yt = torch.as_tensor(y, dtype=torch.float32, device=DEV)
        self.nets = []
        for e in range(self.n_ens):
            torch.manual_seed(self.seed + e)
            net = self._net(Xp.shape[1])
            opt = torch.optim.AdamW(net.parameters(), lr=self.lr, weight_decay=self.wd)
            n = len(yt)
            for ep in range(self.ep):
                net.train()
                perm = torch.randperm(n, device=DEV)
                for i in range(0, n, 256):
                    idx = perm[i:i + 256]
                    opt.zero_grad()
                    loss = torch.nn.functional.binary_cross_entropy_with_logits(net(Xp[idx]).squeeze(1), yt[idx])
                    loss.backward()
                    opt.step()
            net.eval()
            self.nets.append(net)
        return self

    def predict_proba(self, X):
        Xp = torch.as_tensor(self.pre.transform(X), dtype=torch.float32, device=DEV).clamp(-6, 6)
        with torch.no_grad():
            p = torch.stack([torch.sigmoid(n(Xp).squeeze(1)) for n in self.nets]).mean(0).cpu().numpy()
        return np.c_[1 - p, p]


def _mlp(params, seed):
    return MLPClf(hidden=params["hidden"], depth=params["depth"], drop=params["drop"], wd=params["wd"],
                  lr=params["lr"], epochs=params["epochs"], seed=seed)


def _mlp_space(t):
    return dict(hidden=t.suggest_categorical("hidden", [32, 64, 128]), depth=t.suggest_int("depth", 1, 3),
                drop=t.suggest_float("drop", 0.0, 0.5), wd=t.suggest_float("wd", 1e-5, 1e-1, log=True),
                lr=t.suggest_float("lr", 1e-3, 1e-2, log=True), epochs=t.suggest_int("epochs", 30, 90, step=30))


MODELS = {
    "xgb": (_xgb, _xgb_space),
    "lr": (_lr, _lr_space),
    "rf": (_rf, _rf_space),
    "mlp": (_mlp, _mlp_space),
}


def fit_tuned(name, X, y, groups, n_trials=20, seed=0, inner_splits=3):
    """Optuna on grouped inner CV (maximise mean AUPRC). Returns fitted model, best params, inner OOF scores."""
    make, space = MODELS[name]
    folds = inner_cv(y, groups, inner_splits, seed)

    def objective(trial):
        p = space(trial)
        ap = []
        for tr, va in folds:
            m = make(p, seed).fit(X[tr], y[tr])
            ap.append(average_precision_score(y[va], m.predict_proba(X[va])[:, 1]))
        return float(np.mean(ap))

    study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=seed))
    study.optimize(objective, n_trials=n_trials)
    best = study.best_params
    oof = np.zeros(len(y))
    for tr, va in folds:
        oof[va] = make(best, seed).fit(X[tr], y[tr]).predict_proba(X[va])[:, 1]
    final = make(best, seed).fit(X, y)
    return final, best, oof, study.best_value
