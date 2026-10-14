"""Outer split = held-out chromosome groups (identical for every arm); inner = grouped stratified CV."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def chrom_groups(df: pd.DataFrame, k: int = 5, seed: int = 0) -> dict[str, int]:
    """Greedy balanced partition of chromosomes into k groups by (label x type) row counts."""
    key = df.label.astype(str) + "_" + df.svtype
    cnt = pd.crosstab(df.chrom, key).reindex(AUTOSOMES).fillna(0)
    rng = np.random.default_rng(seed)
    order = (cnt.sum(1) + rng.random(len(cnt)) * 1e-3).sort_values(ascending=False).index
    load = np.zeros((k, cnt.shape[1]))
    assign = {}
    for c in order:
        v = cnt.loc[c].to_numpy()
        # place where the squared relative imbalance after adding is smallest
        tot = load.sum(0) + v + 1e-9
        costs = [(((load[g] + v) / tot - 1.0 / k) ** 2).sum() + 0.0 * g for g in range(k)]
        g = int(np.argmin(costs))
        load[g] += v
        assign[c] = g
    return assign


def outer_folds(df: pd.DataFrame, k: int = 5, seed: int = 0):
    m = chrom_groups(df, k, seed)
    fold = df.chrom.map(m).to_numpy()
    for f in range(k):
        yield f, np.where(fold != f)[0], np.where(fold == f)[0]


def inner_cv(y, groups, n_splits=3, seed=0):
    s = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(s.split(np.zeros(len(y)), y, groups))
