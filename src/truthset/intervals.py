"""Interval helpers: reciprocal overlap, clustering, connected components."""
from __future__ import annotations

import numpy as np


def reciprocal_overlap(s1, e1, s2, e2):
    """Vectorised reciprocal overlap fraction (min of the two overlap/length fractions)."""
    ov = np.minimum(e1, e2) - np.maximum(s1, s2)
    ov = np.clip(ov, 0, None)
    l1 = np.maximum(e1 - s1, 1)
    l2 = np.maximum(e2 - s2, 1)
    return np.minimum(ov / l1, ov / l2)


class UnionFind:
    def __init__(self, n: int):
        self.p = np.arange(n)

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def cluster_ro(chrom, start, end, thresh: float, same_key=None) -> np.ndarray:
    """Connected components of intervals linked by reciprocal overlap >= thresh.

    chrom/start/end are equal-length arrays. same_key (optional) is an array whose equal
    values are required for linking (e.g. svtype).  Returns integer component ids.
    """
    n = len(start)
    uf = UnionFind(n)
    chrom = np.asarray(chrom)
    start = np.asarray(start)
    end = np.asarray(end)
    key = None if same_key is None else np.asarray(same_key)
    for c in np.unique(chrom):
        idx = np.where(chrom == c)[0]
        idx = idx[np.argsort(start[idx], kind="stable")]
        s, e = start[idx], end[idx]
        for a in range(len(idx)):
            # candidates must start before e[a]; sorted by start so scan forward until s[b] >= e[a]
            b = a + 1
            while b < len(idx) and s[b] < e[a]:
                if key is None or key[idx[a]] == key[idx[b]]:
                    if reciprocal_overlap(s[a], e[a], s[b], e[b]) >= thresh:
                        uf.union(idx[a], idx[b])
                b += 1
    return np.array([uf.find(i) for i in range(n)])
