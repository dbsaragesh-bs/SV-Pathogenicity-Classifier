"""P1.9: place SVs on pseudo-samples (batches) so windows never overlap/touch inside a batch.

Greedy, seeded, label-balanced:  each batch gets ~equal pathogenic/benign, spread across chromosomes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.utils.common import PROCESSED, load_cfg, stable_seed


def make_batches(df: pd.DataFrame, batch_size: int, min_gap: int, flank: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    df = df.copy()
    df["ws"] = df.start - flank - min_gap // 2
    df["we"] = df.end + flank + min_gap // 2
    # interleave labels so each batch is balanced: shuffle within label, then zip
    pos = df[df.label == 1].sample(frac=1, random_state=int(rng.integers(1e9))).index.tolist()
    neg = df[df.label == 0].sample(frac=1, random_state=int(rng.integers(1e9))).index.tolist()
    order = []
    while pos or neg:
        if pos:
            order.append(pos.pop())
        if neg:
            order.append(neg.pop())
    batches: list[list[int]] = []
    occupied: list[dict] = []  # per batch: chrom -> list of (ws, we)
    cnt: list[list[int]] = []  # per batch label counts
    for i in order:
        r = df.loc[i]
        placed = False
        for b, occ in enumerate(occupied):
            if len(batches[b]) >= batch_size:
                continue
            if any(not (r.we <= s or r.ws >= e) for s, e in occ.get(r.chrom, [])):
                continue
            # keep labels balanced: refuse if this label already exceeds half+1
            if cnt[b][int(r.label)] > batch_size // 2:
                continue
            occ.setdefault(r.chrom, []).append((r.ws, r.we))
            batches[b].append(i)
            cnt[b][int(r.label)] += 1
            placed = True
            break
        if not placed:
            occupied.append({r.chrom: [(r.ws, r.we)]})
            batches.append([i])
            cnt.append([0, 0])
            cnt[-1][int(r.label)] += 1
    out = []
    for b, idxs in enumerate(batches):
        for i in idxs:
            out.append((df.loc[i, "sv_id"], f"B{b:03d}"))
    return pd.DataFrame(out, columns=["sv_id", "batch"])


def main():
    cfg = load_cfg("truthset")
    sim = load_cfg("simulation")
    t = pd.read_parquet(PROCESSED / "truthset_matched.parquet")
    b = make_batches(t, cfg["batch_size"], cfg["batch_min_gap"], sim["flank"], stable_seed("batches", cfg["seed"]))
    out = t.merge(b, on="sv_id")
    out.to_parquet(PROCESSED / "batches.parquet", index=False)
    g = out.groupby("batch").agg(n=("sv_id", "size"), pos=("label", "sum"), mb=("svlen", lambda x: x.sum() / 1e6))
    print(g.describe())


if __name__ == "__main__":
    main()
