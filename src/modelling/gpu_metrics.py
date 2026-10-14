"""GPU-batched weighted AUPRC / AUROC and cluster bootstrap (torch).

One sort of the scores is shared by all B bootstrap replicates; each replicate is a vector of
row weights (multiplicity of the row's cluster), so B x n cumulative sums run on the GPU.
"""
from __future__ import annotations

import numpy as np
import torch


def device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _prep(y, s, dev):
    y = torch.as_tensor(np.asarray(y), dtype=torch.float32, device=dev)
    s = torch.as_tensor(np.asarray(s), dtype=torch.float32, device=dev)
    order = torch.argsort(s, descending=True, stable=True)
    return y[order], s[order], order


def weighted_ap(y_sorted, s_sorted, W):
    """W: [B, n] row weights in sorted order. Average precision with tie handling by score groups."""
    pos = W * y_sorted  # [B,n]
    tp = torch.cumsum(pos, 1)
    fp = torch.cumsum(W * (1 - y_sorted), 1)
    # last index of each tie group
    last = torch.ones_like(s_sorted, dtype=torch.bool)
    last[:-1] = s_sorted[:-1] != s_sorted[1:]
    P = tp[:, -1:].clamp_min(1e-12)
    prec = tp / (tp + fp).clamp_min(1e-12)
    rec = tp / P
    rec_l = rec[:, last]
    prec_l = prec[:, last]
    drec = rec_l - torch.cat([torch.zeros_like(rec_l[:, :1]), rec_l[:, :-1]], 1)
    return (drec * prec_l).sum(1)


def weighted_auc(y_sorted, s_sorted, W):
    pos = W * y_sorted
    neg = W * (1 - y_sorted)
    # AUC = P(score_pos > score_neg) + 0.5 P(tie); process ascending
    # use descending order: for each negative, count positives ranked above (strictly) + half ties
    last = torch.ones_like(s_sorted, dtype=torch.bool)
    last[:-1] = s_sorted[:-1] != s_sorted[1:]
    first = torch.ones_like(last)
    first[1:] = last[:-1]
    tp = torch.cumsum(pos, 1)
    gid = torch.cumsum(first.long(), 0) - 1
    ng = int(gid[-1]) + 1
    B = W.shape[0]
    pos_g = torch.zeros(B, ng, device=W.device).index_add_(1, gid, pos)
    neg_g = torch.zeros(B, ng, device=W.device).index_add_(1, gid, neg)
    pos_above = torch.cumsum(pos_g, 1) - pos_g
    num = (neg_g * (pos_above + 0.5 * pos_g)).sum(1)
    den = (pos.sum(1) * neg.sum(1)).clamp_min(1e-12)
    return num / den


def point_metrics(y, s):
    dev = device()
    ys, ss, _ = _prep(y, s, dev)
    W = torch.ones(1, len(ys), device=dev)
    return {"auprc": float(weighted_ap(ys, ss, W)[0]), "auroc": float(weighted_auc(ys, ss, W)[0])}


def cluster_bootstrap(y, scores: dict[str, np.ndarray], groups, B=2000, seed=0):
    """Paired cluster bootstrap. scores: name -> score array (same rows). Returns dict name -> {auprc: [B], auroc: [B]}."""
    dev = device()
    g = np.asarray(groups)
    uniq, inv = np.unique(g, return_inverse=True)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(uniq), size=(B, len(uniq)))
    counts = np.zeros((B, len(uniq)), dtype=np.float32)
    for b in range(B):
        counts[b] = np.bincount(draws[b], minlength=len(uniq))
    Wrow = torch.as_tensor(counts, device=dev)[:, torch.as_tensor(inv, device=dev)]  # [B,n]
    out = {}
    for name, s in scores.items():
        ys, ss, order = _prep(y, s, dev)
        W = Wrow[:, order]
        out[name] = {"auprc": weighted_ap(ys, ss, W).cpu().numpy(), "auroc": weighted_auc(ys, ss, W).cpu().numpy()}
    return out


def ci(a, lo=2.5, hi=97.5):
    a = np.asarray(a)
    return float(np.percentile(a, lo)), float(np.percentile(a, hi))
