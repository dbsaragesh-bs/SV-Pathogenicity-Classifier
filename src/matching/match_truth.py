"""P4: attach planted truth to called SVs (type-aware reciprocal-overlap matching)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from src.truthset.intervals import reciprocal_overlap

LEVELS = {"loose": 0.3, "std": 0.5, "strict": 0.8}


def match_batch(truth: pd.DataFrame, calls: pd.DataFrame, level: str = "std") -> pd.DataFrame:
    """For every planted SV in `truth` find the best consensus call.

    Returns one row per planted SV: matched (bool), call index, ro, bp errors, n_candidates.
    calls must have chrom,start,end,svtype (consensus coordinates).
    """
    thr = LEVELS[level]
    rows = []
    for _, sv in truth.iterrows():
        r = {"sv_id": sv.sv_id, "matched": False, "call_idx": -1, "ro": 0.0, "bp_err_start": np.nan,
             "bp_err_end": np.nan, "n_candidates": 0}
        if len(calls):
            c = calls[(calls.chrom == sv.chrom) & (calls.svtype == sv.svtype)]
            if len(c):
                ro = reciprocal_overlap(sv.start, sv.end, c.start.values, c.end.values)
                ok = ro >= thr
                r["n_candidates"] = int(ok.sum())
                if ok.any():
                    j = np.argmax(np.where(ok, ro, -1))
                    ci = c.index[j]
                    r.update(matched=True, call_idx=int(ci), ro=float(ro[j]),
                             bp_err_start=float(c.start.values[j] - sv.start),
                             bp_err_end=float(c.end.values[j] - sv.end))
        rows.append(r)
    return pd.DataFrame(rows)


def false_positive_calls(truth: pd.DataFrame, calls: pd.DataFrame, windows: pd.DataFrame, min_ro: float = 0.1):
    """Calls inside a simulated window that overlap no planted SV of the same type (RO >= min_ro)."""
    if calls.empty:
        return calls
    keep = []
    for i, c in calls.iterrows():
        in_win = ((windows.chrom == c.chrom) & (windows.ws <= c.start) & (windows.we >= c.end)).any()
        if not in_win:
            continue
        t = truth[(truth.chrom == c.chrom) & (truth.svtype == c.svtype)]
        ro = reciprocal_overlap(c.start, c.end, t.start.values, t.end.values) if len(t) else np.array([0.0])
        if ro.max(initial=0) < min_ro:
            keep.append(i)
    return calls.loc[keep]
