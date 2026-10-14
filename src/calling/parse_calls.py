"""Parse per-caller VCFs into one call table and merge callers within a technology.

Call record: caller, chrom, start (0-based), end, svtype, svlen, qual, filter, support fields,
ci_start_w, ci_end_w, precise.  Only DEL/DUP are kept (primary scope).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pysam

KEEP = {"DEL", "DUP"}


def _g(rec, key, default=None):
    try:
        v = rec.info.get(key, default)
    except (ValueError, KeyError):
        return default
    return v


def _ci_width(rec, key):
    v = _g(rec, key)
    if not v:
        return 0
    return int(abs(v[1] - v[0]))


def _first(x):
    if isinstance(x, (tuple, list)):
        return x[0] if x else None
    return x


def parse_vcf(path: Path, caller: str) -> pd.DataFrame:
    rows = []
    if not Path(path).exists():
        return pd.DataFrame()
    vf = pysam.VariantFile(str(path))
    sample = list(vf.header.samples)[0] if len(vf.header.samples) else None
    for rec in vf:
        t = _g(rec, "SVTYPE")
        if isinstance(t, tuple):
            t = t[0]
        if t not in KEEP:
            continue
        start = rec.pos - 1
        end = rec.stop
        d = dict(caller=caller, chrom=rec.chrom, start=start, end=end, svtype=t, svlen=end - start,
                 qual=rec.qual if rec.qual is not None else np.nan,
                 filter=";".join(rec.filter.keys()) if len(rec.filter) else "PASS",
                 call_id=rec.id,
                 ci_start_w=_ci_width(rec, "CIPOS"), ci_end_w=_ci_width(rec, "CIEND"),
                 precise=0 if "IMPRECISE" in rec.info else 1)
        fm = rec.samples[sample] if sample else {}
        if caller == "manta":
            pr = fm.get("PR")
            sr = fm.get("SR")
            d.update(pr_ref=pr[0] if pr else np.nan, pr_alt=pr[1] if pr else np.nan,
                     sr_ref=sr[0] if sr else np.nan, sr_alt=sr[1] if sr else np.nan,
                     gq=fm.get("GQ"), homlen=_g(rec, "HOMLEN", 0))
            hl = _g(rec, "HOMLEN")
            d["homlen"] = _first(hl) if hl is not None else 0
        elif caller == "delly":
            d.update(pe=_g(rec, "PE", 0), sr=_g(rec, "SR", 0), mapq=_g(rec, "MAPQ", np.nan),
                     srq=_g(rec, "SRQ", np.nan), ce=_g(rec, "CE", np.nan),
                     dv=fm.get("DV"), dr=fm.get("DR"), rv=fm.get("RV"), rr=fm.get("RR"), gq=fm.get("GQ"))
        elif caller == "sniffles":
            dv, dr = fm.get("DV"), fm.get("DR")
            d.update(support=_first(_g(rec, "SUPPORT")), coverage=_first(_g(rec, "COVERAGE")),
                     af=(dv / (dv + dr)) if (dv is not None and dr is not None and dv + dr > 0) else np.nan,
                     stdev_pos=_g(rec, "STDEV_POS", np.nan), stdev_len=_g(rec, "STDEV_LEN", np.nan),
                     gq=fm.get("GQ"), dv=dv, dr=dr)
            st = _g(rec, "STRAND")
            d["strand_bal"] = _strand_bal(st)
        elif caller == "cutesv":
            d.update(support=_first(_g(rec, "RE")), af=np.nan, gq=fm.get("GQ"),
                     dv=fm.get("DV"), dr=fm.get("DR"))
            d["strand_bal"] = _strand_bal(_g(rec, "STRAND"))
        rows.append(d)
    return pd.DataFrame(rows)


def _strand_bal(st):
    if st is None:
        return np.nan
    if isinstance(st, tuple):
        st = st[0]
    # "+-" style strings from sniffles/cuteSV encode orientation, not read strand balance
    return np.nan


def load_batch_calls(d: Path, tech: str) -> pd.DataFrame:
    """tech in {short, hifi, ont}. Returns all caller records for the batch."""
    parts = []
    if tech == "short":
        parts.append(parse_vcf(d / "short.manta.vcf.gz", "manta"))
        parts.append(parse_vcf(d / "short.delly.vcf.gz", "delly"))
    else:
        parts.append(parse_vcf(d / f"{tech}.sniffles.vcf.gz", "sniffles"))
        parts.append(parse_vcf(d / f"{tech}.cutesv.vcf.gz", "cutesv"))
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def merge_callers(calls: pd.DataFrame, ro: float = 0.5, bp_tol: int = 1000) -> pd.DataFrame:
    """Merge records of different callers into consensus calls.

    Two records link when they have the same SV type and either reciprocal overlap >= ro or both
    breakpoints within bp_tol.  Components are the consensus; coordinates = median over callers.
    Adds n_callers and keeps per-caller columns suffixed by caller (first record per caller).
    """
    if calls.empty:
        return calls
    from src.truthset.intervals import UnionFind, reciprocal_overlap

    calls = calls.reset_index(drop=True)
    uf = UnionFind(len(calls))
    for (c, t), g in calls.groupby(["chrom", "svtype"]):
        g = g.sort_values("start")
        idx = g.index.to_numpy()
        s, e = g.start.to_numpy(), g.end.to_numpy()
        for a in range(len(idx)):
            b = a + 1
            while b < len(idx) and s[b] <= e[a] + bp_tol:
                if reciprocal_overlap(s[a], e[a], s[b], e[b]) >= ro or (abs(s[a] - s[b]) <= bp_tol and abs(e[a] - e[b]) <= bp_tol):
                    uf.union(idx[a], idx[b])
                b += 1
    calls["cid"] = [uf.find(i) for i in range(len(calls))]
    out = []
    for cid, g in calls.groupby("cid"):
        r = {"cid": cid, "chrom": g.chrom.iloc[0], "svtype": g.svtype.iloc[0],
             "start": int(np.median(g.start)), "end": int(np.median(g.end)),
             "n_callers": g.caller.nunique(), "n_records": len(g),
             "any_pass": int((g["filter"] == "PASS").any())}
        r["svlen"] = r["end"] - r["start"]
        for caller, gg in g.groupby("caller"):
            first = gg.sort_values("qual", ascending=False).iloc[0]
            for k, v in first.items():
                if k in ("caller", "chrom", "svtype", "call_id"):
                    continue
                r[f"{caller}_{k}"] = v
            r[f"{caller}_start"] = int(first.start)
            r[f"{caller}_end"] = int(first.end)
        out.append(r)
    return pd.DataFrame(out)
