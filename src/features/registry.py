"""Feature registry: which columns belong to which arm, plus the forbidden list.

Column prefixes (assigned in build_features.py):
  a_ , s_     short-read call geometry / caller fields / alignment evidence   -> Arm A
  a_ann_      annotation recomputed from the SHORT call coordinates           -> Arm A
  b_ , l_     long-read call geometry / alignment evidence / junction         -> Arm B (extra)
  b_ann_      annotation recomputed from the LONG call coordinates            -> Arm B (extra)
  o_ann_ ,o_  annotation / junction from the TRUE planted coordinates         -> Oracle
Absolute depth columns are excluded everywhere (coverage knobs are nuisance variables).
"""
from __future__ import annotations

import re

import pandas as pd

FORBIDDEN_EXACT = {
    "sv_id", "chrom", "start", "end", "svtype", "svlen", "label", "group_id", "source", "stars",
    "coord_exact", "benign_basis", "zygosity", "batch", "replicate", "origin",
}
FORBIDDEN_RE = re.compile(r"(bp_err|_ro$|_n_cand$|^det_|^detected_|dedup|cluster|cov_|depth_in$|depth_flank$)")

NUISANCE_EXCLUDE = re.compile(r"(depth_in$|depth_flank$)")  # absolute depth: nuisance


def _ok(c: str) -> bool:
    return c not in FORBIDDEN_EXACT and not FORBIDDEN_RE.search(c) and not NUISANCE_EXCLUDE.search(c)


def arm_columns(arm: str, columns) -> list[str]:
    cols = [c for c in columns if _ok(c)]
    a = [c for c in cols if c.startswith(("a_", "s_"))]
    b = [c for c in cols if c.startswith(("b_", "l_"))]
    o = [c for c in cols if c.startswith("o_")]
    if arm == "A":
        return a
    if arm == "B":
        return a + b
    if arm == "O":  # oracle: truth-coordinate annotation/junction + shared evidence from both techs
        coord_derived = re.compile(r"^(a|b)_(ann_|svlen_log|mh_)")
        keep = [c for c in a + b if not coord_derived.match(c)]
        return keep + o
    if arm == "OA":  # oracle on short evidence only (upper bound for short-read-only information)
        coord_derived = re.compile(r"^a_(ann_|svlen_log)")
        return [c for c in a if not coord_derived.match(c)] + o
    raise ValueError(arm)


def feature_groups(columns) -> dict[str, list[str]]:
    """Groups for leave-one-group-out ablations (long-read side)."""
    g: dict[str, list[str]] = {}
    for c in columns:
        if not _ok(c):
            continue
        if c.startswith("b_ann_"):
            fam = "long_annotation"
        elif c.startswith("a_ann_"):
            fam = "short_annotation"
        elif c.startswith("o_"):
            fam = "oracle_truth"
        elif c.startswith("l_") and ("depth" in c):
            fam = "long_depth"
        elif c.startswith("l_") and ("bp_disp" in c or "refined" in c):
            fam = "long_breakpoint_precision"
        elif c.startswith("l_"):
            fam = "long_evidence"
        elif c.startswith("b_mh_"):
            fam = "long_junction_microhomology"
        elif c.startswith("b_"):
            fam = "long_call_geometry"
        elif c.startswith("s_") and "depth" in c:
            fam = "short_depth"
        elif c.startswith("s_"):
            fam = "short_evidence"
        elif c.startswith("a_"):
            fam = "short_call_geometry"
        else:
            continue
        g.setdefault(fam, []).append(c)
    return g
