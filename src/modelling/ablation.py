"""P8 experiments beyond the main comparison: ladder, leave-one-family-out, confound controls.

python -m src.modelling.ablation --which ladder,logo,controls --trials 12
Writes reports/ablation/<name>/{oof,metrics}_<tag>.* ; summarised by analysis_extra.py
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.features.registry import arm_columns, feature_groups
from src.modelling.run_experiment import cohort_mask, run
from src.utils.common import PROCESSED, REPORTS

GENE_DOSAGE = re.compile(r"(loeuf|_pli|n_constrained)")


def ladder_cols(cols_all):
    A = arm_columns("A", cols_all)
    B = arm_columns("B", cols_all)
    extra = [c for c in B if c not in A]
    long_geo_ann = [c for c in extra if c.startswith(("b_ann_", "b_")) and not c.startswith("b_mh_")]
    long_geo_ann = [c for c in long_geo_ann if not c.startswith("b_mh_")]
    long_ev = [c for c in extra if c.startswith("l_") and "depth" not in c]
    long_depth = [c for c in extra if c.startswith("l_") and "depth" in c]
    mh = [c for c in extra if c.startswith("b_mh_")]
    coords_only = [c for c in long_geo_ann if c.startswith("b_ann_") or c == "b_svlen_log"]
    return {
        "R0_A": A,
        "R1_A+long_coords_annot": A + coords_only,
        "R2_+long_call_fields": A + long_geo_ann,
        "R3_+long_evidence": A + long_geo_ann + long_ev,
        "R4_+long_depth": A + long_geo_ann + long_ev + long_depth,
        "R5_+junction=B": B,
        "R6_Oracle": arm_columns("O", cols_all),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default=str(PROCESSED / "features_all.parquet"))
    ap.add_argument("--which", default="ladder,logo,controls")
    ap.add_argument("--trials", type=int, default=12)
    ap.add_argument("--model", default="xgb")
    a = ap.parse_args()
    df0 = pd.read_parquet(a.features)
    base = df0[cohort_mask(df0, "both")].reset_index(drop=True)
    cols = list(df0.columns)
    which = a.which.split(",")
    if "ladder" in which:
        lc = ladder_cols(cols)
        run(base, list(lc), [a.model], 5, a.trials, 0, REPORTS / "ablation" / "ladder", "ladder", arm_cols=lc)
    if "logo" in which:
        B = arm_columns("B", cols)
        fam = feature_groups(B)
        lc = {"B_full": B}
        for f, cs in fam.items():
            lc[f"B_minus_{f}"] = [c for c in B if c not in set(cs)]
        run(base, list(lc), [a.model], 5, a.trials, 0, REPORTS / "ablation" / "logo", "logo", arm_cols=lc)
    if "controls" in which:
        out = REPORTS / "ablation" / "controls"
        # (1) without gene-dosage / constraint predictors
        lc = {k: [c for c in arm_columns(k, cols) if not GENE_DOSAGE.search(c)] for k in ("A", "B", "O")}
        run(base, ["A", "B", "O"], [a.model], 5, a.trials, 0, out, "no_gene_dosage", arm_cols=lc)
        # (2) size removed
        lc = {k: [c for c in arm_columns(k, cols) if "svlen" not in c] for k in ("A", "B", "O")}
        run(base, ["A", "B", "O"], [a.model], 5, a.trials, 0, out, "no_size", arm_cols=lc)
        # (3) only SVs that hit an exon (gene context matched across classes)
        sub = base[base.o_ann_n_exons_hit > 0].reset_index(drop=True)
        print("exon-hit subset", len(sub), sub.label.value_counts().to_dict())
        run(sub, ["A", "B", "O"], [a.model], 5, a.trials, 0, out, "exon_hit_only")
        # (4) ClinVar-only labels (removes the gnomAD/ClinVar source difference)
        sub = base[(base.label == 1) | (base.benign_basis == "clinvar")].reset_index(drop=True)
        print("clinvar-only", len(sub), sub.label.value_counts().to_dict())
        run(sub, ["A", "B", "O"], [a.model], 5, a.trials, 0, out, "clinvar_labels_only")
        # (5) unique (non-segdup) breakpoints: removes detection-selection bias
        sub = base[base.o_ann_bp_segdup_max <= 0.5].reset_index(drop=True)
        run(sub, ["A", "B", "O"], [a.model], 5, a.trials, 0, out, "non_segdup_breakpoints")
        # (6) all planted SVs with detection flags (C_all)
        allc = df0.copy()
        allc["a_detected"] = allc.detected_short.astype(float)
        allc["b_detected"] = allc.detected_hifi.astype(float)
        run(allc, ["A", "B", "O"], [a.model], 5, a.trials, 0, out, "C_all_with_detection_flags")


if __name__ == "__main__":
    main()
