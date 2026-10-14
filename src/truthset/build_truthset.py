"""P1: build the labelled truth set, with a full filter waterfall.

Usage: python -m src.truthset.build_truthset
Writes data/processed/truthset_{natural,matched}.parquet and truthset_manifest.json
"""
from __future__ import annotations

import json
import subprocess

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from src.truthset.intervals import cluster_ro, reciprocal_overlap
from src.truthset.parse_sources import parse_clinvar, parse_dbvar_ci, parse_gnomad_sv
from src.utils.common import PROCESSED, RAW, ensure, load_cfg, sha256_file, stable_seed

REF = RAW / "GRCh38_no_alt.fa"


def log_step(wf: list, name: str, df: pd.DataFrame):
    row = {"step": name, "n": int(len(df)),
           "n_path": int((df.label == 1).sum()), "n_benign": int((df.label == 0).sum())}
    wf.append(row)
    print(f"[{name}] n={row['n']} path={row['n_path']} benign={row['n_benign']}", flush=True)


def window_n_frac(df: pd.DataFrame, fa: Fasta, flank: int) -> np.ndarray:
    out = np.ones(len(df))
    for i, (c, s, e) in enumerate(zip(df.chrom, df.start, df.end)):
        ws, we = s - flank, e + flank
        L = len(fa[c])
        if ws < 0 or we > L:
            continue
        # count N in the flanks and at the SV; use a coarse check on 3 segments
        seq = str(fa[c][ws:we]).upper()
        out[i] = seq.count("N") / len(seq)
    return out


def conflict_remove(df: pd.DataFrame, ro: float) -> pd.DataFrame:
    """Remove every record that has >=ro reciprocal overlap with an opposite-label record."""
    keep = np.ones(len(df), bool)
    for (c, t), g in df.groupby(["chrom", "svtype"]):
        pos = g[g.label == 1]
        neg = g[g.label == 0]
        if len(pos) == 0 or len(neg) == 0:
            continue
        ps, pe = pos.start.values[:, None], pos.end.values[:, None]
        ns, ne = neg.start.values[None, :], neg.end.values[None, :]
        m = reciprocal_overlap(ps, pe, ns, ne) >= ro
        keep[df.index.get_indexer(pos.index[m.any(1)])] = False
        keep[df.index.get_indexer(neg.index[m.any(0)])] = False
    return df[keep]


def cap_local_density(df: pd.DataFrame, pad: int, k: int, rng) -> pd.DataFrame:
    """Keep at most k SVs whose padded windows overlap (bounds hotspot dominance and batch count)."""
    order = rng.permutation(len(df))
    acc_ws: dict = {}
    acc_we: dict = {}
    keep = []
    for i in order:
        c, s, e = df.chrom.iat[i], df.start.iat[i] - pad, df.end.iat[i] + pad
        ws, we = acc_ws.setdefault(c, []), acc_we.setdefault(c, [])
        n_ov = sum(1 for a, b in zip(ws, we) if not (b <= s or a >= e))
        if n_ov < k:
            ws.append(s)
            we.append(e)
            keep.append(i)
    return df.iloc[sorted(keep)]


def main():
    cfg = load_cfg("truthset")
    wf: list = []
    rng = np.random.default_rng(cfg["seed"])
    out_dir = ensure(PROCESSED)

    # ---- harmonise
    cv = parse_clinvar(RAW / "variant_summary.txt.gz", cfg)
    wf.append({"step": "clinvar_raw_rows_all_types", "n": cv.attrs["n_raw_rows"]})
    log_step(wf, "clinvar_GRCh38_DELDUP_autosomal_all_labels", cv)
    ci = parse_dbvar_ci(RAW / "nstd102.GRCh38.variant_call.vcf.gz")
    cv = cv.merge(ci, on="variation_id", how="left")
    cv = cv[cv.label >= 0]
    log_step(wf, "clinvar_path_or_benign", cv)
    cv = cv[cv.review.isin(cfg["clinvar_ok_review"])]
    log_step(wf, "clinvar_review>=1star", cv)

    gn = parse_gnomad_sv(RAW / "gnomad.v4.1.sv.sites.vcf.gz", cfg, cfg["gnomad_min_af"])
    log_step(wf, "gnomad_common_PASS_DELDUP", gn)

    df = pd.concat([cv, gn], ignore_index=True)
    df["imprecise"] = df.imprecise.fillna(False).astype(bool)
    df["ci_width"] = df.ci_width.fillna(0)
    log_step(wf, "combined", df)

    # ---- coordinate quality & size
    bad_ci = (df.ci_width > cfg["max_ci_frac"] * df.svlen) & (df.ci_width > 1000)
    df = df[~bad_ci]
    log_step(wf, "drop_wide_CI", df)
    df = df[(df.svlen >= cfg["min_svlen"]) & (df.svlen <= cfg["max_svlen"])]
    log_step(wf, "size_filter", df)
    df = df[df.svtype.isin(cfg["types"])].reset_index(drop=True)

    # ---- gaps / windows
    fa = Fasta(str(REF))
    nfrac = window_n_frac(df, fa, cfg["flank"])
    df = df[nfrac <= cfg["max_window_n_frac"]].reset_index(drop=True)
    log_step(wf, "window_inside_chrom_and_N<=1%", df)

    # ---- conflicts, duplicates
    df = conflict_remove(df, cfg["conflict_ro"]).reset_index(drop=True)
    log_step(wf, "remove_pathogenic_benign_conflicts", df)

    df["dedup_id"] = cluster_ro(df.chrom, df.start, df.end, cfg["dedup_ro"], same_key=df.svtype + df.label.astype(str))
    # representative: highest stars, exact coords, narrowest CI, earliest id
    df = df.sort_values(["dedup_id", "stars", "coord_exact", "ci_width"],
                        ascending=[True, False, False, True])
    df["dedup_size"] = df.groupby("dedup_id").dedup_id.transform("size")
    df = df.drop_duplicates("dedup_id", keep="first").reset_index(drop=True)
    log_step(wf, "dedup_ro80", df)

    df = cap_local_density(df, cfg["locus_pad"], cfg["max_per_locus"], rng).reset_index(drop=True)
    log_step(wf, f"cap_local_density(<= {cfg['max_per_locus']} SVs per overlapping +-{cfg['locus_pad']} window)", df)

    df["group_id"] = cluster_ro(df.chrom, df.start, df.end, cfg["group_ro"])
    df["group_id"] = pd.factorize(df.group_id)[0]
    df["sv_id"] = [f"SV{i:06d}" for i in range(len(df))]

    # a group mixing labels would be a conflict at 50% -> already removed; sanity assert
    n_mixed = int((df.groupby("group_id").label.nunique() > 1).sum())
    wf.append({"step": "groups_with_mixed_labels (DEL/DUP of opposite class overlapping >=50%)", "n": n_mixed})
    print("groups with mixed labels:", n_mixed)

    # ---- size bands, matched sampling (common support)
    df["log10len"] = np.log10(df.svlen)
    w = cfg["size_bin_decades"]
    df["size_bin"] = np.floor(df.log10len / w).astype(int)
    natural = df.copy()
    keep_idx = []
    strata = []
    for (t, b), g in natural.groupby(["svtype", "size_bin"]):
        npos, nneg = int((g.label == 1).sum()), int((g.label == 0).sum())
        n = min(npos, nneg)
        strata.append({"svtype": t, "size_bin": int(b), "lo": float(10 ** (b * w)), "hi": float(10 ** ((b + 1) * w)),
                       "n_path": npos, "n_benign": nneg, "n_kept_each": n})
        if n == 0:
            continue
        for lab in (1, 0):
            gg = g[g.label == lab]
            # prefer higher-quality records, then random
            gg = gg.assign(_r=rng.random(len(gg))).sort_values(["stars", "coord_exact", "_r"], ascending=[False, False, True])
            keep_idx += list(gg.index[:n])
    matched_full = natural.loc[sorted(keep_idx)].reset_index(drop=True)
    log_step(wf, "matched_type_x_sizebin_common_support", matched_full)

    # compute-driven cap: subsample each stratum with the same fraction for both classes
    f = min(1.0, cfg["target_per_class"] / (matched_full.label == 1).sum())
    keep = []
    for (t, b), g in matched_full.groupby(["svtype", "size_bin"]):
        n_each = int((g.label == 1).sum())
        n_keep = min(n_each, max(int(np.ceil(n_each * f)), 3))
        for lab in (1, 0):
            gg = g[g.label == lab]
            gg = gg.assign(_r=rng.random(len(gg))).sort_values(["stars", "coord_exact", "_r"], ascending=[False, False, True])
            keep += list(gg.index[:n_keep])
    matched = matched_full.loc[sorted(keep)].reset_index(drop=True)
    log_step(wf, f"matched_subsampled_to_~{cfg['target_per_class']}_per_class", matched)

    # ---- label-independent genotype
    for d in (natural, matched, matched_full):
        d["zygosity"] = cfg["zygosity"]

    natural.to_parquet(out_dir / "truthset_natural.parquet", index=False)
    matched.to_parquet(out_dir / "truthset_matched.parquet", index=False)
    matched_full.to_parquet(out_dir / "truthset_matched_full.parquet", index=False)
    pd.DataFrame(strata).to_csv(out_dir / "matching_strata.csv", index=False)
    manifest = {
        "build": cfg["build"], "config": cfg, "waterfall": wf,
        "inputs": {p.name: sha256_file(p) for p in
                   [RAW / "variant_summary.txt.gz", RAW / "nstd102.GRCh38.variant_call.vcf.gz"]},
        "matched_source_counts": matched.groupby(["label", "source"]).size().reset_index(name="n").to_dict("records"),
    }
    (out_dir / "truthset_manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    print("wrote truth sets")


if __name__ == "__main__":
    main()
