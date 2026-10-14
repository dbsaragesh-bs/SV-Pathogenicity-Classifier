"""P4+P5 for one batch: call parsing -> truth matching -> evidence + annotation features.

python -m src.features.build_features --batch B000 --replicate 0 [--techs short,hifi] [--tag base]
Output: $DNA_DATA/interim/feat/<tag>/<batch>_r<rep>.parquet  (one row per planted SV)
        $DNA_DATA/interim/feat/<tag>/<batch>_r<rep>_fp.parquet (false-positive calls)
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from pyfaidx import Fasta

from src.calling.parse_calls import load_batch_calls, merge_callers
from src.features.annotation import Annotations
from src.features.evidence import insert_stats, long_features, microhomology, short_features
from src.matching.match_truth import LEVELS, false_positive_calls, match_batch
from src.utils.common import BATCH_FILE, INTERIM, PROCESSED, RAW, ensure, load_cfg

REF = RAW / "GRCh38_no_alt.fa"
_G: dict = {}


def _init():
    _G["ref"] = Fasta(str(REF))
    _G["bam"] = {}


def _bam(path):
    if path not in _G["bam"]:
        _G["bam"][path] = pysam.AlignmentFile(path, "rb")
    return _G["bam"][path]


def geometry_short(c: pd.Series) -> dict:
    g = {"svtype_dup": float(c.svtype == "DUP"), "svlen_log": float(np.log10(max(c.svlen, 1))),
         "n_callers": c.n_callers, "any_pass": c.any_pass,
         "ci_w": float(np.nanmax([c.get("manta_ci_start_w", 0), c.get("manta_ci_end_w", 0),
                                   c.get("delly_ci_start_w", 0), c.get("delly_ci_end_w", 0)]))}
    g["manta_present"] = float(pd.notna(c.get("manta_qual", np.nan)))
    g["delly_present"] = float(pd.notna(c.get("delly_qual", np.nan)))
    for k in ("manta_qual", "manta_pr_alt", "manta_sr_alt", "manta_pr_ref", "manta_sr_ref", "manta_homlen", "manta_precise",
              "delly_qual", "delly_pe", "delly_sr", "delly_mapq", "delly_precise", "delly_dv", "delly_rv"):
        g[k] = c.get(k, np.nan)
    return g


def geometry_long(c: pd.Series) -> dict:
    g = {"svtype_dup": float(c.svtype == "DUP"), "svlen_log": float(np.log10(max(c.svlen, 1))),
         "n_callers": c.n_callers, "any_pass": c.any_pass}
    g["sniffles_present"] = float(pd.notna(c.get("sniffles_qual", np.nan)))
    g["cutesv_present"] = float(pd.notna(c.get("cutesv_qual", np.nan)))
    for k in ("sniffles_qual", "sniffles_support", "sniffles_af", "sniffles_stdev_pos", "sniffles_stdev_len",
              "sniffles_precise", "cutesv_qual", "cutesv_support", "cutesv_precise", "cutesv_dv", "cutesv_dr"):
        g[k] = c.get(k, np.nan)
    return g


def _row(args):
    sv, scall, lcalls, d, techs, ins = args
    ann: Annotations = _G_ANN
    ref = _G["ref"]
    clen = len(ref[sv["chrom"]])
    out = {k: sv[k] for k in ("sv_id", "chrom", "start", "end", "svtype", "svlen", "label", "group_id", "source",
                               "stars", "coord_exact", "benign_basis", "zygosity", "batch")}
    # oracle: truth coordinates
    out.update({f"o_ann_{k}": v for k, v in ann.interval_features(sv["chrom"], sv["start"], sv["end"]).items()})
    out.update({f"o_{k}": v for k, v in microhomology(ref, sv["chrom"], sv["start"], sv["end"]).items()})
    for tech in techs:
        m = scall[tech]
        det = m is not None
        out[f"detected_{tech}"] = int(det)
        if not det:
            continue
        c, mm = m
        out[f"{tech}_bp_err_start"] = mm["bp_err_start"]
        out[f"{tech}_bp_err_end"] = mm["bp_err_end"]
        out[f"{tech}_ro"] = mm["ro"]
        out[f"{tech}_n_cand"] = mm["n_candidates"]
        pre = "a" if tech == "short" else "b"
        geo = geometry_short(c) if tech == "short" else geometry_long(c)
        out.update({f"{pre}_{k}": v for k, v in geo.items()})
        out.update({f"{pre}_ann_{k}": v for k, v in ann.interval_features(c.chrom, c.start, c.end).items()})
        bam = _bam(str(d / ("short.bam" if tech == "short" else f"{tech}.bam")))
        if tech == "short":
            out.update(short_features(bam, c.chrom, int(c.start), int(c.end), ins[0], ins[1], clen))
        else:
            out.update(long_features(bam, c.chrom, int(c.start), int(c.end), clen))
            out.update({f"b_{k}": v for k, v in microhomology(ref, c.chrom, int(c.start), int(c.end)).items()})
    return out


_G_ANN = None


def run(batch, replicate, techs, tag, workers=12, level="std", suffix=""):
    global _G_ANN
    d = INTERIM / "sim" / f"{batch}_r{replicate}{suffix}"
    outdir = ensure(INTERIM / "feat" / tag)
    truth = pd.read_parquet(BATCH_FILE)
    truth = truth[truth.batch == batch].reset_index(drop=True)
    wins = pd.read_csv(d / "windows.bed", sep="\t", header=None, names=["chrom", "ws", "we", "sv_id"])
    ref = Fasta(str(REF))
    if _G_ANN is None:
        _G_ANN = Annotations(ref)
    calls, matches, fp = {}, {}, []
    for tech in techs:
        raw = load_batch_calls(d, tech)
        calls[tech] = merge_callers(raw) if len(raw) else raw
        matches[tech] = match_batch(truth, calls[tech], level)
        # multi-level detection for stringency stability
        for lv in LEVELS:
            mlv = match_batch(truth, calls[tech], lv)
            truth[f"det_{tech}_{lv}"] = mlv.matched.values.astype(int)
        fp.append(false_positive_calls(truth, calls[tech], wins).assign(tech=tech))
    # insert-size stats for short reads
    ins = (400.0, 50.0)
    if "short" in techs:
        sb = pysam.AlignmentFile(str(d / "short.bam"))
        w0 = wins.iloc[0]
        ins = insert_stats(sb, w0.chrom, int(w0.ws), int(w0.we))
    tasks = []
    for i, sv in truth.iterrows():
        sc = {}
        for tech in techs:
            m = matches[tech].iloc[i]
            sc[tech] = (calls[tech].loc[m.call_idx], m.to_dict()) if m.matched else None
        tasks.append((sv.to_dict(), sc, None, d, techs, ins))
    ctx = mp.get_context("fork")
    with ctx.Pool(workers, initializer=_init) as pool:
        rows = pool.map(_row, [(t[0], t[1], t[2], t[3], t[4], t[5]) for t in tasks], chunksize=4)
    df = pd.DataFrame(rows)
    for tech in techs:
        for lv in LEVELS:
            df[f"det_{tech}_{lv}"] = truth[f"det_{tech}_{lv}"].values
    df["replicate"] = replicate
    df.to_parquet(outdir / f"{batch}_r{replicate}{suffix}.parquet", index=False)
    fpdf = pd.concat(fp, ignore_index=True) if fp else pd.DataFrame()
    fpdf.to_parquet(outdir / f"{batch}_r{replicate}{suffix}_fp.parquet", index=False)
    print(f"{batch} r{replicate}: {len(df)} SVs; " + "; ".join(f"{t}: {df[f'detected_{t}'].mean():.2f}" for t in techs)
          + f"; FP calls: {len(fpdf)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True)
    ap.add_argument("--replicate", type=int, default=0)
    ap.add_argument("--techs", default="short,hifi")
    ap.add_argument("--tag", default="base")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    run(a.batch, a.replicate, a.techs.split(","), a.tag, a.workers, suffix=a.suffix)
