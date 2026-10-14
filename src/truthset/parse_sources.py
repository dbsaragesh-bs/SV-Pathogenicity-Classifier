"""Parse raw sources into one harmonised schema (P1.1).

Schema: source, source_id, chrom, start, end, svtype, svlen, label, review, stars,
        coord_exact, ci_width, benign_basis
Coordinates are returned 0-based half-open [start, end) internally (VCF/ClinVar 1-based
inclusive Start..Stop  ->  start = Start-1, end = Stop).
"""
from __future__ import annotations

import gzip
import re
import subprocess

import numpy as np
import pandas as pd

TYPE_MAP = {
    "Deletion": "DEL",
    "copy number loss": "DEL",
    "Duplication": "DUP",
    "copy number gain": "DUP",
    "Tandem duplication": "DUP",
}
STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
}


def parse_clinvar(path, cfg) -> pd.DataFrame:
    usecols = ["Type", "Name", "ClinicalSignificance", "Assembly", "Chromosome", "Start", "Stop",
               "ReviewStatus", "VariationID", "Origin", "GeneSymbol"]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False, dtype=str)
    n0 = len(df)
    df = df[df.Assembly == "GRCh38"]
    df = df[df.Type.isin(TYPE_MAP)]
    df["chrom"] = "chr" + df.Chromosome.astype(str)
    df = df[df.chrom.isin(cfg["chroms"])]
    df["start"] = pd.to_numeric(df.Start, errors="coerce") - 1
    df["end"] = pd.to_numeric(df.Stop, errors="coerce")
    df = df.dropna(subset=["start", "end"])
    df["start"] = df.start.astype(int)
    df["end"] = df.end.astype(int)
    df["svtype"] = df.Type.map(TYPE_MAP)
    df["svlen"] = df.end - df.start
    sig = df.ClinicalSignificance
    df["label"] = np.where(sig.isin(cfg["clinvar_pathogenic"]), 1,
                           np.where(sig.isin(cfg["clinvar_benign"]), 0, -1))
    df["stars"] = df.ReviewStatus.map(STARS).fillna(0).astype(int)
    # precise = HGVS name has no uncertain-bound notation "(?_" / "_?)" / "(" ranges
    df["coord_exact"] = ~df.Name.fillna("").str.contains(r"\?|\(\d", regex=True)
    out = pd.DataFrame({
        "source": "clinvar",
        "source_id": "VCV" + df.VariationID.astype(str),
        "variation_id": df.VariationID.astype(str),
        "chrom": df.chrom, "start": df.start, "end": df.end, "svtype": df.svtype,
        "svlen": df.svlen, "label": df.label, "stars": df.stars, "review": df.ReviewStatus,
        "coord_exact": df.coord_exact, "origin": df.Origin,
        "benign_basis": np.where(df.label == 0, "clinvar", ""),
    })
    out.attrs["n_raw_rows"] = n0
    return out.reset_index(drop=True)


def parse_dbvar_ci(path) -> pd.DataFrame:
    """VariationID -> (imprecise, ci_start_width, ci_end_width) from the ClinVar-in-dbVar VCF."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
            m = re.findall(r"VCV0*(\d+)", info.get("CLNACC", ""))
            if not m:
                continue

            def w(key):
                v = info.get(key)
                if not v:
                    return np.nan
                a, b = (x if x != "." else None for x in (v.split(",") + [None])[:2])
                try:
                    return abs(int(a or 0)) + abs(int(b or 0))
                except ValueError:
                    return np.nan

            rows.append((m[0], "IMPRECISE" in info, w("CIPOS"), w("CIEND")))
    d = pd.DataFrame(rows, columns=["variation_id", "imprecise", "ci_pos", "ci_end"])
    d["ci_width"] = d[["ci_pos", "ci_end"]].max(axis=1)
    return d.groupby("variation_id", as_index=False).agg(
        imprecise=("imprecise", "max"), ci_width=("ci_width", "max"))


def parse_gnomad_sv(vcf, cfg, min_af: float) -> pd.DataFrame:
    """Common (AF>=min_af) PASS biallelic-ish DEL/DUP records from gnomAD-SV via bcftools."""
    expr = (f'FILTER="PASS" && (INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP") && INFO/AF>={min_af}')
    cmd = ["bcftools", "query", "-i", expr, "-f",
           "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\t%INFO/AF\t%INFO/N_BI_GENOS\t%INFO/ALGORITHMS\n",
           str(vcf)]
    txt = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    from io import StringIO
    d = pd.read_csv(StringIO(txt), sep="\t", header=None,
                    names=["chrom", "pos", "end", "svtype", "svlen", "id", "af", "n_geno", "algs"])
    d = d[d.chrom.isin(cfg["chroms"])].copy()
    d["start"] = d.pos.astype(int) - 1
    d["end"] = d.end.astype(int)
    d["svlen"] = d.end - d.start
    out = pd.DataFrame({
        "source": "gnomad_sv", "source_id": d.id, "variation_id": "",
        "chrom": d.chrom, "start": d.start, "end": d.end, "svtype": d.svtype,
        "svlen": d.svlen, "label": 0, "stars": 0, "review": "gnomad_pass_af",
        "coord_exact": True, "origin": "population", "benign_basis": "af",
    })
    return out.reset_index(drop=True)
