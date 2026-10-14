"""P9: real-data realism check on GIAB HG002 (GRCh38).

Steps (all streamed from the GIAB FTP, only the windows around chosen SVs are fetched):
  prepare : choose HG002 benchmark DEL/DUP (v5.0q GRCh38 stvar), fetch windows from the Revio HiFi and
            Illumina 2x250 BAMs, subsample to ~30x, write a pseudo-sample dir R000_r0
  (then)  : python -m src.calling.run_calling --steps call ... ; python -m src.features.build_features
  compare : feature-distribution shift, detection, adversarial validation, model predictions on real SVs

DNA_BATCHES must point at data/processed/real_batches.parquet for the calling/feature steps.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from src.utils.common import INTERIM, PROCESSED, RAW, REPORTS, ensure, load_cfg

GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"
VCF = f"{GIAB}/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz"
HIFI = f"{GIAB}/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam"
ILMN = f"{GIAB}/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
DIR = INTERIM / "sim" / "R000_r0"
FLANK = 50_000


def sh(c):
    r = subprocess.run(["bash", "-o", "pipefail", "-c", c], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(c + "\n" + r.stderr[-800:])
    return r.stdout


def choose(n_per_type=(130, 30), seed=3):
    ensure(INTERIM / "real")
    vcf = INTERIM / "real" / "stvar.vcf.gz"
    if not vcf.exists():
        sh(f"curl -sSfL {VCF} -o {vcf} && curl -sSfL {VCF}.tbi -o {vcf}.tbi && curl -sSfL "
           f"{VCF.replace('.vcf.gz', '.benchmark.bed')} -o {INTERIM / 'real' / 'stvar.bed'}")
    txt = sh(f"bcftools query -f '%CHROM\\t%POS\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%ID\\n' {vcf}")
    from io import StringIO
    d = pd.read_csv(StringIO(txt), sep="\t", header=None, names=["chrom", "pos", "svtype", "svlen", "id"])
    d["svlen"] = pd.to_numeric(d.svlen, errors="coerce").abs()
    d = d.dropna(subset=["svlen"])
    d = d[d.chrom.isin([f"chr{i}" for i in range(1, 23)]) & d.svtype.isin(["DEL", "DUP"]) & d.svlen.between(50, 1_000_000)]
    d["start"] = d.pos  # VCF POS is the base before the event; |error| of 1 bp is irrelevant at this scale
    d["svlen"] = d.svlen.astype(int)
    d["end"] = (d.start + d.svlen).astype(int)
    rng = np.random.default_rng(seed)
    out = []
    for t, n in zip(("DEL", "DUP"), n_per_type):
        g = d[d.svtype == t]
        out.append(g.sample(min(n, len(g)), random_state=int(rng.integers(1e6))))
    d = pd.concat(out)
    # keep windows >= 1 Mb apart and non-overlapping
    d = d.sort_values(["chrom", "start"])
    keep, last = [], {}
    for i, r in d.iterrows():
        if r.chrom in last and r.start - FLANK < last[r.chrom] + 1_000_000:
            continue
        keep.append(i)
        last[r.chrom] = r.end + FLANK
    d = d.loc[keep].reset_index(drop=True)
    d["sv_id"] = [f"RV{i:04d}" for i in range(len(d))]
    d["label"] = 0
    d["group_id"] = range(len(d))
    d["source"] = "giab_hg002"
    d["stars"] = 0
    d["coord_exact"] = True
    d["benign_basis"] = "real"
    d["zygosity"] = "unknown"
    d["batch"] = "R000"
    return d


def fetch(d: pd.DataFrame, cov_target=30.0):
    ensure(DIR)
    ref_chroms = [l.split("\t")[0] for l in open(RAW / "GRCh38_no_alt.fa.fai")]
    regions = [f"{r.chrom}:{max(r.start - FLANK, 1)}-{r.end + FLANK}" for r in d.itertuples()]
    with open(DIR / "windows.bed", "w") as fh:
        for r in d.itertuples():
            fh.write(f"{r.chrom}\t{max(r.start - FLANK, 0)}\t{r.end + FLANK}\t{r.sv_id}\n")
    tmpd = INTERIM / "real" / "parts"
    ensure(tmpd)
    for name, url, cov_full in (("hifi", HIFI, 48.0), ("short", ILMN, None)):
        if (DIR / f"{name}.bam.bai").exists():
            continue
        frac_cmd = ""
        if cov_full is None:  # estimate Illumina coverage from the first window
            reg = regions[0]
            n = int(sh(f"samtools view -c -F 3844 {url} {reg}").strip())
            L = int(reg.split("-")[1]) - int(reg.split(":")[1].split("-")[0])
            cov_full = n * 250 / L
            print("estimated Illumina coverage", cov_full)
        frac = min(1.0, cov_target / cov_full)
        sam = f"-s 7.{int(frac * 10000):04d}" if frac < 1 else ""

        def one(i):
            out = tmpd / f"{name}_{i}.sam"
            for attempt in range(5):
                try:
                    sh(f"samtools view {sam} -h {url} {regions[i]} > {out}")
                    return out
                except RuntimeError:
                    if attempt == 4:
                        raise
            return out

        with ThreadPoolExecutor(4) as ex:
            parts = list(ex.map(one, range(len(regions))))
        hdr = [l for l in open(parts[0]) if l.startswith("@")]
        hdr = [l for l in hdr if not l.startswith("@SQ") or l.split("\t")[1][3:] in ref_chroms]
        # order @SQ lines like the reference
        sq = {l.split("\t")[1][3:]: l for l in hdr if l.startswith("@SQ")}
        other = [l for l in hdr if not l.startswith("@SQ")]
        with open(tmpd / f"{name}_all.sam", "w") as fo:
            fo.writelines([l for l in other if l.startswith("@HD")])
            fo.writelines([sq[c] for c in ref_chroms if c in sq])
            fo.writelines([l for l in other if not l.startswith("@HD")])
            for p in parts:
                with open(p) as fh:
                    fo.writelines(l for l in fh if not l.startswith("@"))
                p.unlink()
        sh(f"samtools sort -@ 8 -o {DIR}/{name}.bam {tmpd}/{name}_all.sam && samtools index {DIR}/{name}.bam && rm {tmpd}/{name}_all.sam")
        print("fetched", name)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare"])
    a = ap.parse_args()
    d = choose()
    d.to_parquet(PROCESSED / "real_batches.parquet", index=False)
    print(d.svtype.value_counts().to_dict(), "SVs;", d.svlen.describe().to_dict())
    fetch(d)
