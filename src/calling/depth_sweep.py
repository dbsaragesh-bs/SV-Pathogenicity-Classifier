"""P8.5 depth sweep: deterministic read-name-hash subsampling of the 30x BAMs, re-calling, re-featurising.

python -m src.calling.depth_sweep --batches B000,B001 --long 5,10,20 --short 15
Dir layout: <batch>_r0_L<d> (long depth d, short at full depth) / <batch>_r0_S<d> (short depth d, long full).
Lower depths are true subsets of higher depths (samtools -s hashes read names).
"""
from __future__ import annotations

import argparse
import json
import os

from src.calling import run_calling as rc
from src.features import build_features as bf
from src.utils.common import INTERIM, load_cfg, stable_seed


def make_variant(batch, rep, kind, depth, cfg):
    src = INTERIM / "sim" / f"{batch}_r{rep}"
    dst = INTERIM / "sim" / f"{batch}_r{rep}_{kind}{depth}"
    dst.mkdir(exist_ok=True)
    man = json.loads((src / "sim_manifest.json").read_text())
    log = dst / "calling.log"
    full = {"L": man["cov_long"], "S": man["cov_short"]}[kind]
    frac = min(1.0, depth / full)
    seed = stable_seed("sweep", batch, kind) % 1000
    for f in ["windows.bed", "sim_manifest.json"]:
        if not (dst / f).exists():
            os.symlink(src / f, dst / f)
    keep_other = {"L": ["short.bam", "short.bam.bai", "short.manta.vcf.gz", "short.manta.vcf.gz.tbi", "short.delly.vcf.gz", "short.delly.vcf.gz.tbi"],
                  "S": ["hifi.bam", "hifi.bam.bai", "hifi.sniffles.vcf.gz", "hifi.sniffles.vcf.gz.tbi", "hifi.cutesv.vcf.gz", "hifi.cutesv.vcf.gz.tbi"]}[kind]
    for f in keep_other:
        if not (dst / f).exists():
            os.symlink(src / f, dst / f)
    tech = "hifi" if kind == "L" else "short"
    bam_in, bam_out = src / f"{tech}.bam", dst / f"{tech}.bam"
    if not bam_out.exists() or not (dst / ("hifi.cutesv.vcf.gz" if kind == "L" else "short.delly.vcf.gz")).exists():
        if frac < 1.0:
            rc.sh(f"samtools view -b -@ 8 -s {seed}.{int(frac * 10000):04d} {bam_in} > {bam_out} && samtools index {bam_out}", log)
        else:
            rc.sh(f"cp {bam_in} {bam_out} && cp {bam_in}.bai {bam_out}.bai", log)
        if kind == "L":
            rc.call_long(dst, "hifi", cfg, log)
        else:
            rc.call_short(dst, cfg, log)
    return f"_{kind}{depth}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batches", required=True)
    ap.add_argument("--long", default="")
    ap.add_argument("--short", default="")
    ap.add_argument("--rep", type=int, default=0)
    ap.add_argument("--s60", action="store_true")
    a = ap.parse_args()
    cfg = load_cfg("calling")
    for b in a.batches.split(","):
        if a.s60:
            suf = short60(b, cfg)
            if not (INTERIM / "feat" / "sweep" / f"{b}_r0{suf}.parquet").exists():
                bf.run(b, a.rep, ["short", "hifi"], "sweep", workers=14, suffix=suf)
            continue
        for kind, ds in (("L", a.long), ("S", a.short)):
            for d in [int(x) for x in ds.split(",") if x]:
                suf = make_variant(b, a.rep, kind, d, cfg)
                if not (INTERIM / "feat" / "sweep" / f"{b}_r{a.rep}{suf}.parquet").exists():
                    bf.run(b, a.rep, ["short", "hifi"], "sweep", workers=14, suffix=suf)



def short60(batch, cfg):
    """60x short arm: merge the r0 30x short BAM with an independent 30x short sim (replicate 1, same haplotypes)."""
    r0 = INTERIM / "sim" / f"{batch}_r0"
    r1 = INTERIM / "sim" / f"{batch}_r1"
    dst = INTERIM / "sim" / f"{batch}_r0_S60"
    dst.mkdir(exist_ok=True)
    log = dst / "calling.log"
    for f in ["windows.bed", "sim_manifest.json", "hifi.bam", "hifi.bam.bai", "hifi.sniffles.vcf.gz", "hifi.sniffles.vcf.gz.tbi",
              "hifi.cutesv.vcf.gz", "hifi.cutesv.vcf.gz.tbi"]:
        if not (dst / f).exists():
            os.symlink(r0 / f, dst / f)
    if not (dst / "short.delly.vcf.gz").exists():
        rc.sh(f"samtools merge -f -@ 8 {dst}/short.bam {r0}/short.bam {r1}/short.bam && samtools index {dst}/short.bam", log)
        rc.call_short(dst, cfg, log)
    return "_S60"


if __name__ == "__main__":
    main()
