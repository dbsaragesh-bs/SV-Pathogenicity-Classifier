"""P3: align a simulated batch to full GRCh38 and call SVs with two callers per technology.

python -m src.calling.run_calling --batch B000 --replicate 0 --techs hifi [--steps align,call]
All tool parameters come from configs/calling.yaml (single source of truth).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

from src.utils.common import INTERIM, RAW, load_cfg

REF = RAW / "GRCh38_no_alt.fa"
MM = os.path.expanduser("~/bin/micromamba")
DELLY = os.path.expanduser("~/mamba/envs/delly/bin/delly")


def sh(cmd: str, log: Path | None = None):
    t = time.time()
    print(f"$ {cmd}", flush=True)
    r = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {cmd}\n{r.stdout}\n{r.stderr}\n")
    if r.returncode != 0:
        raise RuntimeError(f"failed ({r.returncode}): {cmd}\n{r.stderr[-2000:]}")
    return time.time() - t


def align_short(d: Path, cfg, log):
    t = cfg["threads"]
    rg = r"@RG\tID:sim\tSM:sim\tPL:ILLUMINA"
    return sh(f"bwa mem -t {t} -R '{rg}' {REF} {d}/short_1.fq.gz {d}/short_2.fq.gz | "
              f"samtools sort -@ 4 -o {d}/short.bam - && samtools index {d}/short.bam", log)


def align_long(d: Path, tech: str, cfg, log):
    preset = cfg["minimap2_preset"][tech]
    t = cfg["threads"]
    mmi = RAW / f"{preset}.mmi"
    if not mmi.exists():
        sh(f"minimap2 -x {preset} -t {t} -d {mmi} {REF}", log)
    return sh(f"minimap2 -ax {preset} -t {t} -Y -R '@RG\\tID:sim\\tSM:sim' {mmi} {d}/{tech}.fq.gz | "
              f"samtools sort -@ 4 -o {d}/{tech}.bam - && samtools index {d}/{tech}.bam", log)


def call_short(d: Path, cfg, log, skip_manta=False):
    t = cfg["threads"]
    times = {}
    bed = d / "windows.bed"
    sh(f"bgzip -c {bed} > {d}/windows.bed.gz && tabix -f -p bed {d}/windows.bed.gz", log)
    sh(f"rm -rf {d}/manta", log)
    if not skip_manta:
      times["manta"] = sh(
        f"{MM} run -r ~/mamba -n manta configManta.py --bam {d}/short.bam --referenceFasta {REF} "
        f"--runDir {d}/manta --exome --callRegions {d}/windows.bed.gz && "
        f"{MM} run -r ~/mamba -n manta python {d}/manta/runWorkflow.py -m local -j {t}", log)
      sh(f"cp {d}/manta/results/variants/diploidSV.vcf.gz {d}/short.manta.vcf.gz && tabix -f -p vcf {d}/short.manta.vcf.gz", log)
    times["delly"] = sh(
        f"{DELLY} sr -g {REF} -h {t} -o {d}/short.delly.bcf -q {cfg['delly']['min_mapq']} {d}/short.bam && "
        f"bcftools view {d}/short.delly.bcf -Oz -o {d}/short.delly.vcf.gz && tabix -f -p vcf {d}/short.delly.vcf.gz", log)
    return times


def call_long(d: Path, tech: str, cfg, log):
    t = cfg["threads"]
    cs = cfg["cutesv"][tech]
    times = {}
    times["sniffles"] = sh(
        f"sniffles --input {d}/{tech}.bam --vcf {d}/{tech}.sniffles.vcf.gz --reference {REF} --threads {t} "
        f"--minsupport {cfg['sniffles']['minsupport']} --minsvlen {cfg['sniffles']['minsvlen']} "
        f"--mapq {cfg['sniffles']['mapq']} --output-rnames --allow-overwrite", log)
    sh(f"rm -rf {d}/cutesv_wd && mkdir -p {d}/cutesv_wd", log)
    times["cutesv"] = sh(
        f"cuteSV {d}/{tech}.bam {REF} {d}/{tech}.cutesv.vcf {d}/cutesv_wd --threads {t} --genotype "
        f"--min_support {cfg['cutesv']['min_support']} --min_size {cfg['cutesv']['min_size']} --max_size -1 "
        f"--max_cluster_bias_INS {cs['bias_ins']} --diff_ratio_merging_INS {cs['dr_ins']} "
        f"--max_cluster_bias_DEL {cs['bias_del']} --diff_ratio_merging_DEL {cs['dr_del']} --report_readid", log)
    sh(f"bgzip -f {d}/{tech}.cutesv.vcf && tabix -f -p vcf {d}/{tech}.cutesv.vcf.gz", log)
    return times


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True)
    ap.add_argument("--replicate", type=int, default=0)
    ap.add_argument("--techs", default="short,hifi")
    ap.add_argument("--steps", default="align,call")
    ap.add_argument("--rm-fastq", action="store_true")
    ap.add_argument("--skip-manta", action="store_true")
    a = ap.parse_args()
    cfg = load_cfg("calling")
    d = INTERIM / "sim" / f"{a.batch}_r{a.replicate}"
    log = d / "calling.log"
    timing = {}
    for tech in a.techs.split(","):
        if "align" in a.steps:
            timing[f"align_{tech}"] = align_short(d, cfg, log) if tech == "short" else align_long(d, tech, cfg, log)
        if "call" in a.steps:
            timing[f"call_{tech}"] = call_short(d, cfg, log, a.skip_manta) if tech == "short" else call_long(d, tech, cfg, log)
    if a.rm_fastq:
        sh(f"rm -f {d}/*.fq.gz", log)
    prev = {}
    tf = d / "timing.json"
    if tf.exists():
        prev = json.loads(tf.read_text())
    prev.update(timing)
    tf.write_text(json.dumps(prev, indent=1))
    print("done", a.batch, timing)


if __name__ == "__main__":
    main()
