"""Simulate one pseudo-sample (batch): many far-apart windows pooled into one read set.

python -m src.sim.run_batch --batch B000 [--replicate 0]
Inputs : data/processed/batches.parquet (sv_id -> batch), truthset_*.parquet
Outputs: $DNA_DATA/interim/sim/<batch>_r<rep>/{short_1.fq.gz,short_2.fq.gz,hifi.fq.gz,ont.fq.gz,
         windows.bed, small_variants.tsv, sim_manifest.json}
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from src.sim.make_window import build_window
from src.sim.sim_reads import sim_long, sim_short
from src.utils.common import BATCH_FILE, INTERIM, PROCESSED, RAW, ensure, load_cfg, sha256_file, stable_seed

REF = RAW / "GRCh38_no_alt.fa"


def _one_window(args):
    sv, cfg, rep, outdir, cov_s, cov_l, techs = args
    fa = Fasta(str(REF))
    w = build_window(sv, fa, cfg, zygosity=sv.get("zygosity", cfg["zygosity"]))
    haps = [w.hap1, w.hap2]
    res = {"sv_id": sv["sv_id"], "chrom": w.chrom, "ws": w.ws, "we": w.we,
           "hap1_len": len(w.hap1), "hap2_len": len(w.hap2), "n_small": len(w.small)}
    outdir = Path(outdir)
    s_seed = stable_seed(sv["sv_id"], rep, "short")
    st = sim_short(sv["sv_id"], haps, cov_s, s_seed, cfg, outdir / f"{sv['sv_id']}_s1.fq", outdir / f"{sv['sv_id']}_s2.fq")
    res.update({f"short_{k}": v for k, v in st.items()})
    for tech in techs:
        l_seed = stable_seed(sv["sv_id"], rep, tech)
        lt = sim_long(sv["sv_id"], haps, cov_l, l_seed, tech, cfg, outdir / f"{sv['sv_id']}_{tech}.fq")
        res.update({f"{tech}_{k}": v for k, v in lt.items()})
    res["small_variants"] = [(w.chrom, w.ws + p, r, a, g1, g2) for p, r, a, g1, g2 in w.small]
    return res


def concat_gz(files, dest):
    with gzip.open(dest, "wb", compresslevel=3) as out:
        for f in files:
            with open(f, "rb") as fh:
                shutil.copyfileobj(fh, out, 1 << 22)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True)
    ap.add_argument("--replicate", type=int, default=0)
    ap.add_argument("--techs", default="hifi")
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()

    cfg = load_cfg("simulation")
    batches = pd.read_parquet(BATCH_FILE)
    b = batches[batches.batch == a.batch].reset_index(drop=True)
    assert len(b), f"no SVs in batch {a.batch}"
    outdir = ensure(INTERIM / "sim" / f"{a.batch}_r{a.replicate}")
    work = ensure(outdir / "work")
    rng = np.random.default_rng(stable_seed("cov", a.batch, a.replicate))
    cov_s = float(rng.uniform(*cfg["coverage_total"]))
    cov_l = float(rng.uniform(*cfg["coverage_total"]))
    techs = [t for t in a.techs.split(",") if t]
    t0 = time.time()
    tasks = [(r.to_dict(), cfg, a.replicate, str(work), cov_s, cov_l, techs) for _, r in b.iterrows()]
    with ProcessPoolExecutor(a.workers) as ex:
        results = list(ex.map(_one_window, tasks, chunksize=1))
    ids = [r["sv_id"] for r in results]
    concat_gz([work / f"{i}_s1.fq" for i in ids], outdir / "short_1.fq.gz")
    concat_gz([work / f"{i}_s2.fq" for i in ids], outdir / "short_2.fq.gz")
    for t in techs:
        concat_gz([work / f"{i}_{t}.fq" for i in ids], outdir / f"{t}.fq.gz")
    if not a.keep:
        shutil.rmtree(work)
    with open(outdir / "windows.bed", "w") as fh:
        for r in sorted(results, key=lambda r: (r["chrom"], r["ws"])):
            fh.write(f"{r['chrom']}\t{r['ws']}\t{r['we']}\t{r['sv_id']}\n")
    sm = [(c, p, r, al, g1, g2, w["sv_id"]) for w in results for (c, p, r, al, g1, g2) in w["small_variants"]]
    pd.DataFrame(sm, columns=["chrom", "pos0", "ref", "alt", "g1", "g2", "sv_id"]).to_csv(outdir / "small_variants.tsv", sep="\t", index=False)
    for r in results:
        r.pop("small_variants")
    man = {"batch": a.batch, "replicate": a.replicate, "cov_short": cov_s, "cov_long": cov_l, "techs": techs,
           "seconds": time.time() - t0, "windows": results,
           "sim_cfg_hash": sha256_file(Path(__file__).resolve().parents[2] / "configs" / "simulation.yaml")}
    (outdir / "sim_manifest.json").write_text(json.dumps(man, indent=1))
    print(f"batch {a.batch} r{a.replicate}: {len(results)} windows in {man['seconds']:.0f}s "
          f"(cov short {cov_s:.1f}, long {cov_l:.1f})")


if __name__ == "__main__":
    main()
