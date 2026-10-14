"""Read simulation: ART for Illumina, a small in-house simulator for HiFi/ONT.

Everything is seeded; identical inputs give identical bytes.
Limitations (see LEARN.md): long-read error model is uniform (no homopolymer bias), constant
base qualities, no chimeras / sample-prep artefacts, GC bias is a smooth Gaussian efficiency.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

_COMP_TAB = np.zeros(256, np.uint8)
for a, b in zip(b"ACGTNacgtn", b"TGCANtgcan"):
    _COMP_TAB[a] = b
_BASE = np.frombuffer(b"ACGT", np.uint8)
_IDX = np.full(256, 0, np.uint8)
for i, b in enumerate(b"ACGT"):
    _IDX[b] = i


def gc_eff(gc: float, mu: float, sigma: float) -> float:
    return float(np.exp(-((gc - mu) ** 2) / (2 * sigma ** 2)))


def sim_short(window_id: str, haps: list[str], fold_total: float, seed: int, cfg: dict,
              out1: Path, out2: Path) -> dict:
    """Run ART on each haplotype (fold_total/2 each), GC-thin, write plain FASTQ pair."""
    s = cfg["short"]
    tmp = Path(tempfile.mkdtemp(prefix="art_", dir=os.environ.get("DNA_TMP", None)))
    n_pairs = 0
    n_kept = 0
    gc_mu, gc_sig = s["gc_mu"], s["gc_sigma"]
    rng = np.random.default_rng(seed)
    try:
        with open(out1, "w") as f1, open(out2, "w") as f2:
            for h, hap in enumerate(haps):
                fa = tmp / f"h{h}.fa"
                fa.write_text(f">{window_id}_h{h}\n{hap}\n")
                prefix = tmp / f"h{h}_"
                fold = fold_total / len(haps) * s["overshoot"]
                cmd = ["art_illumina", "-ss", s["profile"], "-i", str(fa), "-p", "-l", str(s["read_len"]),
                       "-f", f"{fold:.4f}", "-m", str(s["frag_mean"]), "-s", str(s["frag_sd"]),
                       "-rs", str(seed + h), "-o", str(prefix), "-na"]
                subprocess.run(cmd, check=True, capture_output=True)
                with open(f"{prefix}1.fq") as r1, open(f"{prefix}2.fq") as r2:
                    while True:
                        a = [r1.readline() for _ in range(4)]
                        if not a[0]:
                            break
                        b = [r2.readline() for _ in range(4)]
                        n_pairs += 1
                        seq = a[1] + b[1]
                        gc = (seq.count("G") + seq.count("C")) / max(len(seq) - 2, 1)
                        if rng.random() > gc_eff(gc, gc_mu, gc_sig):
                            continue
                        n_kept += 1
                        # unique read name; strip ART's "/1" "/2" suffix handling for BWA
                        name = f"{window_id}_h{h}_{n_kept}"
                        f1.write(f"@{name}/1\n{a[1]}{a[2]}{a[3]}")
                        f2.write(f"@{name}/2\n{b[1]}{b[2]}{b[3]}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return {"pairs_generated": n_pairs, "pairs_kept": n_kept,
            "bases": n_kept * 2 * s["read_len"]}


def _mutate(seq: np.ndarray, rate: float, fr: tuple, rng) -> np.ndarray:
    n = len(seq)
    k = rng.binomial(n, rate)
    if k == 0:
        return seq
    pos = np.unique(rng.integers(0, n, k))
    p = np.array(fr, float)
    p /= p.sum()
    t = rng.choice(3, size=len(pos), p=p)
    out = seq.copy()
    sp = pos[t == 0]
    if len(sp):
        out[sp] = _BASE[(_IDX[out[sp]] + rng.integers(1, 4, len(sp))) % 4]
    dp = pos[t == 2]
    ip = pos[t == 1]
    if len(ip):
        ins_b = _BASE[rng.integers(0, 4, len(ip))]
    if len(dp):
        keep = np.ones(n, bool)
        keep[dp] = False
        out = out[keep]
        if len(ip):
            ip = ip - np.searchsorted(dp, ip)
    if len(ip):
        out = np.insert(out, ip, ins_b)
    return out


def sim_long(window_id: str, haps: list[str], cov_total: float, seed: int, tech: str, cfg: dict, out: Path) -> dict:
    """Uniform-start long-read simulator. tech in {hifi, ont}."""
    lc = cfg["long"][tech]
    rng = np.random.default_rng(seed)
    n_reads = 0
    bases = 0
    qchar = lc["qchar"]
    with open(out, "w") as fo:
        for h, hap in enumerate(haps):
            arr = np.frombuffer(hap.encode(), np.uint8)
            L = len(arr)
            target = L * cov_total / len(haps)
            got = 0
            while got < target:
                ln = int(np.clip(rng.lognormal(lc["mu"], lc["sigma"]), lc["min_len"], lc["max_len"]))
                if ln >= L:
                    ln = L - 1
                st = int(rng.integers(0, L - ln + 1))
                r = arr[st:st + ln]
                if rng.random() < 0.5:
                    r = _COMP_TAB[r][::-1]
                r = _mutate(r, lc["error_rate"], (lc["sub"], lc["ins"], lc["dele"]), rng)
                n_reads += 1
                got += ln
                bases += len(r)
                s = r.tobytes().decode()
                fo.write(f"@{window_id}_h{h}_{n_reads}\n{s}\n+\n{qchar * len(s)}\n")
    return {"reads": n_reads, "bases": bases}
