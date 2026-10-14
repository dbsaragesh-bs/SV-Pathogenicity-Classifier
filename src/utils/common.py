"""Shared paths, config loading, deterministic seeding."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("DNA_DATA", str(Path.home() / "dna_data")))
RAW = DATA / "raw"
INTERIM = DATA / "interim"
PROCESSED = REPO / "data" / "processed"  # small, committed-size artefacts
CONFIGS = REPO / "configs"
REPORTS = REPO / "reports"

BATCH_FILE = Path(os.environ.get("DNA_BATCHES", str(PROCESSED / "batches.parquet")))
CANON_CHROMS = [f"chr{i}" for i in range(1, 23)]


def load_cfg(name: str) -> dict:
    with open(CONFIGS / f"{name}.yaml") as fh:
        return yaml.safe_load(fh)


def stable_seed(*parts, mod: int = 2**31 - 1) -> int:
    """Deterministic seed from arbitrary parts (never Python's salted hash())."""
    h = hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()
    return int(h[:15], 16) % mod


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while b := fh.read(chunk):
            h.update(b)
    return h.hexdigest()


def ensure(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p
