import json, hashlib, glob
from pathlib import Path
from src.utils.common import sha256_file, PROCESSED, CONFIGS
m = {"dataset": "dataset_v1", "files": {p.name: sha256_file(p) for p in sorted(PROCESSED.glob("*.parquet")) + [PROCESSED/"truthset_manifest.json"]},
     "configs": {p.name: sha256_file(p) for p in sorted(CONFIGS.glob("*.yaml"))}}
(PROCESSED/"dataset_v1_manifest.json").write_text(json.dumps(m, indent=1))
print(len(m["files"]), "files hashed")
