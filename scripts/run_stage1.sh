#!/bin/bash
# Stage 1: simulate reads for every batch, align + call long reads. Resumable.
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
REP=${1:-0}
python - > /tmp/batch_list.txt <<'PY'
import pandas as pd
b = pd.read_parquet("data/processed/batches.parquet")
print("\n".join(sorted(b.batch.unique())))
PY
for B in $(cat /tmp/batch_list.txt); do
  D=$DNA_DATA/interim/sim/${B}_r${REP}
  if [ ! -s $D/sim_manifest.json ]; then
    python -m src.sim.run_batch --batch $B --replicate $REP --techs hifi --workers 14 || { echo "SIM FAIL $B"; continue; }
  fi
  if [ ! -s $D/hifi.cutesv.vcf.gz ]; then
    python -m src.calling.run_calling --batch $B --replicate $REP --techs hifi || { echo "LONG FAIL $B"; continue; }
  fi
  echo "STAGE1 OK $B $(date +%T)"
done
echo STAGE1_DONE
