#!/bin/bash
# Stage 2: short-read align + Manta/Delly, then feature table per batch. Resumable.
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
REP=${1:-0}; TAG=${2:-base}
for B in $(python -c "
import pandas as pd
print(' '.join(sorted(pd.read_parquet('data/processed/batches.parquet').batch.unique())))"); do
  D=$DNA_DATA/interim/sim/${B}_r${REP}
  if [ ! -s $D/short.delly.vcf.gz ]; then
    python -m src.calling.run_calling --batch $B --replicate $REP --techs short --rm-fastq || { echo "SHORT FAIL $B"; continue; }
  fi
  if [ ! -s $DNA_DATA/interim/feat/$TAG/${B}_r${REP}.parquet ]; then
    python -m src.features.build_features --batch $B --replicate $REP --techs short,hifi --tag $TAG --workers 14 || { echo "FEAT FAIL $B"; continue; }
  fi
  echo "STAGE2 OK $B $(date +%T)"
done
echo STAGE2_DONE
