#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
export DNA_BATCHES=/mnt/d/Projects/DNA/data/processed/pilot_batches.parquet
"$@"
