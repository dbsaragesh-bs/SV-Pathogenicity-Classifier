#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
export DNA_BATCHES=/mnt/d/Projects/DNA/data/processed/real_batches.parquet
D=$DNA_DATA/interim/sim/R000_r0
for t in; do
  [ -s $D/$t.raw.bam ] || mv $D/$t.bam $D/$t.raw.bam
  python /mnt/d/Projects/DNA/scripts/filter_bam.py $D/$t.raw.bam $D/$t.bam && samtools index $D/$t.bam
done
rm -rf $D/manta
python -m src.calling.run_calling --batch R000 --replicate 0 --techs short,hifi --steps call --skip-manta 2>&1 | grep -vE '^\$ ' | tail -3
python -m src.features.build_features --batch R000 --replicate 0 --techs short,hifi --tag real --workers 8 2>&1 | grep -v Warn | tail -5
