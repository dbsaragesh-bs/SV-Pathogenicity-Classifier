#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
BS="B000 B001 B002 B003 B004 B005 B006 B007"
for B in $BS; do
  D=$DNA_DATA/interim/sim/${B}_r1
  if [ ! -s $D/short.bam ]; then
    python -m src.sim.run_batch --batch $B --replicate 1 --techs "" --workers 14 && \
    python -m src.calling.run_calling --batch $B --replicate 1 --techs short --steps align --rm-fastq || echo "FAIL $B"
  fi
done
python -m src.calling.depth_sweep --batches B000,B001,B002,B003,B004,B005,B006,B007 --s60 2>&1 | grep -vE '^\$ |Warn'
echo S60_DONE
