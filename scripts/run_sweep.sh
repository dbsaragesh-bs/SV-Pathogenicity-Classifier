#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
nice -n 10 python -m src.calling.depth_sweep --batches B000,B001,B002,B003,B004,B005,B006,B007 --long 5,10,20 --short 15 2>&1 | grep -vE '^\$ |Warn'
echo SWEEP_DONE
