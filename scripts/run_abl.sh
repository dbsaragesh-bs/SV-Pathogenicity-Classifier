#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
python -u -m src.modelling.ablation --which ladder,logo,controls --trials 12 > ~/abl_run.log 2>&1
echo ABL_DONE >> ~/abl_run.log
