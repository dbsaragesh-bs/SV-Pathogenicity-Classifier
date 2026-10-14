#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
python -u -m src.modelling.long_only --trials 20 > /home/dbs/long_only.log 2>&1
