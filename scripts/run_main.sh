#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
python -u -m src.modelling.run_experiment --features data/processed/features_all.parquet --cohort both --arms A,B,O,OA --models xgb,lr,rf,mlp --trials 20 --k 5 --out reports/ablation/main --tag main > ~/main_run.log 2>&1
echo MAIN_DONE >> ~/main_run.log
