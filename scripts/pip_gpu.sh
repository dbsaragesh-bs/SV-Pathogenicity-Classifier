#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
pip install --no-input optuna shap xgboost 2>&1 | tail -3
pip install --no-input torch --index-url https://download.pytorch.org/whl/cu128 2>&1 | tail -3
python - <<'PY'
import torch, xgboost, optuna, shap
print("torch", torch.__version__, "cuda", torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else "")
print("xgboost", xgboost.__version__)
PY
echo PIPDONE
