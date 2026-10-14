#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
python -u -m src.validation.realdata_p9 prepare > ~/p9a.log 2>&1
echo P9A_DONE >> ~/p9a.log
