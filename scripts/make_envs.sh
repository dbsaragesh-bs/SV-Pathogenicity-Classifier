#!/bin/bash
cd ~
./bin/micromamba create -y -r ~/mamba -n manta -f /mnt/d/Projects/DNA/envs/manta.yml > ~/env_manta.log 2>&1 &
./bin/micromamba create -y -r ~/mamba -n svenv -f /mnt/d/Projects/DNA/envs/environment.yml > ~/env_create.log 2>&1
wait
echo ENVS_DONE
