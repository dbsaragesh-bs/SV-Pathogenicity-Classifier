#!/bin/bash
# source me: activates the project env in WSL
export MAMBA_ROOT_PREFIX=$HOME/mamba
export PATH=$HOME/mamba/envs/svenv/bin:$HOME/bin:$PATH
export DNA_DATA=$HOME/dna_data
export PYTHONPATH=/mnt/d/Projects/DNA
export DNA_TMP=$HOME/dna_data/tmp
mkdir -p $DNA_TMP
cd /mnt/d/Projects/DNA
