#!/bin/bash
source /mnt/d/Projects/DNA/scripts/wsl_env.sh
cd ~/dna_data/raw
samtools faidx broad_hg38.fasta
CH="$(for i in $(seq 1 22) X Y M; do printf 'chr%s ' $i; done)"
samtools faidx broad_hg38.fasta $CH > GRCh38_no_alt.fa
samtools faidx GRCh38_no_alt.fa
cut -f1,2 GRCh38_no_alt.fa.fai | head -30
echo "ref ready; indexing minimap2 (hifi) then bwa"
minimap2 -x map-hifi -t 14 -d map-hifi.mmi GRCh38_no_alt.fa
( time bwa index GRCh38_no_alt.fa ) > ~/bwa_index.log 2>&1
echo REFDONE
