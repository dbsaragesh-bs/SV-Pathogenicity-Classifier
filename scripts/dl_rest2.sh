#!/bin/bash
cd ~/dna_data/raw
P=/mnt/d/Projects/DNA/scripts/pdl.sh
bash $P https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v50.annotation.gtf.gz gencode.v50.annotation.gtf.gz 6 > ~/dl_gencode.log 2>&1 &
bash $P https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz rmsk.txt.gz 10 > ~/dl_rmsk.log 2>&1 &
wait
# gnomad after the reference finished
while pgrep -f "pdl.sh .*broad_hg38" >/dev/null; do sleep 5; done
bash $P https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz gnomad.v4.1.sv.sites.vcf.gz 8
curl -sSfLO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
echo ALLDONE2
