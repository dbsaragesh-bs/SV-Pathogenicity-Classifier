#!/bin/bash
cd ~/dna_data/raw
rm -f hg38.analysisSet* GCA*.part
P=/mnt/d/Projects/DNA/scripts/pdl.sh
[ -s broad_hg38.fasta ] || bash $P https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta broad_hg38.fasta 8 &
[ -s gencode.v50.annotation.gtf.gz ] || bash $P https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v50.annotation.gtf.gz gencode.v50.annotation.gtf.gz 6 &
[ -s rmsk.txt.gz ] || bash $P https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz rmsk.txt.gz 10 &
wait
[ -s gnomad.v4.1.sv.sites.vcf.gz ] || bash $P https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz gnomad.v4.1.sv.sites.vcf.gz 8
[ -s gnomad.v4.1.sv.sites.vcf.gz.tbi ] || curl -sSfLO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
echo ALLDONE3
