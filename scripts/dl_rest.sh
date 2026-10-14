#!/bin/bash
R=~/dna_data/raw; cd $R
P=/mnt/d/Projects/DNA/scripts/pdl.sh
[ -s hg38.analysisSet.fa.gz ] || bash $P https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz hg38.analysisSet.fa.gz 12
[ -s gencode.v50.annotation.gtf.gz ] || bash $P https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v50.annotation.gtf.gz gencode.v50.annotation.gtf.gz 6
[ -s rmsk.txt.gz ] || bash $P https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz rmsk.txt.gz 6
[ -s gnomad.v4.1.sv.sites.vcf.gz ] || bash $P https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz gnomad.v4.1.sv.sites.vcf.gz 8
[ -s gnomad.v4.1.sv.sites.vcf.gz.tbi ] || curl -sSfLO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
echo ALLDONE
