#!/bin/bash
set -u
R=~/dna_data/raw; mkdir -p $R; cd $R
GC=$(curl -s https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/ | grep -oE 'gencode\.v[0-9]+\.annotation\.gtf\.gz' | head -1)
echo "gencode file: $GC"
dl() { f=$(basename "$1"); [ -s "$f" ] && { echo "have $f"; return; }; curl -fL --retry 5 --retry-delay 5 -C - -o "$f.part" "$1" && mv "$f.part" "$f" && echo "got $f"; }
dl https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz
dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
dl https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
dl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
dl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
[ -n "$GC" ] && dl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/$GC
dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
dl https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
dl https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi
echo ALLDONE
