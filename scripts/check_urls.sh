#!/bin/bash
while read u; do
  echo "== $u"
  curl -sIL --max-time 30 "$u" | grep -iE "^HTTP|content-length" | tail -2
done <<URLS
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v49.annotation.gtf.gz
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG002-T2TQ100-V1.1_stvar.vcf.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
URLS
