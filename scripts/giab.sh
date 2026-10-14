B=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son
for p in v5.0q v5.0q/GRCh38 v5.0q/GRCh38/SVs CMRG_v1.00/GRCh38; do echo "== $p"; curl -sL $B/$p/ | grep -o 'href="[^"?][^"]*"' | head -15; done
D=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son
for p in PacBio_HiFi-Revio_20231031 PacBio_HiFi-Revio_20231031/GRCh38 NIST_Illumina_2x250bps NIST_Illumina_2x250bps/novoalign_bams; do echo "== $p"; curl -sL $D/$p/ | grep -o 'href="[^"?][^"]*"' | head -12; done
