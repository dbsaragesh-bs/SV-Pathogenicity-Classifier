source /mnt/d/Projects/DNA/scripts/wsl_env.sh
cd ~/dna_data/raw
bcftools view -h gnomad.v4.1.sv.sites.vcf.gz | grep -E "^##INFO=<ID=(SVTYPE|SVLEN|END|AF|AC|AN|N_BI_GENOS|ALGORITHMS|EVIDENCE|CPX|FILTER|MULTIALLELIC|PAR|FREQ)" | cut -c1-160
bcftools view -h gnomad.v4.1.sv.sites.vcf.gz | grep -E "^##FILTER" | head
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%FILTER\t%INFO/AF\t%INFO/ALGORITHMS\n' gnomad.v4.1.sv.sites.vcf.gz | head -5
bcftools query -f '%INFO/SVTYPE\t%FILTER\n' gnomad.v4.1.sv.sites.vcf.gz | sort | uniq -c | sort -k1nr | head -20
