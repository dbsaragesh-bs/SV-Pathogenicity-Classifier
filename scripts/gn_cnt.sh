source /mnt/d/Projects/DNA/scripts/wsl_env.sh
cd ~/dna_data/raw
for af in 0.01 0.005 0.001; do
  echo "AF>=$af"
  bcftools query -i "FILTER=\"PASS\" && (INFO/SVTYPE=\"DEL\"||INFO/SVTYPE=\"DUP\") && INFO/AF>=$af && INFO/SVLEN>=50000 && INFO/SVLEN<=1000000" -f '%INFO/SVTYPE\t%INFO/SVLEN\n' gnomad.v4.1.sv.sites.vcf.gz | awk '{b=($2>=300000)?"300k-1M":($2>=100000?"100-300k":"50-100k"); c[$1" "b]++} END{for(k in c) print k, c[k]}' | sort
done
