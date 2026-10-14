source /mnt/d/Projects/DNA/scripts/wsl_env.sh
D=$HOME/dna_data/interim/sim/R000_r0
head -2 $D/windows.bed
samtools view -H $D/short.bam | grep -c "^@SQ"
samtools view -H $D/short.bam | grep -v "^@SQ" | head -5 | cut -c1-200
samtools view $D/short.bam chr5:210786-312098 | awk '$7!="=" && $7!="*"' | head -3 | cut -f1-9
samtools view -f 1 -F 12 $D/short.bam chr5:210786-312098 | awk '{print $2}' | sort | uniq -c | sort -k1nr | head -8 | tr "\n" " "
samtools view $D/short.bam chr5:210786-312098 | awk '$9>1000000000 || $9<-1000000000' | head -3 | cut -f1-9
