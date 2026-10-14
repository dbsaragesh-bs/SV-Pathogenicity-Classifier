source /mnt/d/Projects/DNA/scripts/wsl_env.sh
for t in samtools bcftools bedtools minimap2 bwa art_illumina delly sniffles cuteSV SURVIVOR whatshap truvari bgzip tabix; do
  if which $t >/dev/null; then echo "$t ok"; else echo "$t MISSING"; fi
done
ls $HOME/mamba/envs/svenv/lib | grep -i boost | head
~/bin/micromamba run -r ~/mamba -n manta configManta.py --version
