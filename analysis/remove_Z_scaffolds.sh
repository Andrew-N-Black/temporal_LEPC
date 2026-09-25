module load biocontainers bcftools
cd /scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs
bcftools view -t ^NW_026294758.1,NW_026294813.1 \
    output.subset.biallelic.snps.vcf.gz -Oz -o output.subset.biallelic.snps.auto.vcf.gz
bcftools index -t output.subset.biallelic.snps.auto.vcf.gz
