#!/bin/bash
# =============================================================================
# remove_Z_scaffolds.sh -- exclude the Z-linked scaffolds from the joint
# biallelic-SNP VCF, producing the autosomal VCF that downstream scripts
# (run_plink.sh, outflank_drift.sh, fst_temporal.py's wrapper) expect at
# output.subset.biallelic.snps.auto.vcf.gz.
#
# The two scaffold IDs below (NW_026294758.1, NW_026294813.1) are the Z
# scaffolds identified by chicken synteny + genotype hemizygosity in
# find_sex_scaffolds.sh / sex_from_vcf.py (see that script's
# results/sex_scaffolds/SEX_CHROMS.txt output). This reference assembly is
# from a male bird, so there is no W scaffold to also exclude.
#
# Run this once, after the joint VCF exists and before any downstream
# autosomal analysis (PCA, FST, OutFLANK, etc.).
#
# Requirements: bcftools (RCAC biocontainers module, loaded below).
# =============================================================================
set -euo pipefail

module load biocontainers bcftools

VCF_DIR="/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs"
IN_VCF="output.subset.biallelic.snps.vcf.gz"
OUT_VCF="output.subset.biallelic.snps.auto.vcf.gz"
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"

cd "$VCF_DIR"
[[ -s "$IN_VCF" ]] || { echo "ERROR: input VCF not found: ${VCF_DIR}/${IN_VCF}" >&2; exit 1; }

echo "[remove_Z_scaffolds] Excluding Z scaffolds (${Z_SCAFFOLDS}) from ${IN_VCF}"
bcftools view -t "^${Z_SCAFFOLDS}" "$IN_VCF" -Oz -o "$OUT_VCF"
bcftools index -t "$OUT_VCF"

echo "[remove_Z_scaffolds] Wrote ${VCF_DIR}/${OUT_VCF}"
echo "[remove_Z_scaffolds] Records: $(bcftools index -n "$OUT_VCF")"
