#!/bin/bash
#SBATCH --job-name=pca_grouse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=pca_%j.log

module load plink2  # adjust to your cluster's module name/version

VCF="/scratch/gautschi/blackan/GROUSE/output_shotgun/variant_calling/normalized/joint_variant_calling/joint_germline.norm.sorted.vcf.gz"
OUTDIR="/scratch/gautschi/blackan/GROUSE/output_shotgun/pca"
PREFIX="joint_germline"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# --- Step 1: Import VCF to PLINK2 binary format, with standard QC filters ---
# Filters: biallelic SNPs only, MAF >= 0.05, missingness < 10%, HWE filter
# (drop --allow-extra-chr / adjust chr-set if your grouse assembly uses non-standard contig names)
plink2 \
    --vcf "$VCF" \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --snps-only just-acgt \
    --max-alleles 2 \
    --maf 0.05 \
    --geno 0.1 \
    --hwe 1e-6 \
    --make-pgen \
    --out "${PREFIX}.qc" \
    --threads 8

# --- Step 2: LD-prune (critical — unpruned PCA is dominated by LD blocks, not population structure) ---
plink2 \
    --pfile "${PREFIX}.qc" \
    --allow-extra-chr \
    --indep-pairwise 50 10 0.1 \
    --out "${PREFIX}.pruned" \
    --threads 8

# --- Step 3: Run PCA on pruned variant set ---
# PLINK2's --pca uses randomized/truncated SVD by default when N is large,
# which is benchmarked equivalent to exact PCA for structure inference but far faster.
# Use 'approx' explicitly if you have hundreds+ samples; drop it for exact PCA on small N.
plink2 \
    --pfile "${PREFIX}.qc" \
    --allow-extra-chr \
    --extract "${PREFIX}.pruned.prune.in" \
    --pca 10 approx \
    --out "${PREFIX}.pca" \
    --threads 8

echo "Done. Eigenvectors: ${PREFIX}.pca.eigenvec, Eigenvalues: ${PREFIX}.pca.eigenval"
