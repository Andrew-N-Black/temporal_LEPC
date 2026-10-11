#!/bin/bash
#SBATCH --job-name=pca_grouse_auto_unrel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=pca_auto_unrel_%j.log
#SBATCH -A dewoody
#SBATCH -p cpu

# =============================================================================
# PLINK2 PCA on autosomes only, without normal_F10.
#   - Z scaffolds are excluded here with --not-chr, whether or not the input
#     VCF was already filtered (so the result doesn't depend on the file name).
#   - normal_F10 is removed (first-degree relative of normal_F21; F21 kept
#     for lower missingness) BEFORE the MAF/missingness/HWE filters, so those
#     are computed on the 19 analysed birds.
#   - Outputs use a new prefix so earlier PCA results are not overwritten.
# =============================================================================
set -euo pipefail
module --force purge
module load biocontainers plink2
# RCAC's xalt injects LD_PRELOAD into containerized tools; a plain
# `unset LD_PRELOAD` doesn't hold. Blank it inside the container instead.
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""

VCF="/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs/output.subset.biallelic.snps.auto.vcf.gz"
OUTDIR="/scratch/gautschi/blackan/GROUSE/old_vs_new/pca"
PREFIX="joint_germline_auto_unrel"
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"   # Z (chicken synteny); male reference, no W
DROP_SAMPLES="normal_F10"                     # space-separated if more than one
THREADS="${SLURM_CPUS_PER_TASK:-8}"

[[ -s "$VCF" ]] || { echo "ERROR: VCF not found: $VCF"; exit 1; }
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Samples to remove (PLINK2 sample-ID file; VCF import sets IID = VCF sample name)
REMOVE="${PREFIX}.remove.txt"
{ echo "#IID"; for s in $DROP_SAMPLES; do echo "$s"; done; } > "$REMOVE"

# --- Step 1: Import VCF, drop Z and F10, standard QC filters -----------------
# --hwe ... keep-fewhet: remove only heterozygote EXCESS (paralog/mapping
# artifacts). Past and Present are pooled here, so real frequency differences
# between periods produce a heterozygote DEFICIT; a two-sided HWE filter
# would preferentially remove exactly the loci that separate the groups.
plink2 \
    --vcf "$VCF" \
    --allow-extra-chr \
    --not-chr "$Z_SCAFFOLDS" \
    --remove "$REMOVE" \
    --set-missing-var-ids @:# \
    --snps-only just-acgt \
    --max-alleles 2 \
    --maf 0.05 \
    --geno 0.1 \
    --hwe 1e-6 keep-fewhet \
    --make-pgen \
    --out "${PREFIX}.qc" \
    --threads "$THREADS"

# Sanity checks: F10 gone, 19 samples, no Z variants
N_SAMPLES=$(grep -vc '^#' "${PREFIX}.qc.psam")
for s in $DROP_SAMPLES; do
    ! awk -v s="$s" '!/^#/ && $1==s {f=1} END {exit !f}' "${PREFIX}.qc.psam" \
        || { echo "ERROR: $s still present"; exit 1; }
done
N_Z=$(awk -v z="$Z_SCAFFOLDS" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1}
      !/^#/ && ($1 in Z) {c++} END {print c+0}' "${PREFIX}.qc.pvar")
echo "Samples after removal: ${N_SAMPLES}; Z variants remaining: ${N_Z}; variants: $(grep -vc '^#' "${PREFIX}.qc.pvar")"
[[ "$N_Z" -eq 0 ]] || { echo "ERROR: Z variants remain"; exit 1; }

# --- Step 2: LD-prune (unpruned PCA is dominated by LD blocks) --------------
plink2 \
    --pfile "${PREFIX}.qc" \
    --allow-extra-chr \
    --indep-pairwise 50 10 0.1 \
    --out "${PREFIX}.pruned" \
    --threads "$THREADS"

# --- Step 3a: Allele frequencies on the pruned set ---------------------------
# n < 50, so PLINK2 needs frequencies supplied explicitly via --read-freq.
plink2 \
    --pfile "${PREFIX}.qc" \
    --allow-extra-chr \
    --extract "${PREFIX}.pruned.prune.in" \
    --freq \
    --out "${PREFIX}.freq" \
    --threads "$THREADS"

# --- Step 3b: PCA on the pruned set (exact PCA at this sample size) ---------
plink2 \
    --pfile "${PREFIX}.qc" \
    --allow-extra-chr \
    --extract "${PREFIX}.pruned.prune.in" \
    --read-freq "${PREFIX}.freq.afreq" \
    --pca 10 \
    --out "${PREFIX}.pca" \
    --threads "$THREADS"

echo "Pruned variants used: $(wc -l < "${PREFIX}.pruned.prune.in")"
echo "Done. Eigenvectors: ${OUTDIR}/${PREFIX}.pca.eigenvec, Eigenvalues: ${OUTDIR}/${PREFIX}.pca.eigenval"
