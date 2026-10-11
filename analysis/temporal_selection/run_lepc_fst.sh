#!/bin/bash
# ============================================================================
# SLURM wrapper: Past vs Present LEPC FST (global + sliding windows)
#   Hudson FST (Bhatia et al. 2013 ratio of averages) + Weir & Cockerham,
#   block-jackknife CIs, sliding windows, Manhattan plot. Uses scikit-allel.
#
# One-time environment setup (login node):
#   module load anaconda
#   conda create -n lepc_py -c conda-forge python=3.11 scikit-allel numpy pandas matplotlib
#
# Input: the biallelic SNP VCF made by run_lepc_outflank_temporal.sh.
# ============================================================================
#SBATCH -J lepc_fst_temporal
#SBATCH -o lepc_fst_temporal_%j.out
#SBATCH -e lepc_fst_temporal_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -A dewoody
#SBATCH -p cpu

set -euo pipefail

module --force purge
module load anaconda
set +u
source activate lepc_py
set -u

cd "$SLURM_SUBMIT_DIR"

# --- Edit these paths/params for your run -----------------------------------
VCF="/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs/output.subset.biallelic.snps.auto.vcf.gz"
POPMAP="/scratch/gautschi/blackan/GROUSE/old_vs_new/popmap_unrelated.txt"           # <sample_id> <Past|Present>, no header
OUTPREFIX="results/fst/old_vs_new_unrelated"
WINDOW=50000                   # window size (bp)
STEP=10000                     # step (bp)
MIN_SNPS=50                    # min SNPs for a window FST
MIN_CALLED=5                   # min genotyped individuals per period per site
JK_BLOCKSIZE=5000000           # jackknife block size (bp); matches the Ne CI
GENERATIONS=5                  # for the drift-implied Ne cross-check
# -----------------------------------------------------------------------------

for f in "$VCF" "$POPMAP" fst_temporal.py; do
    [[ -s "$f" ]] || { echo "ERROR: required file not found or empty: $f"; exit 1; }
done
python -c "import allel" 2>/dev/null || { echo "ERROR: scikit-allel not importable (conda env?)"; exit 1; }

echo "== Job $SLURM_JOB_ID on $(hostname) at $(date) =="
python fst_temporal.py \
    --vcf "$VCF" \
    --popmap "$POPMAP" \
    --outprefix "$OUTPREFIX" \
    --window "$WINDOW" \
    --step "$STEP" \
    --min_snps "$MIN_SNPS" \
    --mincalled "$MIN_CALLED" \
    --jk_blocksize "$JK_BLOCKSIZE" \
    --generations "$GENERATIONS"
echo "== Finished at $(date) =="
