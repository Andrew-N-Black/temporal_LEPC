#!/bin/bash
# ============================================================================
# SLURM wrapper: old vs new LEPC temporal analysis
#   OutFLANK outlier scan + temporal Fk/Ne drift null + site-class comparison
#
# One-time environment setup (run on a login node before submitting):
#   module load anaconda
#   conda create -n lepc_rstats -c conda-forge -c bioconda \
#       r-base r-optparse r-dplyr r-ggplot2 r-vcfr bioconductor-qvalue r-devtools
#   conda activate lepc_rstats
# Rscript -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'
# Rscript -e 'remotes::install_github("whitlock/OutFLANK")'
#
# ============================================================================
#SBATCH -J lepc_outflank_temporal
#SBATCH -o lepc_outflank_temporal_%j.out
#SBATCH -e lepc_outflank_temporal_%j.err
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH -A dewoody
#SBATCH -p cpu

set -euo pipefail

module purge --force
module load biocontainers anaconda
source activate lepc_rstats

cd "$SLURM_SUBMIT_DIR"

# --- Edit these paths/params for your run -----------------------------------
VCF="/scratch/gautschi/blackan/GROUSE/output_shotgun/variant_calling/normalized/joint_variant_calling/joint_germline.norm.sorted.vcf.gz"
POPMAP="popmap.txt"          # <sample_id> <old|new>, no header
SITECLASS="siteclass.txt"    # <CHROM> <POS> <deleterious|neutral>, no header
GENERATIONS=5                # elapsed generations between old and new samples
OUTPREFIX="results/old_vs_new_temporal"
# -----------------------------------------------------------------------------

mkdir -p results

Rscript outflank_temporal_dNe_siteclass.R \
    --vcf "$VCF" \
    --popmap "$POPMAP" \
    --siteclass "$SITECLASS" \
    --generations "$GENERATIONS" \
    --nsim 2000 \
    --freqbins 20 \
    --qthresh 0.05 \
    --outprefix "$OUTPREFIX"
