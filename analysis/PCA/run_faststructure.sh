#!/bin/bash
#SBATCH --job-name=faststructure_auto_unrel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=faststructure_auto_unrel_%j.log
#SBATCH -A fnrdewoody
#SBATCH -p cpu

# =============================================================================
# fastStructure (Raj et al. 2014, v1.0) on the same 19 unrelated birds and
# LD-pruned autosomal SNPs as the PLINK2 PCA (run_plink.sh must run first).
#   K = 2-5, logistic prior; model complexity chosen with chooseK.py.
#   Output prefix "output_results_log" is what R/plot_faststructure.R reads
#   (output_results_log.<K>.meanQ).
# USFWS report, Objective 1, "Population structure" (Figure 6).
#
# NOTE: reconstructed from the reported methods -- the original fastStructure
# commands were not versioned. Check paths/module names against what was run.
# The era-by-cluster Fisher's exact test is done in R (plot_faststructure.R).
# =============================================================================
set -euo pipefail
module --force purge
module load biocontainers plink2 faststructure
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""

PCA_DIR="/scratch/gautschi/blackan/GROUSE/old_vs_new/pca"
PREFIX="joint_germline_auto_unrel"          # from run_plink.sh
OUTDIR="/scratch/gautschi/blackan/GROUSE/old_vs_new/faststructure"
KMIN=2; KMAX=5
THREADS="${SLURM_CPUS_PER_TASK:-4}"

for f in "${PCA_DIR}/${PREFIX}.qc.pgen" "${PCA_DIR}/${PREFIX}.pruned.prune.in"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f (run run_plink.sh first)"; exit 1; }
done
mkdir -p "$OUTDIR"; cd "$OUTDIR"

# fastStructure reads PLINK1 binary (bed/bim/fam): export the pruned set
plink2 --pfile "${PCA_DIR}/${PREFIX}.qc" --allow-extra-chr \
       --extract "${PCA_DIR}/${PREFIX}.pruned.prune.in" \
       --make-bed --out "${PREFIX}.pruned" --threads "$THREADS"

for K in $(seq $KMIN $KMAX); do
    structure.py -K "$K" --input="${PREFIX}.pruned" \
                 --output=output_results_log --prior=logistic --full --seed=100
done

chooseK.py --input=output_results_log | tee chooseK.txt
cut -d' ' -f1-2 "${PREFIX}.pruned.fam" > sample_order.txt
echo "Done: ${OUTDIR}/output_results_log.<K>.meanQ (sample order: sample_order.txt)"
