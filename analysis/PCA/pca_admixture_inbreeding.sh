#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: heterozygosity aggregation + PCA/admixture/inbreeding
# (ANGSD -> pcangsd), whole-genome, all 20 samples.
#
# PCA/admixture here is superseded by run_plink.sh (autosomal, 19-sample
# unrelated -- the one used in the manuscript). This script still owns the
# only source of the pcangsd inbreeding coefficient (--inbreedSamples
# --inbreedSites), a genotype-likelihood-based estimate independent of both
# the bcftools-roh and ROHan fROH pipelines -- not yet re-run autosomal/
# unrelated itself.
#
# USAGE:
#   sbatch pca_admixture_inbreeding.sh
# =============================================================================
#SBATCH --job-name=old_new_pca_inbreed
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

ml biocontainers
ml angsd/0.940
ml pcangsd
# RCAC's xalt accounting hook injects LD_PRELOAD (libxalt_init.so) into
# every command, including containerized ones. singularity forwards it
# into the container by default, and the container's older glibc lacks
# the GLIBC_2.33/2.34 symbols that library needs, so angsd/pcangsd abort
# before running. Blanking it inside the container via these two env vars
# is the reliable fix (a plain `unset LD_PRELOAD` on the host doesn't
# hold — xalt is sticky and re-injects it):
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"
HET_DIR="${PROJECT_DIR}/heterozygosity"

BEAGLE_DIR="${PROJECT_DIR}/beagle"
PCA_DIR="${PROJECT_DIR}/pca"

THREADS=$SLURM_CPUS_PER_TASK
PARALLEL_JOBS=8

mkdir -p logs "$BEAGLE_DIR" "$PCA_DIR"

echo ">>> pca_admixture_inbreeding.sh"
echo ">>> Start time: $(date)"

if [[ ! -f "$FINAL_CRAMLIST" ]]; then
    echo "ERROR: ${FINAL_CRAMLIST} not found. Run 06_downsample_and_finalize.sh first."
    exit 1
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "ERROR: ${REF_FASTA}.fai not found. Run 05_combined_alignment_array.sh's prep step first."
    exit 1
fi

N_SAMPLES=$(wc -l < "$FINAL_CRAMLIST")
MININD=$(awk -v n="$N_SAMPLES" 'BEGIN { printf "%d", (n * 0.75) + 0.5 }')
echo ">>> N samples : ${N_SAMPLES}"
echo ">>> minInd    : ${MININD} (75% of N, matching the original's ~75% stringency)"

# =============================================================================
# STEP 0: Aggregate per-sample heterozygosity results from 08 (each array
# task there wrote its own file to avoid a shared-file race).
# =============================================================================
echo ">>> Step 0: Aggregating heterozygosity results"

HET_SUMMARY="${PROJECT_DIR}/heterozygosity_summary.tsv"
echo -e "sample_id\theterozygosity" > "$HET_SUMMARY"
find "$HET_DIR" -name "*_heterozygosity.txt" -exec cat {} + >> "$HET_SUMMARY" 2>/dev/null || true
echo "  Wrote ${HET_SUMMARY} ($(($(wc -l < "$HET_SUMMARY") - 1)) samples)"

# =============================================================================
# STEP 1: Genotype likelihoods (beagle format), parallelized per chromosome
# Flags match the original beagle.sh exactly (GL model, major/minor, MAF,
# quality, triallelic/SNP filtering),  -minInd is recomputed for N
# =============================================================================
echo ">>> Step 1: ANGSD genotype likelihoods (beagle format)"

CHROM_LIST="${BEAGLE_DIR}/chroms.txt"
cut -f1 "${REF_FASTA}.fai" > "$CHROM_LIST"
N_CHROMS=$(wc -l < "$CHROM_LIST")
echo "  ${N_CHROMS} chromosomes/contigs to process"

BEAGLE_THREADS_PER_JOB=$(( THREADS / PARALLEL_JOBS > 0 ? THREADS / PARALLEL_JOBS : 1 ))

run_beagle_chrom() {
    local chrom="$1"
    local out="${BEAGLE_DIR}/${chrom}"
    if [[ -f "${out}.beagle.gz" ]]; then
        return 0
    fi
    angsd -bam "$FINAL_CRAMLIST" -ref "$REF_FASTA" -r "${chrom}:" \
        -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -minQ 30 \
        -skipTriallelic 1 -SNP_pval 1e-6 -minInd "$MININD" \
        -P "$BEAGLE_THREADS_PER_JOB" -out "$out"
}

export -f angsd
export -f run_beagle_chrom
export FINAL_CRAMLIST REF_FASTA BEAGLE_DIR MININD BEAGLE_THREADS_PER_JOB

xargs -a "$CHROM_LIST" -I{} -P "$PARALLEL_JOBS" bash -c 'run_beagle_chrom "$@"' _ {}

echo ">>> Step 1b: Concatenating per-chromosome beagle files"

FINAL_BEAGLE="${BEAGLE_DIR}/final.beagle.gz"
if [[ ! -f "$FINAL_BEAGLE" ]]; then
    FIRST=1
    > "${BEAGLE_DIR}/final.beagle"
    while IFS= read -r chrom; do
        f="${BEAGLE_DIR}/${chrom}.beagle.gz"
        [[ ! -f "$f" ]] && { echo "  WARNING: missing ${f} — skipping" >&2; continue; }
        if [[ "$FIRST" -eq 1 ]]; then
            zcat "$f" >> "${BEAGLE_DIR}/final.beagle"
            FIRST=0
        else
            zcat "$f" | tail -n +2 >> "${BEAGLE_DIR}/final.beagle"
        fi
    done < "$CHROM_LIST"
    gzip "${BEAGLE_DIR}/final.beagle"
fi

echo "  Final beagle file: ${FINAL_BEAGLE}"

# =============================================================================
# STEP 2: PCA + inbreeding (pcangsd) — flags match pca.sh exactly
# =============================================================================
echo ">>> Step 2: pcangsd"

pcangsd -b "$FINAL_BEAGLE" -o "${PCA_DIR}/final" --threads "$THREADS" --minMaf 0.01 --admix

pcangsd -b "$FINAL_BEAGLE" -o "${PCA_DIR}/final_inbreed" --threads "$THREADS" --minMaf 0.01 \
    --maf_tole 1e-9 --tole 1e-9 --inbreedSamples --inbreedSites \
    --iter 5000 --maf_iter 5000 --inbreed_iter 5000 --inbreed_tole 1e-9

echo "  PCA output      : ${PCA_DIR}/final.cov (+ .admix.Q etc.)"
echo "  Inbreeding output: ${PCA_DIR}/final_inbreed.*"

echo ""
echo ">>> PCA analysis complete."
echo "    Heterozygosity : ${HET_SUMMARY}"
echo "    PCA            : ${PCA_DIR}/final.cov"
echo "    Inbreeding     : ${PCA_DIR}/final_inbreed.*"
echo "    ROH + ROHan    : run bcftools_roh.sh / ROHan\/rohan_array.sh separately"
echo ">>> End time: $(date)"
