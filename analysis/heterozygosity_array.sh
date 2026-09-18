#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: PER-SAMPLE HETEROZYGOSITY (ANGSD + realSFS)
# Step 07 — requires 06_downsample_and_finalize.sh to have completed
# (reads final_cramlist.txt: OLD-cohort CRAMs as-is + depth-matched
# NEW-cohort CRAMs). One array task per sample.
#
# Replicates the original heterozygosity.sh exactly (same ANGSD/realSFS
# flags). Input is CRAM, not BAM — ANGSD reads
# CRAM natively via its htslib backend as long as -ref is provided (already
# required here for the SAF/ancestral-state calculation anyway).
#
# Individual heterozygosity = SFS[1] / (SFS[0] + SFS[1]) from the folded,
# 2-category single-sample SFS — the standard ANGSD single-sample
# heterozygosity estimate.
#
# USAGE:
#    N=$(wc -l < final_cramlist.txt)
#   sbatch --array=0-$((N-1))%20 07_heterozygosity_array.sh
# =============================================================================
#SBATCH --job-name=lepc_het
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
#SBATCH -t 3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

ml biocontainers
ml angsd/0.940
module unload xalt
unset LD_PRELOAD
# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"
HET_DIR="${PROJECT_DIR}/heterozygosity"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$HET_DIR"

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 07_heterozygosity_array.sh"
    exit 1
fi

if [[ ! -f "$FINAL_CRAMLIST" ]]; then
    echo "ERROR: ${FINAL_CRAMLIST} not found."
    echo "Run 06_downsample_and_finalize.sh first."
    exit 1
fi

mapfile -t CRAMS < "$FINAL_CRAMLIST"
CRAM="${CRAMS[$SLURM_ARRAY_TASK_ID]:-}"

if [[ -z "$CRAM" ]]; then
    echo "ERROR: No CRAM at index ${SLURM_ARRAY_TASK_ID} in ${FINAL_CRAMLIST}"
    exit 1
fi

SAMPLE=$(basename "$CRAM")
SAMPLE="${SAMPLE%_filt.cram}"
SAMPLE="${SAMPLE%_ds.cram}"

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> sample: ${SAMPLE}"
echo ">>> CRAM : ${CRAM}"
echo ">>> CPUs : ${THREADS}"
echo ">>> Start: $(date)"

# =============================================================================
# STEP 1: ANGSD site allele frequency likelihood (genome-wide)
# =============================================================================
echo ">>> Step 1: ANGSD -doSaf"

angsd -i "$CRAM" -ref "$REF_FASTA" -anc "$REF_FASTA" \
    -dosaf 1 -minMapQ 30 -GL 1 -P "$THREADS" \
    -out "${HET_DIR}/${SAMPLE}" \
    -doCounts 1 -setMinDepth 3

# =============================================================================
# STEP 2: realSFS — folded single-sample SFS
# =============================================================================
echo ">>> Step 2: realSFS"

SAF_IDX="${HET_DIR}/${SAMPLE}.saf.idx"
SFS_OUT="${HET_DIR}/${SAMPLE}_est.ml"

if [[ ! -f "$SAF_IDX" ]]; then
    echo "ERROR: ANGSD did not produce expected output: ${SAF_IDX}"
    exit 1
fi

realSFS "$SAF_IDX" -P "$THREADS" -fold 1 > "$SFS_OUT"

# =============================================================================
# STEP 3: Heterozygosity = SFS[1] / (SFS[0] + SFS[1])
# =============================================================================
echo ">>> Step 3: Computing heterozygosity"

HET_VALUE=$(awk '{ if (($1+$2) > 0) printf "%.6f", $2/($1+$2); else print "NA" }' "$SFS_OUT")

# Each task writes only its own file — no shared file for parallel array
# tasks to race on. 08_pca_roh.sh (or `cat *_heterozygosity.txt` directly)
# aggregates these once the array finishes.
echo -e "${SAMPLE}\t${HET_VALUE}" > "${HET_DIR}/${SAMPLE}_heterozygosity.txt"

echo "  Heterozygosity: ${HET_VALUE}"
echo ">>> Sample ${SAMPLE} complete — $(date)"
