#!/bin/bash
#SBATCH --job-name=lepc_rohan_arr
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blackan@purdue.edu
# --array is passed at submission time (see USAGE below), not hardcoded
# here, since the task count depends on how many samples remain.

set -euo pipefail

# =============================================================================
# ROHan per-sample analysis (Renaud et al. 2019) -- job-array version of
# run_rohan.sh's Step 1. One SLURM task per remaining sample, so one slow
# sample's runtime no longer eats into another sample's walltime budget
# (the failure mode that stranded F5550 under the original serial loop).
#
# Reads one sample ID per SLURM_ARRAY_TASK_ID from SAMPLE_LIST (one sample
# per line, 1-indexed to match SLURM array indexing). Build/refresh that
# list first, from the project's results/rohan directory:
#
#   cd /scratch/gautschi/blackan/GROUSE/old_vs_new/results/rohan
#   ls *.win50000.hEst.gz | sed -E 's/\.win50000\.hEst\.gz$//' | sort > completed.txt
#   ls /scratch/gautschi/blackan/GROUSE/old_vs_new/crams/*.md.dedup_q20.cram \
#     | xargs -n1 basename | sed -E 's/\.md\.dedup_q20\.cram$//' | sort > all_samples.txt
#   comm -23 all_samples.txt completed.txt > remaining_samples.txt
#
# USAGE:
#   sbatch --array=1-$(wc -l < remaining_samples.txt) 01_rohan_array.sh
#
# Optional: append e.g. "%5" to the --array range (1-10%5) to cap how many
# tasks run concurrently, if your account has a running-jobs limit.
# =============================================================================

# ---------------- CONFIG -- edit these before running ----------------
PROJ=/scratch/gautschi/blackan/GROUSE/old_vs_new
REF=/scratch/gautschi/blackan/GROUSE/old_vs_new/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna
ROHAN_BIN=$PROJ/tools/ROHan/bin/rohan
GSL_PREFIX=$PROJ/tools/gsl

CRAM_DIR=/scratch/gautschi/blackan/GROUSE/old_vs_new/crams
CRAM_SUFFIX=".md.dedup_q20.cram"

OUT_DIR=$PROJ/results/rohan
SAMPLE_LIST=$OUT_DIR/remaining_samples.txt

THREADS=16
ROHMU=1e-5
ROHAN_WINDOW_SIZE=50000

mkdir -p "$OUT_DIR" logs

if [[ ! -f "$SAMPLE_LIST" ]]; then
    echo "ERROR: ${SAMPLE_LIST} not found. Build it first (see header comment)." >&2
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: no sample at line ${SLURM_ARRAY_TASK_ID} of ${SAMPLE_LIST}" >&2
    exit 1
fi
CRAM="${CRAM_DIR}/${SAMPLE}${CRAM_SUFFIX}"

if [ ! -x "$ROHAN_BIN" ]; then
    echo "ERROR: ROHan binary not found at $ROHAN_BIN" >&2
    echo "  Run install_rohan.sh on the login node first." >&2
    exit 1
fi
if [[ ! -f "$CRAM" ]]; then
    echo "ERROR: CRAM not found: ${CRAM}" >&2
    exit 1
fi

echo ">>> Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE}"
echo ">>> Start time: $(date)"

module --force purge
module load gcc/14.1.0
module load biocontainers
module load samtools

# rohan was linked against a custom-built GSL (no 'gsl' module exists on
# Gautschi), so it needs this to find libgsl.so at runtime, not just at
# build time -- see install_rohan.sh.
export LD_LIBRARY_PATH="$GSL_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

OUT_PREFIX=$OUT_DIR/${SAMPLE}.win${ROHAN_WINDOW_SIZE}
if [[ -f "${OUT_PREFIX}.hEst.gz" ]]; then
    echo "  ${SAMPLE}: ${OUT_PREFIX}.hEst.gz already exists -- skipping (safety net; shouldn't happen if SAMPLE_LIST was built correctly)."
    exit 0
fi

# ROHan's most reliably supported input format is BAM. Rather than rely on
# ROHan's own CRAM/reference handling (undocumented in what we could
# verify), convert to a temporary indexed BAM first -- slower, but removes
# any ambiguity about reference resolution.
TMP_BAM=$OUT_DIR/${SAMPLE}.tmp.bam
samtools view -@ "$THREADS" -b -T "$REF" -o "$TMP_BAM" "$CRAM"
samtools index "$TMP_BAM"

"$ROHAN_BIN" \
    -t "$THREADS" \
    --rohmu "$ROHMU" \
    --size "$ROHAN_WINDOW_SIZE" \
    -o "$OUT_PREFIX" \
    "$REF" "$TMP_BAM"

rm -f "$TMP_BAM" "${TMP_BAM}.bai"

echo "  Done: ${OUT_PREFIX}.*"
echo ">>> End time: $(date)"
