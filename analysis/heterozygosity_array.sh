#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: PER-SAMPLE HETEROZYGOSITY (ANGSD + realSFS)
#
# Restricted to autosomes: -rf points ANGSD at every contig in the
# reference EXCEPT the two Z scaffolds (NW_026294758.1, NW_026294813.1;
# see analysis/find_sex_scaffolds.sh). Without this, females (hemizygous
# for Z) show artificially low heterozygosity on Z-linked sites, biasing
# the genome-wide estimate downward for females specifically -- this was
# flagged in AUDIT.md as still needing "autosomal confirmation" and was
# previously whole-genome. Same two scaffold IDs used everywhere else in
# this repo (run_plink.sh, remove_Z_scaffolds.sh).
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
#SBATCH -A fnrdewoody
#SBATCH -t 3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu
#SBATCH -p cpu

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
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"
AUTOSOME_RF="${HET_DIR}/autosomes.rf"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$HET_DIR"

# =============================================================================
# BUILD THE AUTOSOMES-ONLY REGION FILE (ANGSD -rf format: "<contig>:" per
# line = whole contig). Shared across all array tasks, so this is built
# race-safely: each task that finds it missing writes its own temp file,
# then atomically renames it into place with "mv -n" -- a losing task's
# rename becomes a silent no-op rather than clobbering a file another task
# already finished writing, and its temp file is cleaned up either way.
# =============================================================================
if [[ ! -s "$AUTOSOME_RF" ]]; then
    if [[ ! -f "${REF_FASTA}.fai" ]]; then
        echo "ERROR: ${REF_FASTA}.fai not found." >&2
        exit 1
    fi
    # Confirm both excluded scaffolds actually exist in the reference
    # (catches a typo'd/renamed scaffold ID loudly instead of silently
    # excluding nothing).
    for c in ${Z_SCAFFOLDS//,/ }; do
        awk -v c="$c" '$1==c {f=1} END {exit !f}' "${REF_FASTA}.fai" \
            || { echo "ERROR: excluded scaffold '$c' not in ${REF_FASTA}.fai" >&2; exit 1; }
    done
    TMP_RF=$(mktemp "${AUTOSOME_RF}.XXXXXX")
    awk -v z="$Z_SCAFFOLDS" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1}
        !($1 in Z) {print $1":"}' "${REF_FASTA}.fai" > "$TMP_RF"
    if mv -n "$TMP_RF" "$AUTOSOME_RF" 2>/dev/null; then
        :
    else
        rm -f "$TMP_RF"
    fi
fi
N_AUTOSOMES=$(wc -l < "$AUTOSOME_RF")
echo ">>> Autosome region file: ${AUTOSOME_RF} (${N_AUTOSOMES} contigs, excluding ${Z_SCAFFOLDS})"

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
# final_cramlist.txt entries are named <ID>.md.dedup_q20.cram (confirmed
# against processing/final_cramlist.txt and the suffix run_ROHan.sh already
# strips) -- NOT _filt.cram or _ds.cram, which don't match anything in this
# project and previously left SAMPLE as the full, undotted filename. That
# produced garbled sample IDs downstream (e.g. a heterozygosity_summary.tsv
# row for "F10.md.dedup_q20.cram" instead of "F10"), breaking any join
# against popmap.txt or the metadata spreadsheet by sample ID.
SAMPLE="${SAMPLE%.md.dedup_q20.cram}"

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> sample: ${SAMPLE}"
echo ">>> CRAM : ${CRAM}"
echo ">>> CPUs : ${THREADS}"
echo ">>> Start: $(date)"

# =============================================================================
# STEP 1: ANGSD site allele frequency likelihood (autosomes only)
# =============================================================================
echo ">>> Step 1: ANGSD -doSaf (autosomes only, -rf ${AUTOSOME_RF})"

angsd -i "$CRAM" -ref "$REF_FASTA" -anc "$REF_FASTA" \
    -rf "$AUTOSOME_RF" \
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

echo "  Heterozygosity (autosomal): ${HET_VALUE}"
echo ">>> Sample ${SAMPLE} complete — $(date)"
