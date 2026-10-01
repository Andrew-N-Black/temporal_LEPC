#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: BUSCO COMPLETENESS + SEX/ORGANELLE SEQUENCE SPLIT
# Step 03 — requires 02_genome_assembly_array.sh to have completed, i.e. the
# final renamed assemblies must exist in ${PROJECT_DIR}/final as
#   <SPECIES>_<SAMPLE>_<HAP>.pseudo_chr.fasta
#
# n=23 samples x 2 haplotypes = 46 assemblies. One array task per SAMPLE;
# each task handles both haplotypes, matching 02_genome_assembly_array.sh.
#
# Per sample, per haplotype (hap1/hap2):
#   1. Split    — pull chr_W, chr_Z and chr_MT out of *.pseudo_chr.fasta into
#                 their own FASTAs, and write everything that remains to
#                 *.autosome_chr.fasta
#   2. BUSCO    — genome-mode completeness against aves_odb10, run on the FULL
#                 *.pseudo_chr.fasta (not the autosome subset), so the numbers
#                 are comparable to published avian assemblies
#
# Step 1 is seconds; Step 2 is hours. Step 1 runs first so the FASTAs are on
# disk even if BUSCO later fails or is pre-empted.
#
# Both steps are resume-safe: existing outputs are detected and skipped, so a
# failed task can simply be resubmitted. Outputs are timestamp-checked against
# the source assembly (`-nt`), so they regenerate automatically if an assembly
# is ever rebuilt rather than silently keeping stale results.
#
# INPUT MANIFEST: the same assembly_manifest.tsv used by step 02. Only the
# first two columns (sample_id, species) are read here.
#
# USAGE:
#   N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
#   sbatch --array=0-$((N-1))%8 03_busco_autosomes_array.sh
#
#   Rerun a single failed task, e.g. index 7:
#   sbatch --array=7 03_busco_autosomes_array.sh
#
#   NOTE on %8: BUSCO is far lighter than the assembly job, so more tasks fit
#   per node. Raise or lower to suit fair-share.
#
# COLLATING RESULTS AFTERWARDS:
#   The per-haplotype one-line summaries are gathered by:
#     cat ${PROJECT_DIR}/qc/busco/*/short_summary.specific.*.txt \
#       | grep -E '^\s+C:'
#   or, sample-labelled:
#     for f in ${PROJECT_DIR}/qc/busco/*/short_summary.specific.*.txt; do
#         printf '%s\t%s\n' "$(basename "$(dirname "$f")")" \
#             "$(grep -m1 -oE 'C:[0-9.]+%\[S:[0-9.]+%,D:[0-9.]+%\],F:[0-9.]+%,M:[0-9.]+%,n:[0-9]+' "$f")"
#     done
# =============================================================================
#SBATCH --job-name=grouse_busco
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
#SBATCH -t 4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
# 64G: BUSCO genome mode on a ~1.05 Gb assembly with metaeuk is dominated by
# the metaeuk search, which stays well under this at 24 threads. Two
# haplotypes run sequentially within the task, so the peak is per-haplotype,
# not additive. Check `sacct -o MaxRSS` after the first few and trim if the
# real peak is much lower.
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

# Defensive: --export=ALL (the sbatch default) carries the submitting shell's
# Lmod state into the job. Start from a known state, as in step 02.
module unload anaconda 2>/dev/null || true

ml biocontainers
ml busco
ml samtools/1.22.1

# unset LD_PRELOAD: RCAC's XALT usage-tracking library is injected via
# LD_PRELOAD and fails on some nodes (GLIBC_2.33/2.34 mismatch), which can
# kill subshells under `set -e`. It's accounting only — safe to drop.
unset LD_PRELOAD

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
MANIFEST="${SLURM_SUBMIT_DIR}/assembly_manifest.tsv"

# Same PROJECT_DIR as steps 01/02 so the existing assemblies are found.
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
FINAL_DIR="${PROJECT_DIR}/final"
QC_DIR="${PROJECT_DIR}/qc"
BUSCO_DIR="${QC_DIR}/busco"
SPLIT_DIR="${PROJECT_DIR}/final_autosomes"

# ---- Sequences to pull out of the assembly ----------------------------------
# These are the post-rename names produced by Step 7b of 02_genome_assembly_
# array.sh, which come from the NCBI assembly report for GRCg7b: the chicken
# mitochondrion is named "MT" there, so it becomes chr_MT here.
# Any name in this list that is absent from a given assembly is skipped with a
# note rather than treated as an error — chr_W is legitimately absent from a
# clean male assembly, and chr_MT is only present when RagTag actually placed
# a mitochondrial contig.
EXTRACT_SEQS=(chr_W chr_Z chr_MT)

# ---- What goes into *.autosome_chr.fasta ------------------------------------
# true  : autosomes (chr_1 .. chr_33) PLUS the unplaced scaffold_N sequences.
#         Unplaced scaffolds are mostly autosomal sequence that RagTag could
#         not place, so keeping them avoids silently discarding real data.
#         This is the default.
# false : named autosomes only (chr_* minus the extracted set). Use this when
#         you want a strictly chromosome-level autosomal set, e.g. for a
#         window-based scan where unplaced scaffolds would be noise.
#
# Either way nothing is lost: whatever is excluded is written to
# *.excluded_from_autosomes.fasta alongside, so the split is auditable.
INCLUDE_UNPLACED_SCAFFOLDS=true

# ---- BUSCO ------------------------------------------------------------------
# aves_odb10 is the standard avian lineage and is what comparable galliform
# assemblies report against. If the loaded BUSCO is v6+, its default database
# is odb12 — set BUSCO_LINEAGE=aves_odb12 below if you would rather match that,
# but do not mix lineages across samples, the scores are not comparable.
BUSCO_LINEAGE="aves_odb10"
BUSCO_MODE="genome"
BUSCO_DOWNLOAD_PATH="${PROJECT_DIR}/busco_downloads"

THREADS=$SLURM_CPUS_PER_TASK

# Locate BUSCO's completion marker without tripping `set -e`.
#
# BUSCO writes short_summary.specific.<lineage>.<name>.txt into its run
# directory only on success, so that file's presence is the resume signal.
# Finding it needs care: on a first run the run directory does not exist yet,
# `find` then exits 1, `set -o pipefail` carries that through the pipe, and
# under `set -e` a failing command substitution in an assignment kills the
# script outright. That is exactly how the first version of this script died
# — silently, before BUSCO was ever invoked, because 2>/dev/null had also
# swallowed find's error. Guard the directory test first and swallow any
# remaining non-zero status explicitly.
find_busco_summary() {
    local dir="$1"
    [[ -d "$dir" ]] || { printf ''; return 0; }
    find "$dir" -maxdepth 1 -name 'short_summary.specific.*.txt' -print 2>/dev/null \
        | head -n1 || true
}

mkdir -p logs "$BUSCO_DIR" "$SPLIT_DIR" "$BUSCO_DOWNLOAD_PATH"

# Verify the tools actually landed on PATH. `ml busco` after `ml biocontainers`
# is the same pattern quast uses in step 02; if this cluster ever renames the
# module, fail here with a clear message rather than deep inside a BUSCO run.
for BIN in busco samtools; do
    if ! command -v "$BIN" >/dev/null 2>&1; then
        echo "ERROR: '${BIN}' not on PATH after module load."
        echo "       Check: ml spider ${BIN}"
        exit 1
    fi
done
# `busco` is a shell function on this cluster (the biocontainers module wraps
# a container), so `command -v` prints the name rather than a path — that is
# expected, not a broken install. `|| true` because head closes the pipe after
# one line and pipefail would otherwise surface busco's SIGPIPE status.
echo ">>> busco    : $(command -v busco)  [$(busco --version 2>&1 | head -n1 || true)]"
echo ">>> samtools : $(command -v samtools)"

if [[ ! -d "$FINAL_DIR" ]]; then
    echo "ERROR: final assembly directory not found: ${FINAL_DIR}"
    echo "       Run 02_genome_assembly_array.sh first."
    exit 1
fi

# =============================================================================
# SHARED ONE-TIME SETUP — BUSCO lineage download
# Serialized with flock: many array tasks starting at once would otherwise all
# try to download aves_odb10 into the same path and corrupt it. The first task
# to take the lock downloads; the rest wait, then find it present and skip.
# Downloading once up front also lets every BUSCO run use --offline, so no
# task depends on compute-node internet access at run time.
# =============================================================================
exec 9>"${BUSCO_DOWNLOAD_PATH}/.setup.lock"
flock 9

BUSCO_LINEAGE_DIR="${BUSCO_DOWNLOAD_PATH}/lineages/${BUSCO_LINEAGE}"
if [[ ! -s "${BUSCO_LINEAGE_DIR}/dataset.cfg" ]]; then
    echo ">>> Downloading BUSCO lineage ${BUSCO_LINEAGE}"
    busco --download "$BUSCO_LINEAGE" --download_path "$BUSCO_DOWNLOAD_PATH" \
        || echo "  WARNING: busco --download returned non-zero; checking for the dataset anyway"
    if [[ ! -s "${BUSCO_LINEAGE_DIR}/dataset.cfg" ]]; then
        echo "ERROR: BUSCO lineage ${BUSCO_LINEAGE} not present at ${BUSCO_LINEAGE_DIR}"
        echo "       Download it manually, e.g.:"
        echo "         busco --download ${BUSCO_LINEAGE} --download_path ${BUSCO_DOWNLOAD_PATH}"
        echo "       or fetch https://busco-data.ezlab.org/v5/data/lineages/ and untar it there."
        flock -u 9; exec 9>&-
        exit 1
    fi
fi

flock -u 9
exec 9>&-

echo ">>> BUSCO lineage: ${BUSCO_LINEAGE_DIR}"

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 03_busco_autosomes_array.sh"
    exit 1
fi
if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

mapfile -t ROWS < <(grep -v '^#' "$MANIFEST" | tail -n +2 | grep -v '^[[:space:]]*$')
LINE="${ROWS[$SLURM_ARRAY_TASK_ID]:-}"
if [[ -z "$LINE" ]]; then
    echo "ERROR: No manifest row at index ${SLURM_ARRAY_TASK_ID} (${#ROWS[@]} samples in manifest)"
    exit 1
fi

IFS=$'\t' read -r SAMPLE SPECIES _REST <<< "$LINE"
if [[ -z "${SAMPLE:-}" || -z "${SPECIES:-}" ]]; then
    echo "ERROR: Malformed manifest row at index ${SLURM_ARRAY_TASK_ID}: ${LINE}"
    exit 1
fi

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> sample: ${SAMPLE} (${SPECIES})"
echo ">>> Started: $(date)"
echo ">>> Threads: ${THREADS}"

# TMPDIR is deliberately left exactly as SLURM sets it. An earlier version
# pointed it at a per-task subdirectory, which is an unnecessary risk here:
# `busco` runs inside a container that bind-mounts only certain paths, so a
# custom TMPDIR can be invisible to the tool that has to write there. BUSCO
# keeps its working files under --out_path anyway, and each task has its own
# run directory there, so concurrent tasks cannot collide regardless.

# =============================================================================
# PER-HAPLOTYPE PROCESSING
# =============================================================================
for HAP in hap1 hap2; do

    PREFIX="${SPECIES}_${SAMPLE}_${HAP}"
    ASSEMBLY="${FINAL_DIR}/${PREFIX}.pseudo_chr.fasta"

    echo ""
    echo "============================================================"
    echo ">>> ${PREFIX}"
    echo "============================================================"

    if [[ ! -s "$ASSEMBLY" ]]; then
        echo "ERROR: assembly not found: ${ASSEMBLY}"
        echo "       Step 02 must finish for this sample/haplotype first."
        exit 1
    fi

    # -------------------------------------------------------------------------
    # STEP 1 — Split out chr_W / chr_Z / chr_MT, build *.autosome_chr.fasta
    # -------------------------------------------------------------------------
    echo ""
    echo ">>> [1/2] Splitting sex chromosomes and mitochondrion"

    AUTOSOME_FASTA="${SPLIT_DIR}/${PREFIX}.autosome_chr.fasta"
    EXCLUDED_FASTA="${SPLIT_DIR}/${PREFIX}.excluded_from_autosomes.fasta"
    SPLIT_REPORT="${SPLIT_DIR}/${PREFIX}.split_report.tsv"

    if [[ -s "$AUTOSOME_FASTA" && -f "$SPLIT_REPORT" && "$AUTOSOME_FASTA" -nt "$ASSEMBLY" ]]; then
        echo "  already split — skipping"
    else
        # samtools faidx needs the index; Step 02 usually leaves one, but a
        # rebuilt assembly can leave it stale, so refresh when it is older.
        if [[ ! -s "${ASSEMBLY}.fai" || "${ASSEMBLY}.fai" -ot "$ASSEMBLY" ]]; then
            samtools faidx "$ASSEMBLY"
        fi

        # --- pull each requested sequence into its own FASTA ---
        : > "$SPLIT_REPORT"
        printf 'sequence\tstatus\tlength_bp\toutput_file\n' >> "$SPLIT_REPORT"

        for SEQ in "${EXTRACT_SEQS[@]}"; do
            SEQ_LEN=$(awk -F'\t' -v s="$SEQ" '$1 == s {print $2; exit}' "${ASSEMBLY}.fai")
            SEQ_FASTA="${SPLIT_DIR}/${PREFIX}.${SEQ}.fasta"
            if [[ -n "$SEQ_LEN" ]]; then
                samtools faidx "$ASSEMBLY" "$SEQ" > "$SEQ_FASTA"
                printf '%s\tpresent\t%s\t%s\n' "$SEQ" "$SEQ_LEN" "$(basename "$SEQ_FASTA")" >> "$SPLIT_REPORT"
                echo "  ${SEQ}: present (${SEQ_LEN} bp) -> $(basename "$SEQ_FASTA")"
            else
                # Absent is a normal outcome, not a failure: chr_W should be
                # absent from a clean male assembly, and chr_MT only appears
                # when a mitochondrial contig was placed.
                rm -f "$SEQ_FASTA"
                printf '%s\tabsent\tNA\tNA\n' "$SEQ" >> "$SPLIT_REPORT"
                echo "  ${SEQ}: absent from this assembly — nothing to extract"
            fi
        done

        # --- build the keep / exclude lists ---
        KEEP_LIST="${SPLIT_DIR}/${PREFIX}.autosome_keep.txt"
        EXCL_LIST="${SPLIT_DIR}/${PREFIX}.autosome_excluded.txt"
        EXTRACT_CSV=$(IFS=,; echo "${EXTRACT_SEQS[*]}")

        awk -F'\t' -v excl="$EXTRACT_CSV" -v keepscaf="$INCLUDE_UNPLACED_SCAFFOLDS" '
            BEGIN {
                n = split(excl, e, ",")
                for (i = 1; i <= n; i++) drop[e[i]] = 1
            }
            {
                if ($1 in drop)                            next   # chr_W/chr_Z/chr_MT
                if ($1 !~ /^chr_/ && keepscaf != "true")   next   # unplaced scaffolds
                print $1
            }
        ' "${ASSEMBLY}.fai" > "$KEEP_LIST"

        awk -F'\t' -v excl="$EXTRACT_CSV" -v keepscaf="$INCLUDE_UNPLACED_SCAFFOLDS" '
            BEGIN {
                n = split(excl, e, ",")
                for (i = 1; i <= n; i++) drop[e[i]] = 1
            }
            {
                if (($1 in drop) || ($1 !~ /^chr_/ && keepscaf != "true")) print $1
            }
        ' "${ASSEMBLY}.fai" > "$EXCL_LIST"

        if [[ ! -s "$KEEP_LIST" ]]; then
            echo "ERROR: autosome list came out empty for ${PREFIX}"
            echo "       Check the sequence names in ${ASSEMBLY}.fai"
            exit 1
        fi

        # `-r <file>` rather than xargs: the sequence list is long and this is
        # the form verified to work with this samtools build.
        samtools faidx "$ASSEMBLY" -r "$KEEP_LIST" > "$AUTOSOME_FASTA"
        if [[ -s "$EXCL_LIST" ]]; then
            samtools faidx "$ASSEMBLY" -r "$EXCL_LIST" > "$EXCLUDED_FASTA"
        else
            : > "$EXCLUDED_FASTA"
        fi
        samtools faidx "$AUTOSOME_FASTA"

        # --- arithmetic check: nothing gained, nothing lost ---
        N_IN=$(wc -l < "${ASSEMBLY}.fai")
        N_KEEP=$(wc -l < "$KEEP_LIST")
        N_EXCL=$(wc -l < "$EXCL_LIST")
        if (( N_KEEP + N_EXCL != N_IN )); then
            echo "ERROR: sequence counts do not reconcile for ${PREFIX}:"
            echo "       input ${N_IN}, kept ${N_KEEP}, excluded ${N_EXCL}"
            exit 1
        fi

        BP_IN=$(awk -F'\t' '{s+=$2} END {print s+0}' "${ASSEMBLY}.fai")
        BP_KEEP=$(awk -F'\t' '{s+=$2} END {print s+0}' "${AUTOSOME_FASTA}.fai")
        echo "  input      : ${N_IN} sequences, ${BP_IN} bp"
        echo "  autosomes  : ${N_KEEP} sequences, ${BP_KEEP} bp -> $(basename "$AUTOSOME_FASTA")"
        echo "  excluded   : ${N_EXCL} sequences          -> $(basename "$EXCLUDED_FASTA")"
        echo "  unplaced scaffolds kept in autosome file: ${INCLUDE_UNPLACED_SCAFFOLDS}"
    fi

    # -------------------------------------------------------------------------
    # STEP 2 — BUSCO on the full assembly
    # -------------------------------------------------------------------------
    echo ""
    echo ">>> [2/2] BUSCO (${BUSCO_LINEAGE}, ${BUSCO_MODE} mode)"

    BUSCO_OUT_NAME="${PREFIX}"
    BUSCO_RUN_DIR="${BUSCO_DIR}/${BUSCO_OUT_NAME}"
    # BUSCO writes short_summary.specific.<lineage>.<name>.txt on success;
    # its presence is the only reliable "this finished" marker.
    BUSCO_SUMMARY=$(find_busco_summary "$BUSCO_RUN_DIR")

    if [[ -n "$BUSCO_SUMMARY" && -s "$BUSCO_SUMMARY" && "$BUSCO_SUMMARY" -nt "$ASSEMBLY" ]]; then
        echo "  already complete — skipping"
        echo "  $(grep -m1 -E '^[[:space:]]+C:' "$BUSCO_SUMMARY" || true)"
    else
        # A partial run from a killed task leaves a directory BUSCO refuses to
        # write into. Clear it rather than passing -f blindly, so a genuinely
        # complete run above is never overwritten.
        if [[ -d "$BUSCO_RUN_DIR" ]]; then
            echo "  removing incomplete previous run: ${BUSCO_RUN_DIR}"
            rm -rf "$BUSCO_RUN_DIR"
        fi

        busco \
            --in "$ASSEMBLY" \
            --out "$BUSCO_OUT_NAME" \
            --out_path "$BUSCO_DIR" \
            --mode "$BUSCO_MODE" \
            --lineage_dataset "$BUSCO_LINEAGE_DIR" \
            --offline \
            --download_path "$BUSCO_DOWNLOAD_PATH" \
            --cpu "$THREADS"

        BUSCO_SUMMARY=$(find_busco_summary "$BUSCO_RUN_DIR")
        if [[ -z "$BUSCO_SUMMARY" || ! -s "$BUSCO_SUMMARY" ]]; then
            echo "ERROR: BUSCO finished without writing a short summary for ${PREFIX}"
            echo "       Inspect ${BUSCO_RUN_DIR}"
            exit 1
        fi
        echo "  $(grep -m1 -E '^[[:space:]]+C:' "$BUSCO_SUMMARY" || true)"
    fi

done

# =============================================================================
# DONE
# =============================================================================
echo ""
echo "============================================================"
echo ">>> Sample ${SAMPLE} (${SPECIES}) complete: $(date)"
echo "  autosome FASTAs : ${SPLIT_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}.autosome_chr.fasta"
echo "  extracted seqs  : ${SPLIT_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}.chr_{W,Z,MT}.fasta"
echo "  split reports   : ${SPLIT_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}.split_report.tsv"
echo "  BUSCO           : ${BUSCO_DIR}/${SPECIES}_${SAMPLE}_hap{1,2}/"
echo "============================================================"
