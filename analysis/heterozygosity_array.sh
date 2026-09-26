#!/bin/bash
# =============================================================================
# SLURM ARRAY JOB: DE NOVO GENOME ASSEMBLY — hifiasm + yahs + RagTag
# Step 02 — requires 01_download_reference_genomes.sh to have completed first.
# n=23 samples across 3 species (Sharp-tailed Grouse, Lesser Prairie-Chicken,
# Greater Prairie-Chicken). One array task per sample.
#
# Built from 02_test_single_sample.sh after it was validated end-to-end on
# F5545 (LEPC). Per sample, per haplotype (hap1/hap2):
#   1.  HiFiAdapterFilt   — adapter filtering of raw HiFi BAMs (per sample)
#   2.  hifiasm           — Hi-C-phased assembly, +ONT UL when available (per sample)
#   3.  GFA -> FASTA
#   3b. QUAST             — post-hifiasm contig stats
#   4.  bwa mem -5SP      — Hi-C reads -> contigs (MAPQ >= 20)
#   5.  yahs              — Hi-C scaffolding
#   5b. QUAST             — post-yahs scaffold stats
#   6.  juicer pre + juicer_tools — .hic/.assembly for Juicebox curation
#   7.  RagTag            — pseudo-chromosome ordering vs chicken (GRCg7b)
#   7b. Rename            — chr_1..chr_33/chr_Z/chr_W from the NCBI assembly
#                           report; scaffold_N for non-chromosome matches
#   7c. QUAST             — final assembly stats
#   7d. tidk              — telomere (TTAGGG) profiling + plot, named chromosomes
#   8.  Liftoff           — chicken annotation transfer
#   9.  minimap2 dotplot  — independent order/orientation check vs chicken
#   10. Pretext           — Hi-C contact map on the final assembly
#   11. HiFi depth check  — chr_Z/chr_W depth vs autosomes (sex-chromosome QC)
#
# Every step is resume-safe: completed outputs are detected and skipped, so a
# failed task can simply be resubmitted. Steps 7d-11 compare output timestamps
# against the final FASTA (`-nt`), so they rerun automatically if the final
# assembly is ever regenerated, rather than silently keeping stale results.
#
# INPUT MANIFEST (tab-separated, header required, see assembly_manifest.tsv):
#   sample_id  species  hifi_bams  hic_r1  hic_r2  ont_ul
#   - hifi_bams : raw, unaligned PacBio HiFi BAM(s) (*.hifi_reads.bc####.bam),
#                 comma-separated if multiple SMRT cells
#   - hic_r1/r2 : paired-end Hi-C fastq.gz
#   - ont_ul    : ONT ultra-long fastq.gz, or "NA"
#
# USAGE:
#   N=$(grep -v '^#' assembly_manifest.tsv | tail -n +2 | grep -c .)
#   sbatch --array=0-$((N-1))%6 02_genome_assembly_array.sh
#   (%6 caps concurrent tasks — tune to fair-share/node availability)
#
#   Rerun a single failed task, e.g. index 7:
#   sbatch --array=7 02_genome_assembly_array.sh
# =============================================================================
#SBATCH --job-name=grouse_asm
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
# 160G: measured peak RSS across 17 array tasks was 128-143 GB (sacct MaxRSS,
# whole-cgroup, at 128 threads with unpinned samtools sort). With 48 threads
# and `samtools sort -m 1G` below, the sort contribution is capped near 48 GB,
# so the real peak should land comfortably under this. On a 377 GB node this
# fits 2 tasks concurrently. Re-check `sacct -o MaxRSS` after a few samples;
# if the peak drops below ~120 GB, 3 tasks per node becomes viable.
#SBATCH --mem=160G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

# Defensive: if the submitting shell had anaconda loaded, --export=ALL (the
# sbatch default) carries that Lmod state into the job, and liftoff will
# refuse to load later. Start from a known state.
module unload anaconda 2>/dev/null || true

ml biocontainers
ml quast
ml hifiasm
ml bwa
ml samtools/1.22.1
unset LD_PRELOAD
# unset LD_PRELOAD: RCAC's XALT usage-tracking library is injected via
# LD_PRELOAD and fails on some nodes (GLIBC_2.33/2.34 mismatch), which can
# kill subshells under `set -e`. It's accounting only — safe to drop.
#
# ragtag/liftoff/minimap2 are deliberately NOT loaded here. Lmod on this
# cluster refuses to have liftoff and anaconda loaded at the same time, and
# anaconda is still needed below for the conda installs and for Step 1
# (HiFiAdapterFilt needs `conda activate`). They're loaded right after
# Step 1, once anaconda activity for the whole run is finished.
#
# yahs, Pretext, juicer_tools, tidk, HiFiAdapterFilt and the dotplot Python
# env have no modules here — installed via conda/bioconda into shared envs
# under PROJECT_DIR (see below). Tools are invoked by absolute path rather
# than bare names, which also sidesteps the liftoff module's `python3` shell
# function shadowing any later python3 call.

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
MANIFEST="${SLURM_SUBMIT_DIR}/assembly_manifest.tsv"

# Same PROJECT_DIR as the validated test run, so the shared conda envs and
# chicken reference built there are reused rather than rebuilt.
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/grouse_asm"
REF_DIR="${PROJECT_DIR}/ref"

FILT_DIR="${PROJECT_DIR}/hifi_filtered"
ASM_DIR="${PROJECT_DIR}/hifiasm"
SCAFFOLD_DIR="${PROJECT_DIR}/yahs"
RAGTAG_DIR="${PROJECT_DIR}/ragtag"
LIFTOFF_DIR="${PROJECT_DIR}/liftoff"
FINAL_DIR="${PROJECT_DIR}/final"
QC_DIR="${PROJECT_DIR}/qc"

CHICKEN_ASM_DIR="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b"
CHICKEN_FASTA="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.fna"
CHICKEN_GFF="${REF_DIR}/${CHICKEN_ASM_DIR}_genomic.gff"

# NCBI assembly report -> RefSeq accession to chromosome-name map (Step 7b).
CHICKEN_ASSEMBLY_REPORT="${REF_DIR}/${CHICKEN_ASM_DIR}_assembly_report.txt"
CHICKEN_CHR_MAP="${REF_DIR}/${CHICKEN_ASM_DIR}.chr_map.tsv"

HIC_MIN_MAPQ=20

# Minimum length for UNPLACED scaffolds in the final assembly (Step 7b).
# Sequences RagTag assigned to a real chicken chromosome (chr_*) are ALWAYS
# kept regardless of length — several galliform microchromosomes are
# legitimately small (chicken chr33 is only ~2-3 Mb), so a blanket size cut
# would risk discarding real chromosome content. This threshold applies only
# to the unplaced scaffold_N pile, which is largely repeat-derived fragments
# that inflate scaffold counts and clutter downstream QC. Removed sequence is
# not discarded: it is written to <assembly>.unplaced_short.fasta alongside.
MIN_UNPLACED_SCAFFOLD_LEN=50000

CONDA_ENVS_DIR="${PROJECT_DIR}/conda_envs"

YAHS_VERSION="1.2.2"
YAHS_ENV_DIR="${CONDA_ENVS_DIR}/yahs-${YAHS_VERSION}"
YAHS_BIN="${YAHS_ENV_DIR}/bin/yahs"
JUICER_BIN="${YAHS_ENV_DIR}/bin/juicer"   # yahs's own JBAT pre-processor

# juicer_tools via bioconda, pinned. Manual jar downloads proved unreliable
# (the classic S3 mirror 403s) and v3.0.0 dropped the classic `pre` syntax.
JUICER_TOOLS_VERSION="2.20.00"
JUICER_TOOLS_ENV_DIR="${CONDA_ENVS_DIR}/juicertools-${JUICER_TOOLS_VERSION}"

PRETEXT_ENV_DIR="${CONDA_ENVS_DIR}/pretext"
PRETEXTMAP_BIN="${PRETEXT_ENV_DIR}/bin/PretextMap"
PRETEXTSNAPSHOT_BIN="${PRETEXT_ENV_DIR}/bin/PretextSnapshot"

DOTPLOT_ENV_DIR="${CONDA_ENVS_DIR}/dotplot-python"
DOTPLOT_PYTHON_BIN="${DOTPLOT_ENV_DIR}/bin/python3"

# ACTIVATED (not absolute-path invoked) in Step 1 — its bioconda packaging
# relies on conda's activate hook to put its adapter database on PATH.
HIFIADAPTERFILT_ENV_DIR="${CONDA_ENVS_DIR}/hifiadapterfilt"
HIFIADAPTERFILT_SCRIPT="${HIFIADAPTERFILT_ENV_DIR}/bin/hifiadapterfilt.sh"

TIDK_ENV_DIR="${CONDA_ENVS_DIR}/tidk"
TIDK_BIN="${TIDK_ENV_DIR}/bin/tidk"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$CONDA_ENVS_DIR" "$REF_DIR" "$FILT_DIR" "$ASM_DIR" "$SCAFFOLD_DIR" \
    "$RAGTAG_DIR" "$LIFTOFF_DIR" "$FINAL_DIR" "$QC_DIR" \
    "${QC_DIR}/quast" "${QC_DIR}/tidk" "${QC_DIR}/depth_check"

# =============================================================================
# SHARED ONE-TIME SETUP (conda envs + chicken chromosome map)
# Serialized with flock: with many array tasks starting at once, two tasks
# running `conda create` into the same prefix (or writing the same map file)
# simultaneously would corrupt it. The first task to get the lock does the
# work; the rest wait, then find everything already present and skip.
# =============================================================================
exec 9>"${CONDA_ENVS_DIR}/.setup.lock"
flock 9

NEED_ANACONDA=false
[[ ! -x "$YAHS_BIN" || ! -x "$JUICER_BIN" ]]              && NEED_ANACONDA=true
[[ ! -x "$PRETEXTMAP_BIN" || ! -x "$PRETEXTSNAPSHOT_BIN" ]] && NEED_ANACONDA=true
[[ ! -x "$DOTPLOT_PYTHON_BIN" ]]                          && NEED_ANACONDA=true
[[ ! -f "$HIFIADAPTERFILT_SCRIPT" ]]                      && NEED_ANACONDA=true
[[ ! -d "$JUICER_TOOLS_ENV_DIR" ]]                        && NEED_ANACONDA=true
[[ ! -x "$TIDK_BIN" ]]                                    && NEED_ANACONDA=true

if [[ "$NEED_ANACONDA" == true ]]; then
    ml anaconda/2025.12-py313

    if [[ ! -x "$YAHS_BIN" || ! -x "$JUICER_BIN" ]]; then
        echo ">>> Installing yahs v${YAHS_VERSION}"
        conda create --yes --override-channels --prefix "$YAHS_ENV_DIR" -c bioconda -c conda-forge "yahs=${YAHS_VERSION}"
    fi
    if [[ ! -x "$PRETEXTMAP_BIN" || ! -x "$PRETEXTSNAPSHOT_BIN" ]]; then
        echo ">>> Installing PretextMap/PretextSnapshot"
        conda create --yes --override-channels --prefix "$PRETEXT_ENV_DIR" -c bioconda -c conda-forge \
            pretextmap=0.2.4 pretextsnapshot=0.0.7
    fi
    if [[ ! -x "$DOTPLOT_PYTHON_BIN" ]]; then
        echo ">>> Installing Python 3 + matplotlib"
        conda create --yes --override-channels --prefix "$DOTPLOT_ENV_DIR" -c conda-forge python=3.11 matplotlib
    fi
    if [[ ! -f "$HIFIADAPTERFILT_SCRIPT" ]]; then
        echo ">>> Installing HiFiAdapterFilt"
        conda create --yes --override-channels --prefix "$HIFIADAPTERFILT_ENV_DIR" -c bioconda -c conda-forge hifiadapterfilt
    fi
    if [[ ! -d "$JUICER_TOOLS_ENV_DIR" ]]; then
        echo ">>> Installing juicertools v${JUICER_TOOLS_VERSION}"
        conda create --yes --override-channels --prefix "$JUICER_TOOLS_ENV_DIR" -c bioconda -c conda-forge "juicertools=${JUICER_TOOLS_VERSION}"
    fi
    if [[ ! -x "$TIDK_BIN" ]]; then
        echo ">>> Installing tidk"
        conda create --yes --override-channels --prefix "$TIDK_ENV_DIR" -c bioconda -c conda-forge tidk
    fi

    module unload anaconda/2025.12-py313
fi

# Chicken assembly report + chromosome-name map (inside the lock too).
if [[ ! -f "$CHICKEN_ASSEMBLY_REPORT" ]]; then
    echo ">>> Downloading NCBI assembly report for ${CHICKEN_ASM_DIR}"
    CHICKEN_ACCESSION=$(echo "$CHICKEN_ASM_DIR" | cut -d'_' -f1-2)
    ACC_DIGITS=$(echo "$CHICKEN_ACCESSION" | sed -E 's/^GCF_([0-9]+)\..*$/\1/')
    REPORT_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${ACC_DIGITS:0:3}/${ACC_DIGITS:3:3}/${ACC_DIGITS:6:3}/${CHICKEN_ASM_DIR}/${CHICKEN_ASM_DIR}_assembly_report.txt"
    wget -q -O "$CHICKEN_ASSEMBLY_REPORT" "$REPORT_URL" || true
    if [[ ! -s "$CHICKEN_ASSEMBLY_REPORT" ]]; then
        echo "ERROR: Failed to download ${REPORT_URL}"
        echo "Download it manually to: ${CHICKEN_ASSEMBLY_REPORT}"
        rm -f "$CHICKEN_ASSEMBLY_REPORT"
        exit 1
    fi
fi

if [[ ! -s "$CHICKEN_CHR_MAP" ]]; then
    echo ">>> Building chromosome-name map (assembled-molecule rows only)"
    # Report columns: 1 Sequence-Name, 2 Sequence-Role, ..., 7 RefSeq-Accn
    awk -F'\t' '!/^#/ && $2 == "assembled-molecule" { print $7"\t"$1 }' \
        "$CHICKEN_ASSEMBLY_REPORT" > "$CHICKEN_CHR_MAP"
    if [[ ! -s "$CHICKEN_CHR_MAP" ]]; then
        echo "ERROR: Chromosome-name map came out empty — check ${CHICKEN_ASSEMBLY_REPORT}"
        rm -f "$CHICKEN_CHR_MAP"
        exit 1
    fi
fi

flock -u 9
exec 9>&-

# Verify every tool is actually in place.
for BIN in "$YAHS_BIN" "$JUICER_BIN" "$PRETEXTMAP_BIN" "$PRETEXTSNAPSHOT_BIN" "$DOTPLOT_PYTHON_BIN" "$TIDK_BIN"; do
    if [[ ! -x "$BIN" ]]; then
        echo "ERROR: expected tool not found/executable: ${BIN}"
        exit 1
    fi
done
if [[ ! -f "$HIFIADAPTERFILT_SCRIPT" ]]; then
    echo "ERROR: HiFiAdapterFilt not found: ${HIFIADAPTERFILT_SCRIPT}"
    exit 1
fi

# Resolve juicer_tools: bioconda wrapper script if present, else the jar.
if [[ -x "${JUICER_TOOLS_ENV_DIR}/bin/juicer_tools" ]]; then
    JUICER_TOOLS_CMD=("${JUICER_TOOLS_ENV_DIR}/bin/juicer_tools" -Xmx32G)
else
    JUICER_TOOLS_JAR_FOUND=$(find "$JUICER_TOOLS_ENV_DIR" -iname "juicer_tools*.jar" 2>/dev/null | head -n1)
    if [[ -n "$JUICER_TOOLS_JAR_FOUND" ]]; then
        JUICER_TOOLS_CMD=(java -Xmx32G -jar "$JUICER_TOOLS_JAR_FOUND")
    else
        echo "ERROR: no juicer_tools executable or jar under ${JUICER_TOOLS_ENV_DIR}"
        exit 1
    fi
fi

if [[ ! -f "$CHICKEN_FASTA" || ! -f "$CHICKEN_GFF" ]]; then
    echo "ERROR: Chicken reference not found at ${REF_DIR}"
    echo "Run 01_download_reference_genomes.sh first."
    exit 1
fi

# =============================================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# =============================================================================
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit with: sbatch --array=0-N 02_genome_assembly_array.sh"
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

IFS=$'\t' read -r SAMPLE SPECIES HIFI_BAMS HIC_R1 HIC_R2 ONT_UL <<< "$LINE"
ONT_UL="${ONT_UL:-NA}"
ONT_UL="${ONT_UL%$'\r'}"   # tolerate Windows line endings in the manifest

if [[ -z "$SAMPLE" || -z "$SPECIES" || -z "$HIFI_BAMS" || -z "$HIC_R1" || -z "$HIC_R2" ]]; then
    echo "ERROR: Malformed manifest row at index ${SLURM_ARRAY_TASK_ID}: ${LINE}"
    exit 1
fi

echo ">>> Array task ${SLURM_ARRAY_TASK_ID} -> sample: ${SAMPLE} (${SPECIES})"
echo ">>> HiFi BAMs  : ${HIFI_BAMS}"
echo ">>> Hi-C reads : ${HIC_R1} / ${HIC_R2}"
echo ">>> ONT UL     : ${ONT_UL}"
echo ">>> Running on : $(hostname)"
echo ">>> CPUs       : ${THREADS}"
echo ">>> juicer_tools: ${JUICER_TOOLS_CMD[*]}"
echo ">>> Start time : $(date)"

IFS=',' read -ra HIFI_BAM_ARR <<< "$HIFI_BAMS"
for F in "${HIFI_BAM_ARR[@]}" "$HIC_R1" "$HIC_R2"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: Input file not found: ${F}"
        exit 1
    fi
done

HAS_UL=false
if [[ -n "$ONT_UL" && "$ONT_UL" != "NA" ]]; then
    if [[ ! -f "$ONT_UL" ]]; then
        echo "ERROR: ONT UL file not found: ${ONT_UL}"
        exit 1
    fi
    HAS_UL=true
fi

# =============================================================================
# STEP 1: HiFiAdapterFilt (per sample)
#
# IMPORTANT (learned the hard way): hifiadapterfilt.sh sets `outdir=$(pwd)` by
# default and uses TWO different variables internally — it converts the BAM to
# FASTQ/FASTA next to the BAM (its `read_path_str`), but hands BLAST a path
# built from `outdir`. If those two directories differ, BLAST silently finds
# no query file, no blocklist is produced, and the "filtered" output is a
# VERBATIM COPY of the input with nothing removed. That failure is silent:
# the .filt.fastq.gz exists, so any existence-based guard happily skips it.
#
# Fix: symlink the BAM into the output directory and run entirely inside it
# (cd there, -o "."), so read_path_str and outdir are the same directory.
#
# NOTE ON DISK: this step converts each BAM to uncompressed FASTQ *and* FASTA.
# With 130-190 GB BAMs that is roughly 400-600 GB of transient scratch per
# sample, on top of the BAM itself. The script removes the FASTA itself; the
# FASTQ and the symlink are cleaned up below once filtering succeeds.
# =============================================================================
echo ">>> Step 1: Adapter filtering raw HiFi BAMs (HiFiAdapterFilt)"

FILT_SAMPLE_DIR="${FILT_DIR}/${SAMPLE}"
mkdir -p "$FILT_SAMPLE_DIR"

FILT_HIFI_ARR=()
NEED_FILTER=false
for RAW_HIFI_BAM in "${HIFI_BAM_ARR[@]}"; do
    PREFIX=$(basename "$RAW_HIFI_BAM" .bam)
    [[ -f "${FILT_SAMPLE_DIR}/${PREFIX}.filt.fastq.gz" ]] || NEED_FILTER=true
done

if [[ "$NEED_FILTER" == true ]]; then
    ml anaconda/2025.12-py313
    conda activate "$HIFIADAPTERFILT_ENV_DIR"

    # --- BLAST database path shim (bioconda packaging bug) ---
    # hifiadapterfilt.sh line 7 derives its BLAST database location with:
    #   DBpath=$(echo $PATH | sed 's/:/\n/g' | grep "HiFiAdapterFilt/DB" | head -n 1)
    # i.e. it greps PATH for the literal string "HiFiAdapterFilt/DB", which is
    # the upstream GitHub layout. bioconda installs the database at
    # <env>/bin/DB (lowercase package name, "bin" not the package name), so
    # that grep matches NOTHING, DBpath ends up empty, and every blastn call
    # becomes `-db /pacbio_vectors_db` -> "BLAST Database error: No alias or
    # index file found". BLAST writes nothing, the blocklist is empty, and
    # NOTHING is filtered -- silently, while the stats file cheerfully reports
    # "0 adapter contaminated ccs reads (0% of total)".
    #
    # DBpath is assigned with a plain `=`, so exporting it here would just be
    # overwritten. Instead give that grep something to match: a symlink whose
    # path literally contains "HiFiAdapterFilt/DB", appended to PATH.
    HIFIADAPTERFILT_DB_SHIM="${CONDA_ENVS_DIR}/HiFiAdapterFilt"
    mkdir -p "$HIFIADAPTERFILT_DB_SHIM"
    ln -sfn "${HIFIADAPTERFILT_ENV_DIR}/bin/DB" "${HIFIADAPTERFILT_DB_SHIM}/DB"
    export PATH="${PATH}:${HIFIADAPTERFILT_DB_SHIM}/DB"

    if ! echo "$PATH" | tr ':' '\n' | grep -q "HiFiAdapterFilt/DB"; then
        echo "ERROR: BLAST database shim not on PATH — HiFiAdapterFilt would"
        echo "silently filter nothing. Expected: ${HIFIADAPTERFILT_DB_SHIM}/DB"
        exit 1
    fi
    if [[ ! -f "${HIFIADAPTERFILT_DB_SHIM}/DB/pacbio_vectors_db.nin" ]]; then
        echo "ERROR: BLAST database not found via shim:"
        echo "  ${HIFIADAPTERFILT_DB_SHIM}/DB/pacbio_vectors_db.nin"
        exit 1
    fi
fi

for RAW_HIFI_BAM in "${HIFI_BAM_ARR[@]}"; do
    PREFIX=$(basename "$RAW_HIFI_BAM" .bam)
    FILT_FASTQ="${FILT_SAMPLE_DIR}/${PREFIX}.filt.fastq.gz"
    BLOCKLIST="${FILT_SAMPLE_DIR}/${PREFIX}.blocklist"
    STATS="${FILT_SAMPLE_DIR}/${PREFIX}.stats"

    if [[ -f "$FILT_FASTQ" ]]; then
        echo "  ${PREFIX}: already filtered — skipping"
    else
        echo "  Filtering ${PREFIX}"
        ln -sf "$RAW_HIFI_BAM" "${FILT_SAMPLE_DIR}/${PREFIX}.bam"
        # Tee stderr to a file so BLAST failures can be detected: the script
        # backgrounds its blastn calls and never checks their exit status, so
        # a database error otherwise just yields an empty blocklist and a
        # clean-looking "0% contaminated" result.
        FILT_STDERR="${FILT_SAMPLE_DIR}/${PREFIX}.hifiadapterfilt.stderr"
        # BLAST_THREADS is deliberately NOT $THREADS. hifiadapterfilt passes -t
        # straight to `blastn -num_threads`, and at high thread counts BLAST
        # hits "CThread::Run() -- error creating thread" and dies PARTWAY
        # through the search. That leaves a truncated .contaminant.blastout, an
        # empty blocklist, and a stats file reporting a confident but wrong
        # "0% contaminated". BLAST also scales poorly past ~8-16 threads, so
        # capping it costs little. The BAM->FASTQ conversion is single-threaded
        # regardless and dominates this step's runtime.
        BLAST_THREADS=8
        (cd "$FILT_SAMPLE_DIR" && hifiadapterfilt.sh -p "$PREFIX" -o "." -t "$BLAST_THREADS") \
            2> >(tee "$FILT_STDERR" >&2)

        # Any NCBI/BLAST failure means the adapter search did not complete, so
        # the "contaminated reads" count cannot be trusted. Fail loudly rather
        # than assembling from reads that were never fully screened.
        if grep -qi "BLAST Database error\|No alias or index file found\|error creating thread\|NCBI C++ Exception" "$FILT_STDERR"; then
            echo "ERROR: BLAST did not complete — adapter screening is incomplete."
            echo "The read counts in ${PREFIX}.stats are NOT trustworthy."
            echo "See ${FILT_STDERR}"
            exit 1
        fi
    fi

    if [[ ! -f "$FILT_FASTQ" ]]; then
        echo "ERROR: HiFiAdapterFilt did not produce expected output: ${FILT_FASTQ}"
        exit 1
    fi

    # Hard validation: the silent-no-op failure mode above produces a
    # .filt.fastq.gz with NO blocklist and NO stats file. Refuse to continue
    # on an unfiltered copy rather than assembling from it.
    if [[ ! -f "$BLOCKLIST" || ! -s "$STATS" ]]; then
        echo "ERROR: HiFiAdapterFilt produced ${FILT_FASTQ} but no blocklist/stats."
        echo "That means BLAST never ran and NOTHING was filtered — the output is"
        echo "an unfiltered copy of the input. Not continuing."
        echo "Expected: ${BLOCKLIST} and ${STATS}"
        exit 1
    fi
    echo "  --- ${PREFIX} filtering summary ---"
    grep -E "Number of (ccs reads|adapter contaminated|ccs reads retained)" "$STATS" || cat "$STATS"

    # Drop the huge uncompressed intermediates and the BAM symlink.
    rm -f "${FILT_SAMPLE_DIR}/${PREFIX}.fastq" "${FILT_SAMPLE_DIR}/${PREFIX}.fasta" \
          "${FILT_SAMPLE_DIR}/${PREFIX}.fq" "${FILT_SAMPLE_DIR}/${PREFIX}.bam"

    FILT_HIFI_ARR+=("$FILT_FASTQ")
done

if [[ "$NEED_FILTER" == true ]]; then
    conda deactivate
    module unload anaconda/2025.12-py313
fi

echo "  Filtered HiFi reads: ${FILT_HIFI_ARR[*]}"

# All anaconda activity for this run is now finished — safe to load these.
ml ragtag
ml liftoff
ml minimap2

# =============================================================================
# STEP 2: hifiasm (per sample)
# =============================================================================
echo ">>> Step 2: hifiasm assembly"

OUT_PREFIX="${ASM_DIR}/${SAMPLE}"

if [[ -f "${OUT_PREFIX}.hic.hap1.p_ctg.gfa" && -f "${OUT_PREFIX}.hic.hap2.p_ctg.gfa" ]]; then
    echo "  hifiasm hap1/hap2 GFAs already exist — skipping"
else
    HIFIASM_CMD=(hifiasm -o "$OUT_PREFIX" -t "$THREADS" --h1 "$HIC_R1" --h2 "$HIC_R2")
    if [[ "$HAS_UL" == true ]]; then
        echo "  Ultra-long ONT reads detected — adding --ul"
        HIFIASM_CMD+=(--ul "$ONT_UL")
    fi
    HIFIASM_CMD+=("${FILT_HIFI_ARR[@]}")
    echo "  ${HIFIASM_CMD[*]}"
    "${HIFIASM_CMD[@]}"
fi

# =============================================================================
# STEPS 3-11: per haplotype
# =============================================================================
for HAP in hap1 hap2; do
    echo ""
    echo ">>> ===== ${SAMPLE} ${HAP} ====="

    # ---------------------------------------------------------------- Step 3
    echo ">>> Step 3 (${HAP}): GFA -> FASTA"
    GFA="${OUT_PREFIX}.hic.${HAP}.p_ctg.gfa"
    CONTIGS="${ASM_DIR}/${SAMPLE}.${HAP}.contigs.fa"

    if [[ ! -f "$GFA" ]]; then
        echo "ERROR: Expected hifiasm output not found: ${GFA}"
        exit 1
    fi
    if [[ -f "$CONTIGS" && -f "${CONTIGS}.fai" ]]; then
        echo "  contigs FASTA already exists — skipping"
    else
        awk '/^S/{print ">"$2"\n"$3}' "$GFA" > "$CONTIGS"
        samtools faidx "$CONTIGS"
    fi

    # --------------------------------------------------------------- Step 3b
    echo ">>> Step 3b (${HAP}): QUAST — post-hifiasm contigs"
    QUAST_HIFIASM_OUT="${QC_DIR}/quast/${SAMPLE}.${HAP}.post_hifiasm"
    if [[ -f "${QUAST_HIFIASM_OUT}/report.txt" ]]; then
        echo "  already exists — skipping"
    else
        quast.py -o "$QUAST_HIFIASM_OUT" -t "$THREADS" --large "$CONTIGS"
    fi

    # ---------------------------------------------------------------- Step 4
    echo ">>> Step 4 (${HAP}): Align Hi-C reads to contigs (MAPQ >= ${HIC_MIN_MAPQ})"
    if [[ -f "${CONTIGS}.bwt" ]]; then
        echo "  bwa index already exists — skipping"
    else
        bwa index "$CONTIGS"
    fi

    # Stays BAM: yahs reads it directly and has no CRAM support.
    HIC_BAM="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}.hic2contigs.bam"
    if [[ -f "$HIC_BAM" ]]; then
        echo "  Hi-C BAM already exists — skipping"
    else
        bwa mem -5SP -t "$THREADS" "$CONTIGS" "$HIC_R1" "$HIC_R2" \
            | samtools view -@ "$THREADS" -buS -q "$HIC_MIN_MAPQ" - \
            | samtools sort -@ "$THREADS" -m 1G -n -o "${HIC_BAM}.part" -
        mv "${HIC_BAM}.part" "$HIC_BAM"
    fi

    # ---------------------------------------------------------------- Step 5
    echo ">>> Step 5 (${HAP}): yahs Hi-C scaffolding"
    YAHS_PREFIX="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}"
    SCAFFOLDS="${YAHS_PREFIX}_scaffolds_final.fa"
    if [[ -f "$SCAFFOLDS" ]]; then
        echo "  yahs scaffolds already exist — skipping"
    else
        "$YAHS_BIN" "$CONTIGS" "$HIC_BAM" -o "$YAHS_PREFIX"
    fi
    if [[ ! -f "$SCAFFOLDS" ]]; then
        echo "ERROR: yahs did not produce expected output: ${SCAFFOLDS}"
        exit 1
    fi

    # --------------------------------------------------------------- Step 5b
    echo ">>> Step 5b (${HAP}): QUAST — post-yahs scaffolds"
    QUAST_YAHS_OUT="${QC_DIR}/quast/${SAMPLE}.${HAP}.post_yahs"
    if [[ -f "${QUAST_YAHS_OUT}/report.txt" ]]; then
        echo "  already exists — skipping"
    else
        quast.py -o "$QUAST_YAHS_OUT" -t "$THREADS" --large "$SCAFFOLDS"
    fi

    # ---------------------------------------------------------------- Step 6
    echo ">>> Step 6 (${HAP}): Juicebox .hic/.assembly from raw yahs scaffolding"
    JBAT_PREFIX="${SCAFFOLD_DIR}/${SAMPLE}.${HAP}_JBAT"
    if [[ -f "${JBAT_PREFIX}.hic" && -f "${JBAT_PREFIX}.assembly" ]]; then
        echo "  Juicebox files already exist — skipping"
    else
        "$JUICER_BIN" pre -a -o "$JBAT_PREFIX" \
            "${YAHS_PREFIX}.bin" "${YAHS_PREFIX}_scaffolds_final.agp" "${CONTIGS}.fai" \
            > "${JBAT_PREFIX}.log" 2>&1

        # -n skips normalization-vector calculation. juicer_tools 2.x throws a
        # NullPointerException in that stage (AddNorm/NormalizationCalculations)
        # AFTER the .hic body is already written, and JBAT curation does not
        # need those vectors. Failure here is non-fatal: this file is only for
        # optional manual curation, so it must not abort the whole sample.
        "${JUICER_TOOLS_CMD[@]}" pre -n \
            "${JBAT_PREFIX}.txt" "${JBAT_PREFIX}.hic.part" \
            <(grep PRE_C_SIZE "${JBAT_PREFIX}.log" | awk '{print $2" "$3}') \
            || echo "  WARNING: juicer_tools pre failed — continuing without the Juicebox .hic"

        if [[ -s "${JBAT_PREFIX}.hic.part" ]]; then
            mv "${JBAT_PREFIX}.hic.part" "${JBAT_PREFIX}.hic"
        else
            echo "  WARNING: no .hic produced for ${SAMPLE} ${HAP}; JBAT curation unavailable"
        fi
    fi
    echo "  ${JBAT_PREFIX}.hic + ${JBAT_PREFIX}.assembly (open together in Juicebox)"

    # ---------------------------------------------------------------- Step 7
    echo ">>> Step 7 (${HAP}): RagTag vs chicken (GRCg7b)"
    RAGTAG_OUT="${RAGTAG_DIR}/${SAMPLE}.${HAP}"
    PSEUDO_CHR="${RAGTAG_OUT}/ragtag.scaffold.fasta"
    if [[ -f "$PSEUDO_CHR" ]]; then
        echo "  RagTag output already exists — skipping"
    else
        ragtag.py scaffold -o "$RAGTAG_OUT" -t "$THREADS" -u "$CHICKEN_FASTA" "$SCAFFOLDS"
    fi
    if [[ ! -f "$PSEUDO_CHR" ]]; then
        echo "ERROR: RagTag did not produce expected output: ${PSEUDO_CHR}"
        exit 1
    fi

    # --------------------------------------------------------------- Step 7b
    echo ">>> Step 7b (${HAP}): Rename to chicken chromosome names, strip _RagTag"
    FINAL_FASTA="${FINAL_DIR}/${SPECIES}_${SAMPLE}_${HAP}.pseudo_chr.fasta"

    # Guard requires BOTH that renaming happened AND that the length filter ran
    # (the .unplaced_short.fasta companion exists). An assembly produced before
    # the filter was added therefore regenerates rather than being skipped.
    ALREADY_RENAMED=false
    if [[ -f "$FINAL_FASTA" ]] \
        && [[ -f "${FINAL_DIR}/${SPECIES}_${SAMPLE}_${HAP}.unplaced_short.fasta" ]] \
        && grep -qm1 '^>chr_\|^>scaffold_' "$FINAL_FASTA" \
        && ! grep -qm1 '_RagTag' "$FINAL_FASTA"; then
        ALREADY_RENAMED=true
    fi

    if [[ "$ALREADY_RENAMED" == true ]]; then
        echo "  already renamed — skipping"
    else
        samtools faidx "$PSEUDO_CHR"

        RENAME_MAP="${RAGTAG_OUT}/${SAMPLE}.${HAP}.rename_map.tsv"
        awk -F'\t' '$1 ~ /_RagTag$/ {print $1"\t"$2}' "${PSEUDO_CHR}.fai" \
            | sort -k2,2 -nr \
            | awk -F'\t' -v chrmap="$CHICKEN_CHR_MAP" '
                BEGIN {
                    while ((getline line < chrmap) > 0) {
                        split(line, a, "\t")
                        chrname[a[1]] = a[2]
                    }
                }
                {
                    acc = $1
                    sub(/_RagTag$/, "", acc)
                    if (acc in chrname) {
                        print $1"\tchr_"chrname[acc]
                    } else {
                        scafn++
                        print $1"\tscaffold_"scafn
                    }
                }
            ' > "$RENAME_MAP"

        awk -v map="$RENAME_MAP" '
            BEGIN {
                while ((getline line < map) > 0) {
                    split(line, a, "\t")
                    newname[a[1]] = a[2]
                }
            }
            /^>/ {
                split(substr($0, 2), parts, " ")
                name = parts[1]
                if (name in newname) {
                    print ">" newname[name]
                } else {
                    sub(/_RagTag$/, "", name)
                    print ">" name
                }
                next
            }
            { print }
        ' "$PSEUDO_CHR" > "${FINAL_FASTA}.prefilter"
        echo "  Rename map: ${RENAME_MAP}"

        # --- Minimum-length filter on UNPLACED scaffolds only ---
        # Every chr_* sequence is kept regardless of length (microchromosomes
        # are legitimately small). scaffold_N sequences below the threshold go
        # to a separate file rather than being discarded, so nothing is lost
        # and the split is auditable.
        SHORT_FASTA="${FINAL_DIR}/${SPECIES}_${SAMPLE}_${HAP}.unplaced_short.fasta"
        samtools faidx "${FINAL_FASTA}.prefilter"
        awk -F'\t' -v min="$MIN_UNPLACED_SCAFFOLD_LEN" \
            '$1 ~ /^chr_/ || $2 >= min {print $1}' "${FINAL_FASTA}.prefilter.fai" > "${FINAL_FASTA}.keep.txt"
        awk -F'\t' -v min="$MIN_UNPLACED_SCAFFOLD_LEN" \
            '$1 !~ /^chr_/ && $2 < min {print $1}' "${FINAL_FASTA}.prefilter.fai" > "${FINAL_FASTA}.short.txt"

        if [[ -s "${FINAL_FASTA}.keep.txt" ]]; then
            samtools faidx "${FINAL_FASTA}.prefilter" -r "${FINAL_FASTA}.keep.txt" > "$FINAL_FASTA"
        else
            echo "ERROR: length filter would keep nothing — check ${FINAL_FASTA}.prefilter"
            exit 1
        fi
        if [[ -s "${FINAL_FASTA}.short.txt" ]]; then
            samtools faidx "${FINAL_FASTA}.prefilter" -r "${FINAL_FASTA}.short.txt" > "$SHORT_FASTA"
        else
            : > "$SHORT_FASTA"
        fi

        echo "  Kept    : $(grep -c '^>' "$FINAL_FASTA") sequences (all chr_* plus unplaced >= ${MIN_UNPLACED_SCAFFOLD_LEN} bp)"
        echo "  Set aside: $(wc -l < "${FINAL_FASTA}.short.txt") short unplaced scaffolds -> ${SHORT_FASTA}"

        rm -f "${FINAL_FASTA}.prefilter" "${FINAL_FASTA}.prefilter.fai" \
              "${FINAL_FASTA}.keep.txt" "${FINAL_FASTA}.short.txt"

        # Liftoff caches a minimap2 index (<target>.mmi) by path, not content —
        # clear it (and its intermediates) whenever FINAL_FASTA is rewritten.
        rm -f "${FINAL_FASTA}.mmi"
        rm -rf "${LIFTOFF_DIR}/${SAMPLE}.${HAP}_intermediate"
    fi

    samtools faidx "$FINAL_FASTA"
    echo "  Final assembly: ${FINAL_FASTA}"

    # --------------------------------------------------------------- Step 7c
    echo ">>> Step 7c (${HAP}): QUAST — final assembly"
    QUAST_RAGTAG_OUT="${QC_DIR}/quast/${SAMPLE}.${HAP}.post_ragtag"
    if [[ -f "${QUAST_RAGTAG_OUT}/report.txt" ]]; then
        echo "  already exists — skipping"
    else
        quast.py -o "$QUAST_RAGTAG_OUT" -t "$THREADS" --large "$FINAL_FASTA"
    fi

    # --------------------------------------------------------------- Step 7d
    echo ">>> Step 7d (${HAP}): tidk telomere profiling (named chromosomes only)"
    TIDK_PREFIX="${SAMPLE}.${HAP}"
    TIDK_CHR_FASTA="${QC_DIR}/tidk/${TIDK_PREFIX}.chr_only.fasta"
    TIDK_TSV="${QC_DIR}/tidk/${TIDK_PREFIX}_telomeric_repeat_windows.tsv"
    TIDK_PLOT="${QC_DIR}/tidk/${TIDK_PREFIX}.svg"

    if [[ -f "$TIDK_PLOT" && "$TIDK_PLOT" -nt "$FINAL_FASTA" ]]; then
        echo "  tidk output already exists and is current — skipping"
    else
        awk -F'\t' '$1 ~ /^chr_/ {print $1}' "${FINAL_FASTA}.fai" > "${TIDK_CHR_FASTA}.names"
        if [[ ! -s "${TIDK_CHR_FASTA}.names" ]]; then
            echo "  WARNING: no chr_ sequences in ${FINAL_FASTA} — skipping tidk"
            rm -f "${TIDK_CHR_FASTA}.names"
        else
            samtools faidx "$FINAL_FASTA" -r "${TIDK_CHR_FASTA}.names" > "$TIDK_CHR_FASTA"
            rm -f "${TIDK_CHR_FASTA}.names"
            "$TIDK_BIN" search --string TTAGGG --output "$TIDK_PREFIX" --dir "${QC_DIR}/tidk" "$TIDK_CHR_FASTA"
            "$TIDK_BIN" plot --tsv "$TIDK_TSV" --output "${QC_DIR}/tidk/${TIDK_PREFIX}"
            echo "  Telomere windows: ${TIDK_TSV}"
            echo "  Telomere plot   : ${TIDK_PLOT}"
        fi
    fi

    # ---------------------------------------------------------------- Step 8
    echo ">>> Step 8 (${HAP}): Liftoff chicken annotation"
    LIFTOFF_GFF="${FINAL_DIR}/${SPECIES}_${SAMPLE}_${HAP}.liftoff.gff3"
    LIFTOFF_UNMAPPED="${LIFTOFF_DIR}/${SAMPLE}.${HAP}.unmapped_features.txt"
    LIFTOFF_INTERMEDIATE="${LIFTOFF_DIR}/${SAMPLE}.${HAP}_intermediate"

    if [[ -f "$LIFTOFF_GFF" && "$LIFTOFF_GFF" -nt "$FINAL_FASTA" ]]; then
        echo "  Liftoff output already exists and is current — skipping"
    else
        mkdir -p "$LIFTOFF_INTERMEDIATE"
        liftoff \
            -g "$CHICKEN_GFF" \
            -o "$LIFTOFF_GFF" \
            -u "$LIFTOFF_UNMAPPED" \
            -dir "$LIFTOFF_INTERMEDIATE" \
            -p "$THREADS" \
            "$FINAL_FASTA" \
            "$CHICKEN_FASTA"
    fi
    echo "  Liftoff annotation: ${LIFTOFF_GFF}"

    # ---------------------------------------------------------------- Step 9
    echo ">>> Step 9 (${HAP}): Orientation dotplot vs chicken (independent minimap2)"
    PAF="${QC_DIR}/${SAMPLE}.${HAP}.vs_chicken.paf"
    DOTPLOT="${QC_DIR}/${SAMPLE}.${HAP}.orientation_dotplot.png"

    if [[ -f "$DOTPLOT" && "$DOTPLOT" -nt "$FINAL_FASTA" ]]; then
        echo "  dotplot already exists and is current — skipping"
    else
        minimap2 -x asm20 -t "$THREADS" "$CHICKEN_FASTA" "$FINAL_FASTA" > "$PAF"
        "$DOTPLOT_PYTHON_BIN" - "$PAF" "$DOTPLOT" "${SPECIES} ${SAMPLE} ${HAP} vs chicken (GRCg7b)" << 'PYEOF'
import math, sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

paf_path, out_path, title = sys.argv[1:4]

blocks_by_query = defaultdict(list)
matched_bases = defaultdict(lambda: defaultdict(int))

with open(paf_path) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        qname, qstart, qend, strand = f[0], int(f[2]), int(f[3]), f[4]
        tname, tstart, tend = f[5], int(f[7]), int(f[8])
        match_len = int(f[9])
        if match_len < 1000:
            continue
        blocks_by_query[qname].append((qstart, qend, tstart, tend, strand, tname))
        matched_bases[qname][tname] += match_len

if not blocks_by_query:
    print("No alignments above length threshold; skipping dotplot.")
    sys.exit(0)

primary_target = {q: max(tb.items(), key=lambda kv: kv[1])[0] for q, tb in matched_bases.items()}
queries = sorted(blocks_by_query.keys())

ncols = min(6, len(queries))
nrows = math.ceil(len(queries) / ncols)
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), squeeze=False)

for i, qname in enumerate(queries):
    ax = axes[i // ncols][i % ncols]
    tgt = primary_target[qname]
    for qstart, qend, tstart, tend, strand, tname in blocks_by_query[qname]:
        if tname != tgt:
            continue
        color = 'tab:blue' if strand == '+' else 'tab:red'
        ys = (tstart, tend) if strand == '+' else (tend, tstart)
        ax.plot([qstart, qend], ys, color=color, linewidth=1.5)
    ax.set_title(qname, fontsize=7)
    ax.set_xlabel(f"vs {tgt}", fontsize=6)
    ax.tick_params(labelsize=5)

for j in range(len(queries), nrows * ncols):
    axes[j // ncols][j % ncols].axis('off')

fig.suptitle(f"{title}\nblue = + strand (expected orientation), red = - strand (inverted)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(out_path, dpi=150)
print(f"Wrote {out_path}")
PYEOF
    fi
    echo "  Orientation dotplot: ${DOTPLOT}"

    # --------------------------------------------------------------- Step 10
    echo ">>> Step 10 (${HAP}): Pretext Hi-C contact map on final assembly"
    FINAL_HIC_CRAM="${QC_DIR}/${SAMPLE}.${HAP}.hic2final.cram"
    PRETEXT_MAP="${QC_DIR}/${SAMPLE}.${HAP}.pretext"

    if [[ -f "$PRETEXT_MAP" && "$PRETEXT_MAP" -nt "$FINAL_FASTA" ]]; then
        echo "  Pretext map already exists and is current — skipping"
    else
        if [[ ! -f "${FINAL_FASTA}.bwt" || "$FINAL_FASTA" -nt "${FINAL_FASTA}.bwt" ]]; then
            bwa index "$FINAL_FASTA"
        fi
        # CRAM is fine here: only ever read back via `samtools view`.
        bwa mem -5SP -t "$THREADS" "$FINAL_FASTA" "$HIC_R1" "$HIC_R2" \
            | samtools view -@ "$THREADS" -buS -q "$HIC_MIN_MAPQ" - \
            | samtools sort -@ "$THREADS" -m 1G --output-fmt cram --reference "$FINAL_FASTA" -o "$FINAL_HIC_CRAM" -
        samtools index "$FINAL_HIC_CRAM"

        samtools view -h -T "$FINAL_FASTA" "$FINAL_HIC_CRAM" \
            | "$PRETEXTMAP_BIN" -o "$PRETEXT_MAP" --sortby length --sortorder descend --mapq "$HIC_MIN_MAPQ"
        "$PRETEXTSNAPSHOT_BIN" --map "$PRETEXT_MAP" --sequences "=full" \
            --prefix "${SAMPLE}.${HAP}." --folder "$QC_DIR"
    fi
    echo "  Hi-C contact map: ${QC_DIR}/${SAMPLE}.${HAP}.*.png"

    # --------------------------------------------------------------- Step 11
    echo ">>> Step 11 (${HAP}): HiFi depth check — chr_Z/chr_W vs autosomes"
    DEPTH_BAM="${QC_DIR}/depth_check/${SAMPLE}.${HAP}.hifi2final.bam"
    DEPTH_COVERAGE_TSV="${QC_DIR}/depth_check/${SAMPLE}.${HAP}.coverage.tsv"

    if [[ -f "$DEPTH_BAM" && "$DEPTH_BAM" -nt "$FINAL_FASTA" ]]; then
        echo "  HiFi-to-final alignment already exists and is current — skipping"
    else
        minimap2 -ax map-hifi -t "$THREADS" "$FINAL_FASTA" "${FILT_HIFI_ARR[@]}" \
            | samtools sort -@ "$THREADS" -m 1G -o "${DEPTH_BAM}.part" -
        mv "${DEPTH_BAM}.part" "$DEPTH_BAM"
        samtools index "$DEPTH_BAM"
    fi

    samtools coverage "$DEPTH_BAM" > "$DEPTH_COVERAGE_TSV"
    echo "  Full per-sequence coverage table: ${DEPTH_COVERAGE_TSV}"

    echo "  --- chr_Z / chr_W / chr_1-5 ---"
    awk -F'\t' 'NR==1 || $1 ~ /^chr_(Z|W|1|2|3|4|5)$/' "$DEPTH_COVERAGE_TSV" | column -t

    # Baseline = mean depth across all numbered autosomes. Note this is
    # PER-HAPLOTYPE depth, so a real single-copy sequence (e.g. a female's W)
    # is expected near ~1.0x of this baseline, not 0.5x.
    echo "  --- Z/W vs autosome ratio ---"
    awk -F'\t' '
        NR==1 { next }
        $1 ~ /^chr_[0-9]+$/ { auto_sum += $7; auto_n++ }
        $1 == "chr_Z" { z_depth = $7 }
        $1 == "chr_W" { w_depth = $7 }
        END {
            if (auto_n == 0) { print "  No autosomes found to baseline against."; exit }
            auto_mean = auto_sum / auto_n
            printf "  Autosome mean depth (n=%d): %.2fx\n", auto_n, auto_mean
            if (z_depth != "") printf "  chr_Z: %.2fx  (ratio to autosome: %.2f)\n", z_depth, z_depth/auto_mean
            if (w_depth != "") printf "  chr_W: %.2fx  (ratio to autosome: %.2f)\n", w_depth, w_depth/auto_mean
        }
    ' "$DEPTH_COVERAGE_TSV"
done

echo ""
echo ">>> QUAST reports : ${QC_DIR}/quast/"
echo ">>> tidk plots    : ${QC_DIR}/tidk/"
echo ">>> Depth tables  : ${QC_DIR}/depth_check/"
echo ">>> Sample ${SAMPLE} (${SPECIES}) complete — $(date)"
