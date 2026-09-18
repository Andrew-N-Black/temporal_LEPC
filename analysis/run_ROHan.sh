#!/bin/bash
#SBATCH --job-name=lepc_rohan
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

# =============================================================================
# Run ROHan (Renaud et al. 2019) per sample to estimate heterozygosity and
# runs of homozygosity directly from BAM/CRAM -- NOT from a VCF. ROHan
# computes its own genotype likelihoods internally (Bayesian model + HMM),
# so it deliberately works upstream of variant calling rather than on
# already-called genotypes.
#
# Combines the original run_rohan.sh (per-sample ROHan execution) with
# parse_rohan_roh.sh (post-processing the per-segment ROH output into the
# same 100kb-1Mb / >1Mb length classes bcftools roh uses, for direct
# comparison) into a single script -- Step 1 runs ROHan for every sample,
# Step 2 parses all of their outputs once Step 1 is done.
#
# Requires install_rohan.sh to have been run on the login node first.
# =============================================================================

# ---------------- CONFIG -- edit these before running ----------------
PROJ=/scratch/gautschi/blackan/GROUSE/old_vs_new
REF=/scratch/gautschi/blackan/GROUSE/old_vs_new/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna
ROHAN_BIN=$PROJ/tools/ROHan/bin/rohan
GSL_PREFIX=$PROJ/tools/gsl

# Directory + filename pattern for the input alignments. Defaults to the
# naming convention already used elsewhere in this project
# (*.md.dedup_q20.cram) -- change CRAM_DIR/CRAM_PATTERN if your actual
# files for this analysis live somewhere else or are named differently.
CRAM_DIR=/scratch/gautschi/blackan/GROUSE/old_vs_new/crams
CRAM_PATTERN="*.md.dedup_q20.cram"

OUT_DIR=$PROJ/results/rohan
THREADS=16

# ROHan's expected within-ROH heterozygosity rate parameter. The default in
# ROHan's own examples is on the order of 2e-5; adjust if your species'
# expected mutation rate differs substantially. Confirm the flag name is
# still --rohmu with `rohan --help` after building (see install_rohan.sh).
ROHMU=1e-5

# ROHan's HMM operates on the genome tiled into fixed windows (--size,
# confirmed via `rohan --help`); at the 1 Mb default it cannot resolve
# anything in the 100kb-1Mb class at all -- confirmed empirically from the
# first run: ROH segment lengths came back as exact multiples of 1,000,000
# across every sample, and near-zero total ROH, while bcftools roh found
# substantial 100kb-1Mb signal in every one of the same individuals. 50 kb
# gives several windows of resolution within that bin (enough to
# distinguish it from the >1Mb class) without pushing per-window
# segregating-site counts so low that local theta estimates become
# unreliably noisy.
ROHAN_WINDOW_SIZE=50000

mkdir -p "$OUT_DIR" logs

if [ ! -x "$ROHAN_BIN" ]; then
    echo "ERROR: ROHan binary not found at $ROHAN_BIN" >&2
    echo "  Run install_rohan.sh on the login node first." >&2
    exit 1
fi

module --force purge
module load gcc/14.1.0
module load biocontainers
module load samtools

# rohan was linked against a custom-built GSL (no 'gsl' module exists on
# Gautschi), so it needs this to find libgsl.so at runtime, not just at
# build time -- see install_rohan.sh.
export LD_LIBRARY_PATH="$GSL_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# =============================================================================
# STEP 1: Run ROHan per sample
# =============================================================================
echo ">>> Step 1: ROHan per-sample analysis"

for CRAM in "$CRAM_DIR"/$CRAM_PATTERN; do
    SAMPLE=$(basename "$CRAM" | sed -E 's/\.md\.dedup_q20\.cram$//')

    # Window size baked into the output prefix -- an earlier 1Mb-window
    # run's files (${SAMPLE}.hEst etc.) are left alone rather than silently
    # overwritten, since they used a resolution too coarse to see the
    # 100kb-1Mb ROH class at all.
    OUT_PREFIX=$OUT_DIR/${SAMPLE}.win${ROHAN_WINDOW_SIZE}
    if [[ -f "${OUT_PREFIX}.hEst.gz" ]]; then
        echo "  ${SAMPLE}: ${OUT_PREFIX}.hEst.gz already exists -- skipping."
        continue
    fi
    echo "=== Processing $SAMPLE ==="

    # ROHan's most reliably supported input format is BAM. Rather than rely
    # on ROHan's own CRAM/reference handling (undocumented in what we could
    # verify), convert to a temporary indexed BAM first -- slower, but
    # removes any ambiguity about reference resolution.
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

    echo "  Done: ${OUT_PREFIX}.* "
done

echo ""
echo ">>> Step 1 complete. Per-sample outputs in $OUT_DIR/"
echo "Each sample produces (filenames per ROHan's own convention):"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.hEst.gz          -- genome-wide heterozygosity estimate + CI, outside ROH"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.mid.hmmp.gz      -- per-window HMM posterior probabilities"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz   -- called ROH segments (chrom, begin, end, length)"
echo "  <sample>.win${ROHAN_WINDOW_SIZE}.summary.txt      -- genome-wide summary stats"

# =============================================================================
# STEP 2: Parse per-segment ROH output into bcftools-roh-comparable bins
#
# hmmrohl format (confirmed from actual output):
#   #ROH_ID CHROM BEGIN END ROH_LENGTH VALIDATED_SITES
# One row per called ROH segment; a sample with zero ROH has no data rows
# (header only).
# =============================================================================
echo ""
echo ">>> Step 2: Parsing ROH segments into 100kb-1Mb / >1Mb bins"

if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: ${REF}.fai not found." >&2
    exit 1
fi
# Genome length denominator for fROH, summed from the reference .fai --
# matches the convention rohparser.py already uses for the bcftools roh
# side, so the two methods' fROH values are computed the same way.
GENOME_LEN=$(awk '{sum+=$2} END{print sum}' "${REF}.fai")
echo "  Genome length (from .fai): ${GENOME_LEN} bp"

OUT_TSV=$OUT_DIR/rohan_froh_summary.win${ROHAN_WINDOW_SIZE}.tsv
echo -e "sample\tfROH_100kb-1Mb_rohan\tfROH_1Mb_rohan\tfROH_total_rohan\tn_segments_short\tn_segments_long" > "$OUT_TSV"

N_PARSED=0
for HMMROHL in "$OUT_DIR"/*.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz; do
    [[ -e "$HMMROHL" ]] || { echo "ERROR: no *.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz files found in ${OUT_DIR}" >&2; exit 1; }

    SAMPLE=$(basename "$HMMROHL" | sed -E "s/\.win${ROHAN_WINDOW_SIZE}\.mid\.hmmrohl\.gz$//")

    # Bin each segment by length into the same two classes bcftools roh
    # uses: 100kb-1Mb (short/background) and >1Mb (long/recent). Segments
    # below 100kb are excluded from both, matching bcftools roh's own
    # convention (its "total" is the sum of just these two classes, not
    # all called segments regardless of length).
    zcat "$HMMROHL" 2>/dev/null | awk -v genome="$GENOME_LEN" -v sample="$SAMPLE" '
        NR > 1 {
            len = $5
            if (len >= 100000 && len <= 1000000) { short += len; n_short++ }
            else if (len > 1000000) { long += len; n_long++ }
        }
        END {
            froh_short = short / genome
            froh_long  = long / genome
            froh_total = (short + long) / genome
            printf "%s\t%.9f\t%.9f\t%.9f\t%d\t%d\n", sample, froh_short, froh_long, froh_total, n_short+0, n_long+0
        }' >> "$OUT_TSV"

    N_PARSED=$((N_PARSED + 1))
done

echo "  Parsed ${N_PARSED} samples"
echo "  Written: ${OUT_TSV}"
echo ""
column -t "$OUT_TSV" 2>/dev/null || cat "$OUT_TSV"

echo ""
echo ">>> All done."
