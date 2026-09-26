#!/bin/bash
#SBATCH --job-name=lepc_rohan_parse_auto
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -p cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

# =============================================================================
# rohan_parse_autosomal.sh -- autosomal-only re-parse of ROHan's per-segment
# output into the same length bins rohan_parse.sh uses (100kb-1Mb, >1Mb),
# without re-running ROHan itself.
#
# rohan_parse.sh sums ROH segment lengths from every chromosome in each
# sample's *.mid.hmmrohl.gz, Z scaffolds included. Females are hemizygous
# for Z, so Z-linked sites look artificially homozygous -- biasing fROH
# upward for females specifically, and inconsistent with every other
# per-sample statistic in this repo (heterozygosity_array.sh, run_plink.sh,
# roh_analyses.sh/rohparser.py), which are all autosome-restricted.
#
# This is a SEPARATE script from rohan_parse.sh (not an in-place edit), so
# both the whole-genome and autosomal-only summaries stay available side by
# side. Re-parses the SAME already-computed hmmrohl.gz files, so ROHan
# itself does not need to be re-run.
#
# Does NOT go through analysis/rohparser.py: that script parses
# bcftools-roh's "RG" line format (used by roh_analyses.sh, a completely
# different upstream tool from ROHan), not ROHan's hmmrohl format. This
# stays a self-contained awk parser, matching rohan_parse.sh's own
# original style, rather than bolting a second, differently-shaped input
# onto a parser built for the other pipeline's output.
#
# hmmrohl format (confirmed from actual output):
#   #ROH_ID CHROM BEGIN END ROH_LENGTH VALIDATED_SITES
# One row per called ROH segment; a sample with zero ROH has no data rows
# (header only).
#
# USAGE (run any time after rohan_array.sh has produced
# *.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz for every sample):
#   sbatch rohan_parse_autosomal.sh
# =============================================================================

PROJ=/scratch/gautschi/blackan/GROUSE/old_vs_new
REF=/scratch/gautschi/blackan/GROUSE/old_vs_new/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna
OUT_DIR=$PROJ/results/rohan
ROHAN_WINDOW_SIZE=50000
EXPECTED_N_SAMPLES=20
# Space-separated for the awk arrays below; same two scaffolds used
# throughout this repo (run_plink.sh, remove_Z_scaffolds.sh,
# roh_analyses.sh, run_lepc_relatedness.sh).
Z_SCAFFOLDS="NW_026294758.1 NW_026294813.1"

if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: ${REF}.fai not found." >&2
    exit 1
fi

# Confirm both excluded scaffolds actually exist in the reference (catches
# a typo'd/renamed scaffold ID loudly instead of silently excluding
# nothing).
for c in $Z_SCAFFOLDS; do
    awk -v c="$c" '$1==c {f=1} END {exit !f}' "${REF}.fai" \
        || { echo "ERROR: excluded scaffold '$c' not in ${REF}.fai" >&2; exit 1; }
done

# Autosomal genome length denominator: whole .fai length minus the two Z
# scaffold lengths.
GENOME_LEN=$(awk -v z="$Z_SCAFFOLDS" '
    BEGIN { n = split(z, zarr, " "); for (i = 1; i <= n; i++) ZSET[zarr[i]] = 1 }
    !($1 in ZSET) { sum += $2 }
    END { print sum }
' "${REF}.fai")
echo "  Autosomal genome length (from .fai, excluding ${Z_SCAFFOLDS}): ${GENOME_LEN} bp"

OUT_TSV=$OUT_DIR/rohan_froh_summary_autosomal.win${ROHAN_WINDOW_SIZE}.tsv
echo -e "sample\tfROH_100kb-1Mb_rohan\tfROH_1Mb_rohan\tfROH_total_rohan\tn_segments_short\tn_segments_long\tn_segments_excluded_Z" > "$OUT_TSV"

N_PARSED=0
for HMMROHL in "$OUT_DIR"/*.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz; do
    [[ -e "$HMMROHL" ]] || { echo "ERROR: no *.win${ROHAN_WINDOW_SIZE}.mid.hmmrohl.gz files found in ${OUT_DIR}" >&2; exit 1; }

    SAMPLE=$(basename "$HMMROHL" | sed -E "s/\.win${ROHAN_WINDOW_SIZE}\.mid\.hmmrohl\.gz$//")

    # Same two length classes as rohan_parse.sh: 100kb-1Mb (short) and
    # >1Mb (long); segments below 100kb are excluded from both, matching
    # bcftools roh's own convention. Segments on either Z scaffold are
    # dropped before binning (n_z counts them, reported as
    # n_segments_excluded_Z below, for a sanity check against the whole-
    # genome summary), not counted as "excluded for length".
    zcat "$HMMROHL" 2>/dev/null | awk -v genome="$GENOME_LEN" -v sample="$SAMPLE" -v z="$Z_SCAFFOLDS" '
        BEGIN { n = split(z, zarr, " "); for (i = 1; i <= n; i++) ZSET[zarr[i]] = 1 }
        NR > 1 {
            chrom = $2; len = $5
            if (chrom in ZSET) { n_z++; next }
            if (len >= 100000 && len <= 1000000) { short += len; n_short++ }
            else if (len > 1000000) { long += len; n_long++ }
        }
        END {
            froh_short = short / genome
            froh_long  = long / genome
            froh_total = (short + long) / genome
            printf "%s\t%.9f\t%.9f\t%.9f\t%d\t%d\t%d\n", sample, froh_short, froh_long, froh_total, n_short+0, n_long+0, n_z+0
        }' >> "$OUT_TSV"

    N_PARSED=$((N_PARSED + 1))
done

echo "  Parsed ${N_PARSED} samples"
if [[ "$N_PARSED" -ne "$EXPECTED_N_SAMPLES" ]]; then
    echo "  WARNING: expected ${EXPECTED_N_SAMPLES} samples, found ${N_PARSED} -- check ${OUT_DIR} for missing/failed samples." >&2
fi
echo "  Written: ${OUT_TSV}"
echo ""
column -t "$OUT_TSV" 2>/dev/null || cat "$OUT_TSV"

echo ""
echo ">>> All done."
