#!/bin/bash
#SBATCH --job-name=lepc_rohan_parse
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
# STEP 2 (standalone): Parse per-segment ROHan output into bcftools-roh-
# comparable bins, across ALL 20 samples' outputs -- not just the ones the
# array job reprocessed. Depends on 01_rohan_array.sh having finished
# successfully for every task.
#
# hmmrohl format (confirmed from actual output):
#   #ROH_ID CHROM BEGIN END ROH_LENGTH VALIDATED_SITES
# One row per called ROH segment; a sample with zero ROH has no data rows
# (header only).
#
# USAGE (chain after the array job so this only starts once every task
# succeeds -- afterok holds this job if any array task fails):
#   ARRAY_JOBID=$(sbatch --parsable --array=1-$(wc -l < remaining_samples.txt) 01_rohan_array.sh)
#   sbatch --dependency=afterok:${ARRAY_JOBID} 02_rohan_parse.sh
# =============================================================================

PROJ=/scratch/gautschi/blackan/GROUSE/old_vs_new
REF=/scratch/gautschi/blackan/GROUSE/old_vs_new/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna
OUT_DIR=$PROJ/results/rohan
ROHAN_WINDOW_SIZE=50000
EXPECTED_N_SAMPLES=20

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
if [[ "$N_PARSED" -ne "$EXPECTED_N_SAMPLES" ]]; then
    echo "  WARNING: expected ${EXPECTED_N_SAMPLES} samples, found ${N_PARSED} -- check ${OUT_DIR} for missing/failed samples." >&2
fi
echo "  Written: ${OUT_TSV}"
echo ""
column -t "$OUT_TSV" 2>/dev/null || cat "$OUT_TSV"

echo ""
echo ">>> All done."
