#!/bin/bash
#SBATCH --job-name=lepc_roh_parse_auto
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
# roh_parse_autosomal.sh -- self-contained, autosome-only parse of
# bcftools_roh.sh's "RG" output into F(ROH), in the same two length bins
# ROHan/rohan_parse_autosomal.sh uses (100kb-1Mb, >1Mb). Reads the RG lines
# directly and does the binning + F(ROH) arithmetic itself (no external
# parser script, no hardcoded paths), so there's nothing upstream left to
# silently bias the sizes.
#
# Depends only on bcftools_roh.sh having already produced
# ROH_GROUSE_PL_regions.txt -- does not re-run ANGSD or bcftools roh.
#
# RG line format (bcftools roh, tab-separated):
#   RG  Sample  Chromosome  Start  End  Length(bp)  NumMarkers  Quality
#
# Sample-ID note: ANGSD's -dobcf output uses each line of $FINAL_CRAMLIST
# verbatim as the sample identifier -- the full CRAM path, ending in
# ".md.dedup_q20.cram" (confirmed against processing/final_cramlist.txt).
# Strips that suffix, with a bare ".cram" fallback in case the convention
# ever differs, rather than assuming one and risking garbled sample IDs.
#
# Every sample in $FINAL_CRAMLIST gets an output row even with zero called
# ROH regions, seeded before any RG line is read.
#
# USAGE (run any time after bcftools_roh.sh has produced
# ROH_GROUSE_PL_regions.txt):
#   sbatch roh_parse_autosomal.sh
# =============================================================================

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"
ROH_DIR="${PROJECT_DIR}/roh"
ROH_RG_ONLY="${ROH_DIR}/ROH_GROUSE_PL_regions.txt"
OUT_TSV="${ROH_DIR}/roh_froh_summary_autosomal.tsv"

# Same two Z scaffolds used throughout this repo (run_plink.sh,
# remove_Z_scaffolds.sh, heterozygosity_array.sh, ROHan/rohan_parse_autosomal.sh,
# run_lepc_relatedness.sh).
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"
# Minimum mean fwd-bwd phred quality (bcftools roh's RG column 8) to keep a region.
MIN_QUALITY=30

for f in "$REF_FASTA.fai" "$FINAL_CRAMLIST" "$ROH_RG_ONLY"; do
    [[ -f "$f" ]] || { echo "ERROR: required input not found: $f" >&2; exit 1; }
done

# Confirm both excluded scaffolds actually exist in the reference (catches
# a typo'd/renamed scaffold ID loudly instead of silently excluding
# nothing).
for c in ${Z_SCAFFOLDS//,/ }; do
    awk -v c="$c" '$1==c {f=1} END {exit !f}' "${REF_FASTA}.fai" \
        || { echo "ERROR: excluded scaffold '$c' not in ${REF_FASTA}.fai" >&2; exit 1; }
done

# Autosomal genome length denominator: whole .fai length minus the two Z
# scaffold lengths.
GENOME_LEN=$(awk -v z="$Z_SCAFFOLDS" '
    BEGIN { n = split(z, zarr, ","); for (i = 1; i <= n; i++) ZSET[zarr[i]] = 1 }
    !($1 in ZSET) { sum += $2 }
    END { print sum }
' "${REF_FASTA}.fai")
echo "  Autosomal genome length (from .fai, excluding ${Z_SCAFFOLDS}): ${GENOME_LEN} bp"

N_SAMPLES=$(wc -l < "$FINAL_CRAMLIST")

# Clean, expected sample IDs from the cramlist, stripping both the real
# suffix and a bare ".cram" fallback -- same double-strip used below on the
# RG lines' own sample column, so both sides land on identical IDs.
SAMPLE_LIST_CLEAN=$(mktemp)
trap 'rm -f "$SAMPLE_LIST_CLEAN"' EXIT
awk -F'/' '{
    s = $NF
    gsub(/\.md\.dedup_q20\.cram$/, "", s)
    gsub(/\.cram$/, "", s)
    print s
}' "$FINAL_CRAMLIST" > "$SAMPLE_LIST_CLEAN"

echo -e "sample\tfROH_100kb-1Mb\tfROH_1Mb\tfROH_total\tn_segments_short\tn_segments_long\tn_segments_excluded_Z\tn_segments_excluded_qual" > "$OUT_TSV"

# Single pass over the RG-only file (FNR==NR seeds every expected sample
# at zero from SAMPLE_LIST_CLEAN before any RG line is read), binning into
# the same two length classes ROHan/rohan_parse_autosomal.sh uses.
awk -v genome="$GENOME_LEN" -v z="$Z_SCAFFOLDS" -v minq="$MIN_QUALITY" '
    BEGIN { n = split(z, zarr, ","); for (i = 1; i <= n; i++) ZSET[zarr[i]] = 1 }
    FNR == NR { seen[$1] = 1; next }
    {
        sample = $2; chrom = $3; len = $6; qual = $8
        gsub(/^.*\//, "", sample)
        gsub(/\.md\.dedup_q20\.cram$/, "", sample)
        gsub(/\.cram$/, "", sample)
        seen[sample] = 1
        if (chrom in ZSET) { n_z[sample]++; next }
        if (qual < minq)   { n_q[sample]++; next }
        if (len >= 100000 && len <= 1000000) { short[sample] += len; n_short[sample]++ }
        else if (len > 1000000)              { long[sample]  += len; n_long[sample]++ }
    }
    END {
        for (s in seen) {
            sh = short[s] + 0; lg = long[s] + 0
            froh_short = sh / genome
            froh_long  = lg / genome
            froh_total = (sh + lg) / genome
            printf "%s\t%.9f\t%.9f\t%.9f\t%d\t%d\t%d\t%d\n", \
                s, froh_short, froh_long, froh_total, \
                n_short[s] + 0, n_long[s] + 0, n_z[s] + 0, n_q[s] + 0
        }
    }
' "$SAMPLE_LIST_CLEAN" "$ROH_RG_ONLY" | sort -k1,1 >> "$OUT_TSV"

N_PARSED=$(($(wc -l < "$OUT_TSV") - 1))
echo "  Parsed ${N_PARSED} samples (expected ${N_SAMPLES})"
if [[ "$N_PARSED" -ne "$N_SAMPLES" ]]; then
    echo "  WARNING: expected ${N_SAMPLES} samples, found ${N_PARSED} -- check ${ROH_RG_ONLY} and ${FINAL_CRAMLIST} for a sample-ID mismatch." >&2
fi
echo "  Written: ${OUT_TSV}"
echo ""
column -t "$OUT_TSV" 2>/dev/null || cat "$OUT_TSV"

echo ""
echo ">>> All done."
