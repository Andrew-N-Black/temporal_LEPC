#!/bin/bash
# ============================================================================
# make_siteclass.sh — build siteclass.txt (<CHROM> <POS> <deleterious|neutral>)
# from a SnpEff-annotated VCF, optionally requiring GERP RS support.
#
# Classification uses the FIRST (most severe) ANN entry per site:
#   deleterious : Annotation_Impact == HIGH, or missense_variant
#                 (+ GERP RS >= -m if a GERP table is given; sites without a
#                 GERP score are dropped from the deleterious class)
#   neutral     : synonymous_variant
#                 (-n synonymous+intergenic also adds intergenic_region)
#
# Usage:
#   make_siteclass.sh -a snpeff.ann.vcf.gz -o siteclass.txt \
#       [-g gerp_sites.tsv(.gz)] [-m 2] [-n synonymous|synonymous+intergenic]
#
# GERP table: tab-delimited <CHROM> <POS> <RS>, one row per site (header OK,
#   non-matching rows are ignored). CHROM names must match the VCF.
# Requires: bcftools on PATH.
# ============================================================================
set -euo pipefail

ANN_VCF=""; OUT=""; GERP_TSV=""; GERP_MIN=2; NEUTRAL_MODE="synonymous"
while getopts "a:o:g:m:n:" flag; do
  case "$flag" in
    a) ANN_VCF="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    g) GERP_TSV="$OPTARG" ;;
    m) GERP_MIN="$OPTARG" ;;
    n) NEUTRAL_MODE="$OPTARG" ;;
    *) echo "Bad option"; exit 1 ;;
  esac
done
[[ -n "$ANN_VCF" && -n "$OUT" ]] || { echo "Usage: $0 -a ann.vcf.gz -o siteclass.txt [-g gerp.tsv] [-m 2] [-n synonymous]"; exit 1; }
[[ -s "$ANN_VCF" ]] || { echo "ERROR: annotated VCF not found: $ANN_VCF"; exit 1; }
[[ "$NEUTRAL_MODE" == "synonymous" || "$NEUTRAL_MODE" == "synonymous+intergenic" ]] \
  || { echo "ERROR: -n must be synonymous or synonymous+intergenic"; exit 1; }
bcftools view -h "$ANN_VCF" | grep -q '^##INFO=<ID=ANN' \
  || { echo "ERROR: $ANN_VCF has no INFO/ANN field (not SnpEff-annotated?)"; exit 1; }

TMP="${OUT}.snpeff.tmp"
INCLUDE_INTERGENIC=0
[[ "$NEUTRAL_MODE" == "synonymous+intergenic" ]] && INCLUDE_INTERGENIC=1

echo "[siteclass] Classifying sites from SnpEff ANN: $ANN_VCF"
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$ANN_VCF" \
| awk -F'\t' -v OFS='\t' -v ig="$INCLUDE_INTERGENIC" '
    $3 == "." { next }
    {
      split($3, a, ",")          # first ANN entry = most severe (SnpEff sorts)
      split(a[1], f, "|")        # Allele|Annotation|Impact|...
      eff = f[2]; imp = f[3]; cls = ""
      if (imp == "HIGH" || eff ~ /missense_variant/)            cls = "deleterious"
      else if (eff ~ /^synonymous_variant/)                      cls = "neutral"
      else if (ig == 1 && eff ~ /^intergenic_region/)            cls = "neutral"
      if (cls != "") print $1, $2, cls
    }' > "$TMP"

if [[ -n "$GERP_TSV" ]]; then
  [[ -s "$GERP_TSV" ]] || { echo "ERROR: GERP table not found: $GERP_TSV"; exit 1; }
  echo "[siteclass] Requiring GERP RS >= $GERP_MIN for deleterious sites: $GERP_TSV"
  if [[ "$GERP_TSV" == *.gz ]]; then GERP_CAT="zcat"; else GERP_CAT="cat"; fi
  # Load classified sites (small), stream the GERP table, keep RS for matches
  awk -F'\t' -v OFS='\t' -v min="$GERP_MIN" '
      NR == FNR { cls[$1 SUBSEP $2] = $3; next }
      (($1 SUBSEP $2) in cls) { rs[$1 SUBSEP $2] = $3 + 0 }
      END {
        for (k in cls) {
          if (cls[k] == "deleterious" && (!(k in rs) || rs[k] < min)) continue
          split(k, p, SUBSEP); print p[1], p[2], cls[k]
        }
      }' "$TMP" <($GERP_CAT "$GERP_TSV") \
  | sort -k1,1 -k2,2n > "$OUT"
else
  sort -k1,1 -k2,2n -u "$TMP" > "$OUT"
fi
rm -f "$TMP"

echo "[siteclass] Wrote $OUT"
awk '{n[$3]++} END {for (c in n) printf "[siteclass]   %-12s %d sites\n", c, n[c]}' "$OUT"
