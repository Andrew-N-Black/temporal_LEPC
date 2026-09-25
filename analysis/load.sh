#!/bin/bash
# =============================================================================
# LEPC GENETIC LOAD PIPELINE — master script
#
# Generates the SLURM child scripts for:
#   ancestral-allele polarization -> functional annotation (SnpEff) ->
#   site classification -> per-individual genetic load
#
# Run this script to (re)generate scripts/, then submit them in the order
# printed at the end.
#
# Design notes:
#   - Derived alleles are polarized against a chicken-based ancestral
#     sequence (Step 1), not against the reference allele.
#   - Functional impact comes from SnpEff alone. Evolutionary constraint
#     (GERP++ over an 11-taxon galliform Cactus alignment) was built and
#     evaluated, then dropped: that tree totals only ~0.72 substitutions per
#     site, which is too shallow for informative constraint. Per-site RS was
#     effectively binary, and element-level calls (gerpelem) did not enrich
#     missense over synonymous variants. The archived GERP version of this
#     pipeline is kept separately if a deeper alignment is ever built.
#   - Load is reported as total / realized / masked per individual per
#     category, following the standard decomposition.
# =============================================================================
set -euo pipefail

PROJ=/scratch/gautschi/blackan/GROUSE/old_vs_new
OUT=$PROJ/results   # defined early: several config vars below depend on it

REF=$PROJ/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna    # LEPC reference genome
GFF=$PROJ/ref/GCF_026119805.1_pur_lepc_1.0_genomic.gtf    # LEPC gene annotation (GTF)
VCF=/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs/output.subset.biallelic.snps.vcf.gz

# Chicken genome, used in Step 1 to build the ancestral sequence
GALLUS_DIR=$PROJ/ref

# Ancestral sequence in LEPC coordinates (written by Step 1)
ANC=$OUT/ancestral/LEPC_ancestral_from_chicken.fa

THREADS=64

SNPEFF_DB=LEPC_custom                           # Custom SnpEff database name to build
OLD_SAMPLES=$PROJ/sample_lists/old.samples      # Historical/museum-era samples
NEW_SAMPLES=$PROJ/sample_lists/new.samples      # Contemporary samples
ALL_SAMPLES=$PROJ/sample_lists/all.samples      # All 20 samples (old.samples + new.samples)

mkdir -p $PROJ/scripts $PROJ/logs $OUT/{ancestral,snpeff,polarized,load,logs}

# =============================================================================

# =============================================================================
# STEP 0: Download outgroup genome assemblies from NCBI Datasets
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step1_build_ancestral.sh
#!/bin/bash
#SBATCH --job-name=lepc_ancestral
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

module --force purge
module load biocontainers
module load minimap2
module load samtools
module load bcftools
module load htslib
set -euo pipefail

echo "Aligning chicken (Gallus gallus) to LEPC reference: $(date)"

# Auto-detect the chicken FASTA in the pre-downloaded reference directory
# rather than assuming a filename — GALLUS_DIR also holds the LEPC reference
# itself in this setup, so explicitly exclude REF from the candidates
# (otherwise the glob could just as easily match the LEPC genome).
CHICKEN_FASTA=$(find GALLUS_DIR -maxdepth 1 -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" | grep -v "\.gz$" | grep -vxF "REF" | grep -vxF "ANC" | head -n 1)
if [ -z "$CHICKEN_FASTA" ]; then
    CHICKEN_FASTA=$(find GALLUS_DIR -maxdepth 1 \( -iname "*.fa.gz" -o -iname "*.fasta.gz" -o -iname "*.fna.gz" \) | grep -vxF "REF.gz" | head -n 1)
    if [ -n "$CHICKEN_FASTA" ]; then
        echo "Found gzipped chicken FASTA, decompressing: $CHICKEN_FASTA"
        gunzip -k "$CHICKEN_FASTA"
        CHICKEN_FASTA="${CHICKEN_FASTA%.gz}"
    fi
fi
echo "Using chicken reference: $CHICKEN_FASTA"
if [ -z "$CHICKEN_FASTA" ]; then
    echo "ERROR: no chicken FASTA found in GALLUS_DIR — check the path/extension" >&2
    exit 1
fi

# Whole-genome pairwise alignment, chicken -> LEPC coordinates
# asm10 is appropriate for cross-species divergence at this phylogenetic distance;
# switch to asm20 if the alignment rate is low.
minimap2 -ax asm10 -t THREADS REF "$CHICKEN_FASTA" | \
    samtools sort -@ THREADS -o OUT/ancestral/gallus_to_lepc.bam -
samtools index OUT/ancestral/gallus_to_lepc.bam

# Call a consensus base at every LEPC position covered by a unique chicken
# alignment; positions with no alignment or ambiguous (multi-mapping) coverage
# are left as N and dropped downstream at the polarization step.
bcftools mpileup -f REF OUT/ancestral/gallus_to_lepc.bam -Ou | \
    bcftools call -c --ploidy 1 -Oz -o OUT/ancestral/gallus_consensus.vcf.gz
tabix -p vcf OUT/ancestral/gallus_consensus.vcf.gz

bcftools consensus -f REF OUT/ancestral/gallus_consensus.vcf.gz \
    > ANC

samtools faidx ANC

echo "Ancestral FASTA written to: ANC"
echo "Step 1 complete: $(date)"
EOF


# =============================================================================
# STEP 2: SnpEff — build custom database from LEPC annotation, then annotate
#         Produces HIGH/MODERATE/LOW/MODIFIER impact calls per variant
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step2_snpeff_annotate.sh
#!/bin/bash
#SBATCH --job-name=lepc_snpeff
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

module --force purge
module load biocontainers
module load snpeff
module load htslib
module load bcftools
set -euo pipefail

# The module gives you a working "snpEff" command directly (no need to
# manage snpEff.jar/snpEff.config paths by hand), but its install directory
# is shared/read-only, so build the custom genome in our own writable data
# dir and point a local config file at it instead of editing the module's.
DATA_DIR=OUT/snpeff/data
CONFIG_FILE=OUT/snpeff/snpEff.config
mkdir -p ${DATA_DIR}/SNPEFF_DB

cp REF ${DATA_DIR}/SNPEFF_DB/sequences.fa

# Auto-detect gtf vs gff/gff3 from the annotation file's extension -- SnpEff
# needs a different build flag and a different expected filename for each,
# and building a gtf file under the wrong flag silently produces broken or
# empty annotations rather than an obvious error.
ANNOT_FILE=GFF
ANNOT_EXT=$(echo "${ANNOT_FILE##*.}" | tr '[:upper:]' '[:lower:]')
case "$ANNOT_EXT" in
    gtf)
        cp "$ANNOT_FILE" ${DATA_DIR}/SNPEFF_DB/genes.gtf
        BUILD_FLAG="-gtf22"
        ;;
    gff|gff3)
        cp "$ANNOT_FILE" ${DATA_DIR}/SNPEFF_DB/genes.gff
        BUILD_FLAG="-gff3"
        ;;
    *)
        echo "ERROR: unrecognized annotation extension '.${ANNOT_EXT}' -- expected .gtf, .gff, or .gff3" >&2
        exit 1
        ;;
esac
echo "Detected annotation format: .${ANNOT_EXT} -- using SnpEff build flag ${BUILD_FLAG}"

cat > "$CONFIG_FILE" << CFGEOF
data.dir = ${DATA_DIR}
SNPEFF_DB.genome : LEPC_custom
CFGEOF

snpEff build -c "$CONFIG_FILE" $BUILD_FLAG -v SNPEFF_DB -noCheckCds -noCheckProtein

# --- Annotate the joint VCF ---
snpEff -Xmx16g -c "$CONFIG_FILE" -v SNPEFF_DB \
    VCF \
    > OUT/snpeff/lepc_snpeff.vcf

bgzip -f OUT/snpeff/lepc_snpeff.vcf
tabix -p vcf OUT/snpeff/lepc_snpeff.vcf.gz

# Pull out per-impact-class site lists from the ANN field for use downstream
for IMPACT in HIGH MODERATE LOW MODIFIER; do
    bcftools view -i "INFO/ANN[*] ~ '${IMPACT}'" OUT/snpeff/lepc_snpeff.vcf.gz | \
        bcftools query -f '%CHROM\t%POS\n' \
        > OUT/snpeff/sites_${IMPACT}.txt
    echo "  ${IMPACT}: $(wc -l < OUT/snpeff/sites_${IMPACT}.txt) sites"
done

# Amino-acid change per missense variant, for Grantham scoring in Step 5.
# ANN subfields are pipe-separated; subfield 11 is HGVS.p (e.g. p.Ala123Val).
# Only the first (primary) annotation per variant is used.
bcftools view -i "INFO/ANN[*] ~ 'MODERATE'" OUT/snpeff/lepc_snpeff.vcf.gz | \
    bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' | \
    awk -F'\t' 'BEGIN{OFS="\t"} {
        split($3, anns, ",")
        split(anns[1], f, "|")
        if (f[11] != "") print $1, $2, f[11]
    }' > OUT/snpeff/missense_aa_changes.tsv
echo "  amino-acid changes extracted: $(wc -l < OUT/snpeff/missense_aa_changes.tsv)"

echo "Step 2 complete: $(date)"
EOF


# =============================================================================
# STEP 3A: Build the galliform multi-species alignment with Cactus
#          Aligns chicken + all downloaded outgroups to the LEPC reference,
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step4_polarize.sh
#!/bin/bash
#SBATCH --job-name=lepc_polarize
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

module --force purge
module load biocontainers
module load bcftools
module load samtools
module load python3
set -euo pipefail

# Pull the reference/alt alleles for every biallelic SNP in the annotated VCF
bcftools query -f '%CHROM\t%POS\t%RE''F\t%ALT\n' OUT/snpeff/lepc_snpeff.vcf.gz \
    > OUT/polarized/sites_ref_alt.tsv

python3 << 'PYEOF'
import pysam

out_dir = "OUT/polarized"
anc_fa = pysam.FastaFile("ANC")

polarized = []
with open(f"{out_dir}/sites_ref_alt.tsv") as fh:
    for line in fh:
        chrom, pos, ref, alt = line.strip().split("\t")
        pos = int(pos)
        try:
            anc_base = anc_fa.fetch(chrom, pos - 1, pos).upper()
        except (KeyError, ValueError):
            continue

        if anc_base not in ("A", "C", "G", "T"):
            continue  # no confident chicken ortholog at this site

        if anc_base == ref.upper():
            ancestral, derived = ref, alt
        elif anc_base == alt.upper():
            ancestral, derived = alt, ref
        else:
            # Ancestral state matches neither the reference nor alt allele (likely a
            # lineage-specific substitution on the chicken branch, or a
            # third allele) -- exclude from polarized load calculations.
            continue

        polarized.append((chrom, pos, ref, alt, ancestral, derived))

with open(f"{out_dir}/sites_polarized.tsv", "w") as out:
    out.write("chrom\tpos\tref\talt\tancestral\tderived\n")
    for row in polarized:
        out.write("\t".join(map(str, row)) + "\n")

print(f"Polarized {len(polarized)} of the input sites")
PYEOF

echo "Step 4 complete: $(date)"
EOF


# =============================================================================
# STEP 5: Build the deleterious site set
#         Categories: LOF (SnpEff HIGH), MISSENSE (MODERATE), NEUTRAL (LOW)
#         threshold -- intersected with the successfully polarized site set
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step5_deleterious_sites.sh
#!/bin/bash
#SBATCH --job-name=lepc_del_sites
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

module load python3
set -euo pipefail

# =============================================================================
# Classify polarized sites by predicted functional impact (SnpEff only).
#
# No conservation axis: GERP was evaluated and dropped. The 11-taxon
# galliform alignment totals ~0.72 substitutions/site, too shallow for
# informative constraint scores -- per-site RS was effectively binary
# (RS = neutral rate where no substitutions were observed, <= 0 otherwise),
# and element-level calls did not enrich missense over synonymous variants.
# =============================================================================

python3 << 'PYEOF'
import os
import pandas as pd

out = "OUT"

pol = pd.read_csv(f"{out}/polarized/sites_polarized.tsv", sep="\t")
pol["site"] = pol["chrom"].astype(str) + ":" + pol["pos"].astype(str)

def load_site_set(path):
    df = pd.read_csv(path, sep="\t", names=["chrom", "pos"])
    return set(df["chrom"].astype(str) + ":" + df["pos"].astype(str))

high     = load_site_set(f"{out}/snpeff/sites_HIGH.txt")
moderate = load_site_set(f"{out}/snpeff/sites_MODERATE.txt")
low      = load_site_set(f"{out}/snpeff/sites_LOW.txt")

pol["snpeff_HIGH"]     = pol["site"].isin(high)
pol["snpeff_MODERATE"] = pol["site"].isin(moderate)
pol["snpeff_LOW"]      = pol["site"].isin(low)

# Categories carried into per-individual load (Step 6):
#   LOF      = SnpEff HIGH     (stop-gain, frameshift, splice-disrupting)
#   MISSENSE = SnpEff MODERATE (amino-acid changing; severity not ranked)
#   NEUTRAL  = SnpEff LOW      (synonymous; comparison class)
# A site matching more than one class takes the most severe, hence the order.
pol["category"] = "OTHER"
pol.loc[pol["snpeff_LOW"], "category"]      = "NEUTRAL"
pol.loc[pol["snpeff_MODERATE"], "category"] = "MISSENSE"
pol.loc[pol["snpeff_HIGH"], "category"]     = "LOF"

# ---------------------------------------------------------------------------
# Grantham (1974) physicochemical distance for missense variants.
#
# This grades how chemically different the substituted amino acids are
# (composition, polarity, molecular volume). It is NOT a conservation or
# deleteriousness score -- it needs no alignment, so it is unaffected by the
# tree-depth problem that ruled out GERP, but it is correspondingly weaker
# evidence. Report it as a severity gradient, not as a deleterious call.
#
# Distances are computed from Grantham's published formula rather than a
# transcribed matrix; the implementation reproduces published pair values to
# within rounding (e.g. Leu-Ile 4.9 vs 5, Cys-Trp 214.4 vs 215).
# ---------------------------------------------------------------------------
import math

AA_PROPS = {
    "A": (0.00,  8.1,  31.0), "R": (0.65, 10.5, 124.0), "N": (1.33, 11.6,  56.0),
    "D": (1.38, 13.0,  54.0), "C": (2.75,  5.5,  55.0), "Q": (0.89, 10.5,  85.0),
    "E": (0.92, 12.3,  83.0), "G": (0.74,  9.0,   3.0), "H": (0.58, 10.4,  96.0),
    "I": (0.00,  5.2, 111.0), "L": (0.00,  4.9, 111.0), "K": (0.33, 11.3, 119.0),
    "M": (0.00,  5.7, 105.0), "F": (0.00,  5.2, 132.0), "P": (0.39,  8.0,  32.5),
    "S": (1.42,  9.2,  32.0), "T": (0.71,  8.6,  61.0), "W": (0.13,  5.4, 170.0),
    "Y": (0.20,  6.2, 136.0), "V": (0.00,  5.9,  84.0),
}
ALPHA, BETA, GAMMA, RHO = 1.833, 0.1018, 0.000399, 50.723

def grantham(aa1, aa2):
    if aa1 not in AA_PROPS or aa2 not in AA_PROPS:
        return float("nan")
    c1, p1, v1 = AA_PROPS[aa1]
    c2, p2, v2 = AA_PROPS[aa2]
    return RHO * math.sqrt(ALPHA * (c1 - c2) ** 2
                           + BETA * (p1 - p2) ** 2
                           + GAMMA * (v1 - v2) ** 2)

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

import re
HGVS_P = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

aa_path = f"{out}/snpeff/missense_aa_changes.tsv"
gr = {}
n_parsed = n_unparsed = 0
if os.path.exists(aa_path) and os.path.getsize(aa_path) > 0:
    with open(aa_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, pos, hgvs = parts[0], parts[1], parts[2]
            m = HGVS_P.match(hgvs)
            if not m:
                # synonymous (p.Ala123Ala is not emitted as MODERATE),
                # frameshift/extension/unknown notations, etc.
                n_unparsed += 1
                continue
            a1 = THREE_TO_ONE.get(m.group(1))
            a2 = THREE_TO_ONE.get(m.group(3))
            if a1 is None or a2 is None:
                n_unparsed += 1
                continue
            gr[f"{chrom}:{pos}"] = grantham(a1, a2)
            n_parsed += 1
else:
    print("WARNING: no missense_aa_changes.tsv found -- rerun Step 2 to produce it; "
          "Grantham columns will be empty.")

pol["grantham"] = pol["site"].map(gr)

# Grantham's own conventional bins.
def gr_class(d):
    if pd.isna(d):
        return ""
    if d <= 50:
        return "conservative"
    if d <= 100:
        return "moderately_conservative"
    if d <= 150:
        return "moderately_radical"
    return "radical"

pol["grantham_class"] = pol["grantham"].apply(gr_class)

# Finer categories for Step 6: missense split at Grantham 100, the usual
# conservative/radical division. Missense with no parseable amino-acid
# change keeps the unsplit MISSENSE label rather than being guessed at.
pol["category_fine"] = pol["category"]
is_mis = pol["category"] == "MISSENSE"
pol.loc[is_mis & (pol["grantham"] <= 100), "category_fine"] = "MISSENSE_CONSERVATIVE"
pol.loc[is_mis & (pol["grantham"] > 100), "category_fine"]  = "MISSENSE_RADICAL"

pol.to_csv(f"{out}/load/sites_classified.tsv", sep="\t", index=False)

print(f"Amino-acid changes scored: {n_parsed}; unparseable/skipped: {n_unparsed}")
n_mis = int(is_mis.sum())
n_scored = int((is_mis & pol["grantham"].notna()).sum())
print(f"Missense sites: {n_mis}; with a Grantham score: {n_scored} "
      f"({100*n_scored/max(n_mis,1):.1f}%)")
print(pol["category"].value_counts())
print(pol.loc[is_mis, "grantham_class"].value_counts())
PYEOF

echo "Step 5 complete: $(date)"
EOF


# =============================================================================
# STEP 6: Per-individual genetic load (Old vs New)
#         Total load    -- all derived alleles carried (het + hom), any zygosity
#         Realized load -- homozygous-derived only (expressed fitness cost)
#         Masked load   -- heterozygous-derived only (recessive/hidden burden)
#         Computed separately for LOF and DELETERIOUS categories, per sample
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step6_load_per_individual.sh
#!/bin/bash
#SBATCH --job-name=lepc_load
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

module load R
set -euo pipefail

Rscript << 'REOF'
library(vcfR)
library(tidyverse)

out <- "OUT"

old_samples <- readLines("OLD_SAMPLES")
new_samples <- readLines("NEW_SAMPLES")
group_map <- c(setNames(rep("Old", length(old_samples)), old_samples),
               setNames(rep("New", length(new_samples)), new_samples))

sites <- read_tsv(sprintf("%s/load/sites_classified.tsv", out), show_col_types = FALSE)

vcf <- read.vcfR("VCF", verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

vcf_chrom <- getCHROM(vcf)
vcf_pos   <- getPOS(vcf)
vcf_site  <- paste(vcf_chrom, vcf_pos, sep = ":")

# For each site, figure out which allele (0=reference, 1=ALT) is the derived one
site_lookup <- sites %>%
    mutate(site = paste(chrom, pos, sep = ":"),
           derived_is_alt = derived == alt) %>%
    select(site, category, category_fine, derived_is_alt)

# Report BOTH the aggregate MISSENSE class and its Grantham split. A site
# contributes to its coarse category and, if scored, to its fine category,
# so MISSENSE totals stay complete while RADICAL/CONSERVATIVE are also
# available. Rows are duplicated only where the two labels differ.
site_lookup <- bind_rows(
    site_lookup %>% select(site, category, derived_is_alt),
    site_lookup %>%
        filter(category_fine != category) %>%
        mutate(category = category_fine) %>%
        select(site, category, derived_is_alt)
)

results <- list()

for (cat in c("LOF", "MISSENSE", "MISSENSE_RADICAL", "MISSENSE_CONSERVATIVE", "NEUTRAL")) {
    cat_sites <- site_lookup %>% filter(category == cat)
    if (nrow(cat_sites) == 0) next

    idx <- match(cat_sites$site, vcf_site)
    idx_valid <- !is.na(idx)
    cat_sites <- cat_sites[idx_valid, ]
    idx <- idx[idx_valid]

    gt_sub <- gt[idx, , drop = FALSE]

    for (ind in colnames(gt_sub)) {
        calls <- gt_sub[, ind]
        called <- !is.na(calls) & calls != "./." & calls != ".|."

        # Normalize genotype separators and pull the two alleles per site
        alleles <- strsplit(gsub("\\|", "/", calls[called]), "/")
        derived_flag <- cat_sites$derived_is_alt[called]

        n_derived <- mapply(function(al, der_is_alt) {
            target <- if (der_is_alt) "1" else "0"
            sum(al == target)
        }, alleles, derived_flag)

        n_sites   <- length(n_derived)
        n_hom_der <- sum(n_derived == 2)
        n_het_der <- sum(n_derived == 1)

        total_load    <- sum(n_derived) / (2 * n_sites)
        realized_load <- n_hom_der / n_sites
        masked_load   <- n_het_der / n_sites

        results[[length(results) + 1]] <- data.frame(
            individual     = ind,
            group          = group_map[[ind]],
            category       = cat,
            n_sites        = n_sites,
            n_het_derived  = n_het_der,
            n_hom_derived  = n_hom_der,
            total_load     = total_load,
            realized_load  = realized_load,
            masked_load    = masked_load
        )
    }
}

load_df <- bind_rows(results)
write_tsv(load_df, sprintf("%s/load/genetic_load_per_individual.tsv", out))

cat("\n=== Per-individual load (head) ===\n")
print(head(load_df, 20))

# Group-level summary and Old vs New comparison per category
summary_df <- load_df %>%
    group_by(group, category) %>%
    summarise(
        n_individuals   = n(),
        mean_total      = mean(total_load, na.rm = TRUE),
        sd_total        = sd(total_load, na.rm = TRUE),
        mean_realized   = mean(realized_load, na.rm = TRUE),
        sd_realized     = sd(realized_load, na.rm = TRUE),
        mean_masked     = mean(masked_load, na.rm = TRUE),
        sd_masked       = sd(masked_load, na.rm = TRUE),
        .groups = "drop"
    )
write_tsv(summary_df, sprintf("%s/load/genetic_load_summary.tsv", out))

cat("\n=== Group summary ===\n")
print(summary_df)

# Old vs New Wilcoxon test per category/metric (small n, non-parametric;
# swap for a paired test if Old/New individuals are matched pairs)
test_results <- list()
for (cat in unique(load_df$category)) {
    d <- load_df %>% filter(category == cat)
    for (metric in c("total_load", "realized_load", "masked_load")) {
        old_vals <- d[[metric]][d$group == "Old"]
        new_vals <- d[[metric]][d$group == "New"]
        if (length(old_vals) > 1 && length(new_vals) > 1) {
            wt <- wilcox.test(new_vals, old_vals)
            test_results[[length(test_results) + 1]] <- data.frame(
                category = cat, metric = metric,
                mean_old = mean(old_vals), mean_new = mean(new_vals),
                W = wt$statistic, p_value = wt$p.value
            )
        }
    }
}
test_df <- bind_rows(test_results)
write_tsv(test_df, sprintf("%s/load/old_vs_new_tests.tsv", out))

cat("\n=== Old vs New tests ===\n")
print(test_df)

cat("\nStep 6 complete.\n")
REOF

echo "Step 6 complete: $(date)"
EOF


# =============================================================================
# Substitute placeholder variables into all scripts
# =============================================================================
for SCRIPT in step1 step2 step4 step5 step6; do
    FILE=$PROJ/scripts/${SCRIPT}*.sh
    sed -i \
        -e "s|REF.fai|${REF}.fai|g" \
        -e "s|REF|$REF|g" \
        -e "s|GFF|$GFF|g" \
        -e "s|VCF|$VCF|g" \
        -e "s|GALLUS_DIR|$GALLUS_DIR|g" \
        -e "s|ANC|$ANC|g" \
        -e "s|SNPEFF_DB|$SNPEFF_DB|g" \
        -e "s|OLD_SAMPLES|$OLD_SAMPLES|g" \
        -e "s|NEW_SAMPLES|$NEW_SAMPLES|g" \
        -e "s|ALL_SAMPLES|$ALL_SAMPLES|g" \
        -e "s|OUT|$OUT|g" \
        -e "s|PROJ|$PROJ|g" \
        -e "s|THREADS|$THREADS|g" \
        $FILE
    chmod +x $FILE
done

echo ""
echo "============================================================"
echo " All scripts generated in $PROJ/scripts/"
echo " Suggested submission order:"
echo ""
echo "   JOB1=\$(sbatch --parsable scripts/step1_build_ancestral.sh)"
echo "   JOB2=\$(sbatch --parsable scripts/step2_snpeff_annotate.sh)"
echo "   JOB4=\$(sbatch --parsable --dependency=afterok:\$JOB1,afterok:\$JOB2 scripts/step4_polarize.sh)"
echo "   JOB5=\$(sbatch --parsable --dependency=afterok:\$JOB4 scripts/step5_deleterious_sites.sh)"
echo "   JOB6=\$(sbatch --parsable --dependency=afterok:\$JOB5 scripts/step6_load_per_individual.sh)"
echo ""
echo " Output files:"
echo "   $OUT/load/sites_classified.tsv               — per-site category (LOF/MISSENSE/NEUTRAL)"
echo "   $OUT/load/genetic_load_per_individual.tsv     — per-individual total/realized/masked load"
echo "   $OUT/load/genetic_load_summary.tsv            — Old vs New group summary"
echo "   $OUT/load/old_vs_new_tests.tsv                — Wilcoxon tests, Old vs New per category/metric"
echo "============================================================"
echo ""
echo "============================================================"
