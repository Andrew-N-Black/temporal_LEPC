#!/bin/bash
# ============================================================================
# SLURM wrapper: Past vs Present LEPC temporal analysis
#   (0a) Input checks (fail in seconds, not hours)
#   (0b) Build siteclass.txt from SnpEff (+ optional GERP) if not present
#   (0c) bcftools prefilter: popmap samples, biallelic SNPs, polymorphic,
#        F_MISSING < threshold (no MAF filter, so rare deleterious alleles stay)
#   (1)  R: OutFLANK + temporal Fk/Ne + drift null + site-class comparison
#        (checkpointed: reruns reuse the genotype matrix and OutFLANK results)
#
# One-time environment setup (run on a login node before submitting):
#   module load anaconda
#   conda create -n lepc_rstats -c conda-forge -c bioconda \
#       r-base r-optparse r-dplyr r-ggplot2 r-vcfr bioconductor-qvalue r-devtools
#   conda activate lepc_rstats
#   Rscript -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'
#   Rscript -e 'remotes::install_github("whitlock/OutFLANK")'
#
# bcftools comes from the RCAC biocontainers module loaded below.
# ============================================================================
#SBATCH -J lepc_outflank_temporal
#SBATCH -o lepc_outflank_temporal_%j.out
#SBATCH -e lepc_outflank_temporal_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH -A dewoody
#SBATCH -p cpu

set -euo pipefail

module --force purge
module load biocontainers bcftools
module load anaconda
# conda's activate scripts reference unset variables, so relax -u briefly
set +u
source activate lepc_rstats
set -u

cd "$SLURM_SUBMIT_DIR"

# --- Edit these paths/params for your run -----------------------------------
RAW_VCF="/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs/output.subset.biallelic.snps.auto.vcf.gz"
FILT_VCF="${RAW_VCF%.vcf.gz}.temporal_filt.vcf.gz"
POPMAP="/scratch/gautschi/blackan/GROUSE/old_vs_new/popmap_unrelated.txt"
# Site classes: use SITECLASS if it exists; otherwise build it from SNPEFF_VCF.
# If neither is available, the site-class step is skipped (with a warning).
SITECLASS="siteclass.txt"    # <CHROM> <POS> <deleterious|neutral>, no header
SNPEFF_VCF="/scratch/gautschi/blackan/GROUSE/old_vs_new/results/snpeff/lepc_snpeff.format.vcf.gz"                # SnpEff-annotated VCF (INFO/ANN), e.g. .../snpeff/old_vs_new.ann.vcf.gz
GERP_TSV=""                  # optional per-site <CHROM> <POS> <RS>; "" = SnpEff only
GERP_MIN=2                   # min GERP RS for a site to count as deleterious
NEUTRAL_MODE="synonymous"    # synonymous | synonymous+intergenic

GENERATIONS=5                # elapsed generations between Past and Present samples
MAX_MISSING=0.2              # bcftools: drop sites with >= this fraction missing
MIN_CALLED=5                 # R: min called individuals in EACH period
MAX_FIT_LOCI=100000          # R: random loci used to fit the OutFLANK null
RUN_OUTFLANK=TRUE            # R: TRUE/FALSE (result is reused from checkpoint on reruns)
NE_MINFREQ=0.05              # R: exclude rare-allele loci from the Ne estimate
JK_BLOCKSIZE=5000000         # R: jackknife block size (bp) for the Ne CI
OUTPREFIX="results/old_vs_new_temporal_unrelated"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
# -----------------------------------------------------------------------------

mkdir -p "$(dirname "$OUTPREFIX")"
echo "== Job $SLURM_JOB_ID on $(hostname) at $(date) =="

# --- Step 0a: input checks ---------------------------------------------------
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found on PATH"; exit 1; }
command -v Rscript  >/dev/null || { echo "ERROR: Rscript not found (conda env?)"; exit 1; }
for f in "$RAW_VCF" "$POPMAP" outflank_temporal_dNe_siteclass.R; do
    [[ -s "$f" ]] || { echo "ERROR: required file not found or empty: $f"; exit 1; }
done

SAMPLES="${OUTPREFIX}_popmap_samples.txt"
awk 'NF>=1 {print $1}' "$POPMAP" | tr -d '\r' > "$SAMPLES"
MISSING_IN_VCF=$(comm -23 <(sort "$SAMPLES") <(bcftools query -l "$RAW_VCF" | sort) || true)
if [[ -n "$MISSING_IN_VCF" ]]; then
    echo "WARNING: popmap samples not found in VCF:"
    echo "$MISSING_IN_VCF"
fi

# --- Step 0b: site classes ---------------------------------------------------
SITECLASS_ARGS=()
if [[ -s "$SITECLASS" ]]; then
    echo "[siteclass] Using existing $SITECLASS"
    SITECLASS_ARGS=(--siteclass "$SITECLASS")
elif [[ -n "$SNPEFF_VCF" ]]; then
    GERP_ARGS=()
    [[ -n "$GERP_TSV" ]] && GERP_ARGS=(-g "$GERP_TSV" -m "$GERP_MIN")
    bash make_siteclass.sh -a "$SNPEFF_VCF" -o "$SITECLASS" -n "$NEUTRAL_MODE" "${GERP_ARGS[@]}"
    SITECLASS_ARGS=(--siteclass "$SITECLASS")
else
    echo "WARNING: $SITECLASS not found and SNPEFF_VCF not set; site-class step will be skipped."
fi

if [[ ${#SITECLASS_ARGS[@]} -gt 0 ]]; then
    # CHROM naming check: site-class contigs should appear in the VCF header
    SC_CHROM=$(awk 'NR==1{print $1; exit}' "$SITECLASS")
    if ! bcftools view -h "$RAW_VCF" | grep -q "^##contig=<ID=${SC_CHROM},"; then
        echo "WARNING: contig '$SC_CHROM' from $SITECLASS not in VCF header; CHROM naming may differ."
    fi
fi

# --- Step 0c: bcftools prefilter (skipped if an up-to-date output exists) ----
if [[ ! -s "$FILT_VCF" || "$RAW_VCF" -nt "$FILT_VCF" ]]; then
    echo "[prefilter] $(date): filtering $RAW_VCF"
    bcftools view --threads "$THREADS" -a -S "$SAMPLES" --force-samples -Ou "$RAW_VCF" \
      | bcftools view --threads "$THREADS" \
            -m2 -M2 -v snps -c 1:minor \
            -i "F_MISSING<${MAX_MISSING}" \
            -Oz -o "$FILT_VCF"
    bcftools index -t --threads "$THREADS" -f "$FILT_VCF"
else
    echo "[prefilter] Using existing filtered VCF: $FILT_VCF"
fi
echo "[prefilter] Filtered records: $(bcftools index -n "$FILT_VCF")"

# --- Step 1: R analysis ------------------------------------------------------
echo "[R] $(date): starting analysis"
Rscript outflank_temporal_dNe_siteclass.R \
    --vcf "$FILT_VCF" \
    --popmap "$POPMAP" \
    "${SITECLASS_ARGS[@]}" \
    --generations "$GENERATIONS" \
    --mincalled "$MIN_CALLED" \
    --maxfitloci "$MAX_FIT_LOCI" \
    --run_outflank "$RUN_OUTFLANK" \
    --ne_minfreq "$NE_MINFREQ" \
    --jk_blocksize "$JK_BLOCKSIZE" \
    --nsim 2000 \
    --nperm 2000 \
    --freqbins 20 \
    --qthresh 0.05 \
    --outprefix "$OUTPREFIX"

echo "== Finished at $(date) =="
