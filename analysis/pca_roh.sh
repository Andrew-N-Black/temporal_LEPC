#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: PCA + RUNS OF HOMOZYGOSITY (combined cohort)
# Step 08 — requires 06_downsample_and_finalize.sh (final_cramlist.txt) and,
# for the heterozygosity aggregation step, 07_heterozygosity_array.sh to
# have completed for all samples. Input is CRAM (ANGSD's -bam flag accepts
# a list of CRAM paths the same way it does BAM — it reads them through
# the same htslib backend regardless of format).
#
# Replicates beagle.sh + pca.sh + ROH.sh + rohparser.py from
# https://github.com/Andrew-N-Black/LEPC-popgen, with two deliberate,
# explicitly-requested deviations from the original:
#   - Whole-genome analysis, not the original's 100kb-window reference
#     subset (no windowing scheme to replicate/fabricate).
#   - ANGSD -doGlf2 (beagle) is parallelized per chromosome instead of per
#     100kb window, since we're not subsetting the genome — this still
#     covers 100% of the genome, just chunked for tractability.
#
# Everything else (ANGSD/pcangsd/bcftools flags) matches the original
# scripts exactly, verified against their actual source rather than
# guessed — see the flag comments at each step below.
#
# -minInd is set to round(0.75 x N), matching the ~75% stringency the
# original used (-minInd 348 of their N=~464) rather than a hardcoded
# number, since N here depends on how many NEW samples end up sequenced.
#
# USAGE:
#   sbatch 08_pca_roh.sh
# =============================================================================
#SBATCH --job-name=old.new_pca_roh
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

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
set -euo pipefail

ml biocontainers
ml bcftools
ml angsd/0.940
ml pcangsd
ml htslib
# RCAC's xalt accounting hook injects LD_PRELOAD (libxalt_init.so) into
# every command, including containerized ones. singularity forwards it
# into the container by default, and the container's older glibc lacks
# the GLIBC_2.33/2.34 symbols that library needs, so angsd/pcangsd abort
# before running. Blanking it inside the container via these two env vars
# is the reliable fix (a plain `unset LD_PRELOAD` on the host doesn't
# hold — xalt is sticky and re-injects it):
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"
HET_DIR="${PROJECT_DIR}/heterozygosity"

BEAGLE_DIR="${PROJECT_DIR}/beagle"
PCA_DIR="${PROJECT_DIR}/pca"
ROH_DIR="${PROJECT_DIR}/roh"

# rohparser.py — vendored verbatim from the original repo rather than
# reimplemented, so ROH size-class/FROH logic matches exactly. Its one
# hardcoded path (a .fai file, for total genome length) is patched below
# to point at our reference instead of the original's (different cluster).
ROHPARSER_URL="https://raw.githubusercontent.com/Andrew-N-Black/LEPC-popgen/main/analysis/rohparser.py"
ROHPARSER="${ROH_DIR}/rohparser.py"
ROHPARSER_ORIG_FAI="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna.fai"

THREADS=$SLURM_CPUS_PER_TASK
ROH_PARALLEL_JOBS=8

mkdir -p logs "$BEAGLE_DIR" "$PCA_DIR" "$ROH_DIR"

echo ">>> 08_pca_roh.sh"
echo ">>> Start time: $(date)"

if [[ ! -f "$FINAL_CRAMLIST" ]]; then
    echo "ERROR: ${FINAL_CRAMLIST} not found. Run 06_downsample_and_finalize.sh first."
    exit 1
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "ERROR: ${REF_FASTA}.fai not found. Run 05_combined_alignment_array.sh's prep step first."
    exit 1
fi

N_SAMPLES=$(wc -l < "$FINAL_CRAMLIST")
MININD=$(awk -v n="$N_SAMPLES" 'BEGIN { printf "%d", (n * 0.75) + 0.5 }')
echo ">>> N samples : ${N_SAMPLES}"
echo ">>> minInd    : ${MININD} (75% of N, matching the original's ~75% stringency)"

# =============================================================================
# STEP 0: Aggregate per-sample heterozygosity results from 08 (each array
# task there wrote its own file to avoid a shared-file race).
# =============================================================================
echo ">>> Step 0: Aggregating heterozygosity results"

HET_SUMMARY="${PROJECT_DIR}/heterozygosity_summary.tsv"
echo -e "sample_id\theterozygosity" > "$HET_SUMMARY"
find "$HET_DIR" -name "*_heterozygosity.txt" -exec cat {} + >> "$HET_SUMMARY" 2>/dev/null || true
echo "  Wrote ${HET_SUMMARY} ($(($(wc -l < "$HET_SUMMARY") - 1)) samples)"

# =============================================================================
# STEP 1: Genotype likelihoods (beagle format), parallelized per chromosome
# Flags match the original beagle.sh exactly (GL model, major/minor, MAF,
# quality, triallelic/SNP filtering) — only -minInd is recomputed for N,
# and region scope is per-chromosome instead of per-100kb-window.
# =============================================================================
echo ">>> Step 1: ANGSD genotype likelihoods (beagle format)"

CHROM_LIST="${BEAGLE_DIR}/chroms.txt"
cut -f1 "${REF_FASTA}.fai" > "$CHROM_LIST"
N_CHROMS=$(wc -l < "$CHROM_LIST")
echo "  ${N_CHROMS} chromosomes/contigs to process"

BEAGLE_THREADS_PER_JOB=$(( THREADS / ROH_PARALLEL_JOBS > 0 ? THREADS / ROH_PARALLEL_JOBS : 1 ))

run_beagle_chrom() {
    local chrom="$1"
    local out="${BEAGLE_DIR}/${chrom}"
    if [[ -f "${out}.beagle.gz" ]]; then
        return 0
    fi
    angsd -bam "$FINAL_CRAMLIST" -ref "$REF_FASTA" -r "${chrom}:" \
        -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -minQ 30 \
        -skipTriallelic 1 -SNP_pval 1e-6 -minInd "$MININD" \
        -P "$BEAGLE_THREADS_PER_JOB" -out "$out"
}
# angsd is a bash function (from `ml angsd/0.940`'s Lmod setup, wrapping the
# singularity call), not a real binary on PATH. Functions don't propagate
# into the fresh bash process `xargs ... bash -c` spawns unless each one is
# individually exported with `export -f` — exporting run_beagle_chrom alone
# is not enough, since its body calls angsd, which was never exported.
export -f angsd
export -f run_beagle_chrom
export FINAL_CRAMLIST REF_FASTA BEAGLE_DIR MININD BEAGLE_THREADS_PER_JOB

xargs -a "$CHROM_LIST" -I{} -P "$ROH_PARALLEL_JOBS" bash -c 'run_beagle_chrom "$@"' _ {}

echo ">>> Step 1b: Concatenating per-chromosome beagle files"

FINAL_BEAGLE="${BEAGLE_DIR}/final.beagle.gz"
if [[ ! -f "$FINAL_BEAGLE" ]]; then
    FIRST=1
    > "${BEAGLE_DIR}/final.beagle"
    while IFS= read -r chrom; do
        f="${BEAGLE_DIR}/${chrom}.beagle.gz"
        [[ ! -f "$f" ]] && { echo "  WARNING: missing ${f} — skipping" >&2; continue; }
        if [[ "$FIRST" -eq 1 ]]; then
            zcat "$f" >> "${BEAGLE_DIR}/final.beagle"
            FIRST=0
        else
            zcat "$f" | tail -n +2 >> "${BEAGLE_DIR}/final.beagle"
        fi
    done < "$CHROM_LIST"
    gzip "${BEAGLE_DIR}/final.beagle"
fi

echo "  Final beagle file: ${FINAL_BEAGLE}"

# =============================================================================
# STEP 2: PCA + inbreeding (pcangsd) — flags match pca.sh exactly
# =============================================================================
echo ">>> Step 2: pcangsd"

pcangsd -b "$FINAL_BEAGLE" -o "${PCA_DIR}/final" --threads "$THREADS" --minMaf 0.01 --admix

pcangsd -b "$FINAL_BEAGLE" -o "${PCA_DIR}/final_inbreed" --threads "$THREADS" --minMaf 0.01 \
    --maf_tole 1e-9 --tole 1e-9 --inbreedSamples --inbreedSites \
    --iter 5000 --maf_iter 5000 --inbreed_iter 5000 --inbreed_tole 1e-9

echo "  PCA output      : ${PCA_DIR}/final.cov (+ .admix.Q etc.)"
echo "  Inbreeding output: ${PCA_DIR}/final_inbreed.*"

# =============================================================================
# STEP 3: ANGSD genome-wide variant calling -> BCF (flags match ROH.sh
# exactly, extracted directly from its source — no -doGeno needed)
# =============================================================================
echo ">>> Step 3: ANGSD variant calling (BCF output)"

JOINT_OUT="${ROH_DIR}/joint"
JOINT_BCF="${JOINT_OUT}.bcf"

if [[ ! -f "$JOINT_BCF" ]]; then
    angsd -bam "$FINAL_CRAMLIST" -ref "$REF_FASTA" \
        -GL 1 -dobcf 1 -dopost 1 -domajorminor 1 -domaf 1 \
        -minQ 30 -SNP_pval 1e-6 -P "$THREADS" -out "$JOINT_OUT"
fi

if [[ ! -f "$JOINT_BCF" ]]; then
    echo "ERROR: ANGSD did not produce expected output: ${JOINT_BCF}"
    exit 1
fi

# =============================================================================
# STEP 4: Allele frequency file for bcftools roh
# =============================================================================
echo ">>> Step 4: Building allele-frequency file"

FREQS="${ROH_DIR}/freqs.tab.gz"
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%AF\n' "$JOINT_BCF" | bgzip -c > "$FREQS"
tabix -s1 -b2 -e2 "$FREQS"

# =============================================================================
# STEP 5: bcftools roh (flags match ROH.sh exactly)
# =============================================================================
echo ">>> Step 5: bcftools roh"

ROH_RAW="${ROH_DIR}/ROH_GROUSE_PLraw.txt"
bcftools roh --AF-file "$FREQS" --output "$ROH_RAW" --threads "$THREADS" "$JOINT_BCF"

echo "  Raw ROH output: ${ROH_RAW}"

# =============================================================================
# STEP 6: Per-sample ROH parsing with rohparser.py (vendored from the
# original repo, patched to use our reference's .fai for genome length)
# =============================================================================
echo ">>> Step 6: Per-sample ROH parsing"

if [[ ! -f "$ROHPARSER" ]]; then
    echo ">>> Downloading rohparser.py"
    wget -q -O "$ROHPARSER" "$ROHPARSER_URL"
    sed -i "s|${ROHPARSER_ORIG_FAI}|${REF_FASTA}.fai|g" "$ROHPARSER"
fi

# Split the joint RG lines out per sample (bcftools roh RG columns:
# RG, sample, chrom, start, end, length, n_markers, quality — matches
# rohparser.py's expected field[5]=length, field[7]=quality).
while IFS= read -r cram; do
    sample=$(basename "$cram")
    sample="${sample%_filt.cram}"
    sample="${sample%_ds.cram}"
    grep "^RG" "$ROH_RAW" | awk -v s="$sample" '$2==s' > "${ROH_DIR}/${sample}ROH.txt"
done < "$FINAL_CRAMLIST"

run_rohparser() {
    local roh_file="$1"
    python3 "$ROHPARSER" "$roh_file" > "${roh_file}_results.txt"
}
export -f run_rohparser
export ROHPARSER

find "$ROH_DIR" -name "*ROH.txt" | xargs -I{} -P "$ROH_PARALLEL_JOBS" bash -c 'run_rohparser "$@"' _ {}

N_ROH_RESULTS=$(find "$ROH_DIR" -name "*ROH.txt_results.txt" | wc -l)
echo "  Parsed ROH results for ${N_ROH_RESULTS} samples"

echo ""
echo ">>> PCA + ROH analysis complete."
echo "    Heterozygosity : ${HET_SUMMARY}"
echo "    PCA            : ${PCA_DIR}/final.cov"
echo "    Inbreeding     : ${PCA_DIR}/final_inbreed.*"
echo "    Raw ROH        : ${ROH_RAW}"
echo "    Per-sample ROH : ${ROH_DIR}/*ROH.txt_results.txt"
echo ">>> End time: $(date)"
