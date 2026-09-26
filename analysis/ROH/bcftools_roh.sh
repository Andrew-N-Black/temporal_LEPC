#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: bcftools roh (ANGSD variant calling -> bcftools roh)
#
# Produces ROH_GROUSE_PL_regions.txt (the raw "RG" region lines, all
# samples mixed together) -- pair this with roh_parse_autosomal.sh for the
# actual per-sample F(ROH) summary (autosomal-only, same length bins as
# ROHan's parser, no rohparser.py dependency). Replicates ROH.sh from
# https://github.com/Andrew-N-Black/LEPC-popgen.
#
# For the independent ROHan (Renaud et al. 2019) cross-check -- estimated
# directly from CRAM via its own genotype-likelihood model rather than
# from these called genotypes -- see ROHan/rohan_array.sh +
# ROHan/rohan_parse_autosomal.sh instead; that pipeline no longer lives in
# this script (previously Step 7 here, now standalone and array-parallel).
#
# USAGE:
#   sbatch bcftools_roh.sh
# =============================================================================
#SBATCH --job-name=old_new_bcftools_roh
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -t 3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250G
#SBATCH -p cpu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=blackan@purdue.edu

set -euo pipefail

PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"
ROH_DIR="${PROJECT_DIR}/roh"
THREADS=$SLURM_CPUS_PER_TASK

mkdir -p logs "$ROH_DIR"

echo ">>> bcftools_roh.sh"
echo ">>> Start time: $(date)"

if [[ ! -f "$FINAL_CRAMLIST" ]]; then
    echo "ERROR: ${FINAL_CRAMLIST} not found." >&2
    exit 1
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "ERROR: ${REF_FASTA}.fai not found." >&2
    exit 1
fi

N_SAMPLES=$(wc -l < "$FINAL_CRAMLIST")
echo ">>> N samples : ${N_SAMPLES}"

# =============================================================================
# ANGSD/bcftools ENVIRONMENT
# =============================================================================
ml biocontainers
ml bcftools
ml angsd/0.940
ml htslib
# RCAC's xalt accounting hook injects LD_PRELOAD (libxalt_init.so) into
# every command, including containerized ones. singularity forwards it
# into the container by default, and the container's older glibc lacks
# the GLIBC_2.33/2.34 symbols that library needs, so angsd aborts before
# running. Blanking it inside the container via these two env vars is the
# reliable fix (a plain `unset LD_PRELOAD` on the host doesn't hold --
# xalt is sticky and re-injects it):
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""

# =============================================================================
# STEP 1: ANGSD genome-wide variant calling -> BCF (flags match ROH.sh
# exactly, extracted directly from its source — no -doGeno needed)
# =============================================================================
echo ">>> Step 1: ANGSD variant calling (BCF output)"

JOINT_OUT="${ROH_DIR}/joint"
JOINT_BCF="${JOINT_OUT}.bcf"

if [[ ! -f "$JOINT_BCF" ]]; then
    angsd -bam "$FINAL_CRAMLIST" -ref "$REF_FASTA" \
        -GL 1 -dobcf 1 -dopost 1 -domajorminor 1 -domaf 1 \
        -minQ 30 -SNP_pval 1e-6 -P "$THREADS" -out "$JOINT_OUT"
else
    echo "  ${JOINT_BCF} already exists -- skipping ANGSD call. Delete it first for a clean rerun."
fi

if [[ ! -f "$JOINT_BCF" ]]; then
    echo "ERROR: ANGSD did not produce expected output: ${JOINT_BCF}"
    exit 1
fi

# =============================================================================
# STEP 2: Allele frequency file for bcftools roh
# =============================================================================
echo ">>> Step 2: Building allele-frequency file"

FREQS="${ROH_DIR}/freqs.tab.gz"
if [[ ! -f "$FREQS" ]]; then
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%AF\n' "$JOINT_BCF" | bgzip -c > "$FREQS"
    tabix -s1 -b2 -e2 "$FREQS"
else
    echo "  ${FREQS} already exists -- skipping."
fi

# =============================================================================
# STEP 3: bcftools roh (flags match ROH.sh exactly), streamed through grep
# so the multi-GB per-site ST output never touches disk -- only the RG
# (called-region) lines, which is all downstream parsing actually needs.
# =============================================================================
echo ">>> Step 3: bcftools roh"

ROH_RG_ONLY="${ROH_DIR}/ROH_GROUSE_PL_regions.txt"
bcftools roh --AF-file "$FREQS" --threads "$THREADS" "$JOINT_BCF" \
    | grep "^RG" > "$ROH_RG_ONLY"

echo "  RG (called-region) lines: ${ROH_RG_ONLY}"
echo "  $(wc -l < "$ROH_RG_ONLY") regions called across all samples"

echo ""
echo ">>> bcftools roh complete."
echo "    Joint BCF          : ${JOINT_BCF}"
echo "    RG region lines    : ${ROH_RG_ONLY}"
echo "    Next               : run roh_parse_autosomal.sh for per-sample F(ROH)"
echo ">>> End time: $(date)"
