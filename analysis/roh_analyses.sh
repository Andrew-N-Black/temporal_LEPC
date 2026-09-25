#!/bin/bash
# =============================================================================
# SLURM JOB SUBMISSION: ROH + ROHan (combined cohort)
# Step 08b — RESUMES AFTER PCA/ADMIXTURE. This is a trimmed copy of
# 08_pca_roh.sh: Steps 0-2 (heterozygosity aggregation, beagle genotype
# likelihoods, pcangsd PCA/admixture/inbreeding) are assumed already
# complete and have been removed. Requires 06_downsample_and_finalize.sh
# (final_cramlist.txt) to have completed. Input is CRAM (ANGSD's -bam flag
# accepts a list of CRAM paths the same way it does BAM — it reads them
# through the same htslib backend regardless of format).
#
# Replicates ROH.sh + rohparser.py from
# https://github.com/Andrew-N-Black/LEPC-popgen, plus ROHan (Renaud et al.
# 2019) as a cross-check estimated directly from BAM/CRAM via its own
# genotype-likelihood model rather than from called genotypes. ROHan runs
# last since its module environment (a clean gcc/samtools setup) is
# incompatible with the angsd/bcftools modules the earlier steps need, and
# there's no reason to juggle both at once when nothing later in the
# script needs the earlier modules again. install_rohan.sh must have been
# run on the login node first.
#
# USAGE:
#   sbatch 08b_roh_resume.sh
# =============================================================================
#SBATCH --job-name=old.new_roh_resume
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

# =============================================================================
# USER-DEFINED VARIABLES
# =============================================================================
# PROJECT_DIR: this script is for the old_vs_new temporal cohort (n=20) --
# same project as every other script in this repo. (A previously committed
# version of this file had PROJECT_DIR pointing at the separate GROUSE/nexus
# 433-sample project by mistake, left over from adapting this script for
# that cohort; fixed here. ROHAN_BIN/GSL_PREFIX below correctly point at
# old_vs_new regardless, since those tools are built once and shared.)
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
FINAL_CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"

ROH_DIR="${PROJECT_DIR}/roh"

# rohparser.py — a corrected, locally vendored copy (analysis/rohparser.py,
# same directory as this script), not downloaded at runtime. The original
# Andrew-N-Black/LEPC-popgen version had a hardcoded absolute
# path_to_directory pointing at GROUSE/nexus (a different project) that a
# runtime sed patch never touched (only its ref_index_file path was
# patched), and divided every sample's own ROH counts/lengths by a
# hardcoded cohort-size constant left over from that project (506 samples)
# before computing F(ROH) -- silently deflating every F(ROH) value by that
# factor. The vendored copy here takes --fai and --exclude-scaffolds as
# CLI arguments instead of hardcoding either, and drops the cohort-size
# division entirely (an individual's F(ROH) has no cohort-size term). See
# that file's own header for the full list of changes.
ROHPARSER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/rohparser.py"
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"

THREADS=$SLURM_CPUS_PER_TASK
ROH_PARALLEL_JOBS=8

# ROHan-specific config (Step 7, run last -- see module note there)
ROHAN_BIN="/scratch/gautschi/blackan/GROUSE/old_vs_new/tools/ROHan/bin/rohan"
GSL_PREFIX="/scratch/gautschi/blackan/GROUSE/old_vs_new/tools/gsl"
ROHAN_OUT_DIR="${PROJECT_DIR}/results/rohan"
ROHAN_THREADS=16
# ROHan's expected within-ROH heterozygosity rate parameter. ROHan's own
# examples use something on the order of 2e-5; adjust if your species'
# expected mutation rate differs substantially. Confirm the flag name is
# still --rohmu with `rohan --help` if ROHan is ever rebuilt/updated.
ROHMU=2e-5

mkdir -p logs "$ROH_DIR" "$ROHAN_OUT_DIR"

echo ">>> 08b_roh_resume.sh (resuming after PCA/admixture -- runs Steps 3-7)"
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
echo ">>> N samples : ${N_SAMPLES}"

# =============================================================================
# ANGSD/bcftools ENVIRONMENT (Steps 3-5 below)
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
else
    echo "  ${JOINT_BCF} already exists -- skipping ANGSD call. Delete it first for a clean rerun."
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
if [[ ! -f "$FREQS" ]]; then
    bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%AF\n' "$JOINT_BCF" | bgzip -c > "$FREQS"
    tabix -s1 -b2 -e2 "$FREQS"
else
    echo "  ${FREQS} already exists -- skipping."
fi

# =============================================================================
# STEP 5: bcftools roh (flags match ROH.sh exactly), streamed through grep
# so the multi-GB per-site ST output never touches disk -- only the RG
# (called-region) lines, which is all downstream parsing actually needs.
# =============================================================================
echo ">>> Step 5: bcftools roh"

ROH_RG_ONLY="${ROH_DIR}/ROH_GROUSE_PL_regions.txt"
bcftools roh --AF-file "$FREQS" --threads "$THREADS" "$JOINT_BCF" \
    | grep "^RG" > "$ROH_RG_ONLY"

echo "  RG (called-region) lines: ${ROH_RG_ONLY}"
echo "  $(wc -l < "$ROH_RG_ONLY") regions called across all samples"

# =============================================================================
# STEP 6: Per-sample ROH parsing with rohparser.py (locally vendored, see
# CONFIG above), restricted to autosomes (excludes ${Z_SCAFFOLDS})
# =============================================================================
echo ">>> Step 6: Per-sample ROH parsing (autosomal F(ROH), excluding ${Z_SCAFFOLDS})"

if [[ ! -f "$ROHPARSER" ]]; then
    echo "ERROR: ${ROHPARSER} not found -- expected analysis/rohparser.py next to this script." >&2
    exit 1
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "ERROR: ${REF_FASTA}.fai not found." >&2
    exit 1
fi

# Single pass over the RG-only file, splitting by sample. NOTE: bcftools
# roh's sample column here is whatever the BCF's own header used as the
# sample name -- and ANGSD's -dobcf output uses the full BAM/CRAM file path
# as that identifier, not a bare sample ID (confirmed from the actual RG
# lines: column 2 is a full "/scratch/.../crams/F10.cram" path). Using
# that path verbatim as a filename produces a broken, doubled path, so
# derive a clean sample ID from its basename instead, stripping this
# project's actual CRAM suffix (plain .cram -- confirmed against
# final_cramlist.txt; NOT .md.dedup_q20.cram as an earlier version of
# this comment assumed).
# Clear out per-sample split files (and their parsed results) from any
# previous run before regenerating -- including ones left over under the
# old, buggy .md.dedup_q20.cram naming, which would otherwise coexist
# with the correctly-named files below and double-count samples.
rm -f "$ROH_DIR"/*ROH.txt "$ROH_DIR"/*ROH.txt_results.txt

# Touch a placeholder <sample>ROH.txt for every sample in the cramlist
# FIRST, before the awk split below. awk's `print >` only ever opens a
# file for samples that actually appear in the RG-only data, so a sample
# with zero called ROH regions would otherwise get no file at all. An
# empty input is handled gracefully by rohparser.py (confirmed: it just
# reports all-zero counts, no error) -- but only if the file exists.
# awk (next) truncates+overwrites these placeholders for samples that do
# have RG data; samples with none are left as empty (zero-ROH) files.
while IFS= read -r CRAM; do
    s=$(basename "$CRAM")
    s="${s%.cram}"
    : > "${ROH_DIR}/${s}ROH.txt"
done < "$FINAL_CRAMLIST"

awk -v dir="$ROH_DIR" '
{
    n = split($2, parts, "/")
    sample = parts[n]
    gsub(/\.cram$/, "", sample)
    print > (dir"/"sample"ROH.txt")
}' "$ROH_RG_ONLY"

N_SAMPLE_FILES=$(find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt" | wc -l)
N_ZERO_ROH=$(find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt" -empty | wc -l)
echo "  ${N_SAMPLE_FILES} per-sample files present (expected ${N_SAMPLES}), of which ${N_ZERO_ROH} have zero called ROH regions"
if [[ "$N_SAMPLE_FILES" -ne "$N_SAMPLES" ]]; then
    echo "  WARNING: per-sample file count doesn't match N_SAMPLES -- unexpected" >&2
    echo "  now that every sample gets a placeholder file regardless of ROH" >&2
    echo "  calls, so this points at a real problem, not just zero calls." >&2
    echo "  Missing sample(s):" >&2
    MISSING_SAMPLES=$(comm -23 \
        <(awk -F'/' '{s=$NF; sub(/\.cram$/,"",s); print s}' "$FINAL_CRAMLIST" | sort -u) \
        <(find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt" -printf "%f\n" | sed 's/ROH\.txt$//' | sort -u))
    echo "$MISSING_SAMPLES" | sed 's/^/    /' >&2
    echo "  Of those, present in the joint BCF at all (in BCF = zero ROH" >&2
    echo "  regions called, not in BCF = dropped before bcftools roh):" >&2
    # Diagnostic only -- reload bcftools defensively in case module state
    # has changed by this point in the job; isolated in a subshell with
    # set +e so nothing here can ever abort the real pipeline under
    # set -euo pipefail.
    (
        set +e
        ml biocontainers >/dev/null 2>&1
        ml bcftools >/dev/null 2>&1
        export SINGULARITYENV_LD_PRELOAD=""
        export APPTAINERENV_LD_PRELOAD=""
        if ! command -v bcftools >/dev/null 2>&1; then
            echo "    (bcftools unavailable here -- run 'bcftools query -l ${JOINT_BCF}' manually)" >&2
        else
            BCF_SAMPLES=$(bcftools query -l "$JOINT_BCF" 2>/dev/null)
            while IFS= read -r s; do
                [[ -z "$s" ]] && continue
                if grep -qF -- "${s}.cram" <<< "$BCF_SAMPLES"; then
                    echo "    ${s}: in joint BCF -- zero ROH regions called" >&2
                else
                    echo "    ${s}: NOT in joint BCF -- dropped before bcftools roh" >&2
                fi
            done <<< "$MISSING_SAMPLES"
        fi
    ) || true
fi

run_rohparser() {
    # rohparser.py builds its own input path internally from a hardcoded
    # directory + a bare filename (matching its documented usage: `cd` into
    # the ROH directory, then `python ROHparser.py SAMPLEROH.txt`). Passing
    # it a full path instead -- as a naive `find`-based invocation would --
    # makes it concatenate a doubled, nonexistent path and fail. So: cd into
    # the file's directory and pass only the basename.
    local roh_file="$1"
    local bn
    bn=$(basename "$roh_file")
    (cd "$(dirname "$roh_file")" && python3 "$ROHPARSER" "$bn" \
        --fai "${REF_FASTA}.fai" \
        --exclude-scaffolds "$Z_SCAFFOLDS") > "${roh_file}_results.txt"
}
export -f run_rohparser
export ROHPARSER REF_FASTA Z_SCAFFOLDS

find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt" \
    | xargs -I{} -P "$ROH_PARALLEL_JOBS" bash -c 'run_rohparser "$@"' _ {}

N_ROH_RESULTS=$(find "$ROH_DIR" -maxdepth 1 -name "*ROH.txt_results.txt" | wc -l)
echo "  Parsed ANGSD/bcftools ROH results for ${N_ROH_RESULTS} samples"

# =============================================================================
# STEP 7: ROHan (Renaud et al. 2019) -- per-sample heterozygosity/ROH
# estimated directly from BAM/CRAM via its own genotype-likelihood model,
# as a cross-check against the ANGSD/bcftools estimates above. Does NOT
# take a VCF/BCF as input by design (that's the point of the tool -- it
# stays upstream of hard-called genotypes). Requires install_rohan.sh to
# have been run on the login node first.
#
# Runs LAST and does its own `module --force purge`: its build needs a
# plain gcc/samtools environment incompatible with the angsd/bcftools
# modules loaded above, and nothing after this point needs those modules
# again.
# =============================================================================
echo ">>> Step 7: ROHan per-sample analysis"

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

while IFS= read -r CRAM; do
    SAMPLE=$(basename "$CRAM" | sed -E 's/\.cram$//')
    OUT_PREFIX="${ROHAN_OUT_DIR}/${SAMPLE}"
    if [[ -f "${OUT_PREFIX}.hEst" ]]; then
        echo "  ${SAMPLE}: ${OUT_PREFIX}.hEst already exists -- skipping."
        continue
    fi
    echo "=== ROHan: $SAMPLE ==="

    # ROHan's most reliably supported input format is BAM. Rather than rely
    # on ROHan's own CRAM/reference handling (undocumented in what we could
    # verify), convert to a temporary indexed BAM first -- slower, but
    # removes any ambiguity about reference resolution.
    TMP_BAM="${ROHAN_OUT_DIR}/${SAMPLE}.tmp.bam"
    samtools view -@ "$ROHAN_THREADS" -b -T "$REF_FASTA" -o "$TMP_BAM" "$CRAM"
    samtools index "$TMP_BAM"

    "$ROHAN_BIN" \
        -t "$ROHAN_THREADS" \
        --rohmu "$ROHMU" \
        -o "$OUT_PREFIX" \
        "$REF_FASTA" "$TMP_BAM"

    rm -f "$TMP_BAM" "${TMP_BAM}.bai"
    echo "  Done: ${OUT_PREFIX}.*"
done < "$FINAL_CRAMLIST"

echo ""
echo ">>> ROH + ROHan analysis complete."
echo "    Joint BCF           : ${JOINT_BCF}"
echo "    ANGSD/bcftools ROH  : ${ROH_DIR}/*ROH.txt_results.txt"
echo "    ROHan (per sample)  : ${ROHAN_OUT_DIR}/<sample>.hEst, .mid.hmmp, .mid.ROH"
echo ">>> End time: $(date)"
