#!/bin/bash
#SBATCH --job-name=ne_temporal_input
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH -A dewoody
#SBATCH -p cpu
#SBATCH -t 0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# =============================================================================
# Build a NeEstimator (temporal method) GENEPOP input from the autosomal,
# unrelated-sample LEPC VCF.
#
#   1. keep the 19 unrelated samples (normal_F10 dropped) and biallelic SNPs
#   2. require >= 5 genotyped birds in EACH era (AN >= 10 per group)
#   3. drop sites with >= 20% missing overall and minor-allele freq < 0.05
#      (same thresholds as the Fk / Ne analysis in outflank_drift.sh)
#   4. thin to one random SNP per THIN window (NeEstimator cannot take
#      millions of loci; ~1 SNP / 10 kb gives ~100k loci on ~1.05 Gb)
#   5. write GENEPOP: Past samples as Pop 1 (generation 0), Present as Pop 2
#
# Copy the .gen file to your Mac and run NeEstimator there with
# info_temporal.txt (see that file).
# =============================================================================
set -euo pipefail

BASE=/scratch/gautschi/blackan/GROUSE/old_vs_new
VCF=$BASE/vcfs/output.subset.biallelic.snps.auto.vcf.gz   # autosomal (Z removed)
POPMAP=$BASE/popmap_unrelated.txt                           # 19 samples, F10 excluded
OUTDIR=$BASE/results/neestimator
THIN=10kb             # window for thinning: 1 SNP kept per window
MAF=0.05
MAX_MISSING=0.2
MIN_AN_PER_ERA=10     # = 5 diploid birds genotyped per era
SEED=42
THREADS=${SLURM_CPUS_PER_TASK:-4}

# bcftools on Gautschi comes from biocontainers (same setup as load.sh Step 6).
module --force purge
module load biocontainers
module load bcftools
# xalt injects LD_PRELOAD into containerised commands; blank it or bcftools
# aborts with "GLIBC_2.33/2.34 not found (required by libxalt_init.so)".
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""
command -v bcftools >/dev/null || { echo "bcftools not found" >&2; exit 1; }

mkdir -p "$OUTDIR" logs
cd "$OUTDIR"
cut -f1 "$POPMAP" > samples_unrelated.txt
PREFIX=lepc_auto_unrel

FILT=${PREFIX}.filt.vcf.gz
THINNED=${PREFIX}.thin${THIN}.vcf.gz

# Steps 1-4 are skipped if their outputs already exist (delete them to redo).
if [[ ! -s "$FILT" ]]; then
    echo "[1-3] subset samples, per-era call rate, missingness, MAF: $(date)"
    bcftools view --threads "$THREADS" -S samples_unrelated.txt -a -Ou "$VCF" \
      | bcftools view --threads "$THREADS" -m2 -M2 -v snps -Ou \
      | bcftools +fill-tags -Ou -- -S "$POPMAP" -t AN,AC,MAF \
      | bcftools view --threads "$THREADS" \
            -i "INFO/AN_Past>=${MIN_AN_PER_ERA} && INFO/AN_Present>=${MIN_AN_PER_ERA} && F_MISSING<${MAX_MISSING} && INFO/MAF>=${MAF}" \
            -Oz -o "$FILT"
    bcftools index -t -f "$FILT"
fi
echo "    sites after filtering: $(bcftools index -n "$FILT")"

if [[ ! -s "$THINNED" ]]; then
    echo "[4] thin: 1 random SNP per ${THIN}: $(date)"
    bcftools +prune -n 1 -w "$THIN" -N rand --random-seed "$SEED" \
        "$FILT" -Oz -o "$THINNED"
    bcftools index -t -f "$THINNED"
fi
echo "    sites after thinning: $(bcftools index -n "$THINNED")"

echo "[5] write GENEPOP: $(date)"
# bcftools is a module shell function here, invisible to Python's subprocess,
# so the genotype dump is done in this shell and handed to the converter.
bcftools query -l "$THINNED" > ${PREFIX}.thin${THIN}.samples.txt
bcftools query -f '%CHROM:%POS[\t%GT]\n' "$THINNED" > ${PREFIX}.thin${THIN}.geno.tsv
python3 "${SLURM_SUBMIT_DIR:-$PWD}/vcf2genepop.py" \
    --samples ${PREFIX}.thin${THIN}.samples.txt \
    --geno ${PREFIX}.thin${THIN}.geno.tsv \
    --popmap "$POPMAP" \
    --out ${PREFIX}.thin${THIN}.gen --order Past,Present
rm -f ${PREFIX}.thin${THIN}.geno.tsv

echo "Done. Copy to your Mac with:"
echo "  scp $(whoami)@gautschi.rcac.purdue.edu:$OUTDIR/${PREFIX}.thin${THIN}.gen ~/NeEstimator2.X/"
