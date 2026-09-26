#!/bin/bash
# ============================================================================
# SLURM: identify Z (and W) scaffolds in the scaffold-level LEPC reference
# (pur_lepc_1.0, GCF_026119805.1), then confirm with the VCF and call sex.
#
#   Step 1  Synteny: align LEPC scaffolds to the chicken chromosome assembly
#           (GRCg7b, GCF_016699485.2) with minimap2 -x asm20. Each scaffold is
#           assigned the chicken chromosome with the most aligned bases.
#   Step 2  Genotype check: sex_from_vcf.py. Females (ZW) are hemizygous for Z,
#           so true Z scaffolds show ~no heterozygosity in females. Confirms the
#           synteny list, finds Z scaffolds synteny missed, and calls each
#           bird's sex.
#
# Output: $OUTDIR/SEX_CHROMS.txt, a comma-separated list ready to paste into
#         SEX_CHROMS in run_lepc_relatedness.sh.
#
# Requirements:
#   - minimap2 from the RCAC biocontainers module (loaded below)
#   - lepc_py conda env (scikit-allel) for the genotype check
#   - LEPC and chicken references already in REF_DIR. The chicken FASTA is
#     auto-detected there (GRCg7b / bGalGal1 / Gallus in the name); set
#     CHICKEN_FNA explicitly if detection picks the wrong file.
# ============================================================================
#SBATCH -J lepc_sex_scaffolds
#SBATCH -o lepc_sex_scaffolds_%j.out
#SBATCH -e lepc_sex_scaffolds_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH -A dewoody
#SBATCH -p cpu

set -euo pipefail

module --force purge
module load biocontainers minimap2
module load anaconda
cd "$SLURM_SUBMIT_DIR"

# --- Edit these paths/params for your run -----------------------------------
REF_DIR="/scratch/gautschi/blackan/GROUSE/old_vs_new/ref"
LEPC_REF="$REF_DIR/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
CHICKEN_FNA=""          # leave "" to auto-detect in REF_DIR, or give the full path
VCF="/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs/output.subset.biallelic.snps.vcf.gz"
POPMAP="/scratch/gautschi/blackan/GROUSE/old_vs_new/popmap.txt"
MIN_MAPQ=20             # PAF mapping quality filter
MIN_FRAC=0.6            # scaffold's best chicken chromosome must hold >= this fraction of aligned bp
MIN_ALIGNED=50000       # and at least this many aligned bp
OUTDIR="results/sex_scaffolds"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
# -----------------------------------------------------------------------------

# Auto-detect the chicken reference in REF_DIR
if [[ -z "$CHICKEN_FNA" ]]; then
    mapfile -t hits < <(find "$REF_DIR" -maxdepth 2 -type f \
        \( -iname '*GRCg7b*' -o -iname '*bGalGal1*' -o -iname '*GRCg6a*' -o -iname '*gallus*' \) \
        \( -name '*.fna' -o -name '*.fna.gz' -o -name '*.fa' -o -name '*.fa.gz' \
           -o -name '*.fasta' -o -name '*.fasta.gz' \) | sort)
    if [[ ${#hits[@]} -eq 1 ]]; then
        CHICKEN_FNA="${hits[0]}"
    else
        echo "ERROR: found ${#hits[@]} candidate chicken FASTAs in $REF_DIR; set CHICKEN_FNA explicitly."
        printf '  %s\n' "${hits[@]}"
        echo "FASTA files in $REF_DIR:"
        find "$REF_DIR" -maxdepth 2 -type f \( -name '*.fna*' -o -name '*.fa' -o -name '*.fa.gz' -o -name '*.fasta*' \) | sed 's/^/  /'
        exit 1
    fi
fi
echo "LEPC reference:    $LEPC_REF"
echo "Chicken reference: $CHICKEN_FNA"

mkdir -p "$OUTDIR"
echo "== Job $SLURM_JOB_ID on $(hostname) at $(date) =="
for f in "$LEPC_REF" "$CHICKEN_FNA" "$VCF" "$POPMAP" sex_from_vcf.py; do
    [[ -s "$f" ]] || { echo "ERROR: required file not found: $f"; exit 1; }
done

# --- Step 1: synteny to chicken -------------------------------------------------
type minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 not available (module load biocontainers minimap2)"; exit 1; }

CHR_FA="$OUTDIR/chicken_chromosomes.fa"
CHR_MAP="$OUTDIR/chicken_acc2chr.tsv"
if [[ ! -s "$CHR_FA" || "$CHICKEN_FNA" -nt "$CHR_FA" ]]; then
    echo "[1] Extracting chicken assembled chromosomes"
    # Accession -> chromosome name. Handles NCBI headers ("... chromosome Z, ...")
    # and Ensembl/UCSC-style names (">Z", ">chrZ"). Unplaced/unlocalized and
    # mitochondrial sequences are dropped.
    zcat -f "$CHICKEN_FNA" | awk '/^>/ {
            acc = substr($1, 2); name = ""
            if ($0 ~ /mitochondri/ || $0 ~ /unlocalized|unplaced/) next
            if (match($0, /chromosome [^ ,]+/)) name = substr($0, RSTART + 11, RLENGTH - 11)
            else { t = acc; sub(/^chr/, "", t); if (t ~ /^([0-9]+|Z|W)$/) name = t }
            if (name != "" && name != "MT" && name != "Un") print acc "\t" name }' > "$CHR_MAP"
    zcat -f "$CHICKEN_FNA" | awk 'NR==FNR {keep[$1]=1; next}
            /^>/ {k = (substr($1, 2) in keep)} k' "$CHR_MAP" - > "$CHR_FA"
fi
[[ -s "$CHR_MAP" ]] || { echo "ERROR: no chromosome-level sequences recognized in $CHICKEN_FNA"; exit 1; }
awk '$2=="Z"' "$CHR_MAP" | grep -q . || { echo "ERROR: no chicken Z found in $CHICKEN_FNA headers"; exit 1; }
echo "[1] Chicken chromosomes: $(wc -l < "$CHR_MAP") (Z = $(awk '$2=="Z"{print $1}' "$CHR_MAP"), W = $(awk '$2=="W"{print $1}' "$CHR_MAP"))"

PAF="$OUTDIR/lepc_vs_chicken.paf"
if [[ ! -s "$PAF" ]]; then
    echo "[1] $(date): minimap2 -x asm20 (LEPC scaffolds -> chicken chromosomes)"
    minimap2 -x asm20 -t "$THREADS" --secondary=no "$CHR_FA" "$LEPC_REF" > "$PAF"
fi

# Per-scaffold assignment: aligned matching bases per chicken chromosome
ASSIGN="$OUTDIR/lepc_scaffold_to_chicken.tsv"
awk -F'\t' -v OFS='\t' -v mq="$MIN_MAPQ" -v minf="$MIN_FRAC" -v mina="$MIN_ALIGNED" '
    NR == FNR { chr[$1] = $2; next }
    ($12 >= mq) && /tp:A:P/ {
        c = ($6 in chr) ? chr[$6] : $6
        bp[$1, c] += $10; tot[$1] += $10; qlen[$1] = $2; seen[$1, c] = 1; chrs[c] = 1
    }
    END {
        print "scaffold", "length", "aligned_bp", "best_chicken_chr", "best_frac", "assignment"
        for (q in tot) {
            best = ""; bmax = 0
            for (c in chrs) if (((q, c) in seen) && bp[q, c] > bmax) { bmax = bp[q, c]; best = c }
            frac = bmax / tot[q]
            a = "unassigned"
            if (tot[q] >= mina && frac >= minf) a = (best == "Z" ? "Z" : (best == "W" ? "W" : "autosome"))
            print q, qlen[q], tot[q], best, sprintf("%.3f", frac), a
        }
    }' "$CHR_MAP" "$PAF" | { IFS= read -r h; echo "$h"; sort -t$'\t' -k2,2nr; } > "$ASSIGN"

awk 'NR>1{print $6}' "$ASSIGN" | sort | uniq -c | awk '{printf "[1]   %-10s %d scaffolds\n", $2, $1}'

awk -F'\t' 'NR>1 && $6=="Z" {print $1}' "$ASSIGN" > "$OUTDIR/Z_synteny.txt"
awk -F'\t' 'NR>1 && $6=="W" {print $1}' "$ASSIGN" > "$OUTDIR/W_synteny.txt"
echo "[1] Synteny Z scaffolds ($(wc -l < "$OUTDIR/Z_synteny.txt")):"
awk -F'\t' 'NR==1 || $6=="Z" || $6=="W"' "$ASSIGN" | column -t

[[ -s "$OUTDIR/Z_synteny.txt" ]] || { echo "ERROR: no Z scaffolds found by synteny; check the PAF"; exit 1; }

# --- Step 2: confirm with genotypes, call sex -----------------------------------
# Clean environment for Python: drop biocontainers (minimap2 is done), keep
# threads and matplotlib's cache local to the job.
module --force purge
module load anaconda
set +u; source activate lepc_py; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
export MPLCONFIGDIR="${TMPDIR:-/tmp}/mpl_${SLURM_JOB_ID:-$$}"; mkdir -p "$MPLCONFIGDIR"
echo "[2] $(date): hemizygosity check in VCF"
W_EXCL=$(paste -sd, "$OUTDIR/W_synteny.txt")
python -X faulthandler sex_from_vcf.py \
    --vcf "$VCF" \
    --popmap "$POPMAP" \
    --z_scaffolds "$OUTDIR/Z_synteny.txt" \
    --exclude "$W_EXCL" \
    --outprefix "$OUTDIR/vcf"
set +u; conda deactivate; set -u

# --- Final list: VCF-confirmed Z if available, else synteny Z; plus W ----------
if [[ -s "$OUTDIR/vcf_Z_confirmed.txt" ]]; then
    Z_FINAL="$OUTDIR/vcf_Z_confirmed.txt"
    echo "[final] Using VCF-confirmed Z list, plus any synteny W"
else
    Z_FINAL="$OUTDIR/Z_synteny.txt"
    echo "[final] No female/male contrast available; using synteny Z list, plus any synteny W"
fi
cat "$Z_FINAL" "$OUTDIR/W_synteny.txt" | awk 'NF' | sort -u | paste -sd, > "$OUTDIR/SEX_CHROMS.txt"
echo "[final] SEX_CHROMS=\"$(cat "$OUTDIR/SEX_CHROMS.txt")\""
echo "== Finished at $(date) =="
