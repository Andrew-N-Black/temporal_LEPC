#!/bin/bash
# ============================================================================
# SLURM wrapper: relatedness among Past and Present LEPC samples
#
#   Arm A (VCF, fast): KING-robust kinship, R0, R1 from called genotypes
#          (allele-frequency-free; Manichaikul et al. 2010; Waples et al. 2019)
#          -> relatedness_vcf.py, lepc_py conda env
#
#   Arm B (CRAMs, genotype likelihoods): ANGSD genotype likelihoods, then
#          NgsRelate v2 (Hanghoej et al. 2019), which reports rab, theta,
#          KING, R0, R1 without calling genotypes. Preferred for low-coverage
#          or historical samples, where called genotypes under-call
#          heterozygotes. Run separately for Past, Present, and All samples
#          (NgsRelate's rab/theta use allele frequencies from the run's samples).
#          Transitions are excluded (-noTrans 1) to limit post-mortem damage.
#
# Sex chromosomes are excluded in both arms (hemizygous Z in ZW females biases
# heterozygosity-based statistics). Set SEX_CHROMS to your assembly's names.
#
# One-time environment setup (login node):
#   module load anaconda
#   conda create -y -n lepc_py  -c conda-forge python=3.11 scikit-allel numpy pandas matplotlib
#   conda create -y -n lepc_ngs -c conda-forge -c bioconda angsd ngsrelate samtools
# ============================================================================
#SBATCH -J lepc_relatedness
#SBATCH -o lepc_relatedness_%j.out
#SBATCH -e lepc_relatedness_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH -A dewoody
#SBATCH -p cpu

set -euo pipefail

module --force purge
module load anaconda
cd "$SLURM_SUBMIT_DIR"

# --- Edit these paths/params for your run -----------------------------------
RUN_VCF=TRUE                  # Arm A
RUN_NGSRELATE=FALSE            # Arm B (needs CRAMs and REF)

VCF="/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs/output.subset.biallelic.snps.vcf.gz"
POPMAP="popmap.txt"           # <sample_id> <Past|Present>, no header
# Z scaffolds from find_sex_scaffolds.sh (chicken synteny); reference bird is male, so no W
SEX_CHROMS="NW_026294758.1,NW_026294813.1"

REF="/scratch/gautschi/blackan/GROUSE/old_vs_new/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
# CRAMs: give CRAMS_TSV (<sample_id> <path>) directly, or leave it missing and set
# CRAM_DIR: the wrapper then looks for <sample_id minus ID_STRIP><CRAM_SUFFIX>,
# e.g. normal_F340 -> F340.md.dedup_q20.cram
CRAMS_TSV="crams.txt"
CRAM_DIR="/scratch/gautschi/blackan/GROUSE/old_vs_new/crams"
ID_STRIP="normal_"
CRAM_SUFFIX=".md.dedup_q20.cram"

N_CONTIGS=5                   # N longest autosomal scaffolds for ANGSD (~400 Mb; plenty for NgsRelate)
MIN_IND_FRAC=0.8              # ANGSD: site must have data in this fraction of samples
MIN_MAF=0.05                  # ANGSD: minor allele frequency for NgsRelate sites
MIN_MAPQ=20                   # matches the q20 filtering already applied to the CRAMs
ONLY_PROPER_PAIRS=1           # set 0 if any samples have collapsed/merged (single-end) reads

OUTDIR="results/relatedness"
THREADS="${SLURM_CPUS_PER_TASK:-1}"
# -----------------------------------------------------------------------------

mkdir -p "$OUTDIR"
echo "== Job $SLURM_JOB_ID on $(hostname) at $(date) =="
[[ -s "$POPMAP" ]] || { echo "ERROR: popmap not found: $POPMAP"; exit 1; }

# --- Arm A: VCF, called genotypes --------------------------------------------
if [[ "$RUN_VCF" == "TRUE" ]]; then
    for f in "$VCF" relatedness_vcf.py; do
        [[ -s "$f" ]] || { echo "ERROR: required file not found: $f"; exit 1; }
    done
    set +u; source activate lepc_py; set -u
    echo "[A] $(date): KING/R0/R1 from VCF"
    python relatedness_vcf.py \
        --vcf "$VCF" \
        --popmap "$POPMAP" \
        --exclude_chroms "$SEX_CHROMS" \
        --outprefix "$OUTDIR/vcf"
    set +u; conda deactivate; set -u
fi

# --- Arm B: CRAMs, genotype likelihoods --------------------------------------
if [[ "$RUN_NGSRELATE" == "TRUE" ]]; then
    for f in "$REF" "${REF}.fai"; do
        [[ -s "$f" ]] || { echo "ERROR: required file not found: $f"; exit 1; }
    done
    # Build the sample -> CRAM table from CRAM_DIR if it doesn't exist yet
    if [[ ! -s "$CRAMS_TSV" ]]; then
        [[ -d "$CRAM_DIR" ]] || { echo "ERROR: no $CRAMS_TSV and CRAM_DIR not found: $CRAM_DIR"; exit 1; }
        echo "[B] Building $CRAMS_TSV from $CRAM_DIR"
        : > "$CRAMS_TSV"
        while read -r sid _; do
            sid="${sid%$'\r'}"; [[ -n "$sid" ]] || continue
            cram="$CRAM_DIR/${sid#"$ID_STRIP"}${CRAM_SUFFIX}"
            [[ -s "$cram" ]] || { echo "ERROR: CRAM not found for $sid: $cram"; rm -f "$CRAMS_TSV"; exit 1; }
            printf '%s\t%s\n' "$sid" "$cram" >> "$CRAMS_TSV"
        done < "$POPMAP"
    fi
    echo "[B] CRAMs:"; column -t "$CRAMS_TSV"
    set +u; source activate lepc_ngs; set -u
    command -v angsd >/dev/null     || { echo "ERROR: angsd not found (lepc_ngs env?)"; exit 1; }
    command -v ngsRelate >/dev/null || { echo "ERROR: ngsRelate not found (lepc_ngs env?)"; exit 1; }

    # Region file: N longest contigs, excluding sex chromosomes
    REGIONS="$OUTDIR/angsd_regions.txt"
    awk -v sex="$SEX_CHROMS" 'BEGIN{n=split(sex,s,","); for(i=1;i<=n;i++) ex[s[i]]=1}
         !($1 in ex) {print $1"\t"$2}' "${REF}.fai" \
      | sort -k2,2nr | head -n "$N_CONTIGS" | awk '{print $1}' > "$REGIONS"
    echo "[B] Using $(wc -l < "$REGIONS") contigs for ANGSD: $(paste -sd, "$REGIONS")"

    # Per-group CRAM lists and ID files, in matching order
    make_lists () {   # $1 = group name, $2 = awk condition on popmap label
        local grp="$1" cond="$2"
        awk -v OFS='\t' 'NR==FNR {gsub(/\r/,""); lab=tolower($2); pop[$1]=lab; next}
             { gsub(/\r/,""); if ($1 in pop) print $1, $2, pop[$1] }' "$POPMAP" "$CRAMS_TSV" \
          | awk -F'\t' "$cond" > "$OUTDIR/${grp}.map"
        cut -f2 "$OUTDIR/${grp}.map" > "$OUTDIR/${grp}.bamlist"
        cut -f1 "$OUTDIR/${grp}.map" > "$OUTDIR/${grp}.ids"
        local n; n=$(wc -l < "$OUTDIR/${grp}.bamlist")
        [[ "$n" -ge 2 ]] || { echo "ERROR: $grp has $n CRAMs matched to popmap"; exit 1; }
        while read -r c; do [[ -s "$c" ]] || { echo "ERROR: CRAM not found: $c"; exit 1; }; done \
          < "$OUTDIR/${grp}.bamlist"
        echo "[B] $grp: $n samples"
    }
    make_lists past    '$3=="past" || $3=="old"'
    make_lists present '$3=="present" || $3=="new"'
    make_lists all     '1'

    for grp in past present all; do
        n=$(wc -l < "$OUTDIR/${grp}.bamlist")
        min_ind=$(awk -v n="$n" -v f="$MIN_IND_FRAC" 'BEGIN{m=int(n*f+0.5); print (m<2?2:m)}')
        echo "[B] $(date): ANGSD $grp (n=$n, minInd=$min_ind)"
        angsd -b "$OUTDIR/${grp}.bamlist" -ref "$REF" -rf "$REGIONS" \
            -out "$OUTDIR/${grp}" -P "$THREADS" \
            -gl 2 -doGlf 3 -doMajorMinor 1 -doMaf 1 \
            -SNP_pval 1e-6 -minMaf "$MIN_MAF" -noTrans 1 \
            -minMapQ "$MIN_MAPQ" -minQ 20 -C 50 -baq 1 \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs "$ONLY_PROPER_PAIRS" \
            -minInd "$min_ind" -setMinDepthInd 1

        # Allele frequencies for NgsRelate: the knownEM column, found by header name
        zcat "$OUTDIR/${grp}.mafs.gz" \
          | awk 'NR==1 {for (i=1;i<=NF;i++) if ($i=="knownEM") c=i; if (!c) exit 1; next} {print $c}' \
          > "$OUTDIR/${grp}.freq"
        echo "[B] $grp: $(wc -l < "$OUTDIR/${grp}.freq") SNPs"

        echo "[B] $(date): NgsRelate $grp"
        ngsRelate -g "$OUTDIR/${grp}.glf.gz" -n "$n" -f "$OUTDIR/${grp}.freq" \
            -z "$OUTDIR/${grp}.ids" -p "$THREADS" -O "$OUTDIR/${grp}.ngsrelate.res"

        # Compact summary, sorted by KING, with degree calls
        awk -v OFS='\t' '
            NR==1 {for (i=1;i<=NF;i++) h[$i]=i
                   print "id1","id2","nSites","rab","theta","KING","R0","R1","degree"; next}
            { k=$h["KING"]
              d = (k>0.354)?"duplicate/identical":(k>0.177)?"1st degree":(k>0.0884)?"2nd degree":(k>0.0442)?"3rd degree":"unrelated"
              print $h["ida"],$h["idb"],$h["nSites"],$h["rab"],$h["theta"],k,$h["R0"],$h["R1"],d }' \
          "$OUTDIR/${grp}.ngsrelate.res" \
          | { IFS= read -r header; echo "$header"; sort -t$'\t' -k6,6gr; } \
          > "$OUTDIR/${grp}.ngsrelate.summary.tsv"
        echo "[B] $grp related pairs (KING > 0.0442):"
        awk -F'\t' 'NR==1 || $9!="unrelated"' "$OUTDIR/${grp}.ngsrelate.summary.tsv" | column -t
    done
    set +u; conda deactivate; set -u
fi

echo "== Finished at $(date) =="
