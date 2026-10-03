#!/bin/bash
# =============================================================================
# angsd_fst_temporal.sh -- ANGSD pairwise FST, genome-wide and in autosomal sliding windows
# Temporal NM LEPC: Past (2019, n=9) vs Present (2026, n=10)
#
# Stages (chained with SLURM dependencies by "submit"):
#   sites  [SITES_MODE=snps only] array over region chunks: SNP discovery on
#          ALL samples together (-SNP_pval), so every group is scored at the
#          same polymorphic sites
#   saf    array over (group x chunk): angsd -doSaf, autosomes only
#   merge  array over groups: realSFS cat -> one autosomal SAF per group
#   fst    array over group pairs: folded 2D-SFS prior -> realSFS fst index
#          (-whichFst 1 = Hudson, Bhatia et al. 2013) -> genome-wide stats
#          and sliding windows
#
# Each task skips work whose output already exists, and writes to a .tmp
# name first, so rerunning "submit" after a failure only redoes what is
# missing and a killed task never leaves a half-written file behind.
#
# USAGE (login node, from the directory holding this script):
#   bash angsd_fst_temporal.sh check          # build/verify group lists + chunks, no jobs
#   bash angsd_fst_temporal.sh submit         # submit all stages with dependencies
#   bash angsd_fst_temporal.sh submit fst     # resubmit from a later stage only
#
# OUTPUT (${OUT}):
#   <A>_<B>.fst.global.txt                genome-wide FST (unweighted, weighted)
#   <A>_<B>.fst.windows_<WIN>_<STEP>.txt  region chr midPos Nsites fst
#   <A>_<B>.2dsfs.ml, <A>_<B>.fst.idx/.gz per-site files (for reruns/plots)
# =============================================================================
#SBATCH --job-name=angsd_fst_temporal
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A fnrdewoody
#SBATCH -p cpu
#SBATCH -t 2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=blackan@purdue.edu
set -euo pipefail

# =============================================================================
# USER SETTINGS
# =============================================================================
# Temporal: NM LEPC 2019 (Past, n=9; F10 excluded) vs 2026 (Present, n=10).
# n is small, so all autosomal sites (variable and invariant) are used.
PROJECT_DIR="${CLUSTER_SCRATCH}/GROUSE/old_vs_new"
REF_FASTA="${PROJECT_DIR}/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna"
CRAMLIST="${PROJECT_DIR}/final_cramlist.txt"       # 20 CRAMs; F10 is dropped
POPMAP="${PROJECT_DIR}/popmap_unrelated.txt"       # normal_<ID><TAB>Past|Present
RUN_NAME="temporal_past_present"
JOB="fst_temporal"
GROUPS_ORDER=(Past Present)
declare -A EXPECTED_N=([Past]=9 [Present]=10)     # stop if the lists differ
SITES_MODE=all          # all = every autosomal site; snps = jointly called SNPs only
N_CHUNKS=40
# 2D-SFS prior from <=100 M sites (keeps realSFS memory bounded on ~1 Gb)
SFS_EXTRA="-nSites 100000000 -maxIter 200"
# per-stage resources: cpus mem time
SITES_RES=(8 32G 1-00:00:00)
SAF_RES=(8 32G 1-00:00:00)
FST_RES=(16 96G 2-00:00:00)

# ---- shared ANGSD filters (same read filters as the heterozygosity jobs) ---
GL=1                     # samtools GL model
MINMAPQ=30
MINQ=30
MIN_IND_FRAC=0.5         # keep a site if >= this fraction of the group has reads
SNP_PVAL=1e-6            # SITES_MODE=snps: SNP discovery threshold
WIN=50000                # sliding window (bp)
STEP=10000               # window step (bp)
Z_SCAFFOLDS="NW_026294758.1,NW_026294813.1"
MAX_PARALLEL=50          # max simultaneous array tasks per stage

OUT="${PROJECT_DIR}/fst/${RUN_NAME}"
AUTOSOME_RF="${OUT}/autosomes.rf"
THREADS=${SLURM_CPUS_PER_TASK:-4}
mkdir -p "$OUT"/{lists,chunks,sites,saf} logs

# =============================================================================
# HELPERS
# =============================================================================
# Sample ID from a CRAM path: basename up to the first dot
# (F17.md.dedup_q20.cram -> F17).
cram_id() { local b; b=$(basename "$1"); echo "${b%%.*}"; }
# Sample ID from a popmap: drop a leading "normal_" (VCF-style names).
pop_id()  { echo "${1#normal_}"; }

pairs() {   # every unordered pair of GROUPS_ORDER, "A B" per line
    local i j
    for ((i=0; i<${#GROUPS_ORDER[@]}; i++)); do
        for ((j=i+1; j<${#GROUPS_ORDER[@]}; j++)); do
            echo "${GROUPS_ORDER[$i]} ${GROUPS_ORDER[$j]}"
        done
    done
}

in_groups() { local g; for g in "${GROUPS_ORDER[@]}"; do [[ "$g" == "$1" ]] && return 0; done; return 1; }

# Per-group CRAM lists from CRAMLIST + POPMAP, checked against EXPECTED_N.
build_lists() {
    [[ -s "$CRAMLIST" ]] || { echo "ERROR: CRAMLIST not found: $CRAMLIST" >&2; exit 1; }
    [[ -s "$POPMAP"   ]] || { echo "ERROR: POPMAP not found: $POPMAP" >&2; exit 1; }
    declare -A GROUP_OF=()
    local id grp rest cram
    # Column 1 = sample ID, everything after the first tab (or first run of
    # spaces) = group label. Labels can contain spaces, so a hybrid written
    # "STGR / GRPC" stays one label, matches no group and is left out.
    local line
    while IFS= read -r line; do
        line="${line%%$'\r'}"
        [[ -z "${line// }" || "${line:0:1}" == "#" ]] && continue
        if [[ "$line" == *$'\t'* ]]; then id="${line%%$'\t'*}"; grp="${line#*$'\t'}"
        else read -r id grp <<< "$line"; fi
        grp="$(sed -E 's/^[[:space:]]+|[[:space:]]+$//g' <<< "$grp")"
        GROUP_OF["$(pop_id "$id")"]="$grp"
    done < "$POPMAP"

    local -a no_pop=() other=()
    for g in "${GROUPS_ORDER[@]}"; do : > "${OUT}/lists/${g}.cramlist.tmp"; done
    while read -r cram; do
        [[ -z "$cram" ]] && continue
        id=$(cram_id "$cram")
        grp="${GROUP_OF[$id]:-}"
        if   [[ -z "$grp" ]];       then no_pop+=("$id")
        elif ! in_groups "$grp";    then other+=("$id:$grp")
        else echo "$cram" >> "${OUT}/lists/${grp}.cramlist.tmp"
        fi
    done < <(tr -d '\r' < "$CRAMLIST")

    (( ${#no_pop[@]} )) && echo "NOTE: ${#no_pop[@]} CRAM(s) not in popmap, left out: ${no_pop[*]:0:20}"
    (( ${#other[@]}  )) && echo "NOTE: ${#other[@]} CRAM(s) in a group not analysed, left out: ${other[*]:0:20}"
    local n bad=0
    for g in "${GROUPS_ORDER[@]}"; do
        n=$(wc -l < "${OUT}/lists/${g}.cramlist.tmp")
        if [[ -n "${EXPECTED_N[$g]:-}" && "$n" != "${EXPECTED_N[$g]}" ]]; then
            echo "ERROR: group $g has $n CRAMs, expected ${EXPECTED_N[$g]}" >&2; bad=1
        fi
        (( n > 1 )) || { echo "ERROR: group $g has $n CRAMs" >&2; bad=1; }
    done
    (( bad == 0 )) || exit 1
    for g in "${GROUPS_ORDER[@]}"; do
        mv -f "${OUT}/lists/${g}.cramlist.tmp" "${OUT}/lists/${g}.cramlist"
    done
    cat $(printf "${OUT}/lists/%s.cramlist " "${GROUPS_ORDER[@]}") > "${OUT}/lists/ALL.cramlist"
    while read -r cram; do
        [[ -s "$cram" ]] || { echo "ERROR: CRAM missing on disk: $cram" >&2; exit 1; }
        [[ -s "${cram}.crai" || -s "${cram%.cram}.crai" ]] || { echo "ERROR: no index for $cram" >&2; exit 1; }
    done < "${OUT}/lists/ALL.cramlist"
}

# Autosomal contigs (every contig in the .fai except the Z scaffolds), split
# into N_CHUNKS region files balanced by total length (largest first).
build_chunks() {
    [[ -s "${OUT}/chunks/chunk_000.rf" ]] && return 0
    [[ -f "${REF_FASTA}.fai" ]] || { echo "ERROR: ${REF_FASTA}.fai not found" >&2; exit 1; }
    local c
    for c in ${Z_SCAFFOLDS//,/ }; do
        awk -v c="$c" '$1==c {f=1} END {exit !f}' "${REF_FASTA}.fai" \
            || { echo "ERROR: Z scaffold $c not in ${REF_FASTA}.fai" >&2; exit 1; }
    done
    awk -v z="$Z_SCAFFOLDS" 'BEGIN{n=split(z,a,","); for(i=1;i<=n;i++) Z[a[i]]=1}
        !($1 in Z) {print $1"\t"$2}' "${REF_FASTA}.fai" | sort -k2,2nr > "${OUT}/autosomes.len"
    awk '{print $1":"}' "${OUT}/autosomes.len" > "$AUTOSOME_RF"
    awk -v n="$N_CHUNKS" -v d="${OUT}/chunks" '
        { b=0; for(i=1;i<n;i++) if (t[i] < t[b]) b=i
          t[b]+=$2; f=sprintf("%s/chunk_%03d.rf", d, b); print $1":" >> f; close(f) }' "${OUT}/autosomes.len"
    echo "  autosomes: $(wc -l < "$AUTOSOME_RF") contigs, $(awk '{s+=$2} END{printf "%.0f", s/1e6}' "${OUT}/autosomes.len") Mb, $(ls "${OUT}"/chunks/chunk_*.rf | wc -l) chunks"
}

load_angsd() {
    module --force purge 2>/dev/null || true
    ml biocontainers
    ml angsd/0.940
    # xalt's LD_PRELOAD breaks the container (GLIBC_2.33/2.34 errors)
    unset LD_PRELOAD || true
    export SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""
    # Read CRAM reference sequences from a local cache instead of the EBI
    # server (the slow hts-ref fetches seen in the ROH job). Build it once:
    #   seq_cache_populate.pl -root ${PROJECT_DIR}/ref/hts-cache ${REF_FASTA}
    local cache="${PROJECT_DIR}/ref/hts-cache"
    if [[ -d "$cache" ]]; then
        export REF_CACHE="${cache}/%2s/%2s/%s" REF_PATH="${cache}/%2s/%2s/%s"
        export SINGULARITYENV_REF_CACHE="$REF_CACHE" SINGULARITYENV_REF_PATH="$REF_PATH"
        export APPTAINERENV_REF_CACHE="$REF_CACHE" APPTAINERENV_REF_PATH="$REF_PATH"
    fi
}

min_ind() { awk -v n="$1" -v f="$MIN_IND_FRAC" 'BEGIN{m=int(n*f+0.5); print (m<1?1:m)}'; }
n_chunks_real() { ls "${OUT}"/chunks/chunk_*.rf | wc -l; }

STAGE="${1:-}"

# =============================================================================
# check / submit (login node)
# =============================================================================
if [[ "$STAGE" == "check" || "$STAGE" == "submit" ]]; then
    build_lists
    build_chunks
    G=${#GROUPS_ORDER[@]}; NC=$(n_chunks_real); NP=$(pairs | wc -l)
    for g in "${GROUPS_ORDER[@]}"; do
        n=$(wc -l < "${OUT}/lists/${g}.cramlist")
        echo "  $g: $n CRAMs (minInd $(min_ind "$n"))"
    done
    echo "  sites: ${SITES_MODE}   pairs: $(pairs | paste -sd ';')"
    echo "  output: $OUT"
    [[ "$STAGE" == "check" ]] && exit 0

    FROM="${2:-$([[ "$SITES_MODE" == snps ]] && echo sites || echo saf)}"
    SELF="$(readlink -f "$0")"
    DEP=""; started=0
    sub() {  # sub <stage> <array> <cpus> <mem> <time>
        local id
        id=$(sbatch --parsable ${DEP:+--dependency=afterok:$DEP} --job-name="${JOB}_$1" \
             --array="$2" --cpus-per-task="$3" --mem="$4" -t "$5" "$SELF" "$1")
        echo "  $1: job $id (array $2)"; DEP="$id"
    }
    for s in sites saf merge fst; do
        [[ "$s" == "$FROM" ]] && started=1
        (( started )) || continue
        case $s in
          sites) [[ "$SITES_MODE" == snps ]] && sub sites "0-$((NC-1))%${MAX_PARALLEL}" "${SITES_RES[@]}" ;;
          saf)   sub saf   "0-$((G*NC-1))%${MAX_PARALLEL}" "${SAF_RES[@]}" ;;
          merge) sub merge "0-$((G-1))" 4 32G 1-00:00:00 ;;
          fst)   sub fst   "0-$((NP-1))" "${FST_RES[@]}" ;;
        esac
    done
    (( started )) || { echo "ERROR: unknown stage '$FROM' (sites|saf|merge|fst)" >&2; exit 1; }
    exit 0
fi

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || { echo "Run with: bash $0 submit   (or: bash $0 check)" >&2; exit 1; }
TASK=$SLURM_ARRAY_TASK_ID
echo ">>> stage=$STAGE task=$TASK host=$(hostname) start=$(date)"

# =============================================================================
# sites (task = chunk): SNPs called on all samples together
# =============================================================================
if [[ "$STAGE" == "sites" ]]; then
    CH=$(printf "%03d" "$TASK"); RF="${OUT}/chunks/chunk_${CH}.rf"
    S="${OUT}/sites/chunk_${CH}"
    if [[ -s "${S}.sites.idx" ]]; then echo "${S}.sites exists -- skipping"; exit 0; fi
    N=$(wc -l < "${OUT}/lists/ALL.cramlist")
    load_angsd
    angsd -bam "${OUT}/lists/ALL.cramlist" -ref "$REF_FASTA" -rf "$RF" \
        -GL "$GL" -doMajorMinor 1 -doMaf 1 -SNP_pval "$SNP_PVAL" \
        -minMapQ "$MINMAPQ" -minQ "$MINQ" -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
        -minInd "$(min_ind "$N")" -P "$THREADS" -out "${S}.tmp"
    zcat "${S}.tmp.mafs.gz" | awk 'NR>1 {print $1"\t"$2}' > "${S}.sites.tmp"
    mv -f "${S}.sites.tmp" "${S}.sites"
    angsd sites index "${S}.sites"
    rm -f "${S}.tmp.mafs.gz"; mv -f "${S}.tmp.arg" "${S}.arg" 2>/dev/null || true
    echo ">>> chunk $CH: $(wc -l < "${S}.sites") SNPs  $(date)"
    exit 0
fi

# =============================================================================
# saf (task = group x chunk)
# =============================================================================
if [[ "$STAGE" == "saf" ]]; then
    NC=$(n_chunks_real)
    GRP="${GROUPS_ORDER[$((TASK / NC))]}"; CH=$(printf "%03d" $((TASK % NC)))
    RF="${OUT}/chunks/chunk_${CH}.rf"; LIST="${OUT}/lists/${GRP}.cramlist"
    PFX="${OUT}/saf/${GRP}.chunk_${CH}"
    if [[ -s "${PFX}.saf.idx" ]]; then echo "${PFX}.saf.idx exists -- skipping"; exit 0; fi
    SITES_ARG=()
    if [[ "$SITES_MODE" == snps ]]; then
        S="${OUT}/sites/chunk_${CH}.sites"
        [[ -s "${S}.idx" ]] || { echo "ERROR: $S not indexed (sites stage incomplete?)" >&2; exit 1; }
        SITES_ARG=(-sites "$S")
    fi
    N=$(wc -l < "$LIST")
    echo ">>> saf: group=$GRP n=$N minInd=$(min_ind "$N") chunk=$CH"
    load_angsd
    angsd -bam "$LIST" -ref "$REF_FASTA" -anc "$REF_FASTA" -rf "$RF" "${SITES_ARG[@]}" \
        -doSaf 1 -GL "$GL" \
        -minMapQ "$MINMAPQ" -minQ "$MINQ" -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
        -minInd "$(min_ind "$N")" -P "$THREADS" -out "${PFX}.tmp"
    # realSFS locates .saf.gz/.saf.pos.gz from the .saf.idx prefix, so the
    # three files can be renamed together; idx last marks the chunk complete.
    mv -f "${PFX}.tmp.saf.gz" "${PFX}.saf.gz"
    mv -f "${PFX}.tmp.saf.pos.gz" "${PFX}.saf.pos.gz"
    mv -f "${PFX}.tmp.arg" "${PFX}.arg"
    mv -f "${PFX}.tmp.saf.idx" "${PFX}.saf.idx"
    echo ">>> done $(date)"
    exit 0
fi

# =============================================================================
# merge (task = group): concatenate chunk SAFs
# =============================================================================
if [[ "$STAGE" == "merge" ]]; then
    GRP="${GROUPS_ORDER[$TASK]}"; M="${OUT}/${GRP}"
    if [[ -s "${M}.saf.idx" ]]; then echo "${M}.saf.idx exists -- skipping"; exit 0; fi
    mapfile -t PARTS < <(ls "${OUT}"/saf/"${GRP}".chunk_*.saf.idx 2>/dev/null | sort)
    NC=$(n_chunks_real)
    (( ${#PARTS[@]} == NC )) || { echo "ERROR: $GRP has ${#PARTS[@]} of $NC chunk SAFs" >&2; exit 1; }
    load_angsd
    realSFS cat "${PARTS[@]}" -outnames "${M}.tmp"
    mv -f "${M}.tmp.saf.gz" "${M}.saf.gz"; mv -f "${M}.tmp.saf.pos.gz" "${M}.saf.pos.gz"
    mv -f "${M}.tmp.saf.idx" "${M}.saf.idx"
    echo ">>> merged $GRP ($(date))"
    exit 0
fi

# =============================================================================
# fst (task = pair)
# =============================================================================
if [[ "$STAGE" == "fst" ]]; then
    read -r A B < <(pairs | sed -n "$((TASK+1))p")
    P="${OUT}/${A}_${B}"
    for g in "$A" "$B"; do
        [[ -s "${OUT}/${g}.saf.idx" ]] || { echo "ERROR: ${OUT}/${g}.saf.idx missing" >&2; exit 1; }
    done
    load_angsd
    # 1) folded 2D-SFS: the prior for the per-site FST estimates.
    #    (The reference is used as "ancestral", so the SFS is folded.)
    if [[ ! -s "${P}.2dsfs.ml" ]]; then
        echo ">>> 2D-SFS $A x $B  $(date)"
        realSFS "${OUT}/${A}.saf.idx" "${OUT}/${B}.saf.idx" -fold 1 -P "$THREADS" ${SFS_EXTRA} \
            > "${P}.2dsfs.ml.tmp"
        mv -f "${P}.2dsfs.ml.tmp" "${P}.2dsfs.ml"
    fi
    # 2) per-site numerators/denominators, Hudson FST (Bhatia et al. 2013)
    if [[ ! -s "${P}.fst.idx" ]]; then
        echo ">>> fst index $A x $B  $(date)"
        realSFS fst index "${OUT}/${A}.saf.idx" "${OUT}/${B}.saf.idx" \
            -sfs "${P}.2dsfs.ml" -fold 1 -whichFst 1 -P "$THREADS" -fstout "${P}.tmp"
        mv -f "${P}.tmp.fst.gz" "${P}.fst.gz"; mv -f "${P}.tmp.fst.idx" "${P}.fst.idx"
    fi
    # 3) genome-wide (ratio of averages = "weighted") and sliding windows
    realSFS fst stats "${P}.fst.idx" > "${P}.fst.global.txt" 2> "${P}.fst.global.log"
    W="${P}.fst.windows_${WIN}_${STEP}.txt"
    realSFS fst stats2 "${P}.fst.idx" -win "$WIN" -step "$STEP" -type 2 > "${W}.tmp"
    mv -f "${W}.tmp" "$W"
    echo ">>> $A x $B  FST unweighted / weighted: $(tail -n1 "${P}.fst.global.txt")"
    echo ">>> $(($(wc -l < "$W") - 1)) windows -> $W"
    exit 0
fi

echo "Unknown stage '$STAGE' (use: bash $0 check | submit [stage])" >&2
exit 1
