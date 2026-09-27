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
VCF=/scratch/gautschi/blackan/GROUSE/old_vs_new/vcfs/analysis.autosomes.unrel.vcf.gz

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
module load bedtools
set -euo pipefail

echo "Aligning chicken (Gallus gallus) to LEPC reference: $(date)"

# Auto-detect the chicken FASTA in the pre-downloaded reference directory
# rather than assuming a filename — GALLUS_DIR also holds the LEPC reference
# itself in this setup, so explicitly exclude the LEPC reference from the candidates
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
CHICKEN_BAM=OUT/ancestral/gallus_to_lepc.bam

# The alignment dominates the runtime of this step and its output is reusable,
# so skip it when a complete indexed BAM is already present. This makes
# re-running Step 1 purely to rebuild the consensus/mask cheap.
# `samtools quickcheck` is what makes reuse safe: it verifies the header and the
# bgzf EOF block, so a BAM left half-written by a cancelled or timed-out job is
# realigned rather than silently reused. A plain -s test would accept it.
if [ -s "$CHICKEN_BAM" ] && [ -s "${CHICKEN_BAM}.bai" ] \
        && samtools quickcheck -q "$CHICKEN_BAM"; then
    echo "Reusing existing chicken alignment: $CHICKEN_BAM"
else
    if [ -s "$CHICKEN_BAM" ]; then
        echo "Existing $CHICKEN_BAM failed quickcheck (truncated/incomplete) -- realigning."
    fi
    minimap2 -ax asm10 -t THREADS REF "$CHICKEN_FASTA" | \
        samtools sort -@ THREADS -o "$CHICKEN_BAM" -
    samtools index "$CHICKEN_BAM"
fi

# Call a base at every LEPC position the chicken alignment covers. Without -v,
# `bcftools call -c` emits a record for EVERY pileup position, not just variant
# ones, so the position set of that file is exactly the set of LEPC coordinates
# for which chicken data exists.
CONS_CALLS=OUT/ancestral/gallus_consensus.vcf.gz

# Also reusable, and also expensive. Note that `bcftools index --nrecords` is
# NOT a valid integrity test here: it reads the index, not the data, so a stale
# index reports a full record count for a half-written file. The bgzf EOF
# marker is the real O(1) truncation test -- htslib appends a fixed 28-byte
# empty block to every complete bgzf file, so its absence means the writing job
# was cut short.
BGZF_EOF_HEX="1f8b08040000000000ff0600424302001b0003000000000000000000"
CALLS_OK=0
if [ -s "$CONS_CALLS" ] && [ -s "${CONS_CALLS}.tbi" ]; then
    if [ "$(tail -c 28 "$CONS_CALLS" | od -An -tx1 -v | tr -d ' \n')" = "$BGZF_EOF_HEX" ]; then
        CALLS_OK=1
    else
        echo "Existing $CONS_CALLS lacks its bgzf EOF marker (truncated) -- recalling."
    fi
fi
if [ "$CALLS_OK" -eq 1 ]; then
    echo "Reusing existing chicken consensus calls: $CONS_CALLS"
    echo "  records: $(bcftools index --nrecords "$CONS_CALLS")"
else
    bcftools mpileup -f REF "$CHICKEN_BAM" -Ou | \
        bcftools call -c --ploidy 1 -Oz -o "$CONS_CALLS"
    tabix -p vcf "$CONS_CALLS"
fi

# ---------------------------------------------------------------------------
# CRITICAL: `bcftools consensus` only EDITS the reference where that call set
# has a record. Every position the chicken alignment never reached is emitted
# unchanged -- i.e. as the LEPC reference base -- which is indistinguishable
# downstream from a confident "ancestral equals the reference base" call.
# Left unmasked, that
# silently reimposes the naive assumption that the reference allele is
# ancestral across the majority of the genome (chicken covers only ~40% of the
# LEPC assembly), inflating derived-allele counts and therefore every load
# metric. An unmasked run of this step produced only 2 "no ancestral call"
# sites out of 13.2M -- implausible, and the tell that this was happening.
#
# Fix: derive the complement of the called positions and hard-mask it to N.
# Because the mask is built from that same call set, the invariant is
# exact -- every non-N base in the ancestral FASTA is backed by a chicken
# record. Step 4 already skips N ancestral calls, so this is what makes the
# polarization honest.
# ---------------------------------------------------------------------------
GENOME_SIZES=OUT/ancestral/lepc_genome_sizes.txt
COVERED_BED=OUT/ancestral/gallus_covered.bed
INDEL_PAD_BED=OUT/ancestral/gallus_indel_pad.bed
RELIABLE_BED=OUT/ancestral/gallus_reliable.bed
UNCOV_BED=OUT/ancestral/gallus_uncovered.bed
INDEL_RAW_BED=OUT/ancestral/gallus_indel_raw.bed
CONS_SNV_CALLS=OUT/ancestral/gallus_consensus.snvonly.vcf.gz
INDEL_PAD_BP=10

cut -f1,2 REF.fai > "$GENOME_SIZES"

# ---------------------------------------------------------------------------
# SECOND CRITICAL ISSUE: `bcftools consensus` APPLIES indels, which shifts every
# downstream coordinate on that contig. An ancestral FASTA built from a call set
# containing indels is therefore NOT in reference coordinates: Step 4 looks up
# LEPC position POS and silently gets the base from a different position, with
# the error growing along the contig. On this dataset that shifted 22 contigs by
# up to +611 bp, including most of the macrochromosomes -- so the majority of
# "polarized" calls were reading the wrong base entirely.
#
# Fix: apply SNVs only, which are 1->1 substitutions and preserve length exactly.
# Indels themselves are not usable as ancestral states anyway, and the alignment
# immediately around them is ambiguous, so those positions get masked instead.
# ---------------------------------------------------------------------------
bcftools view -V indels -Oz -o "$CONS_SNV_CALLS" "$CONS_CALLS"
tabix -p vcf "$CONS_SNV_CALLS"

# Extract indel spans once and reuse. Two deliberate choices here:
#   * `view -v indels` rather than a -i 'TYPE="indel"' filter expression. The
#     filter form is rejected by some bcftools builds, which parse `indel` as an
#     undefined INFO tag ("the tag \"indel\" is not defined in the VCF header");
#     the -v flag is portable across versions.
#   * %END rather than %POS, so the interval spans the entire REF allele. A
#     multi-base deletion must be masked across its full length, not just at its
#     first base.
# The record count is taken from this file so the call set is not streamed a
# second time just to count indels -- it holds hundreds of millions of records.
bcftools view -v indels "$CONS_CALLS" -Ou \
    | bcftools query -f '%CHROM\t%POS0\t%END\n' > "$INDEL_RAW_BED"
INDEL_COUNT=$(wc -l < "$INDEL_RAW_BED")

echo "Consensus calls: $(bcftools index --nrecords "$CONS_CALLS") total, \
$(bcftools index --nrecords "$CONS_SNV_CALLS") substitutions applied, \
$INDEL_COUNT indels excluded"

# Coverage comes from ALL called positions -- that is where chicken data exists.
bcftools query -f '%CHROM\t%POS0\t%POS\n' "$CONS_CALLS" \
    | bedtools merge -i - > "$COVERED_BED"

# Positions within INDEL_PAD_BP of a called indel are alignment-ambiguous, so
# their ancestral state is not trustworthy even though coverage exists.
awk -v PAD="$INDEL_PAD_BP" 'BEGIN { OFS = "\t" }
    { start = $2 - PAD; if (start < 0) start = 0; print $1, start, $3 + PAD }' \
    "$INDEL_RAW_BED" | bedtools sort -i - | bedtools merge -i - > "$INDEL_PAD_BED"

bedtools subtract -a "$COVERED_BED" -b "$INDEL_PAD_BED" > "$RELIABLE_BED"
# complement emits fully-uncovered contigs in their entirety, so scaffolds with
# no chicken alignment at all are masked end to end rather than skipped.
bedtools complement -i "$RELIABLE_BED" -g "$GENOME_SIZES" > "$UNCOV_BED"

bcftools consensus -f REF --mask "$UNCOV_BED" --mask-with N \
    "$CONS_SNV_CALLS" \
    > ANC

samtools faidx ANC

# Assert the ancestral FASTA is in exact reference coordinates. This is the
# check whose absence let the indel-shift bug above corrupt a whole run
# undetected: every downstream output looked entirely normal.
CTG_TOTAL=$(wc -l < REF.fai)
CTG_JOINED=$(join -j1 <(cut -f1,2 REF.fai | sort -k1,1) \
                      <(cut -f1,2 ANC.fai | sort -k1,1) | wc -l)
LEN_DIFFS=$(join -j1 <(cut -f1,2 REF.fai | sort -k1,1) \
                     <(cut -f1,2 ANC.fai | sort -k1,1) | awk '$2 != $3' | wc -l)
echo "Coordinate check: $CTG_JOINED / $CTG_TOTAL contigs present, $LEN_DIFFS length mismatch(es)"
if [ "$CTG_JOINED" -ne "$CTG_TOTAL" ] || [ "$LEN_DIFFS" -ne 0 ]; then
    echo "ERROR: ancestral FASTA is NOT in reference coordinates." >&2
    join -j1 <(cut -f1,2 REF.fai | sort -k1,1) <(cut -f1,2 ANC.fai | sort -k1,1) \
        | awk '$2 != $3 { printf "       %s: ref=%s ancestral=%s\n", $1, $2, $3 }' >&2
    echo "       Indels were applied, so every position after the first indel" >&2
    echo "       on these contigs is shifted and polarization would read the" >&2
    echo "       wrong base. Do NOT use this file." >&2
    exit 1
fi

# Guard: the ancestral FASTA MUST come out substantially masked. Chicken and
# LEPC are ~35 My diverged and chicken aligns to a minority of this assembly,
# so a low N fraction means masking silently failed. Fail loudly here rather
# than emit plausible-looking load estimates from an unmasked ancestral
# sequence -- that failure mode is invisible in every downstream output.
read MASKED_BASES TOTAL_BASES < <(awk '!/^>/ {
        TOTAL  += length($0)
        MASKED += gsub(/[Nn]/, "")
    } END { print MASKED+0, TOTAL+0 }' ANC)
MASK_PCT=$(awk -v m="$MASKED_BASES" -v t="$TOTAL_BASES" \
    'BEGIN { if (t > 0) printf "%.2f", 100*m/t; else print "0" }')

echo "Ancestral FASTA masked: $MASKED_BASES / $TOTAL_BASES bp = ${MASK_PCT}% N"
if awk -v p="$MASK_PCT" 'BEGIN { exit (p < 20) ? 0 : 1 }'; then
    echo "ERROR: only ${MASK_PCT}% of the ancestral FASTA is masked to N." >&2
    echo "       Expected roughly 50-65% (chicken covers a minority of the" >&2
    echo "       LEPC genome). Masking did not take effect -- do NOT use this" >&2
    echo "       file for polarization; derived-allele counts would be biased." >&2
    exit 1
fi

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

# --- Annotate the joint call set ---
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
set -euo pipefail

# There is no "python3" module on Gautschi (Lmod reports it as unknown), but
# python3 is on PATH by default and that is what this step uses. It runs on
# the standard library alone -- no pysam, no pandas -- because the system
# interpreter has no third-party packages installed.
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found on PATH" >&2; exit 1; }

# Pull the reference/alt alleles for every biallelic SNP in the annotated file
bcftools query -f '%CHROM\t%POS\t%RE''F\t%ALT\n' OUT/snpeff/lepc_snpeff.vcf.gz \
    > OUT/polarized/sites_ref_alt.tsv

python3 << 'PYEOF'
import os
import sys

out_dir = "OUT/polarized"
anc_path = "ANC"


class IndexedFasta:
    """Random access to a samtools-indexed FASTA, using only the standard library.

    pysam is not installed in the system python3 here and there is no python
    module to load, so the .fai index is read directly. It carries everything
    needed: each contig's length, the byte offset where its sequence starts,
    the number of bases per line, and the number of bytes per line (which
    includes the line terminator). A base's byte position is then just
    arithmetic, so this is an exact stand-in for pysam's fetch().
    """

    def __init__(self, path):
        fai_path = path + ".fai"
        if not os.path.exists(fai_path):
            sys.exit("ERROR: %s not found -- run samtools faidx on the "
                     "ancestral FASTA first (Step 1 does this)." % fai_path)
        self.handle = open(path, "rb")
        self.index = {}
        with open(fai_path) as fai:
            for line in fai:
                fields = line.split()
                if len(fields) < 5:
                    continue
                name = fields[0]
                self.index[name] = (
                    int(fields[1]),   # contig length in bases
                    int(fields[2]),   # byte offset of first base
                    int(fields[3]),   # bases per line
                    int(fields[4]),   # bytes per line (bases + newline)
                )

    def base_at(self, contig, position):
        """Return the single base at a 1-based position, or '' if unavailable."""
        record = self.index.get(contig)
        if record is None:
            return ""
        length, offset, bases_per_line, bytes_per_line = record
        if position < 1 or position > length:
            return ""
        zero_based = position - 1
        byte = offset + (zero_based // bases_per_line) * bytes_per_line \
            + (zero_based % bases_per_line)
        self.handle.seek(byte)
        return self.handle.read(1).decode("ascii", "replace").upper()


anc = IndexedFasta(anc_path)

n_in = n_kept = 0
n_skip_not_snv = n_skip_no_call = n_skip_mismatch = 0

# Streamed rather than accumulated in memory: the input has millions of sites.
with open("%s/sites_ref_alt.tsv" % out_dir) as fh, \
        open("%s/sites_polarized.tsv" % out_dir, "w") as sink:
    sink.write("chrom\tpos\tref\talt\tancestral\tderived\n")
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            continue
        chrom, pos_text, ref, alt = fields[0], fields[1], fields[2], fields[3]
        n_in += 1

        # Polarization is only meaningful for biallelic single-nucleotide
        # sites. The input is already filtered to those, so anything
        # else here is unexpected -- count it rather than silently coercing.
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            n_skip_not_snv += 1
            continue

        anc_base = anc.base_at(chrom, int(pos_text))
        if anc_base not in ("A", "C", "G", "T"):
            n_skip_no_call += 1   # no confident chicken ortholog at this site
            continue

        if anc_base == ref.upper():
            ancestral, derived = ref, alt
        elif anc_base == alt.upper():
            ancestral, derived = alt, ref
        else:
            # Ancestral state matches neither the reference nor alt allele
            # (likely a lineage-specific substitution on the chicken branch,
            # or a third allele) -- exclude from polarized load calculations.
            n_skip_mismatch += 1
            continue

        sink.write("\t".join((chrom, pos_text, ref, alt, ancestral, derived)) + "\n")
        n_kept += 1

print("Sites read:                    %d" % n_in)
print("  skipped, not a biallelic SNV: %d" % n_skip_not_snv)
print("  skipped, no ancestral call:   %d" % n_skip_no_call)
print("  skipped, ancestral matches neither allele: %d" % n_skip_mismatch)
print("Polarized: %d (%.1f%% of input sites)"
      % (n_kept, 100.0 * n_kept / max(n_in, 1)))
if n_kept == 0:
    sys.exit("ERROR: no sites could be polarized -- check that the ancestral "
             "FASTA contig names match the VCF's.")

# Guard against an UNMASKED ancestral FASTA. Chicken aligns to only a minority
# of the LEPC assembly, so a correctly masked ancestral sequence must yield a
# large "no ancestral call" fraction. If that count is near zero, the ancestral
# FASTA is almost certainly carrying the LEPC reference base at unaligned
# positions, and every such site is being polarized as ancestral-equals-
# reference by
# construction. That produces a full, healthy-looking output table whose
# derived-allele counts are systematically inflated -- so fail here instead.
frac_no_call = 100.0 * n_skip_no_call / max(n_in, 1)
print("  (no-ancestral-call fraction: %.1f%%)" % frac_no_call)
if frac_no_call < 5.0:
    sys.exit(
        "ERROR: only %.2f%% of sites had no ancestral call. Expected a large\n"
        "       fraction, because chicken covers a minority of the LEPC genome.\n"
        "       This means the ancestral FASTA is NOT masked: unaligned\n"
        "       positions still hold the LEPC reference base and are being\n"
        "       polarized as ancestral-equals-reference, biasing load metrics.\n"
        "       Re-run Step 1 (the alignment is reused; only the consensus and\n"
        "       mask are rebuilt) and confirm it reports a masked %% of ~50-65."
        % frac_no_call)
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

set -euo pipefail

# There is no "python3" module on Gautschi; python3 is on PATH by default.
# This step deliberately uses the standard library only (no pandas), since
# the system interpreter has no third-party packages installed.
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found on PATH" >&2; exit 1; }

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
import csv
import math
import os
import re
import sys
from collections import Counter

out = "OUT"


def load_site_set(path):
    """Read a two-column (contig, position) list into a set of 'contig:pos' keys."""
    keys = set()
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                keys.add(fields[0] + ":" + fields[1])
    return keys


high     = load_site_set("%s/snpeff/sites_HIGH.txt" % out)
moderate = load_site_set("%s/snpeff/sites_MODERATE.txt" % out)
low      = load_site_set("%s/snpeff/sites_LOW.txt" % out)

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
        return None
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

HGVS_P = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

aa_path = "%s/snpeff/missense_aa_changes.tsv" % out
distances = {}
n_parsed = n_unparsed = 0
if os.path.exists(aa_path) and os.path.getsize(aa_path) > 0:
    with open(aa_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            contig, position, hgvs = fields[0], fields[1], fields[2]
            matched = HGVS_P.match(hgvs)
            if not matched:
                # frameshift / extension / unknown notations, etc.
                n_unparsed += 1
                continue
            a1 = THREE_TO_ONE.get(matched.group(1))
            a2 = THREE_TO_ONE.get(matched.group(3))
            if a1 is None or a2 is None:
                n_unparsed += 1
                continue
            distance = grantham(a1, a2)
            if distance is None:
                n_unparsed += 1
                continue
            distances[contig + ":" + position] = distance
            n_parsed += 1
else:
    print("WARNING: no missense_aa_changes.tsv found -- rerun Step 2 to "
          "produce it; Grantham columns will be empty.")


def grantham_class(distance):
    """Grantham's own conventional bins."""
    if distance is None:
        return ""
    if distance <= 50:
        return "conservative"
    if distance <= 100:
        return "moderately_conservative"
    if distance <= 150:
        return "moderately_radical"
    return "radical"


source_path = "%s/polarized/sites_polarized.tsv" % out
sink_path = "%s/load/sites_classified.tsv" % out

category_counts = Counter()
class_counts = Counter()
n_missense = n_missense_scored = 0

# Streamed row by row rather than loaded into a dataframe: this file has
# millions of rows and the classification is purely per-site.
with open(source_path, newline="") as src, open(sink_path, "w", newline="") as sink:
    reader = csv.DictReader(src, delimiter="\t")
    if reader.fieldnames is None:
        sys.exit("ERROR: %s is empty -- did Step 4 finish?" % source_path)
    added = ["site", "snpeff_HIGH", "snpeff_MODERATE", "snpeff_LOW",
             "category", "grantham", "grantham_class", "category_fine"]
    writer = csv.DictWriter(sink, fieldnames=list(reader.fieldnames) + added,
                            delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
    writer.writeheader()

    for row in reader:
        site = "%s:%s" % (row["chrom"], row["pos"])
        is_high = site in high
        is_moderate = site in moderate
        is_low = site in low

        # Categories carried into per-individual load (Step 6):
        #   LOF      = SnpEff HIGH     (stop-gain, frameshift, splice-disrupting)
        #   MISSENSE = SnpEff MODERATE (amino-acid changing; severity not ranked)
        #   NEUTRAL  = SnpEff LOW      (synonymous; comparison class)
        # A site in more than one class takes the most severe, hence the order.
        category = "OTHER"
        if is_low:
            category = "NEUTRAL"
        if is_moderate:
            category = "MISSENSE"
        if is_high:
            category = "LOF"

        distance = distances.get(site)
        category_fine = category
        if category == "MISSENSE":
            n_missense += 1
            if distance is not None:
                n_missense_scored += 1
                # Split at Grantham 100, the usual conservative/radical
                # division. Missense with no parseable amino-acid change
                # keeps the unsplit MISSENSE label rather than being guessed at.
                category_fine = ("MISSENSE_CONSERVATIVE" if distance <= 100
                                 else "MISSENSE_RADICAL")
            class_counts[grantham_class(distance)] += 1

        category_counts[category] += 1

        row["site"] = site
        row["snpeff_HIGH"] = is_high
        row["snpeff_MODERATE"] = is_moderate
        row["snpeff_LOW"] = is_low
        row["category"] = category
        row["grantham"] = "" if distance is None else "%.6f" % distance
        row["grantham_class"] = grantham_class(distance)
        row["category_fine"] = category_fine
        writer.writerow(row)

print("Amino-acid changes scored: %d; unparseable/skipped: %d"
      % (n_parsed, n_unparsed))
print("Missense sites: %d; with a Grantham score: %d (%.1f%%)"
      % (n_missense, n_missense_scored,
         100.0 * n_missense_scored / max(n_missense, 1)))
print("\nSite categories:")
for name, count in category_counts.most_common():
    print("  %-10s %d" % (name, count))
print("\nGrantham class (missense only):")
for name, count in class_counts.most_common():
    print("  %-25s %d" % (name if name else "(unscored)", count))
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

module --force purge
module load biocontainers
module load bcftools
# The R module is lowercase "r" and requires a version (bare "module load R"
# fails, and `ml r` is Lmod shorthand for `module reset`, not a load).
module load r/4.4.1
set -euo pipefail

# xalt injects LD_PRELOAD into containerised commands; blank it or bcftools
# aborts on a glibc symbol mismatch inside the container.
export SINGULARITYENV_LD_PRELOAD=""
export APPTAINERENV_LD_PRELOAD=""

command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found after loading r/4.4.1" >&2; exit 1; }
Rscript -e 'q(status = as.integer(!requireNamespace("tidyverse", quietly = TRUE)))' || {
    echo "ERROR: the tidyverse R package is not available to r/4.4.1." >&2
    exit 1
}

# ---------------------------------------------------------------------------
# Genotypes are pulled with bcftools rather than vcfR, which is not installed
# for this R build. This is also the lighter option: vcfR loads an entire
# variant file into memory, whereas only the classified sites are needed here
# -- a few hundred thousand rows rather than millions.
# ---------------------------------------------------------------------------
SITES=OUT/load/sites_classified.tsv
REGIONS=OUT/load/classified_regions.tsv
GT_TABLE=OUT/load/genotypes_at_classified.tsv
SAMPLE_LIST=OUT/load/vcf_sample_order.txt

[ -s "$SITES" ] || { echo "ERROR: $SITES missing or empty -- did Step 5 finish?" >&2; exit 1; }

# Sites in the load categories only (OTHER is never tallied), sorted, as a
# two-column CHROM/POS regions file for bcftools -R.
awk -F'\t' '
    NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    $col["category"] != "OTHER" { print $col["chrom"] "\t" $col["pos"] }
' "$SITES" | sort -u -k1,1 -k2,2n > "$REGIONS"
echo "  Classified sites to genotype: $(wc -l < "$REGIONS")"

bcftools query -l VCF > "$SAMPLE_LIST"

# `bcftools query -R` uses random access and therefore requires a tabix/csi
# index. The analysis call set is not necessarily indexed, so build one rather
# than failing here. If indexing is not possible (read-only filesystem, say),
# fall back to -T, which takes the same targets file but streams the whole call
# set instead of seeking -- slower, but it needs no index and gives identical
# output. The longest LEPC scaffold is ~109 Mb, well under tbi's 512 Mb limit.
if [ -s "VCF.tbi" ] || [ -s "VCF.csi" ]; then
    QUERY_MODE=-R
else
    echo "  No index found for the analysis call set -- creating one."
    if bcftools index -t -f VCF 2>/dev/null || bcftools index -f VCF 2>/dev/null; then
        QUERY_MODE=-R
        echo "  Index created."
    else
        QUERY_MODE=-T
        echo "  Could not create an index -- falling back to a streaming query."
    fi
fi

bcftools query "$QUERY_MODE" "$REGIONS" -f '%CHROM\t%POS[\t%GT]\n' VCF > "$GT_TABLE"

[ -s "$GT_TABLE" ] || { echo "ERROR: bcftools returned no genotypes -- check that the" >&2; \
    echo "  contig names in $SITES match the variant file's." >&2; exit 1; }
echo "  Genotype rows returned: $(wc -l < "$GT_TABLE")"

# Rscript does not read a script from stdin the way `python3 -` does -- given no
# file argument it just prints its usage and exits. Write the analysis to a file
# and run that. This is also more convenient operationally: the R stage can be
# re-run on its own without repeating the expensive genotype query above.
R_SCRIPT=OUT/load/compute_load.R
cat << 'REOF' > "$R_SCRIPT"
library(tidyverse)

out <- "OUT"

old_samples <- readLines("OLD_SAMPLES")
new_samples <- readLines("NEW_SAMPLES")
group_map <- c(setNames(rep("Old", length(old_samples)), old_samples),
               setNames(rep("New", length(new_samples)), new_samples))

sites <- read_tsv(sprintf("%s/load/sites_classified.tsv", out), show_col_types = FALSE)

# Genotype matrix from the bcftools dump written above: columns are CHROM,
# POS, then one GT column per sample in the order bcftools reports them.
samples <- readLines(sprintf("%s/load/vcf_sample_order.txt", out))
gt_raw <- read_tsv(sprintf("%s/load/genotypes_at_classified.tsv", out),
                   col_names = c("chrom", "pos", samples),
                   col_types = cols(.default = col_character()),
                   na = character(), progress = FALSE)

vcf_site <- paste(gt_raw$chrom, gt_raw$pos, sep = ":")
gt <- as.matrix(gt_raw[, samples, drop = FALSE])
rownames(gt) <- vcf_site
cat(sprintf("Loaded genotypes: %d sites x %d samples\n", nrow(gt), ncol(gt)))

# The call set may carry samples that are not in the Old/New lists -- the
# excluded related individual, for example. Looping over them would fail on
# group_map[[ind]] with "subscript out of bounds", and it would fail deep in
# the loop, after the expensive genotype query. Reconcile up front instead.
analysis_samples <- intersect(colnames(gt), names(group_map))
unknown_samples  <- setdiff(colnames(gt), names(group_map))
missing_samples  <- setdiff(names(group_map), colnames(gt))
if (length(unknown_samples) > 0) {
    cat(sprintf("NOTE: %d sample(s) in the call set are not in the Old/New lists and are skipped: %s\n",
                length(unknown_samples), paste(unknown_samples, collapse = ", ")))
}
if (length(missing_samples) > 0) {
    cat(sprintf("WARNING: %d listed sample(s) are absent from the call set: %s\n",
                length(missing_samples), paste(missing_samples, collapse = ", ")))
}
if (length(analysis_samples) == 0) {
    stop("No call-set sample matches the Old/New lists -- check that the sample lists use the same names as the call set.")
}
cat(sprintf("Analysing %d samples (%d Old, %d New)\n", length(analysis_samples),
            sum(group_map[analysis_samples] == "Old"),
            sum(group_map[analysis_samples] == "New")))

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

    for (ind in analysis_samples) {
        calls <- gt_sub[, ind]
        # Require BOTH alleles to be called. A half-call such as "./1" would
        # otherwise pass and be scored as a diploid genotype with the missing
        # allele silently treated as reference, biasing derived counts down
        # while still counting the site in the denominator.
        called <- !is.na(calls) & !grepl("\\.", calls)

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

Rscript "$R_SCRIPT"

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
