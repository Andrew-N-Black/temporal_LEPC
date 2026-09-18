#!/bin/bash
# =============================================================================
# Lesser Prairie-Chicken (LEPC) Genetic Load Pipeline
# Steps 1–6: Ancestral polarization (chicken) -> SnpEff annotation ->
#            GERP++ constraint scoring -> deleterious site classification ->
#            per-individual genetic load (Old vs New samples, n=20 total)
#
# Built off the same structure as the C. tularosa load/Rxy SLURM pipeline
# (master script -> child SLURM scripts -> sed substitution -> dependency
# chain), adapted here for:
#   - 20x hard-called genotypes (GATK4, via nf-core/sarek) instead of ANGSD GLs
#   - SnpEff instead of VEP for functional annotation
#   - GERP++ constraint scores (chicken + other galliforms) instead of SIFT
#   - Chicken (Gallus gallus) as the ancestral-allele outgroup
#   - Temporal design: OLD vs NEW sample sets instead of spatial ESUs
#
# Assumes:
#   - A joint-genotyped, hard-filtered biallelic SNP VCF already exists
#     (e.g. sarek's GATK4_GENOTYPEGVCFS output, post hard-filtering)
#   - A chicken (Gallus gallus) reference genome for ancestral polarization
#   - A multi-species galliform alignment (MAF, e.g. from Cactus) for GERP++
#   - 20 individuals total, split into two sample lists: OLD and NEW
#
# Adjust all paths in the CONFIG section before running.
# Submit each step individually, or chain with --dependency=afterok:<jobid>
# =============================================================================

# =============================================================================
# CONFIG — edit these paths before running
# =============================================================================
PROJ=/scratch/gautschi/blackan/GROUSE/old_vs_new
REF=$PROJ/ref/GCF_026119805.1_pur_lepc_1.0_genomic.fna                       # LEPC reference genome
GFF=$PROJ/ref/GCF_026119805.1_pur_lepc_1.0_genomic.gtf              # LEPC gene annotation — .gtf, .gff, or .gff3 all work (Step 2 auto-detects from the extension)
VCF=/scratch/gautschi/blackan/GROUSE/out_new_sarek/variant_calling/normalized/joint_variant_calling/joint_germline.norm.sorted.vcf.gz       # Hard-filtered, biallelic SNP VCF (all 20 samples)
OUT=$PROJ/results   # defined early since several config vars below (ANC, GERP_MSA) depend on it

# Ancestral outgroup (chicken) — needed for allele polarization
GALLUS_DIR=/scratch/gautschi/blackan/GROUSE/old_vs_new/ref   # Chicken reference already downloaded here
# ANC is defined further down (after OUT) — it must NOT live in GALLUS_DIR,
# since Step 1/3a's chicken-FASTA auto-detection scans that directory and a
# leftover/previous ANC output sitting there would get mistaken for the
# chicken genome on a re-run.

# GERP++ multi-species alignment (galliform MSA, built in Steps 0/3a below)
# LEPC is the alignment reference; chicken + the outgroups below are aligned to it
GERP_MSA=$OUT/gerp/msa/galliform_alignment.maf     # Multi-species alignment, LEPC-referenced (built in Step 3a)
GERP_TREE=$PROJ/gerp/galliform.nwk                 # Newick guide tree for the alignment taxa — EDIT BRANCH LENGTHS (see note below)

# Outgroup genomes for the GERP alignment, pulled from NCBI Datasets.
# Recommended subset from your 15 available galliform assemblies: prioritizes
# phylogenetic spread over raw count, and mostly chromosome-level assemblies.
#   - Centrocercus urophasianus (greater sage-grouse) — closest relative here,
#     same subfamily (Tetraoninae) as LEPC. Scaffold-level (1502 scaffolds)
#     but the phylogenetic proximity outweighs that.
#   - Tetrao urogallus (western capercaillie) — Tetraoninae, chromosome-level
#   - Lagopus muta (rock ptarmigan) — Tetraoninae, chromosome-level
#   - Meleagris gallopavo (turkey) — Phasianidae, mid-distance, chromosome-level
#   - Coturnix japonica (Japanese quail) — Phasianidae, mid-distance, chromosome-level
#   - Numida meleagris (helmeted guineafowl) — Numididae, more distant, chromosome-level
# Chicken (already downloaded at GALLUS_DIR) is added to this alignment too.
# Skipped: the other 9 galliform assemblies are either redundant with one of
# the above (e.g. Excalfactoria/Lagopus leucura duplicate a clade already
# covered) or heavily fragmented (Phasianus colchicus, Callipepla,
# Bambusicola, Odontophorus, Penelope, Alectura — thousands to hundreds of
# thousands of scaffolds), which tends to hurt alignment quality more than
# the extra taxon helps GERP's power.
#
# Added on request — three taxa spanning much deeper divergence than the
# galliform set above, which widens the range gerpcol's rate model sees:
#   - Tympanuchus cupido pinnatus (greater prairie-chicken) — essentially as
#     close a relative as exists (same genus as LEPC). GenBank only, no
#     RefSeq assembly currently exists for this species. Very little
#     independent constraint signal on its own at this shallow a divergence,
#     but useful as an internal check on the alignment/polarization pipeline
#     since LEPC and GRPC hybridize in the wild (Northern DPS).
#   - Anas platyrhynchos (mallard) — Anseriformes, the other half of
#     Galloanserae; deep outgroup relative to all the galliforms above.
#     Older mallard RefSeq assemblies (GCF_015476345.1, GCF_003850225.1) are
#     now RefSeq-suppressed, so this uses the current T2T reference instead.
#   - Taeniopygia guttata (zebra finch) — Passeriformes (songbirds), the
#     deepest split in this set (outside Galloanserae entirely). Anchors the
#     far end of the divergence range for the neutral-rate model; a common
#     choice as a distant reference point in avian GERP/phyloP work.
# Edit this array to add/remove taxa — format: "Name|Accession"
OUTGROUP_ACCESSIONS=(
    "Centrocercus_urophasianus|GCF_019232065.1"
    "Tetrao_urogallus|GCF_951394365.1"
    "Lagopus_muta|GCF_023343835.1"
    "Meleagris_gallopavo|GCF_905368555.1"
    "Coturnix_japonica|GCF_054131305.1"
    "Numida_meleagris|GCF_002078875.1"
    "Tympanuchus_cupido|GCA_001870855.1"
    "Anas_platyrhynchos|GCF_047663525.1"
    "Taeniopygia_guttata|GCF_003957565.2"
)
OUTGROUP_DIR=$PROJ/ref/outgroups

# Cactus (multi-species whole-genome alignment) — HPC clusters typically run
# this via container rather than a module; point at whichever is available
CACTUS_SIF=$PROJ/containers/cactus.sif

# Conda env for GERP++ / PHAST, built in Step 3-prep (see note above Step 3b)
CONDA_ENV=lepc_gerp

# SnpEff
SNPEFF_DB=LEPC_custom                           # Custom SnpEff database name to build

ANC=$OUT/ancestral/LEPC_ancestral_from_chicken.fa    # Built in Step 1 (chicken-derived ancestral FASTA in LEPC coordinates) — deliberately NOT in GALLUS_DIR, see note above

# Sample lists — one sample ID per line, matching VCF sample names exactly
OLD_SAMPLES=$PROJ/sample_lists/old.samples      # Historical/museum-era samples
NEW_SAMPLES=$PROJ/sample_lists/new.samples      # Contemporary samples
ALL_SAMPLES=$PROJ/sample_lists/all.samples      # All 20 samples (old.samples + new.samples)

# GERP threshold — sites with RS score above this are treated as "constrained"
# and eligible for the deleterious set alongside SnpEff HIGH/MODERATE impact
GERP_RS_THRESHOLD=2

THREADS=64

mkdir -p $PROJ/scripts $PROJ/logs $PROJ/gerp $OUT/{ancestral,snpeff,gerp,gerp/msa,polarized,load,logs} $OUTGROUP_DIR

# =============================================================================
# NOTE on GERP_TREE branch lengths
# Cactus needs a guide tree with (rough) branch lengths, and gerpcol uses the
# same tree for its neutral-rate model. The topology below reflects standard
# avian phylogeny, but the branch lengths are PLACEHOLDERS — replace them
# with divergence-time-scaled values (e.g. from TimeTree.org) before running
# Step 3a/3b in earnest. Rough topology used here:
#   (((((LEPC,Tympanuchus_cupido),Centrocercus_urophasianus),
#      (Tetrao_urogallus,Lagopus_muta)),
#     ((Meleagris_gallopavo,Gallus_gallus),Coturnix_japonica),Numida_meleagris),
#    Anas_platyrhynchos,Taeniopygia_guttata);
# Anas_platyrhynchos (Anseriformes) and Taeniopygia_guttata (Passeriformes)
# both sit outside Galliformes, so they attach at the root relative to
# everything else here — Taeniopygia is the more distant of the two.
# =============================================================================


# =============================================================================
# STEP 0: Download outgroup genome assemblies from NCBI Datasets
#         Pulls each accession in OUTGROUP_ACCESSIONS via the `datasets` CLI
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step0_download_outgroups.sh
#!/bin/bash
#SBATCH --job-name=lepc_dl_outgroups
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
module load ncbi-datasets/16.10.3
set -euo pipefail

mkdir -p OUTGROUP_DIR
cd OUTGROUP_DIR

OUTGROUPS_STR
EOF
# Inject the OUTGROUP_ACCESSIONS array as literal bash into the child script
# (arrays don't survive the sed substitution pass cleanly, so this is built
# separately and spliced in before substitution below)
{
    echo "TAXA_LIST=("
    for entry in "${OUTGROUP_ACCESSIONS[@]}"; do
        echo "    \"$entry\""
    done
    echo ")"
    cat << 'EOF'
for entry in "${TAXA_LIST[@]}"; do
    name="${entry%%|*}"
    accession="${entry##*|}"

    if [ -f "${name}.fa" ]; then
        echo "  ${name} already present, skipping download"
        continue
    fi

    echo "Downloading ${name} (${accession})..."
    datasets download genome accession "$accession" --include genome --filename "${accession}.zip"
    unzip -o "${accession}.zip" -d "${accession}_unzipped" > /dev/null

    FASTA=$(find "${accession}_unzipped" -iname "*.fna" | head -n 1)
    if [ -z "$FASTA" ]; then
        echo "  ERROR: no FASTA found for ${accession}" >&2
        continue
    fi
    cp "$FASTA" "${name}.fa"
    samtools faidx "${name}.fa"
    rm -rf "${accession}.zip" "${accession}_unzipped"

    echo "  ${name}: $(grep -c '^>' ${name}.fa) sequences downloaded"
done

echo "Step 0 complete: $(date)"
EOF
} > $PROJ/scripts/.step0_body.sh
sed -i "/^OUTGROUPS_STR$/r $PROJ/scripts/.step0_body.sh" $PROJ/scripts/step0_download_outgroups.sh
sed -i "/^OUTGROUPS_STR$/d" $PROJ/scripts/step0_download_outgroups.sh
rm -f $PROJ/scripts/.step0_body.sh


# =============================================================================
# STEP 1: Build the ancestral FASTA from chicken alignment
#         Aligns chicken genome to the LEPC reference and calls a consensus
#         base at each LEPC position where an unambiguous chicken ortholog
#         exists. This is the polarization reference for Step 4.
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

echo "Step 2 complete: $(date)"
EOF


# =============================================================================
# STEP 3A: Build the galliform multi-species alignment with Cactus
#          Aligns chicken + all downloaded outgroups to the LEPC reference,
#          producing a HAL file, then exports a LEPC-referenced MAF for GERP++
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step3a_build_msa.sh
#!/bin/bash
#SBATCH --job-name=lepc_cactus
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
module load cactus/2.6.5
set -euo pipefail

# Auto-detect the chicken FASTA (same logic as Step 1 — excludes REF since
# GALLUS_DIR also holds the LEPC reference itself in this setup)
CHICKEN_FASTA=$(find GALLUS_DIR -maxdepth 1 -iname "*.fa" -o -iname "*.fasta" -o -iname "*.fna" | grep -v "\.gz$" | grep -vxF "REF" | grep -vxF "ANC" | head -n 1)
if [ -z "$CHICKEN_FASTA" ]; then
    echo "ERROR: no chicken FASTA found in GALLUS_DIR" >&2
    exit 1
fi

# Build the Cactus seqFile: one line per taxon, "Name path/to/genome.fa",
# with LEPC included as a normal leaf (Cactus produces a reference-free HAL;
# hal2maf below is what fixes LEPC as the MAF reference row).
SEQFILE=OUT/gerp/msa/seqfile.txt
echo "LEPC REF" > "$SEQFILE"
echo "Gallus_gallus $CHICKEN_FASTA" >> "$SEQFILE"
for f in OUTGROUP_DIR/*.fa; do
    taxon=$(basename "$f" .fa)
    echo "${taxon} ${f}" >> "$SEQFILE"
done

echo "Cactus seqFile contents:"
cat "$SEQFILE"

# GERP_TREE (edited with real branch lengths per the note in CONFIG) goes at
# the top of the seqFile as the guide tree, per Cactus convention
if [ ! -s GERP_TREE ]; then
    echo "ERROR: GERP_TREE (GERP_TREE) is missing or empty." >&2
    echo "  Create it with a real Newick guide tree (topology + branch lengths," >&2
    echo "  e.g. from TimeTree.org) for all taxa in the seqFile above before running this step." >&2
    exit 1
fi
cat GERP_TREE "$SEQFILE" > OUT/gerp/msa/seqfile_with_tree.txt

# Cactus alignment can take the better part of a day -- skip it entirely if
# galliform.hal is already sitting there from a prior successful run (e.g.
# this run is only here to retry the hal2maf export). Delete the .hal
# yourself first if you actually want to rebuild the alignment from scratch.
if [ -s OUT/gerp/msa/galliform.hal ]; then
    echo "galliform.hal already exists -- skipping Cactus alignment, going straight to the hal2maf export"
else
    # Toil (Cactus's job manager) refuses to reuse an existing jobstore from a
    # prior attempt without --restart, which this pipeline doesn't set up for --
    # clear it so each run starts clean. Toil creates this directory itself
    # during initialization and errors if it already exists, so only remove it
    # here -- do NOT mkdir it back, or Cactus will hit the same error again.
    rm -rf OUT/gerp/msa/jobstore

    cactus OUT/gerp/msa/jobstore OUT/gerp/msa/seqfile_with_tree.txt OUT/gerp/msa/galliform.hal \
        --maxCores THREADS
fi

# Export a MAF with LEPC fixed as the reference genome/row. hal2maf isn't
# listed among the cactus module's documented commands on this cluster, but
# the HAL toolkit is normally bundled inside the same Cactus install -- try
# it directly first. If it's not exposed via the module's PATH, fall back to
# reusing the exact same container image the cactus module itself wraps,
# discovered from the module's own wrapper script rather than a separately
# downloaded/built .sif (Apptainer/Singularity is system-wide on this
# cluster, so no module load is needed for it).
if command -v hal2maf >/dev/null 2>&1; then
    hal2maf OUT/gerp/msa/galliform.hal GERP_MSA \
        --refGenome LEPC --noAncestors --noDupes
else
    echo "hal2maf not found via the cactus module -- looking up the module's container image" >&2

    # Every command below is guarded with `|| true` / `2>/dev/null || true` so
    # nothing here can trip `set -e` and kill the script silently -- a bare,
    # unguarded pipeline did exactly that on a prior attempt (no error, no
    # success message, job just stopped). Progress is echoed after each
    # attempt so if this still comes up empty, the log shows exactly which
    # methods were tried and found nothing, rather than dying with no trace.

    echo "  [1/4] checking 'module show cactus/2.6.5' output..." >&2
    MODULE_SHOW_TEXT=$(module show cactus/2.6.5 2>&1 || true)
    IMAGE_PATH=$(echo "$MODULE_SHOW_TEXT" | grep -oE '(/[^[:space:]"'"'"']+\.(sif|simg|sqsh))' 2>/dev/null | head -n 1 || true)
    echo "  module show result: ${IMAGE_PATH:-<none found>}" >&2

    if [ -z "$IMAGE_PATH" ]; then
        echo "  [2/4] checking the cactus wrapper script content..." >&2
        CACTUS_PATH=$(command -v cactus 2>/dev/null || true)
        if [ -n "$CACTUS_PATH" ] && [ -f "$CACTUS_PATH" ]; then
            IMAGE_PATH=$(grep -oE '(/[^[:space:]"'"'"']+\.(sif|simg|sqsh))' "$CACTUS_PATH" 2>/dev/null | head -n 1 || true)
        fi
        echo "  wrapper script result: ${IMAGE_PATH:-<none found>} (wrapper path: ${CACTUS_PATH:-<not found>})" >&2
    fi

    if [ -z "$IMAGE_PATH" ]; then
        echo "  [3/4] checking 'type cactus' in case it's a shell function..." >&2
        TYPE_CACTUS_TEXT=$(type cactus 2>/dev/null || true)
        IMAGE_PATH=$(echo "$TYPE_CACTUS_TEXT" | grep -oE '(/[^[:space:]"'"'"']+\.(sif|simg|sqsh))' 2>/dev/null | head -n 1 || true)
        echo "  type cactus result: ${IMAGE_PATH:-<none found>}" >&2
    fi

    if [ -z "$IMAGE_PATH" ]; then
        echo "  [4/4] searching /apps and /opt for a cactus image (up to 60s)..." >&2
        IMAGE_PATH=$(timeout 60 find /apps /opt -maxdepth 6 -iname "*cactus*.sif" 2>/dev/null | head -n 1 || true)
        echo "  filesystem search result: ${IMAGE_PATH:-<none found>}" >&2
    fi

    if [ -z "$IMAGE_PATH" ] || [ ! -f "$IMAGE_PATH" ]; then
        echo "" >&2
        echo "ERROR: could not automatically determine the cactus container image path." >&2
        echo "  All four automated discovery methods failed (see attempts logged above)." >&2
        echo "  Diagnostic dump follows -- please share this output:" >&2
        echo "  --- module show cactus/2.6.5 ---" >&2
        echo "$MODULE_SHOW_TEXT" >&2
        echo "  --- type cactus ---" >&2
        echo "$TYPE_CACTUS_TEXT" >&2
        echo "  --- env | grep -i cactus ---" >&2
        (env | grep -i cactus || true) >&2
        echo "  --- env | grep -i sif ---" >&2
        (env | grep -i sif || true) >&2
        echo "  --- which cactus / command -v cactus ---" >&2
        echo "${CACTUS_PATH:-<not found>}" >&2
        exit 1
    fi
    echo "Found container image: $IMAGE_PATH"
    singularity exec --bind PROJ:PROJ "$IMAGE_PATH" \
        hal2maf OUT/gerp/msa/galliform.hal GERP_MSA \
        --refGenome LEPC --noAncestors --noDupes
fi

echo "MAF written to: GERP_MSA"
echo "Step 3a complete: $(date)"
EOF


# =============================================================================
# STEP 3-PREP: Build a conda env for GERP++ and PHAST
#              Neither has a module on this cluster, but both are on bioconda,
#              so build a small dedicated env once instead. bedtools/samtools
#              are included here too so Step 3b only needs one environment.
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step3prep_setup_gerp_env.sh
#!/bin/bash
#SBATCH --job-name=lepc_gerp_env
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

module load anaconda
set -euo pipefail

# GERP++ (gerpcol) and PHAST (msa_view, used in Step 3b to convert MAF blocks
# to the FASTA alignment gerpcol expects) aren't available as modules here,
# but both are on bioconda. Build one small env for them (plus bedtools/
# samtools, so Step 3b only needs to activate one environment).
if conda env list | grep -q "^CONDA_ENV "; then
    echo "Conda env 'CONDA_ENV' already exists, skipping creation"
else
    conda create -y -n CONDA_ENV -c bioconda -c conda-forge gerp phast bedtools samtools
fi

echo "Step 3-prep complete: $(date)"
EOF


# =============================================================================
# STEP 3B: GERP++ — per-site constraint (RS) scores from the galliform MSA
#         Requires gerpcol (GERP++) and a MAF alignment with LEPC as reference
# =============================================================================
cat << 'EOF' > $PROJ/scripts/step3b_gerp_scores.sh
#!/bin/bash
#SBATCH --job-name=lepc_gerp
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

module load anaconda
set -euo pipefail
# A gerpcol segfault ("core dumped") writes a full process memory image to
# disk, and on this cluster that appears to have exhausted disk space/quota
# -- every command after the first crash (all subsequent msa_view calls)
# failed identically afterward, not just the one bad chromosome. Disable
# core dumps entirely so one crash can never cascade into a full-run
# failure this way again.
ulimit -c 0
source activate CONDA_ENV

# GERP++'s gerpcol takes a per-chromosome FASTA-format alignment, but Cactus/
# hal2maf produces a genome-wide MAF with one block per aligned segment. Split
# the MAF by LEPC-reference chromosome, then use PHAST's msa_view to stitch
# each chromosome's blocks into a single gapped FASTA alignment (unaligned/
# missing taxon positions become gaps, which gerpcol treats as missing data).

mkdir -p OUT/gerp/per_chrom
python3 << 'PYEOF'
import os

maf_path = "GERP_MSA"
out_dir = "OUT/gerp/per_chrom"
os.makedirs(out_dir, exist_ok=True)

handles = {}

def get_handle(chrom):
    if chrom not in handles:
        handles[chrom] = open(f"{out_dir}/{chrom}.maf", "w")
    return handles[chrom]

block = []
block_chrom = None
with open(maf_path) as fh:
    for line in fh:
        if line.startswith("a"):
            if block and block_chrom:
                get_handle(block_chrom).writelines(block)
                get_handle(block_chrom).write("\n")
            block = [line]
            block_chrom = None
        elif line.startswith("s"):
            block.append(line)
            fields = line.split()
            # MAF sequence name is "<taxon>.<chrom>" -- Cactus/hal2maf names
            # the LEPC row "LEPC.<chrom>"
            if fields[1].startswith("LEPC."):
                block_chrom = fields[1].split(".", 1)[1]
        elif line.strip() == "":
            if block and block_chrom:
                get_handle(block_chrom).writelines(block)
                get_handle(block_chrom).write("\n")
            block = []
            block_chrom = None
        else:
            block.append(line)

    if block and block_chrom:
        get_handle(block_chrom).writelines(block)
        get_handle(block_chrom).write("\n")

chrom_files = list(handles.keys())
for h in handles.values():
    h.close()

# msa_view requires blocks sorted by reference-sequence position within each
# file (it warns and silently drops/mishandles out-of-order blocks otherwise,
# which corrupts its internal position tracking against --refseq and produces
# spurious "character does not match" errors downstream). Cactus/hal2maf emit
# blocks in HAL graph traversal order, not reference-coordinate order, so
# re-read each per-chromosome file and sort its blocks by the LEPC row's
# start position before handing off to msa_view. MAF start coordinates are
# always given in the frame of the block's own strand; convert '-' strand
# starts to the equivalent '+' strand position so sorting reflects true
# genomic order regardless of which strand a given block aligned on.
#
# Sorting alone isn't sufficient: hal2maf --noDupes still leaves some
# overlapping reference blocks in practice (confirmed -- the "invalid
# idx_offset" error recurs even against a freshly rebuilt --noDupes MAF).
# msa_view needs genuinely non-overlapping (tiled) reference coverage, so
# after sorting, greedily drop any block whose reference range overlaps the
# previously kept block's range. This can discard real alignment (typically
# at duplicated/repetitive regions, the same places prone to producing
# pathological data elsewhere in this pipeline), which is an acceptable
# trade for GERP scoring -- better to lose ambiguous columns there than to
# have gerpcol choke on them entirely.
for chrom in chrom_files:
    path = f"{out_dir}/{chrom}.maf"
    with open(path) as fh:
        content = fh.read()

    raw_blocks = [b for b in content.split("\n\n") if b.strip()]

    def ref_range(block_text):
        for bline in block_text.splitlines():
            if bline.startswith("s"):
                bfields = bline.split()
                if bfields[1].startswith("LEPC."):
                    start = int(bfields[2])
                    size = int(bfields[3])
                    strand = bfields[4]
                    src_size = int(bfields[5])
                    if strand == "-":
                        eff_start = src_size - (start + size)
                    else:
                        eff_start = start
                    return (eff_start, eff_start + size)
        return (0, 0)  # block had no LEPC row (shouldn't happen -- these were the split key)

    ranged_blocks = [(ref_range(b), b) for b in raw_blocks]
    ranged_blocks.sort(key=lambda rb: rb[0])

    kept_blocks = []
    prev_end = -1
    n_dropped = 0
    for (rstart, rend), block_text in ranged_blocks:
        if rstart < prev_end:
            n_dropped += 1
            continue
        kept_blocks.append(block_text)
        prev_end = rend

    if n_dropped:
        print(f"  {chrom}: dropped {n_dropped} overlapping block(s) of {len(raw_blocks)} to produce non-overlapping reference coverage")

    with open(path, "w") as fh:
        fh.write("\n\n".join(kept_blocks) + "\n")

print(f"Split MAF into {len(chrom_files)} per-chromosome files in {out_dir}")
PYEOF

CHROM_LIST=$(cut -f1 REF.fai | sort -u)
FAILED_LOG=OUT/gerp/per_chrom/FAILED_CHROMS.txt
> "$FAILED_LOG"

for CHR in $CHROM_LIST; do
    CHR_MAF=OUT/gerp/per_chrom/${CHR}.maf
    MFA=OUT/gerp/per_chrom/${CHR}.mfa
    CHR_FASTA=OUT/gerp/per_chrom/${CHR}.ref.fa

    if [ ! -s "$CHR_MAF" ]; then
        echo "  Skipping ${CHR}: no alignment blocks found in GERP_MSA"
        continue
    fi

    # msa_view --refseq was being handed the entire 206-contig genome FASTA
    # on every call, and appears to resolve the reference sequence by some
    # positional/sequential logic rather than a clean per-call name lookup
    # -- confirmed by the exact failure pattern: each chromosome's error
    # reported the immediately PRECEDING scaffold's length, cascading
    # through the whole genome in file order. Extracting just this one
    # contig into its own single-sequence FASTA removes any possibility of
    # it resolving to the wrong entry, since there's now only one sequence
    # in the file for it to find.
    if ! samtools faidx REF "$CHR" > "$CHR_FASTA" 2>OUT/gerp/per_chrom/${CHR}.faidx.err; then
        echo "  FAILED (samtools faidx): ${CHR} -- see ${CHR}.faidx.err" >&2
        echo -e "${CHR}\tsamtools_faidx" >> "$FAILED_LOG"
        continue
    fi

    # A single chromosome's alignment failing (bad character, degenerate
    # data, gerpcol segfaulting on an edge case, etc.) shouldn't take down
    # the whole 200+ chromosome run -- catch failures at each step, log
    # which chromosome and which step, and move on. FAILED_CHROMS.txt at
    # the end gives a complete, one-run picture of exactly what needs
    # investigating, rather than discovering problem chromosomes one at a
    # time across repeated full reruns.
    if ! msa_view "$CHR_MAF" --in-format MAF --out-format FASTA --refseq "$CHR_FASTA" \
        > "$MFA" 2>OUT/gerp/per_chrom/${CHR}.msa_view.err; then
        echo "  FAILED (msa_view): ${CHR} -- see ${CHR}.msa_view.err" >&2
        echo -e "${CHR}\tmsa_view" >> "$FAILED_LOG"
        continue
    fi

    if [ ! -s "$MFA" ]; then
        echo "  Skipping ${CHR}: msa_view produced an empty alignment"
        continue
    fi

    # GERP++'s bundled parser rejects any character outside its expected
    # alphabet (throws BIO::E_InvalidCharacterEx and aborts). Uppercasing
    # alone (soft-masked lowercase repeats) did NOT fix this -- confirmed
    # by rerunning and hitting the identical crash -- so something else in
    # the alphabet is the culprit, most likely IUPAC ambiguity codes
    # (R/Y/S/W/K/M/B/D/H/V, used for ambiguous base calls) which GERP++
    # doesn't accept even though they're valid FASTA. Rather than guess at
    # which exact character(s), sanitize comprehensively: uppercase, then
    # replace anything that isn't A/C/G/T/N/- with N (GERP++'s own
    # missing-data character) so no unexpected symbol can reach gerpcol.
    awk '/^>/{print; next} {line=toupper($0); gsub(/[^ACGTN-]/, "N", line); print line}' "$MFA" > "${MFA}.clean" \
        && mv "${MFA}.clean" "$MFA"

    if ! gerpcol -t GERP_TREE -f "$MFA" -e Gallus_gallus \
        -a > OUT/gerp/per_chrom/${CHR}.rates 2>OUT/gerp/per_chrom/${CHR}.gerpcol.err; then
        echo "  FAILED (gerpcol): ${CHR} -- see ${CHR}.gerpcol.err" >&2
        echo -e "${CHR}\tgerpcol" >> "$FAILED_LOG"
        continue
    fi

    # gerpcol output columns: neutral_rate, RS_score, one row per alignment
    # column of the reference (LEPC) sequence, in order. This must be a
    # proper BED interval (chrom, 0-based start, end) for tabix -p bed --
    # RS_score is NOT a valid "end" coordinate (it's a float, often
    # negative), which is what caused tabix's "end < begin" indexing
    # failure: score was landing in the end-coordinate column.
    awk -v chr="$CHR" 'BEGIN{OFS="\t"; pos=0}
        {pos++; print chr, pos-1, pos, $2}' \
        OUT/gerp/per_chrom/${CHR}.rates \
        >> OUT/gerp/lepc_gerp_scores.bed

    echo "  ${CHR}: GERP scores written"
done

N_FAILED=$(wc -l < "$FAILED_LOG")
if [ "$N_FAILED" -gt 0 ]; then
    echo ""
    echo "  WARNING: ${N_FAILED} chromosome(s) failed and were skipped -- see $FAILED_LOG"
    echo "  and the per-chromosome .msa_view.err / .gerpcol.err files in OUT/gerp/per_chrom/"
    echo "  for details. GERP scores for all other chromosomes were still computed."
fi

# Sort and bgzip/tabix the genome-wide GERP track for fast lookup
sort -k1,1 -k2,2n OUT/gerp/lepc_gerp_scores.bed | \
    bgzip > OUT/gerp/lepc_gerp_scores.bed.gz
tabix -p bed OUT/gerp/lepc_gerp_scores.bed.gz

# Flag sites above the constraint threshold -- score is now column 4
# (chrom, start, end, RS_score), not column 3 as it was before the BED fix.
# Use the "end" column (3) for the extracted site list, not "start" (2):
# end equals the original 1-based position, matching the 1-based convention
# every other site-list file here uses (SnpEff's sites_HIGH.txt etc., via
# bcftools query %POS) -- using the 0-based start instead would silently
# shift every GERP-based site by one position relative to those.
zcat OUT/gerp/lepc_gerp_scores.bed.gz | \
    awk -v thr=GERP_RS_THRESHOLD '$4 > thr {print $1"\t"$3}' \
    > OUT/gerp/sites_GERP_constrained.txt

echo "  Constrained sites (RS > GERP_RS_THRESHOLD): $(wc -l < OUT/gerp/sites_GERP_constrained.txt)"
echo "Step 3 complete: $(date)"
EOF


# =============================================================================
# STEP 4: Polarize alleles against the chicken-derived ancestral FASTA
#         For every biallelic SNP: whichever allele matches ANC is ancestral,
#         the other is derived. Sites where ANC is N/missing are dropped.
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
#         Deleterious = SnpEff HIGH, or SnpEff MODERATE, or GERP RS above
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

python3 << 'PYEOF'
import pandas as pd

out = "OUT"

pol = pd.read_csv(f"{out}/polarized/sites_polarized.tsv", sep="\t")
pol["site"] = pol["chrom"].astype(str) + ":" + pol["pos"].astype(str)

def load_site_set(path):
    df = pd.read_csv(path, sep="\t", names=["chrom", "pos"])
    return set(df["chrom"].astype(str) + ":" + df["pos"].astype(str))

high      = load_site_set(f"{out}/snpeff/sites_HIGH.txt")
moderate  = load_site_set(f"{out}/snpeff/sites_MODERATE.txt")
low       = load_site_set(f"{out}/snpeff/sites_LOW.txt")
gerp_hit  = load_site_set(f"{out}/gerp/sites_GERP_constrained.txt")

pol["snpeff_HIGH"]     = pol["site"].isin(high)
pol["snpeff_MODERATE"] = pol["site"].isin(moderate)
pol["snpeff_LOW"]      = pol["site"].isin(low)
pol["gerp_constrained"] = pol["site"].isin(gerp_hit)

# Category definitions used for per-individual load in Step 6:
#   LOF        = SnpEff HIGH                          (most severe)
#   DELETERIOUS = SnpEff MODERATE AND GERP-constrained (missense + conserved)
#   NEUTRAL    = SnpEff LOW, used as a standardization / comparison class
pol["category"] = "OTHER"
pol.loc[pol["snpeff_LOW"], "category"] = "NEUTRAL"
pol.loc[pol["snpeff_MODERATE"] & pol["gerp_constrained"], "category"] = "DELETERIOUS"
pol.loc[pol["snpeff_HIGH"], "category"] = "LOF"

pol.to_csv(f"{out}/load/sites_classified.tsv", sep="\t", index=False)

print(pol["category"].value_counts())
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
    select(site, category, derived_is_alt)

results <- list()

for (cat in c("LOF", "DELETERIOUS", "NEUTRAL")) {
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
for SCRIPT in step0 step1 step2 step3a step3prep step3b step4 step5 step6; do
    FILE=$PROJ/scripts/${SCRIPT}*.sh
    sed -i \
        -e "s|REF.fai|${REF}.fai|g" \
        -e "s|REF|$REF|g" \
        -e "s|GFF|$GFF|g" \
        -e "s|VCF|$VCF|g" \
        -e "s|GALLUS_DIR|$GALLUS_DIR|g" \
        -e "s|OUTGROUP_DIR|$OUTGROUP_DIR|g" \
        -e "s|CACTUS_SIF|$CACTUS_SIF|g" \
        -e "s|CONDA_ENV|$CONDA_ENV|g" \
        -e "s|ANC|$ANC|g" \
        -e "s|GERP_MSA|$GERP_MSA|g" \
        -e "s|GERP_TREE|$GERP_TREE|g" \
        -e "s|GERP_RS_THRESHOLD|$GERP_RS_THRESHOLD|g" \
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
echo "   JOB0=\$(sbatch --parsable scripts/step0_download_outgroups.sh)"
echo "   JOB1=\$(sbatch --parsable scripts/step1_build_ancestral.sh)          # can run in parallel with JOB0"
echo "   JOB2=\$(sbatch --parsable scripts/step2_snpeff_annotate.sh)          # can run in parallel with JOB0/JOB1"
echo "   JOBPREP=\$(sbatch --parsable scripts/step3prep_setup_gerp_env.sh)    # one-time conda env build, can run anytime"
echo "   JOB3A=\$(sbatch --parsable --dependency=afterok:\$JOB0 scripts/step3a_build_msa.sh)"
echo "   JOB3B=\$(sbatch --parsable --dependency=afterok:\$JOB3A,afterok:\$JOBPREP scripts/step3b_gerp_scores.sh)"
echo "   JOB4=\$(sbatch --parsable --dependency=afterok:\$JOB1,afterok:\$JOB2 scripts/step4_polarize.sh)"
echo "   JOB5=\$(sbatch --parsable --dependency=afterok:\$JOB3B,afterok:\$JOB4 scripts/step5_deleterious_sites.sh)"
echo "   JOB6=\$(sbatch --parsable --dependency=afterok:\$JOB5 scripts/step6_load_per_individual.sh)"
echo ""
echo " Output files:"
echo "   $OUT/load/sites_classified.tsv               — per-site category (LOF/DELETERIOUS/NEUTRAL)"
echo "   $OUT/load/genetic_load_per_individual.tsv     — per-individual total/realized/masked load"
echo "   $OUT/load/genetic_load_summary.tsv            — Old vs New group summary"
echo "   $OUT/load/old_vs_new_tests.tsv                — Wilcoxon tests, Old vs New per category/metric"
echo "============================================================"
echo ""
echo " Before running: edit GERP_TREE (galliform.nwk) with real branch lengths"
echo " (e.g. from TimeTree.org) — Cactus and gerpcol both use it directly."
echo "============================================================"
