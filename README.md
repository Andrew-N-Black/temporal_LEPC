# temporal_LEPC

Temporal genomics of the Lesser Prairie-Chicken (*Tympanuchus
pallidicinctus*): comparing 20 individuals sampled in two eras ("Past" /
2018-ish and "Present" / 2026-ish, 10 per era) from a New Mexico
population, to test whether genomic erosion (heterozygosity, realized
inbreeding, genetic load, allele-frequency drift) increased over the
interval.

Reference assembly: `GCF_026119805.1` (pur_lepc_1.0). This is a male-bird
assembly, so there is no W scaffold -- only Z needs excluding for
autosomal analyses. The two Z scaffolds (identified by chicken synteny +
genotype hemizygosity, see `analysis/sex_and_relatedness/`) are
`NW_026294758.1` and `NW_026294813.1`.

All SLURM jobs run on Purdue RCAC's Gautschi cluster, account `dewoody`
(one exception -- see AUDIT.md item 3).

**Before running anything here, read AUDIT.md.** It lists what's still
open, most importantly: `load.sh` isn't yet on the autosomal/unrelated
dataset, and four `R/` plotting scripts read the pre-relatedness-filter
metadata.

## Two sample sets

| | n | Used for |
|---|---|---|
| **All** (`processing/final_cramlist.txt`, `processing/popmap.txt`) | 20 | Sex calling, raw QC |
| **Autosomal, unrelated** (`processing/final_cramlist_unrel.txt`, `processing/popmap_unrelated.txt`) | 19 | Every downstream population-genetic result |

`normal_F10` is removed from the 19-sample set as a first-degree relative
of `normal_F21` (KING-robust kinship = 0.261, from
`analysis/sex_and_relatedness/relatedness_vcf.py`; more consistent with
full siblings than parent-offspring, based on R0/IBS0). F21 was kept over
F10 for lower missingness.

CRAM-list sample IDs are bare (`F10`); VCF/popmap sample IDs carry a
`normal_` prefix (`normal_F10`) from the joint-genotyping sample sheet.
Both refer to the same 20 (or 19) birds.

## Repository layout

```
processing/               Upstream: dedup/QC filtering, sample-ID
                           manifests (cramlists, popmaps -- both sample
                           sets), install_rohan.sh
nf-core/                   Purdue RCAC nf-core/sarek config (mapping ->
                           GATK4 joint calling)
analysis/
  sex_and_relatedness/     Z-scaffold identification, autosomal VCF,
                           relatedness/kinship (both the VCF-based and
                           genotype-likelihood arms)
  PCA/                     Per-sample heterozygosity, PCA/admixture/
                           inbreeding (pcangsd), the autosomal+unrelated
                           PCA (PLINK2)
  ROH/                     ANGSD + bcftools roh, and its parser
    ROHan/                 ROHan (independent genotype-likelihood
                           cross-check) and its parser
  load/                    Genetic load pipeline
  temporal_selection/      Site classification, OutFLANK/Ne, FST
R/                         Local (non-SLURM) plotting scripts, run on a
                           downloaded results spreadsheet + genotype
                           outputs
```

Scripts call sibling scripts in the same folder by bare relative filename
(e.g. `temporal_selection/outflank_drift.sh` calls
`make_siteclass.sh` and `outflank_temporal_dNe_siteclass.R` assuming
they're in the same directory, and several scripts `cd
"$SLURM_SUBMIT_DIR"` for the same reason) -- keep each folder's contents
together if you reorganize further, and submit jobs from within the
folder the script lives in.

## Pipeline order

Numbers below are logical order, not filenames.

1. **Alignment + joint calling** -- `nf-core/` configs, run through
   nf-core/sarek (`mapping` -> GATK4 HaplotypeCaller/GenotypeGVCFs).
2. **CRAM QC filtering** -- `processing/cram_filtering.sh` (dedup + MAPQ20
   + proper-pair filtering).
3. **Sex calling / Z-scaffold identification** --
   `sex_and_relatedness/find_sex_scaffolds.sh` (chicken synteny via
   minimap2) + `sex_from_vcf.py` (genotype hemizygosity confirmation).
4. **Autosomal VCF** -- `sex_and_relatedness/remove_Z_scaffolds.sh`.
5. **Relatedness / kinship** --
   `sex_and_relatedness/run_lepc_relatedness.sh`: Arm A (VCF, KING-robust
   via `relatedness_vcf.py`) and Arm B (ANGSD genotype likelihoods +
   NgsRelate). Produces the F10-first-degree-relative conclusion and
   `processing/popmap_unrelated.txt`.
6. **Per-sample heterozygosity** -- `PCA/heterozygosity_array.sh` (ANGSD
   `-doSaf` + `realSFS`, one SLURM array task per sample, all 20 --
   heterozygosity doesn't depend on relatedness the way PCA/FST/Ne do).
7. **PCA + admixture + inbreeding** --
   `PCA/pca_admixture_inbreeding.sh` (ANGSD beagle -> pcangsd;
   whole-genome, all 20 -- PCA/admixture output superseded by the next
   item, but still the only source for the pcangsd inbreeding
   coefficient) and `PCA/run_plink.sh` (PLINK2, autosomal, 19-sample
   unrelated -- the PCA used in the manuscript).
8. **ROH** -- `ROH/bcftools_roh.sh` (ANGSD -> bcftools roh) with
   `ROH/roh_parse_autosomal.sh`, and `ROH/ROHan/rohan_array.sh` (ROHan, an
   independent genotype-likelihood cross-check that doesn't need called
   genotypes) with `ROH/ROHan/rohan_parse_autosomal.sh`. Both parsers
   re-parse already-computed output -- neither ANGSD/bcftools roh nor
   ROHan itself needs to be re-run.
9. **Genetic load** -- `load/load.sh`. Current methodology (SnpEff +
   Grantham, no GERP) but not yet on the autosomal/unrelated dataset --
   see AUDIT.md item 1.
10. **Site classification** -- `temporal_selection/make_siteclass.sh`
    (SnpEff ANN -> deleterious/neutral sites, optional GERP support).
11. **Temporal selection / drift** --
    `temporal_selection/outflank_drift.sh` +
    `outflank_temporal_dNe_siteclass.R` (temporal Fk / Ne, OutFLANK
    outlier scan, site-class Fk comparison).
12. **FST** -- `temporal_selection/run_lepc_fst.sh` (wrapper) +
    `fst_temporal.py` (Hudson + Weir & Cockerham, sliding windows,
    block-jackknife CIs).
13. **Plotting** -- `R/*.R`, run locally against a downloaded results
    spreadsheet. Three of seven are on the current (autosomal, unrelated)
    data; four are not yet -- see AUDIT.md item 2.

## Environment setup

Two conda environments, built once on a login node:

```bash
module load anaconda

# lepc_py -- scikit-allel-based Python scripts (fst_temporal.py,
# sex_from_vcf.py, relatedness_vcf.py)
conda create -n lepc_py -c conda-forge python=3.11 scikit-allel numpy pandas matplotlib

# lepc_rstats -- R side (OutFLANK + friends)
conda create -n lepc_rstats -c conda-forge -c bioconda \
    r-base r-optparse r-dplyr r-ggplot2 r-vcfr bioconductor-qvalue r-devtools
conda activate lepc_rstats
Rscript -e 'install.packages("remotes", repos = "https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("whitlock/OutFLANK")'
```

ROHan is built from source once via `processing/install_rohan.sh` (needs a
from-source GSL build too, handled by the same script -- Gautschi has no
`gsl` module).

RCAC's `xalt` accounting hook injects `LD_PRELOAD` into every job,
including containerized ones, which breaks `angsd`/`pcangsd`/`plink2`
against the biocontainers' older glibc. Scripts that load `biocontainers`
blank it via `SINGULARITYENV_LD_PRELOAD=""` / `APPTAINERENV_LD_PRELOAD=""`
-- a plain `unset LD_PRELOAD` does not hold, since `xalt` re-injects it.

## Known gaps

See AUDIT.md for the full list. Short version: `load.sh` needs pointing at
the autosomal/unrelated dataset; four `R/` plotting scripts need the same
metadata update three of their siblings already got; one script's SLURM
account name needs reconciling with the rest of the repo.
