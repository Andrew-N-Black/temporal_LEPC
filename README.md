# temporal_LEPC

Temporal genomics of the Lesser Prairie-Chicken (*Tympanuchus
pallidicinctus*): comparing 20 individuals sampled in two eras ("Past" /
2018-ish and "Present" / 2026-ish, 10 per era) from a New Mexico
population, to test whether genomic erosion (heterozygosity, realized
inbreeding, genetic load, allele-frequency drift) increased over the
interval.

Reference assembly: `GCF_026119805.1` (pur_lepc_1.0). This is a male-bird
assembly, so there is no W scaffold -- only Z needs excluding for autosomal
analyses. The two Z scaffolds (identified by chicken synteny + genotype
hemizygosity, see "Sex and relatedness" below) are `NW_026294758.1` and
`NW_026294813.1`.

All SLURM jobs run on Purdue RCAC's Gautschi cluster, account `dewoody`.

**Before running anything here, read AUDIT.md.** A full script review
turned up several bugs (now fixed) and -- more importantly -- several gaps
that are flagged but deliberately *not* silently patched, including the
genetic load script (`analysis/load.sh`) being an outdated draft. Numbers
this repo produces as currently checked out will not all match the ones
reported in the manuscript until those gaps are closed.

## Two sample sets

| | n | Used for |
|---|---|---|
| **All** (`processing/final_cramlist.txt`) | 20 | Sex calling, raw QC |
| **Autosomal, unrelated** (`processing/final_cramlist_unrel.txt`, `popmap_unrelated.txt`) | 19 | Every downstream population-genetic result |

`normal_F10` is removed from the 19-sample set as a first-degree relative
of `normal_F21` (KING-robust kinship = 0.261; more consistent with full
siblings than parent-offspring, based on R0/IBS0). F21 was kept over F10
for lower missingness. See AUDIT.md, section 3 -- the kinship script that
produced this conclusion isn't in this repo yet.

CRAM-list sample IDs are bare (`F10`); VCF/popmap sample IDs carry a
`normal_` prefix (`normal_F10`) from the joint-genotyping sample sheet.
Both refer to the same 20 (or 19) birds.

## Repository layout

```
processing/     Upstream: dedup/QC filtering, the two sample-ID manifests
nf-core/        Purdue RCAC nf-core/sarek config (mapping -> GATK4 calling)
analysis/       Everything past joint genotyping: sex, relatedness-consuming
                filters, PCA/ROH/heterozygosity, genetic load, temporal
                selection/drift, FST
R/              Local (non-SLURM) plotting scripts, run on a downloaded
                results spreadsheet + genotype outputs
```

Scripts inside `analysis/` call each other by bare relative filename (e.g.
`outflank_drift.sh` calls `outflank_temporal_dNe_siteclass.R` and
`make_siteclass.sh` assuming they're in the same directory) -- keep them
together if you reorganize further.

## Pipeline order

Numbers below are logical order, not filenames (nothing was renamed, to
avoid breaking the relative-path calls above).

1. **Alignment + joint calling** -- `nf-core/` configs, run through
   nf-core/sarek (`mapping` -> GATK4 HaplotypeCaller/GenotypeGVCFs).
2. **CRAM QC filtering** -- `processing/cram_filtering.sh`
   (dedup + MAPQ20 + proper-pair filtering).
3. **Sex calling / Z-scaffold identification** --
   `analysis/find_sex_scaffolds.sh` (chicken synteny via minimap2) +
   `analysis/sex_from_vcf.py` (genotype hemizygosity confirmation).
   Produces `SEX_CHROMS.txt`.
4. **Autosomal VCF** -- `analysis/remove_Z_scaffolds.sh`, using the two Z
   scaffold IDs from step 3.
5. **Relatedness / kinship** -- *(script not yet in this repo, see
   AUDIT.md section 3)*. Produces the F10-first-degree-relative
   conclusion and `popmap_unrelated.txt`.
6. **Per-sample heterozygosity** -- `analysis/heterozygosity_array.sh`
   (ANGSD `-doSaf` + `realSFS`, one SLURM array task per sample).
7. **PCA + admixture + inbreeding** -- `analysis/pca_roh.sh` (ANGSD
   beagle -> pcangsd; whole-genome, all 20) and
   `analysis/run_plink.sh` (PLINK2, autosomal, 19-sample unrelated --
   this is the one used in the manuscript).
8. **ROH** -- `analysis/roh_analyses.sh` (ANGSD -> bcftools roh ->
   vendored `rohparser.py`) and `analysis/rohan_array.sh` (formerly
   `run_ROHan.sh` -- ROHan, a genotype-likelihood cross-check that doesn't
   need called genotypes), parsed by `analysis/rohan_parse.sh`
   (whole-genome) and `analysis/rohan_parse_autosomal.sh` (autosomal-only
   -- the one to use; re-parses the same ROHan output, no need to re-run
   it).
9. **Genetic load** -- `analysis/load.sh`. **Stale -- see AUDIT.md
   section 1 before using.**
10. **Site classification** -- `analysis/make_siteclass.sh` (SnpEff
    ANN -> deleterious/neutral sites, optional GERP support).
11. **Temporal selection / drift** -- `analysis/outflank_drift.sh` +
    `R/outflank_temporal_dNe_siteclass.R` (temporal Fk / Ne, OutFLANK
    outlier scan, site-class Fk comparison).
12. **FST** -- `analysis/run_lepc_fst.sh` (wrapper) +
    `analysis/fst_temporal.py` (Hudson + Weir & Cockerham, sliding
    windows, block-jackknife CIs).
13. **Plotting** -- `R/*.R`, run locally against a downloaded results
    spreadsheet. Three of seven are on the current (autosomal, unrelated)
    data; four are not yet -- see AUDIT.md section 2.

## Environment setup

Two conda environments, built once on a login node:

```bash
module load anaconda

# lepc_py -- scikit-allel-based Python scripts (fst_temporal.py,
# sex_from_vcf.py, the not-yet-committed relatedness script)
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
against the biocontainers' older glibc. Every script that loads
`biocontainers` blanks it via `SINGULARITYENV_LD_PRELOAD=""` /
`APPTAINERENV_LD_PRELOAD=""` -- a plain `unset LD_PRELOAD` does not hold,
since `xalt` re-injects it.

## Known gaps

See AUDIT.md for the full list with reasoning. Short version:

- `analysis/load.sh` needs replacing with the final GERP-free,
  Grantham-scored, autosomal/unrelated pipeline.
- Four `R/` plotting scripts need repointing at the unrelated metadata
  workbook (once its column schema is confirmed).
- The relatedness/kinship script and `popmap.txt` /
  `popmap_unrelated.txt` aren't committed yet.
- Z-scaffold exclusion is now applied everywhere it needs to be
  (`heterozygosity_array.sh`, `roh_analyses.sh`/`rohparser.py`, and
  `rohan_parse_autosomal.sh`). The `rohparser.py` fix also turned up an
  unrelated bug that was silently deflating every per-sample F(ROH) by
  ~506x; see AUDIT.md section 5.
