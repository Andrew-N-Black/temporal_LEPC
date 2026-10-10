# temporal_LEPC

Temporal genomics of the Lesser Prairie-Chicken (*Tympanuchus
pallidicinctus*) in New Mexico: 20 birds resequenced to ~20x, 10 sampled in
2019 ("Past") and 10 in 2026 ("Present"), testing whether heterozygosity,
inbreeding (f_ROH), genetic load and allele frequencies changed over
~2.6 generations.

This repository is the code for **Objective 1** of the USFWS report
*Grouse genomics, October 2026* ("Temporal comparison of New Mexico
Lesser Prairie-Chicken"). Each script header names the report section,
table or figure it produces. Sequence data: NCBI BioProject PRJNA1513026
(embargoed).

- Reference: LEPC `pur_lepc_1.0` (GCF_026119805.1). Male assembly, so no W;
  the two Z-linked scaffolds (`NW_026294758.1`, `NW_026294813.1`; 74.2 Mb)
  are excluded from every analysis. Autosomal length used as the f_ROH
  denominator: ~920 Mb.
- Cluster: Purdue RCAC Gautschi, SLURM account `fnrdewoody`.

## Sample sets

| Set | n | Files | Used for (report) |
|---|---|---|---|
| All | 20 (10 + 10) | `processing/final_cramlist.txt`, `processing/popmap.txt` | Sex, relatedness, sequencing QC, heterozygosity, f_ROH, ROHan (Table 1, S1; Figs 2-4, S1-S3) |
| Unrelated | 19 (9 + 10) | `processing/final_cramlist_unrel.txt`, `processing/popmap_unrelated.txt` | Genetic load, F_ST, N_e, OutFLANK/drift null, PCA, fastStructure (Table 2, S2; Figs 5-6, S4) |

`normal_F10` is dropped from the unrelated set as a first-degree relative of
`normal_F21` (KING-robust kinship 0.261, R1 0.699); F21 was kept for its
lower missingness. CRAM lists use bare IDs (`F10`); VCF/popmap IDs carry a
`normal_` prefix from the sarek sample sheet.

## Pipeline (report section -> scripts)

| Step | Report (Objective 1) | Scripts | Key settings |
|---|---|---|---|
| 1. Mapping + joint calling | Sampling and sequencing | `nf-core/` (nf-core/sarek) | fastp, bwa-mem2, MarkDuplicates, GATK4 HaplotypeCaller -> GenotypeGVCFs, hard filters |
| 2. CRAM filtering | (input to ANGSD/ROHan) | `processing/cram_filtering.sh` | dedup, MAPQ >= 20, proper pairs |
| 3. Z scaffolds + sex | Sex chromosomes, sex and relatedness | `analysis/sex_and_relatedness/find_sex_scaffolds.sh`, `sex_from_vcf.py`, `remove_Z_scaffolds.sh` | minimap2 `-x asm20` to chicken GRCg7b; Z:autosome heterozygosity ratio |
| 4. Relatedness | same | `run_lepc_relatedness.sh`, `relatedness_vcf.py` | KING-robust, R0, R1 (+ NgsRelate check) |
| 5. Heterozygosity | Heterozygosity and inbreeding | `analysis/PCA/heterozygosity_array.sh` | ANGSD 0.940 `-dosaf 1 -GL 1 -minMapQ 30 -doCounts 1 -setMinDepth 3`, autosomes (`-rf`), folded SFS (`realSFS -fold 1`); all 20 |
| 6. ROH / f_ROH | same | `analysis/ROH/bcftools_roh.sh`, `roh_parse_autosomal.sh` | ANGSD `-GL 1 -doBcf 1 -doPost 1 -minQ 30 -SNP_pval 1e-6` -> `bcftools roh --AF-file`; ROH with quality >= 30; Z excluded; classes 100 kb-1 Mb and > 1 Mb (ROH < 100 kb not counted) |
| 7. ROHan cross-check | same (Figs S2-S3; `R/angsd_vs_rohan_heterozygosity.R`, `R/bcftools_vs_rohan.R`) | `analysis/ROH/ROHan/rohan_array.sh`, `rohan_parse_autosomal.sh` | Z excluded |
| 8. Genetic load | Genetic load (Table S2, Fig S4) | `analysis/load/load.sh`, `R/load.R`, `R/plot_load.R` | chicken-polarized ancestral alleles (sites within 10 bp of an indel masked), SnpEff HIGH / MODERATE / LOW; total, realized, masked load; deleterious:neutral ratios; 19 birds |
| 9. F_ST | Temporal differentiation | `analysis/temporal_selection/run_lepc_fst.sh` + `fst_temporal.py` | Hudson and W&C; >= 5 called birds per era; 5-Mb block jackknife (363 blocks) |
| 10. N_e, outliers | Temporal differentiation, N_e and selection | `analysis/temporal_selection/outflank_drift.sh` (+ `make_siteclass.sh`, `outflank_temporal_dNe_siteclass.R`) | t = 2.63 generations; F_k on loci with mean frequency 0.05-0.95; OutFLANK and Wright-Fisher drift null, q < 0.05 |
| 11. NeEstimator | same | `analysis/temporal_selection/neestimator/` | 1 SNP per 10 kb (~90,500 SNPs); Pollak, Jorde-Ryman, Nei-Tajima |
| 12. PCA | Population structure (Fig 5) | `analysis/PCA/run_plink.sh`, `R/plot_pca_plink.R`, `R/PC_correlations.R` | PLINK2; biallelic SNPs, MAF >= 0.05, missingness <= 0.1, HWE 1e-6; LD-pruned (`--indep-pairwise 50 10 0.1`); 19 birds |
| 13. fastStructure | same (Fig 6) | `analysis/PCA/run_faststructure.sh`, `R/plot_faststructure.R` | same pruned SNPs; K = 2-5, logistic prior, chooseK.py; Fisher's exact test of cluster vs era |
| 14. Figures | Results | `R/*.R` (run locally) | read `temporal_metadata.xlsx` (20 birds) or `temporal_metadata_unrelated.xlsx` (19) |

Not used in the report (kept for reference): `analysis/PCA/pca_admixture_inbreeding.sh`
(PCAngsd, whole genome, all 20) and `analysis/fst/angsd_fst_temporal.sh`
(genotype-likelihood F_ST cross-check).

Scripts call sibling scripts by bare filename (e.g. `outflank_drift.sh` ->
`make_siteclass.sh`, `outflank_temporal_dNe_siteclass.R`); submit each job
from the folder it lives in.

## Layout

```
nf-core/                    sarek config, params, sample sheet
processing/                 CRAM filtering, sample lists/popmaps, ROHan install
analysis/
  sex_and_relatedness/      Z scaffolds, sex, autosomal VCF, kinship
  PCA/                      heterozygosity, PLINK2 PCA, fastStructure
  ROH/                      ANGSD + bcftools roh and parser; ROHan/
  load/                     genetic load pipeline
  temporal_selection/       F_ST, F_k/N_e, OutFLANK, NeEstimator
  fst/                      ANGSD F_ST cross-check (not reported)
R/                          plotting / statistics for report figures
```

## Environment

```bash
module load anaconda
conda create -n lepc_py -c conda-forge python=3.11 scikit-allel numpy pandas matplotlib
conda create -n lepc_rstats -c conda-forge -c bioconda \
    r-base r-optparse r-dplyr r-ggplot2 r-vcfr bioconductor-qvalue r-devtools
conda activate lepc_rstats
Rscript -e 'remotes::install_github("whitlock/OutFLANK")'
```

ROHan is built once with `processing/install_rohan.sh` (Gautschi has no
`gsl` module). ANGSD, PLINK2, PCAngsd and fastStructure come from RCAC
`biocontainers`. RCAC's `xalt` injects `LD_PRELOAD` into containers; scripts
blank it with `SINGULARITYENV_LD_PRELOAD="" APPTAINERENV_LD_PRELOAD=""`.

See `AUDIT.md` for the change log and remaining open items.
