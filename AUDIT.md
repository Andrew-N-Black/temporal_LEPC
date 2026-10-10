# Repository audit

Started 2026-09-25 as a full read-through of every script, checking
whether the autosomal (Z-scaffold exclusion) and first-degree-relative
(normal_F10 removal) changes are actually implemented where they need to
be. Updated across several follow-up passes (ROHan parsing, ROH parsing,
directory reorganization). This file keeps only what's still relevant:
resolved items are condensed to one line; open items keep full detail.

## Concordance pass against the USFWS report (2026-10-10)

Checked every script against Objective 1 of the October 2026 USFWS report.
Changes made:

- **SLURM account** standardized to `fnrdewoody` in every script (the three
  most recently edited scripts already used it).
- **`outflank_drift.sh`**: `GENERATIONS` 5 -> 2.63 (7 yr / 2.66 yr per
  generation), the value the report uses for every N_e.
- **fastStructure**: no script was versioned. Added
  `analysis/PCA/run_faststructure.sh` (reconstructed from the reported
  methods: PLINK2-pruned SNPs, K = 2-5, logistic prior, chooseK.py) --
  check against what was actually run. `R/plot_finestructure.R` renamed
  `R/plot_faststructure.R` and given the era x cluster Fisher's test.
- **R scripts**: `map.R` and `plot_roh.R` had an unterminated string in the
  `read_xlsx()` path (would not parse); fixed. Stale "CHECK BEFORE USE"
  headers in `QA_plots.R`, `map.R`, `plot_heterozygosity.R`, `plot_roh.R`
  replaced: they now read `temporal_metadata.xlsx`, and using all 20 birds
  for per-individual statistics matches the report.
- **`load.sh`**: already reads `analysis.autosomes.unrel.vcf.gz`; sample-list
  comments corrected to 9 + 10 = 19 birds (old item 1 closed, see below).
- **`analysis/busco/`** removed: assembly BUSCO belongs to Objective 3
  (`prairie_grouse`, steps 07-09).
- Headers of `pca_admixture_inbreeding.sh` (PCAngsd) and
  `fst/angsd_fst_temporal.sh` now state they are not used in the report.
- README rewritten: 2019/2026 eras, step-by-step map to report sections
  and settings.

Report text corrected to match the code (report v55): heterozygosity has no
`-minQ` filter (ANGSD default) and uses `-setMinDepth 3`; f_ROH settings
(ANGSD `-doBcf`, quality >= 30, ROH < 100 kb not counted); PCA is PLINK2 on
LD-pruned SNPs (not PCAngsd).

## Open items

1. **`processing/cram_filtering.sh`** reads the reference from
   `GROUSE/grouse_asm/ref/...`; every other script uses
   `GROUSE/old_vs_new/ref/...`. Same assembly either way; confirm.
2. Upstream scripts named in comments (`05_combined_alignment_array.sh`,
   `06_downsample_and_finalize.sh`) are not in this repo (analysis-stage
   code only).

## Incident: heterozygosity_array.sh was overwritten with unrelated content

Between the last two passes, `analysis/heterozygosity_array.sh` (via a
GitHub web edit, commit `70ca5b1`) was replaced with an entirely different
script -- a genome-assembly pipeline (hifiasm/yahs/RagTag, job name
`grouse_asm`, "n=23 samples across 3 species") from a different project.
Every other file in the repo was individually checked against its own
header/content for the same kind of mismatch; none found. Restored from
the last known-good commit (`e3d2031`), preserving the legitimate
`fnrdewoody` account-name edit made in between (account now `fnrdewoody` everywhere). If you
still need that assembly script, it wasn't otherwise saved by this
restore -- worth checking its source project's own repo.

## Resolved (condensed)

- **File-content swaps/corruption:** `fst_temporal.py` <-> `run_lepc_fst.sh`
  (fixed, verified with `py_compile`/`bash -n`); `heterozygosity_array.sh`
  overwritten with unrelated content (see incident above, restored).
- **Wrong project path:** former `roh_analyses.sh`'s `PROJECT_DIR` pointed
  at `GROUSE/nexus` (fixed; script since trimmed and renamed, see below).
- **CRAM-suffix mismatches** (three separate instances of the same class
  of bug -- code assumed the wrong suffix for this project's actual CRAM
  names, `<ID>.md.dedup_q20.cram`): `heterozygosity_array.sh` (fixed);
  the old `roh_analyses.sh` Step 6 comment/logic (moot -- that code was
  removed; `roh_parse_autosomal.sh` strips the correct suffix, with a
  defensive fallback).
- **`rohparser.py`'s cohort-size bug:** divided every sample's own F(ROH)
  by a hardcoded `num_sam=506` left over from a different project's copy
  -- deflating every value by ~506x. `rohparser.py` has since been
  deleted entirely (see "Removed" below); if any F(ROH) from that old
  path was ever reported, it should be re-derived.
- **Z-scaffold exclusion**, previously missing from three per-sample
  summary steps, is now applied in all three: `heterozygosity_array.sh`
  (ANGSD `-rf`), `ROH/roh_parse_autosomal.sh`, and
  `ROH/ROHan/rohan_parse_autosomal.sh`.
- **Missing files** (relatedness script, popmap manifests) flagged in an
  earlier pass have since been added directly: `relatedness_vcf.py`,
  `run_lepc_relatedness.sh`, `popmap.txt`, `popmap_unrelated.txt`.
- **`load.sh` was the pre-revision (GERP-based) draft** -- since replaced
  with the current SnpEff-impact methodology (now on the autosomal, unrelated
  dataset).
- **`load.sh` dataset**: now on the autosomal, unrelated VCF; cluster `old.samples` confirmed n = 9 (F10 absent).
- **`outflank_drift.sh`'s `FILT_VCF`** produced a redundant filename
  (`...auto.biallelic.snps.AUTO.vcf.gz`); simplified.

## Removed

- `analysis/rohparser.py` -- superseded by `ROH/roh_parse_autosomal.sh`'s
  self-contained parser (no external script, no hardcoded paths, no
  cohort-size bug).
- `analysis/rohan_parse.sh` -- whole-genome (Z-included) ROHan parser,
  superseded by `ROH/ROHan/rohan_parse_autosomal.sh`.
- Former `roh_analyses.sh` Steps 6 (rohparser.py-based parsing) and 7
  (serial ROHan) -- both superseded (Step 6 by `roh_parse_autosomal.sh`
  reading Step 5's output directly; Step 7 by the standalone, array-
  parallel `ROHan/rohan_array.sh`). The script itself was trimmed to just
  the remaining steps and renamed `ROH/bcftools_roh.sh`.

## Renamed / reorganized

See README.md's "Repository layout" for the current directory structure.
Renames: `roh_analyses.sh` -> `ROH/bcftools_roh.sh` (trimmed, see above);
`pca_roh.sh` -> `PCA/pca_admixture_inbreeding.sh` (header/USAGE no longer
matched its own filename; clarified that its PCA output is superseded by
`PCA/run_plink.sh` while its pcangsd inbreeding coefficient is still the
only source for that statistic). `popmap.txt`/`popmap_unrelated.txt`
moved from `analysis/` to `processing/`; the four scripts that read them
(`find_sex_scaffolds.sh`, `run_lepc_relatedness.sh`, `outflank_drift.sh`,
`run_lepc_fst.sh`) didn't use a `$PROJECT_DIR`-style variable for this,
so their `POPMAP` was changed from a bare, `$SLURM_SUBMIT_DIR`-relative
filename to the same kind of absolute path they already use for their
other inputs -- this doesn't assume the deployed cluster layout mirrors
this repo's folders, only that the file lives on scratch where it always
did.

## Confirmed accurate

`run_plink.sh`, `find_sex_scaffolds.sh`/`sex_from_vcf.py`,
`remove_Z_scaffolds.sh`, `outflank_drift.sh`, `run_lepc_fst.sh`/
`fst_temporal.py`, `run_lepc_relatedness.sh`/`relatedness_vcf.py` all
correctly implement the autosomal and/or first-degree-relative filtering
they need.

## Not reviewed line-by-line

`R/outflank_temporal_dNe_siteclass.R` (562 lines), `analysis/load/load.sh`
beyond its header/config, and the nf-core
configuration files -- read for structure and cross-checked against
caller variable names, not verified statement-by-statement.
