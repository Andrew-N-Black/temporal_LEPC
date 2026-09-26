# Repository audit

Started 2026-09-25 as a full read-through of every script, checking
whether the autosomal (Z-scaffold exclusion) and first-degree-relative
(normal_F10 removal) changes are actually implemented where they need to
be. Updated across several follow-up passes (ROHan parsing, ROH parsing,
directory reorganization). This file keeps only what's still relevant:
resolved items are condensed to one line; open items keep full detail.

## Open items

### 1. `analysis/load/load.sh` is on the whole-genome, all-20-sample dataset

The genetic load pipeline's *methodology* is current and correct (SnpEff
impact + Grantham distance, GERP++ evaluated and deliberately dropped —
see the script's own header). But it isn't on the same dataset as every
other analysis in this repo: `VCF` has no `.auto` in its name (unlike
`remove_Z_scaffolds.sh`'s output and everything downstream of it), there's
no Z-scaffold exclusion logic anywhere in the file, and `ALL_SAMPLES` is
explicitly documented as "All 20 samples (old.samples + new.samples)" --
normal_F10 included. Not modified here: this is a 674-line script that
*generates* the actual child SLURM scripts, and getting a mid-pipeline
edit wrong in something this load-bearing would be worse than flagging it.
**Needs your input on where to point it** (the `.auto` VCF, and an
unrelated sample list built from `processing/final_cramlist_unrel.txt`).

### 2. Four R plotting scripts still read the pre-relatedness-filter metadata

`R/plot_heterozygosity.R`, `R/plot_roh.R`, `R/map.R`, and `R/QA_plots.R`
read `old_new_heterozygosity.xlsx` (20 samples, F10 included).
`R/plot_pca.R`, `R/plot_pca_plink.R`, and `R/PC_correlations.R` already
read `old_new_heterozygosity_unrel.xlsx` (19, F10 removed) instead. Still
not repointed automatically -- the `_unrel` workbook's column schema isn't
available to confirm `plot_roh.R`'s `fROH_*`, `QA_plots.R`'s
`DOC`/`properlyPaired`, or `map.R`'s `GPS` columns exist there. Each file
has a header comment flagging this.

### 3. SLURM account name: `dewoody` vs `fnrdewoody`

`analysis/PCA/heterozygosity_array.sh` uses `-A fnrdewoody` (changed
directly on GitHub). Every other `#SBATCH -A` line in the repo (17 of
them, across 13 files) still says `dewoody`. Not changed elsewhere without
confirming which is actually correct -- if `fnrdewoody` is the real
account, every other job in this repo would currently fail to submit.

### 4. Two smaller things worth a second look

- `processing/cram_filtering.sh` references the reference genome at
  `GROUSE/grouse_asm/ref/...`, while every other script uses a copy at
  `GROUSE/old_vs_new/ref/...`. May be an intentional shared reference, or
  a leftover from a different project.
- Upstream numbered pipeline scripts referenced in comments throughout
  (`05_combined_alignment_array.sh`, `06_downsample_and_finalize.sh`)
  aren't in this repo. Likely fine if only analysis-stage scripts were
  meant to be versioned here -- worth confirming that's deliberate.

## Incident: heterozygosity_array.sh was overwritten with unrelated content

Between the last two passes, `analysis/heterozygosity_array.sh` (via a
GitHub web edit, commit `70ca5b1`) was replaced with an entirely different
script -- a genome-assembly pipeline (hifiasm/yahs/RagTag, job name
`grouse_asm`, "n=23 samples across 3 species") from a different project.
Every other file in the repo was individually checked against its own
header/content for the same kind of mismatch; none found. Restored from
the last known-good commit (`e3d2031`), preserving the legitimate
`fnrdewoody` account-name edit made in between (see item 3 above). If you
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
  with the current SnpEff+Grantham methodology (see Open item 1 above for
  what's still outstanding on the *dataset* it points to).
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
beyond its header/config (674 lines, see Open item 1), and the nf-core
configuration files -- read for structure and cross-checked against
caller variable names, not verified statement-by-statement.
