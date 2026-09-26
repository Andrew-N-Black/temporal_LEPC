# Repository audit — 2026-09-25

A full read-through of every script in this repo, focused on whether the
autosomal (Z-scaffold exclusion) and first-degree-relative (normal_F10
removal) changes are actually implemented where they need to be. This file
records what was found and fixed, and — just as important — what was found
and **not** fixed, with the reasoning either way.

Fixes with clear, verifiable justification were made directly (see commit
history). Anything that would have required guessing at scientific
methodology or reconstructing missing code from memory/description alone
was left alone and is flagged below instead — silently "fixing" a load or
relatedness calculation incorrectly would be worse than leaving a visible
gap.

## Fixed in this pass

| File | Problem | Fix |
|---|---|---|
| `analysis/fst_temporal.py` + `analysis/run_lepc_fst.sh` | Contents were swapped: `fst_temporal.py` held the bash SLURM wrapper, `run_lepc_fst.sh` held the 315-line Python FST implementation. The wrapper's own `python fst_temporal.py` call would have tried to run itself (bash) as Python and crashed immediately. | Swapped back. Verified `fst_temporal.py` parses as Python (`py_compile`) and `run_lepc_fst.sh` passes `bash -n`. |
| `analysis/roh_analyses.sh` | `PROJECT_DIR` pointed at `GROUSE/nexus` — a separate, unrelated 433-sample project — while every other variable in the file (`ROHAN_BIN`, `GSL_PREFIX`, the job name `old.new_roh_resume`) assumed `GROUSE/old_vs_new`. | `PROJECT_DIR` corrected to `GROUSE/old_vs_new`. |
| `analysis/heterozygosity_array.sh` | Sample-ID suffix stripping assumed `_filt.cram` / `_ds.cram`. Neither matches this project's actual CRAM names (`<ID>.md.dedup_q20.cram`), so `SAMPLE` was left as the full filename — e.g. output files named `F10.md.dedup_q20.cram_heterozygosity.txt` instead of `F10_heterozygosity.txt`, breaking any downstream join by sample ID. | Suffix corrected to `.md.dedup_q20.cram`. |
| `analysis/remove_Z_scaffolds.sh` | No shebang (so `./remove_Z_scaffolds.sh` would fail), no `set -euo pipefail`, no input-file check, no header comment — the only script in the repo without any of these. | Added all four; same `bcftools` logic and same Z scaffold IDs as before, unchanged. |
| `processing/final_cramlist_unrel.txt` | Missing entirely. `R/plot_pca.R`, `R/plot_pca_plink.R`, and `R/PC_correlations.R` already reference it as a required input. | Added: `final_cramlist.txt` with only the `F10.md.dedup_q20.cram` line removed (20 → 19 lines), verified by `diff`. |

## Flagged, not fixed — needs your input

### 1. `analysis/load.sh` is the first (pre-revision) draft — highest priority

This is a single, never-updated commit. Two rounds of revisions made in
chat since then are **not** reflected in the committed file:

- **GERP++ was dropped.** The 11-taxon galliform alignment totals only
  ~0.72 substitutions/site — too shallow for informative per-site
  constraint scores — and `gerpelem` element calls showed no missense
  enrichment over synonymous variants (22.1% vs 28.1%). The final
  classification uses SnpEff impact alone (LOF=HIGH, MISSENSE=MODERATE,
  NEUTRAL=LOW/synonymous) plus Grantham (1974) physicochemical distance as
  a severity gradient within missense variants (conservative ≤100, radical
  >100). None of that is in the committed `load.sh` — it still runs the
  full Cactus/GERP++ alignment-and-scoring steps (Steps 0, 3-PREP, 3A, 3B).
- **Not on the autosomal, unrelated dataset.** The VCF path hardcoded in
  `load.sh` (`out_new_sarek/.../joint_germline.norm.sorted.vcf.gz`) is a
  different, older VCF than the one every other script in this repo now
  uses, and it was never run through the Z-scaffold exclusion or the
  normal_F10 removal.

The final pipeline (5 steps: 1 ancestral sequence from chicken; 2 SnpEff +
`missense_aa_changes.tsv` extraction; 4 polarization; 5 classification with
Grantham scoring; 6 per-individual load with Wilcoxon tests) exists only in
prior chat sessions as generated SLURM scripts — it was never pushed to
this repo. I did not attempt to reconstruct ~1000 lines of load-
classification code from a description; getting a detail wrong here would
silently corrupt a number that ends up in the manuscript. **Recommended
next step:** regenerate those five scripts (I can do this directly if you
want — I have the design from prior sessions, just not byte-exact code) and
commit them in place of `load.sh`, pointed at the autosomal + unrelated
VCF.

### 2. Four R plotting scripts read the pre-relatedness-filter metadata

`R/plot_heterozygosity.R`, `R/plot_roh.R`, `R/map.R`, and `R/QA_plots.R`
all read `old_new_heterozygosity.xlsx` (all 20 samples, including F10).
`R/plot_pca.R`, `R/plot_pca_plink.R`, and `R/PC_correlations.R` were
already updated to read `old_new_heterozygosity_unrel.xlsx` (19 samples,
F10 removed) instead. The four were not, and were left alone rather than
repointed blindly — I don't have the `_unrel` workbook's column schema, and
if it was built as a trimmed PCA-only sheet rather than a full superset, an
automatic path swap could silently break `plot_roh.R`'s `fROH_*` columns,
`QA_plots.R`'s `DOC`/`properlyPaired` columns, or `map.R`'s `GPS` column
rather than failing loudly. Each file now has a header comment flagging
this. **Recommended next step:** confirm what's in
`old_new_heterozygosity_unrel.xlsx`, then either point these four at it too
or add the missing columns to it.

### 3. Relatedness/kinship script is missing entirely

`analysis/find_sex_scaffolds.sh`'s own header comment says its output is
"ready to paste into `SEX_CHROMS` in `run_lepc_relatedness.sh`" — but no
`run_lepc_relatedness.sh` exists in this repo, and neither does the
KING-robust kinship script it would call (built in chat as
`relatedness_vcf.py`: pairwise KING/R0/R1 from a VCF, no allele-frequency
dependence, robust to population structure — this is what actually
identified the F10–F21 pair as first-degree). Every downstream script that
consumes its conclusion (`run_plink.sh`'s `DROP_SAMPLES="normal_F10"`,
`popmap_unrelated.txt`) is present and correct, but the analysis that
*produced* that conclusion isn't committed anywhere. Recommend adding it
so the F10 exclusion is reproducible from this repo alone, not just
asserted.

### 4. Two sample-manifest files are referenced but not present

`popmap.txt` (used by `find_sex_scaffolds.sh`) and `popmap_unrelated.txt`
(used by `outflank_drift.sh` and `run_lepc_fst.sh`'s wrapper) — `<sample>
<Past|Present>` mappings for all 20 and 19 samples respectively — aren't in
the repo. I did not fabricate these: getting even one sample's era wrong
would be a real, silent scientific error, and I don't have a verified
source for all 20 assignments. Only `processing/final_cramlist.txt` (CRAM
paths, no era info) is present.

### 5. Z-scaffold exclusion is missing from three per-sample summary steps

**Update:** all three below are now fixed (as of the most recent pass).
Unlike `run_plink.sh`, `remove_Z_scaffolds.sh`, and `outflank_drift.sh`
(all correctly autosome-restricted from the start), three scripts computed
their statistic **genome-wide**, including the Z scaffolds:

- ~~`analysis/heterozygosity_array.sh` — ANGSD `-doSaf` had no region
  restriction.~~ **Fixed.** Now builds a race-safe, shared ANGSD `-rf`
  region file excluding the two Z scaffolds and passes it to `-doSaf`.
- ~~`analysis/roh_analyses.sh` Step 6 / `rohparser.py` — the fROH
  denominator was computed from the whole reference `.fai`.~~ **Fixed**,
  and while fetching and actually reading the vendored `rohparser.py`
  source to fix this, found it also had a second, unrelated, more severe
  bug: it divided every sample's own ROH region counts/lengths by a
  hardcoded `num_sam = 506` (left over from a different, 433-sample
  cohort's copy of this script) before computing F(ROH) — an individual's
  F(ROH) has no cohort-size term in its definition, so this was silently
  deflating every reported per-sample F(ROH) by ~506x. `rohparser.py` is
  now vendored locally (`analysis/rohparser.py`, no more runtime download
  + sed patch) with both issues fixed, plus `--exclude-scaffolds` support;
  verified against mock data with a hand-computed expected answer before
  committing. **If any F(ROH) numbers from the old (ANGSD/bcftools-roh,
  not ROHan) path have already been reported anywhere, they should be
  re-derived — the ~506x deflation bug predates this session.**
- ~~`rohan_array.sh`/`rohan_parse.sh` (formerly `run_ROHan.sh`) — the
  length-class parsing summed ROH segment lengths from every chromosome in
  `.mid.hmmrohl.gz`, Z included.~~ **Fixed via a new script,
  `analysis/rohan_parse_autosomal.sh`**, rather than editing
  `rohan_parse.sh` in place, so the whole-genome and autosomal-only
  summaries both stay available for comparison. Re-parses the same
  already-computed `hmmrohl.gz` files (ROHan itself doesn't need
  re-running), excludes the two Z scaffolds from both the ROH sum and the
  genome-length denominator, uses the same two length bins as
  `rohan_parse.sh` (100kb-1Mb, >1Mb), and does not depend on
  `rohparser.py` (that script parses bcftools-roh's differently-shaped
  "RG" line format, not ROHan's `hmmrohl` format — kept as a separate,
  self-contained awk parser instead). Verified against mock `hmmrohl`
  data + a mock `.fai` with a hand-computed expected answer (0.1 / 0.8 /
  0.9) before committing — matched exactly, including a zero-ROH sample
  and the Z-segment exclusion counter.

This matters because females are hemizygous for Z, so Z-linked sites look
artificially homozygous — inflating fROH and deflating heterozygosity for
females specifically, in a sex-biased way, if left in. This is exactly
what a prior manuscript-editing session flagged as still needing
"autosomal confirmation" for heterozygosity and fROH.

### 6. Two smaller things worth a second look

- `processing/cram_filtering.sh` references the reference genome at
  `GROUSE/grouse_asm/ref/...`, while every other script in this repo uses
  a copy at `GROUSE/old_vs_new/ref/...`. This may be intentional (a
  shared, centrally-stored reference), or it may be a leftover from a
  different project — worth confirming.
- `analysis/outflank_drift.sh`'s `FILT_VCF` is derived as
  `${RAW_VCF%.vcf.gz}.biallelic.snps.AUTO.vcf.gz` where `RAW_VCF` is
  already named `output.subset.biallelic.snps.auto.vcf.gz` — this produces
  a working but redundant filename
  (`...auto.biallelic.snps.AUTO.vcf.gz`). Purely cosmetic; left alone to
  avoid touching a working, already-correct script.

## Confirmed accurate — implements both the autosomal and first-degree-relative changes correctly

- `analysis/run_plink.sh` — excludes both known Z scaffolds via
  `--not-chr`, drops `normal_F10` before the QC filters, and includes its
  own sanity checks (F10 gone, 19 samples, zero Z variants remaining).
- `analysis/find_sex_scaffolds.sh` + `analysis/sex_from_vcf.py` — the
  synteny + hemizygosity pipeline that identified the two Z scaffold IDs
  used throughout the rest of the repo.
- `analysis/remove_Z_scaffolds.sh` (now hardened, see above).
- `analysis/outflank_drift.sh` — reads the autosomal VCF and
  `popmap_unrelated.txt`.
- `analysis/run_lepc_fst.sh` (the bash wrapper, post-swap-fix) — reads the
  autosomal VCF and `popmap_unrelated.txt`.
- `analysis/fst_temporal.py` (the Python implementation, post-swap-fix) —
  no Z-specific logic of its own, correctly inherits autosomal filtering
  from whatever VCF it's pointed at.

## Not reviewed line-by-line

`R/outflank_temporal_dNe_siteclass.R` (562 lines) and the nf-core
configuration files (`nf-core/nextflow.config`, `nf.params`,
`samplesheet.csv`) were read for structure and cross-checked against the
variable names their callers pass in, but not verified statement-by-
statement. Nothing autosome- or relatedness-related stood out in either.
