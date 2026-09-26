#!/usr/bin/env Rscript
# ============================================================================
# Past vs Present LEPC temporal genomics
#   (1) Temporal Fk / drift-based Ne estimate  (Nei & Tajima 1981; Waples 1989)
#   (2) OutFLANK genome-wide FST outlier scan  (Whitlock & Lotterhos 2015)
#   (3) Site-class comparison of Fk: deleterious vs neutral sites
#
# Usage:
#   Rscript outflank_temporal_dNe_siteclass.R \
#       --vcf old_vs_new.biallelic.snps.vcf.gz \
#       --popmap popmap.txt \
#       --siteclass siteclass.txt \
#       --generations 5 \
#       --outprefix results/old_vs_new
#
# Input VCF: prefiltered to biallelic SNPs (see the SLURM wrapper).
#
# Input file formats (no header, whitespace/tab-delimited):
#   popmap.txt    : <sample_id>  <Past|Present>   (old/new also accepted)
#   siteclass.txt : <CHROM>  <POS>  <deleterious|neutral>  (make_siteclass.sh)
#
# Checkpoints: the genotype matrix (step 1) and OutFLANK results (step 2) are
#   saved as <outprefix>_ckpt*.rds and reused on reruns when the VCF, popmap
#   and relevant options are unchanged. Use --no_reuse to force recompute.
#
# Ne estimate (step 3):
#   - Sampling correction is applied per locus using the number of genotyped
#     individuals at that locus (missing data increase sampling variance).
#   - Loci with mean ALT frequency < --ne_minfreq (or > 1 - ne_minfreq) are
#     excluded from Ne, because rare alleles inflate Fk (Waples 1989).
#   - A block jackknife (--jk_blocksize bp blocks) gives an approximate 95% CI,
#     since linked SNPs are not independent.
#
# Drift null (step 3): simulated separately for each combination of starting
#   frequency bin and per-locus sample sizes, so missing data don't inflate the
#   outlier count.
#
# Site-class test (step 4): Fk depends strongly on allele frequency, and
#   deleterious alleles are skewed rare, so an unstratified comparison is
#   confounded. The primary test compares within-frequency-stratum percentile
#   ranks of Fk and permutes class labels within strata. The unstratified
#   Wilcoxon test is reported for reference only.
# ============================================================================

suppressMessages({
  library(optparse)
  library(vcfR)
  library(OutFLANK)
  library(qvalue)
  library(dplyr)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# 0. CLI args and input checks
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--vcf", type = "character", help = "Input VCF (can be .vcf.gz)"),
  make_option("--popmap", type = "character", help = "Sample -> Past/Present map"),
  make_option("--siteclass", type = "character", default = NULL,
              help = "Optional CHROM POS class(deleterious/neutral) table"),
  make_option("--generations", type = "numeric", default = 1,
              help = "Generations elapsed between past and present samples [default %default]"),
  make_option("--nsim", type = "numeric", default = 2000,
              help = "WF drift simulations per null group [default %default]"),
  make_option("--freqbins", type = "numeric", default = 20,
              help = "Starting-frequency bins for the drift null and strata [default %default]"),
  make_option("--qthresh", type = "numeric", default = 0.05,
              help = "q-value threshold for outliers [default %default]"),
  make_option("--mincalled", type = "numeric", default = 5,
              help = "Min called individuals required in EACH period [default %default]"),
  make_option("--maxfitloci", type = "numeric", default = 100000,
              help = "Random loci used to fit the OutFLANK null [default %default]"),
  make_option("--chunksize", type = "numeric", default = 250000,
              help = "Loci per chunk when computing FST [default %default]"),
  make_option("--run_outflank", type = "logical", default = TRUE,
              help = "Run the OutFLANK scan (TRUE/FALSE) [default %default]"),
  make_option("--ne_minfreq", type = "numeric", default = 0.05,
              help = "Exclude loci with mean ALT freq outside [x, 1-x] from Ne [default %default]"),
  make_option("--jk_blocksize", type = "numeric", default = 5e6,
              help = "Block size (bp) for the jackknife CI on Ne [default %default]"),
  make_option("--nperm", type = "numeric", default = 2000,
              help = "Permutations for the stratified site-class test [default %default]"),
  make_option("--no_reuse", action = "store_true", default = FALSE,
              help = "Ignore checkpoint files and recompute everything"),
  make_option("--outprefix", type = "character", default = "outflank_temporal",
              help = "Output file prefix [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$vcf) || is.null(opt$popmap)) {
  stop("Must supply --vcf and --popmap. Use --help for usage.")
}
# Fail fast on missing inputs, before hours of computation
for (f in c(opt$vcf, opt$popmap, opt$siteclass)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}

out_dir <- dirname(opt$outprefix)
if (out_dir != ".") dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("== Past vs Present LEPC temporal analysis ==\n")
cat("VCF:         ", opt$vcf, "\n")
cat("Popmap:      ", opt$popmap, "\n")
cat("Siteclass:   ", ifelse(is.null(opt$siteclass), "(none provided)", opt$siteclass), "\n")
cat("Generations: ", opt$generations, "\n")
cat("Min called:  ", opt$mincalled, "per period\n")
cat("OutFLANK:    ", opt$run_outflank, "\n\n")

file_sig <- function(f) {
  fi <- file.info(f)
  paste(normalizePath(f), fi$size, as.numeric(fi$mtime))
}
load_ckpt <- function(path, key) {
  if (opt$no_reuse || !file.exists(path)) return(NULL)
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(obj) || !identical(obj$key, key)) return(NULL)
  cat("    Reusing checkpoint:", path, "\n")
  obj
}
quiet <- function(expr) {            # swallow OutFLANK's progress printing
  invisible(capture.output(val <- expr))
  val
}

# ---------------------------------------------------------------------------
# 1. Load VCF and build genotype dosage matrix (0/1/2 = alt allele count)
# ---------------------------------------------------------------------------
cat("[1/5] Reading VCF and building genotype matrix...\n")
ckpt1_path <- paste0(opt$outprefix, "_ckpt1_genotypes.rds")
ckpt1_key  <- file_sig(opt$vcf)
ck1 <- load_ckpt(ckpt1_path, ckpt1_key)

if (!is.null(ck1)) {
  geno_mat   <- ck1$geno_mat
  locus_info <- ck1$locus_info
  rm(ck1); invisible(gc())
} else {
  vcf <- read.vcfR(opt$vcf, verbose = FALSE)
  gt  <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  gt  <- gsub("|", "/", gt, fixed = TRUE)

  # Only biallelic 0/1 genotypes are coded; anything else becomes NA.
  geno_mat <- matrix(NA_integer_, nrow = nrow(gt), ncol = ncol(gt),
                     dimnames = list(NULL, colnames(gt)))
  geno_mat[which(gt == "0/0")] <- 0L
  geno_mat[which(gt == "0/1" | gt == "1/0")] <- 1L
  geno_mat[which(gt == "1/1")] <- 2L
  rm(gt); invisible(gc())

  locus_info <- data.frame(
    CHROM = vcf@fix[, "CHROM"],
    POS   = as.numeric(vcf@fix[, "POS"]),
    stringsAsFactors = FALSE
  )
  rm(vcf); invisible(gc())

  saveRDS(list(key = ckpt1_key, geno_mat = geno_mat, locus_info = locus_info),
          ckpt1_path, compress = FALSE)
  cat("    Saved checkpoint:", ckpt1_path, "\n")
}
locus_id <- paste(locus_info$CHROM, locus_info$POS, sep = "_")

# --- Popmap: accept Past/Present (or old/new), map internally to old/new ---
popmap <- read.table(opt$popmap, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("sample", "pop"), strip.white = TRUE)
popmap$sample <- trimws(popmap$sample)
popmap$pop    <- tolower(trimws(popmap$pop))

pop_lookup <- c(past = "old", present = "new", old = "old", new = "new")
bad_labels <- setdiff(unique(popmap$pop), names(pop_lookup))
if (length(bad_labels) > 0) {
  stop(sprintf("Unrecognized population label(s) in popmap: %s (expected Past/Present)",
               paste(bad_labels, collapse = ", ")))
}
popmap$pop <- unname(pop_lookup[popmap$pop])

if (anyDuplicated(popmap$sample)) {
  stop("Duplicate sample IDs in popmap: ",
       paste(unique(popmap$sample[duplicated(popmap$sample)]), collapse = ", "))
}

common_samples <- intersect(colnames(geno_mat), popmap$sample)
if (length(common_samples) < nrow(popmap)) {
  missing <- setdiff(popmap$sample, colnames(geno_mat))
  warning(sprintf("Only %d of %d popmap samples found in VCF. Missing: %s",
                   length(common_samples), nrow(popmap),
                   paste(missing, collapse = ", ")))
}
if (length(common_samples) == 0) {
  stop("No popmap sample IDs match VCF sample names. Check naming (e.g. 'normal_F92').")
}
geno_mat <- geno_mat[, common_samples, drop = FALSE]
pop_vec  <- popmap$pop[match(common_samples, popmap$sample)]
if (!all(c("old", "new") %in% pop_vec)) {
  stop("Need at least one Past and one Present sample present in the VCF.")
}

# --- Locus filter: enough calls in both periods, polymorphic overall ---
n_old_called <- rowSums(!is.na(geno_mat[, pop_vec == "old", drop = FALSE]))
n_new_called <- rowSums(!is.na(geno_mat[, pop_vec == "new", drop = FALSE]))
alt_count    <- rowSums(geno_mat, na.rm = TRUE)
allele_total <- 2 * (n_old_called + n_new_called)

keep <- n_old_called >= opt$mincalled &
        n_new_called >= opt$mincalled &
        alt_count > 0 & alt_count < allele_total

cat(sprintf("    Locus filter: kept %d of %d loci (>= %d called per period, polymorphic, biallelic)\n",
            sum(keep), length(keep), opt$mincalled))
geno_mat   <- geno_mat[keep, , drop = FALSE]
locus_id   <- locus_id[keep]
locus_info <- locus_info[keep, , drop = FALSE]
rm(n_old_called, n_new_called, alt_count, allele_total, keep); invisible(gc())
if (nrow(geno_mat) == 0) stop("No loci left after filtering.")

cat(sprintf("    %d loci x %d samples (%d past, %d present)\n",
            nrow(geno_mat), ncol(geno_mat),
            sum(pop_vec == "old"), sum(pop_vec == "new")))

# ---------------------------------------------------------------------------
# 2. OutFLANK genome-wide FST outlier scan
# ---------------------------------------------------------------------------
outflank_table <- NULL
if (isTRUE(opt$run_outflank)) {
  cat("[2/5] Running OutFLANK FST outlier scan...\n")
  ckpt2_path <- paste0(opt$outprefix, "_ckpt2_outflank.rds")
  ckpt2_key  <- paste(ckpt1_key, unname(tools::md5sum(opt$popmap)),
                      opt$mincalled, opt$maxfitloci, opt$qthresh)
  ck2 <- load_ckpt(ckpt2_path, ckpt2_key)

  if (!is.null(ck2)) {
    of_result      <- ck2$of_result
    outflank_table <- ck2$outflank_table
    n_fit          <- ck2$n_fit
    rm(ck2)
  } else {
    chunks <- split(seq_len(nrow(geno_mat)),
                    ceiling(seq_len(nrow(geno_mat)) / opt$chunksize))
    fst_df <- do.call(rbind, lapply(seq_along(chunks), function(k) {
      ix <- chunks[[k]]
      if (k %% 5 == 1 || k == length(chunks)) {
        cat(sprintf("    FST chunk %d/%d  [%s]\n", k, length(chunks), format(Sys.time(), "%H:%M:%S")))
      }
      sm <- t(geno_mat[ix, , drop = FALSE])
      sm[is.na(sm)] <- 9L
      mode(sm) <- "integer"
      quiet(MakeDiploidFSTMat(SNPmat = sm, locusNames = locus_id[ix], popNames = pop_vec))
    }))
    fst_df_clean <- fst_df[!is.na(fst_df$FST) & !is.na(fst_df$FSTNoCorr) & fst_df$He > 0, ]
    rm(fst_df); invisible(gc())

    set.seed(42)
    n_fit   <- min(opt$maxfitloci, nrow(fst_df_clean))
    fit_set <- fst_df_clean[sort(sample.int(nrow(fst_df_clean), n_fit)), ]

    of_result <- OutFLANK(fit_set,
                           LeftTrimFraction = 0.05,
                           RightTrimFraction = 0.05,
                           Hmin = 0.1,
                           NumberOfSamples = 2,
                           qthreshold = opt$qthresh)

    outflank_table <- pOutlierFinderChiSqNoCorr(fst_df_clean,
                                                 Fstbar = of_result$FSTNoCorrbar,
                                                 dfInferred = of_result$dfInferred,
                                                 qthreshold = opt$qthresh,
                                                 Hmin = 0.1)
    outflank_table$OutlierFlag <- outflank_table$OutlierFlag & !is.na(outflank_table$qvalues)
    rm(fst_df_clean); invisible(gc())

    saveRDS(list(key = ckpt2_key, of_result = of_result,
                 outflank_table = outflank_table, n_fit = n_fit), ckpt2_path)
    cat("    Saved checkpoint:", ckpt2_path, "\n")
  }

  cat(sprintf("    Fitted null: FSTNoCorrbar = %.5f, dfInferred = %.3f\n",
              of_result$FSTNoCorrbar, of_result$dfInferred))
  n_outliers <- sum(outflank_table$OutlierFlag, na.rm = TRUE)
  cat(sprintf("    Null fit on %d loci; %d of %d loci are outliers at q < %.3f\n",
              n_fit, n_outliers, nrow(outflank_table), opt$qthresh))

  write.csv(outflank_table[outflank_table$OutlierFlag %in% TRUE, ],
            paste0(opt$outprefix, "_outflank_outliers.csv"), row.names = FALSE)

  p_of <- OutFLANKResultsPlotter(of_result, withOutliers = TRUE,
                                  NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005,
                                  Zoom = FALSE, RightZoomFraction = 0.05,
                                  titletext = "OutFLANK: past vs present LEPC (fit set)")
  ggsave(paste0(opt$outprefix, "_outflank_fst_fit.png"), plot = p_of,
         width = 7, height = 5, dpi = 300)
} else {
  cat("[2/5] Skipping OutFLANK (--run_outflank FALSE)\n")
}

# ---------------------------------------------------------------------------
# 3. Temporal Fk, Ne (Nei & Tajima 1981; Waples 1989), drift null
# ---------------------------------------------------------------------------
cat("[3/5] Estimating Fk, Ne, and the drift null...\n")

old_cols <- pop_vec == "old"
n_past    <- rowSums(!is.na(geno_mat[, old_cols,  drop = FALSE]))
n_present <- rowSums(!is.na(geno_mat[, !old_cols, drop = FALSE]))
p_old <- rowSums(geno_mat[, old_cols,  drop = FALSE], na.rm = TRUE) / (2 * n_past)
p_new <- rowSums(geno_mat[, !old_cols, drop = FALSE], na.rm = TRUE) / (2 * n_present)

freq_df <- data.frame(locus = locus_id, CHROM = locus_info$CHROM, POS = locus_info$POS,
                      n_past = n_past, n_present = n_present,
                      p_old = p_old, p_new = p_new, stringsAsFactors = FALSE)
rm(n_past, n_present, p_old, p_new, geno_mat); invisible(gc())

freq_df <- freq_df[complete.cases(freq_df) &
                    !(freq_df$p_old == freq_df$p_new & freq_df$p_old %in% c(0, 1)), ]
freq_df$Fk <- with(freq_df,
  (p_old - p_new)^2 / (((p_old + p_new) / 2) - p_old * p_new)
)
freq_df <- freq_df[is.finite(freq_df$Fk), ]

# Per-locus sampling correction using actual genotyped individuals
freq_df$Fk_corr <- freq_df$Fk - 1 / (2 * freq_df$n_past) - 1 / (2 * freq_df$n_present)

t_gen   <- opt$generations
ne_from <- function(fc) ifelse(fc > 0, t_gen / (2 * fc), Inf)

p_mean <- (freq_df$p_old + freq_df$p_new) / 2
ne_use <- p_mean >= opt$ne_minfreq & p_mean <= 1 - opt$ne_minfreq

Fk_all      <- mean(freq_df$Fk)
Fcorr_all   <- mean(freq_df$Fk_corr)
Fcorr_ne    <- mean(freq_df$Fk_corr[ne_use])
Ne_all_loci <- ne_from(Fcorr_all)
Ne_est      <- ne_from(Fcorr_ne)

# Block jackknife CI (linked SNPs are not independent)
jk_blocks <- paste(freq_df$CHROM[ne_use],
                   floor(freq_df$POS[ne_use] / opt$jk_blocksize), sep = ":")
x_ne  <- freq_df$Fk_corr[ne_use]
b_sum <- tapply(x_ne, jk_blocks, sum)
b_n   <- tapply(x_ne, jk_blocks, length)
n_blk <- length(b_sum)
if (n_blk >= 2) {
  loo   <- (sum(x_ne) - b_sum) / (length(x_ne) - b_n)
  jk_se <- sqrt((n_blk - 1) / n_blk * sum((loo - mean(loo))^2))
} else {
  jk_se <- NA_real_
}
F_lo  <- Fcorr_ne - 1.96 * jk_se
F_hi  <- Fcorr_ne + 1.96 * jk_se
Ne_lo <- ne_from(F_hi)
Ne_hi <- ne_from(F_lo)
rm(jk_blocks, x_ne, b_sum, b_n, p_mean); invisible(gc())

fmt_ne <- function(x) ifelse(is.finite(x), sprintf("%.1f", x), "Inf")
cat(sprintf("    Mean Fk (all %d loci) = %.5f; sampling-corrected = %.5f -> Ne = %s\n",
            nrow(freq_df), Fk_all, Fcorr_all, fmt_ne(Ne_all_loci)))
cat(sprintf("    Ne loci (mean freq in [%.2f, %.2f]): %d; corrected Fk = %.5f\n",
            opt$ne_minfreq, 1 - opt$ne_minfreq, sum(ne_use), Fcorr_ne))
cat(sprintf("    Ne = %s  (95%% jackknife CI %s - %s; %d blocks of %.0f bp)\n",
            fmt_ne(Ne_est), fmt_ne(Ne_lo), fmt_ne(Ne_hi), n_blk, opt$jk_blocksize))

Ne_for_sim <- if (is.finite(Ne_est) && Ne_est > 1) Ne_est else 500
if (!(is.finite(Ne_est) && Ne_est > 1)) {
  cat("    NOTE: Ne undefined/very large; using Ne = 500 for the drift null.\n")
}

# Drift null, simulated per (starting-frequency bin, n_past, n_present) group
cat("    Simulating Wright-Fisher drift null...\n")
simulate_null_Fk <- function(p0, Ne, t_gen, S0_genes, St_genes, nsim) {
  p_drift <- rep(p0, nsim)
  for (g in seq_len(max(1, round(t_gen)))) {
    p_drift <- rbinom(nsim, size = round(2 * Ne), prob = p_drift) / round(2 * Ne)
  }
  p0_samp <- rbinom(nsim, size = S0_genes, prob = p0) / S0_genes
  pt_samp <- rbinom(nsim, size = St_genes, prob = p_drift) / St_genes
  Fk_null <- (p0_samp - pt_samp)^2 / (((p0_samp + pt_samp) / 2) - p0_samp * pt_samp)
  Fk_null[is.finite(Fk_null)]
}

breaks  <- seq(0, 1, length.out = opt$freqbins + 1)
bin_mid <- pmin(pmax((head(breaks, -1) + tail(breaks, -1)) / 2, 0.01), 0.99)
freq_df$bin <- cut(freq_df$p_old, breaks = breaks, include.lowest = TRUE)
bin_idx <- as.integer(freq_df$bin)

grp_idx <- split(seq_len(nrow(freq_df)),
                 paste(bin_idx, freq_df$n_past, freq_df$n_present))
set.seed(7)
p_null <- rep(NA_real_, nrow(freq_df))
for (g in names(grp_idx)) {
  ix <- grp_idx[[g]]
  parts <- as.integer(strsplit(g, " ", fixed = TRUE)[[1]])
  nullv <- sort(simulate_null_Fk(bin_mid[parts[1]], Ne_for_sim, t_gen,
                                 2 * parts[2], 2 * parts[3], opt$nsim))
  if (length(nullv) == 0) next
  n_ge <- length(nullv) - findInterval(freq_df$Fk[ix], nullv, left.open = TRUE)
  p_null[ix] <- (n_ge + 1) / (length(nullv) + 1)
}
rm(grp_idx, bin_idx); invisible(gc())
freq_df$p_drift_null <- p_null
freq_df$q_drift_null <- p.adjust(p_null, method = "BH")
rm(p_null); invisible(gc())

n_drift_outliers <- sum(freq_df$q_drift_null < opt$qthresh, na.rm = TRUE)
cat(sprintf("    %d loci deviate from the simulated drift null at q < %.3f\n",
            n_drift_outliers, opt$qthresh))

write.csv(freq_df[, c("locus", "CHROM", "POS", "n_past", "n_present", "p_old", "p_new",
                      "Fk", "Fk_corr", "p_drift_null", "q_drift_null")],
          paste0(opt$outprefix, "_temporal_Fk_driftnull.csv"), row.names = FALSE)

sink(paste0(opt$outprefix, "_Ne_summary.txt"))
cat("Temporal Ne / Fk summary (Nei & Tajima 1981; Waples 1989)\n")
cat("===========================================================\n")
cat(sprintf("Past individuals:     %d\n", sum(pop_vec == "old")))
cat(sprintf("Present individuals:  %d\n", sum(pop_vec == "new")))
cat(sprintf("Elapsed generations:  %.2f\n", t_gen))
cat(sprintf("Sampling correction:  per locus, from genotyped individuals\n\n"))
cat(sprintf("All loci (n = %d)\n", nrow(freq_df)))
cat(sprintf("  mean Fk = %.5f; corrected = %.5f; Ne = %s\n\n",
            Fk_all, Fcorr_all, fmt_ne(Ne_all_loci)))
cat(sprintf("Ne loci, mean freq in [%.2f, %.2f] (n = %d)  <- primary estimate\n",
            opt$ne_minfreq, 1 - opt$ne_minfreq, sum(ne_use)))
cat(sprintf("  corrected Fk = %.5f (jackknife SE %.5f, %d blocks of %.0f bp)\n",
            Fcorr_ne, jk_se, n_blk, opt$jk_blocksize))
cat(sprintf("  Ne = %s (95%% CI %s - %s)\n\n", fmt_ne(Ne_est), fmt_ne(Ne_lo), fmt_ne(Ne_hi)))
cat(sprintf("Ne used for drift null: %.1f\n", Ne_for_sim))
sink()

# ---------------------------------------------------------------------------
# 4. Site-class comparison (deleterious vs neutral)
# ---------------------------------------------------------------------------
if (!is.null(opt$siteclass)) {
  cat("[4/5] Comparing allele frequency change between site classes...\n")
  siteclass <- read.table(opt$siteclass, header = FALSE, stringsAsFactors = FALSE,
                           col.names = c("CHROM", "POS", "class"))
  siteclass$class <- tolower(trimws(siteclass$class))
  siteclass <- siteclass[siteclass$class %in% c("deleterious", "neutral"), ]
  siteclass$locus <- paste(siteclass$CHROM, siteclass$POS, sep = "_")
  siteclass <- siteclass[!duplicated(siteclass$locus), ]

  m  <- match(siteclass$locus, freq_df$locus)
  ok <- !is.na(m)
  merged <- freq_df[m[ok], c("locus", "CHROM", "POS", "p_old", "p_new", "Fk",
                             "bin", "p_drift_null", "q_drift_null")]
  merged$class <- siteclass$class[ok]
  cat(sprintf("    %d of %d site-class entries matched analysed loci\n",
              nrow(merged), nrow(siteclass)))
  if (nrow(merged) == 0) {
    cat("    WARNING: no matches. Check that CHROM names in siteclass.txt match the VCF.\n")
  }

  if (length(unique(merged$class)) == 2) {
    is_del <- merged$class == "deleterious"

    # Frequency-stratified test: percentile of Fk within starting-frequency bin
    merged$Fk_pct <- ave(merged$Fk, merged$bin,
                         FUN = function(x) (rank(x) - 0.5) / length(x))
    obs_stat <- mean(merged$Fk_pct[is_del]) - mean(merged$Fk_pct[!is_del])
    strata <- split(seq_len(nrow(merged)), merged$bin, drop = TRUE)
    set.seed(1)
    perm_stat <- replicate(opt$nperm, {
      lab <- is_del
      for (ix in strata) if (length(ix) > 1) lab[ix] <- lab[ix][sample.int(length(ix))]
      mean(merged$Fk_pct[lab]) - mean(merged$Fk_pct[!lab])
    })
    perm_p <- (sum(abs(perm_stat) >= abs(obs_stat)) + 1) / (opt$nperm + 1)

    # Unstratified test, reference only (confounded by allele frequency)
    wt <- wilcox.test(Fk ~ class, data = merged)

    merged$drift_outlier <- merged$q_drift_null < opt$qthresh
    merged$outflank_outlier <- if (!is.null(outflank_table)) {
      merged$locus %in% outflank_table$LocusName[outflank_table$OutlierFlag %in% TRUE]
    } else NA
    # Signed change in ALT frequency. ALT is not necessarily the derived allele;
    # polarize with the chicken outgroup before interpreting direction.
    merged$delta_p_alt <- merged$p_new - merged$p_old

    class_summary <- merged %>%
      group_by(class) %>%
      summarise(n = n(),
                mean_p_past = mean(p_old),
                median_Fk = median(Fk),
                mean_Fk_pct_within_freq_bin = mean(Fk_pct),
                mean_delta_p_alt = mean(delta_p_alt),
                pct_drift_outlier = 100 * mean(drift_outlier, na.rm = TRUE),
                pct_outflank_outlier = 100 * mean(outflank_outlier, na.rm = TRUE),
                .groups = "drop")

    cat(sprintf("    Deleterious: n=%d, median Fk=%.5f, mean within-bin pct=%.4f\n",
                sum(is_del), median(merged$Fk[is_del]), mean(merged$Fk_pct[is_del])))
    cat(sprintf("    Neutral:     n=%d, median Fk=%.5f, mean within-bin pct=%.4f\n",
                sum(!is_del), median(merged$Fk[!is_del]), mean(merged$Fk_pct[!is_del])))
    cat(sprintf("    Frequency-stratified permutation p (%d perms) = %.4g  <- primary\n",
                opt$nperm, perm_p))
    cat(sprintf("    Unstratified Wilcoxon p = %.4g (confounded by frequency; reference only)\n",
                wt$p.value))

    write.csv(class_summary, paste0(opt$outprefix, "_siteclass_summary.csv"), row.names = FALSE)
    write.csv(merged, paste0(opt$outprefix, "_siteclass_merged.csv"), row.names = FALSE)

    p_class <- ggplot(merged, aes(x = bin, y = Fk, fill = class)) +
      geom_boxplot(outlier.alpha = 0.2, position = position_dodge(preserve = "single")) +
      scale_y_continuous(trans = "sqrt") +
      labs(title = "Temporal allele frequency change (Fk) by site class",
           subtitle = "Stratified by Past ALT frequency bin",
           x = "Past ALT frequency bin", y = expression(F[k] ~ "(sqrt scale)"), fill = NULL) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
    ggsave(paste0(opt$outprefix, "_siteclass_Fk_by_freqbin.png"), plot = p_class,
           width = 10, height = 5, dpi = 300)

    sink(paste0(opt$outprefix, "_siteclass_test.txt"))
    cat("Site-class comparison of temporal Fk (deleterious vs neutral)\n")
    cat("================================================================\n")
    cat("PRIMARY: frequency-stratified permutation test\n")
    cat(sprintf("  Statistic: mean within-bin percentile of Fk (deleterious - neutral) = %.5f\n", obs_stat))
    cat(sprintf("  Strata: %d Past-frequency bins; %d permutations within strata\n",
                length(strata), opt$nperm))
    cat(sprintf("  p = %.4g\n\n", perm_p))
    cat("REFERENCE ONLY: unstratified Wilcoxon rank-sum (confounded by frequency)\n")
    print(wt)
    cat("\nPer-class summary:\n")
    print(as.data.frame(class_summary))
    sink()
  } else {
    cat("    Skipping: need both 'deleterious' and 'neutral' classes among analysed loci.\n")
  }
} else {
  cat("[4/5] No --siteclass file provided; skipping site-class comparison.\n")
}

# ---------------------------------------------------------------------------
# 5. Combined outlier table (loci flagged by either method)
# ---------------------------------------------------------------------------
cat("[5/5] Writing combined outlier table...\n")

freq_df$drift_outlier <- freq_df$q_drift_null < opt$qthresh & !is.na(freq_df$q_drift_null)
if (!is.null(outflank_table)) {
  mo <- match(freq_df$locus, outflank_table$LocusName)
  freq_df$FST <- outflank_table$FST[mo]
  freq_df$outflank_q <- outflank_table$qvalues[mo]
  freq_df$outflank_outlier <- !is.na(mo) & (outflank_table$OutlierFlag[mo] %in% TRUE)
  rm(mo)
} else {
  freq_df$FST <- NA_real_
  freq_df$outflank_q <- NA_real_
  freq_df$outflank_outlier <- FALSE
}
freq_df$both_methods_outlier <- freq_df$drift_outlier & freq_df$outflank_outlier

flagged <- freq_df[freq_df$drift_outlier | freq_df$outflank_outlier,
                   c("locus", "CHROM", "POS", "n_past", "n_present", "p_old", "p_new",
                     "Fk", "p_drift_null", "q_drift_null", "FST", "outflank_q",
                     "drift_outlier", "outflank_outlier", "both_methods_outlier")]
names(flagged)[names(flagged) == "p_old"] <- "p_past"
names(flagged)[names(flagged) == "p_new"] <- "p_present"
write.csv(flagged, paste0(opt$outprefix, "_combined_outliers.csv"), row.names = FALSE)

cat(sprintf("\nDone. %d loci flagged by OutFLANK, %d by drift-null, %d by both.\n",
            sum(freq_df$outflank_outlier), sum(freq_df$drift_outlier),
            sum(freq_df$both_methods_outlier)))
cat(sprintf("Outputs written with prefix: %s\n", opt$outprefix))
