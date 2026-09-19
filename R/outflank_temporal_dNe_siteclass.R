#!/usr/bin/env Rscript
# ============================================================================
# Old vs New LEPC temporal genomics
#   (1) Temporal Fk / drift-based Ne estimate  (Nei & Tajima 1981; Waples 1989)
#   (2) OutFLANK genome-wide FST outlier scan  (Whitlock & Lotterhos 2015)
#   (3) Site-class comparison of Fk: deleterious vs neutral sites
#
# Usage:
#   Rscript outflank_temporal_dNe_siteclass.R \
#       --vcf old_vs_new.filtered.vcf.gz \
#       --popmap popmap.txt \
#       --siteclass siteclass.txt \
#       --generations 5 \
#       --outprefix results/old_vs_new
#
# Input file formats (no header, whitespace/tab-delimited):
#   popmap.txt    : <sample_id>  <old|new>
#   siteclass.txt : <CHROM>  <POS>  <deleterious|neutral>
#                   (derive this from your SnpEff/GERP annotation output —
#                   e.g. GERP score > threshold or SnpEff HIGH/MODERATE impact
#                   = "deleterious"; synonymous/intergenic = "neutral")
#
# R package requirements:
#   install.packages(c("optparse","dplyr","ggplot2"))
#   install.packages("vcfR")
#   BiocManager::install("qvalue")
#   devtools::install_github("whitlock/OutFLANK")   # not on CRAN
#
# Notes on the drift null (step 1):
#   Rather than relying on OutFLANK's own FST outlier model alone, this script
#   also builds an independent Wright-Fisher drift null: for each locus's
#   starting (old) frequency, it simulates allele frequency drift over the
#   specified number of generations under the estimated genome-wide Ne, adds
#   binomial sampling noise matching your actual sample sizes, and asks how
#   often *neutral* drift alone would produce a Fk as extreme as observed.
#   This gives you a result that isn't dependent on OutFLANK's FST-based
#   assumptions, useful for cross-checking.
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
# 0. CLI args
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--vcf", type = "character", help = "Input VCF (can be .vcf.gz)"),
  make_option("--popmap", type = "character", help = "Sample -> old/new map"),
  make_option("--siteclass", type = "character", default = NULL,
              help = "Optional CHROM POS class(deleterious/neutral) table"),
  make_option("--generations", type = "numeric", default = 1,
              help = "Generations elapsed between old and new samples [default %default]"),
  make_option("--nsim", type = "numeric", default = 2000,
              help = "WF drift simulations per starting-frequency bin [default %default]"),
  make_option("--freqbins", type = "numeric", default = 20,
              help = "Number of starting-frequency bins for the null [default %default]"),
  make_option("--qthresh", type = "numeric", default = 0.05,
              help = "q-value threshold for outliers [default %default]"),
  make_option("--outprefix", type = "character", default = "outflank_temporal",
              help = "Output file prefix [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$vcf) || is.null(opt$popmap)) {
  stop("Must supply --vcf and --popmap. Use --help for usage.")
}

out_dir <- dirname(opt$outprefix)
if (out_dir != ".") dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("== Old vs New LEPC temporal analysis ==\n")
cat("VCF:         ", opt$vcf, "\n")
cat("Popmap:      ", opt$popmap, "\n")
cat("Siteclass:   ", ifelse(is.null(opt$siteclass), "(none provided)", opt$siteclass), "\n")
cat("Generations: ", opt$generations, "\n\n")

# ---------------------------------------------------------------------------
# 1. Load VCF and build genotype dosage matrix (0/1/2 = alt allele count)
# ---------------------------------------------------------------------------
cat("[1/5] Reading VCF and building genotype matrix...\n")
vcf <- read.vcfR(opt$vcf, verbose = FALSE)

gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

dosage_from_gt <- function(g) {
  g <- gsub("\\|", "/", g)
  sapply(strsplit(g, "/"), function(alleles) {
    if (length(alleles) != 2 || any(alleles == ".")) return(NA_integer_)
    sum(as.integer(alleles))
  })
}
geno_mat <- t(apply(gt, 1, dosage_from_gt))   # loci x samples
colnames(geno_mat) <- colnames(gt)

locus_info <- data.frame(
  CHROM = vcf@fix[, "CHROM"],
  POS   = as.numeric(vcf@fix[, "POS"]),
  stringsAsFactors = FALSE
)
locus_id <- paste(locus_info$CHROM, locus_info$POS, sep = "_")
rownames(geno_mat) <- locus_id

popmap <- read.table(opt$popmap, header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("sample", "pop"))
popmap$pop <- tolower(popmap$pop)
stopifnot(all(popmap$pop %in% c("old", "new")))

common_samples <- intersect(colnames(geno_mat), popmap$sample)
if (length(common_samples) < nrow(popmap)) {
  warning(sprintf("Only %d of %d popmap samples found in VCF.",
                   length(common_samples), nrow(popmap)))
}
geno_mat <- geno_mat[, common_samples]
pop_vec  <- popmap$pop[match(common_samples, popmap$sample)]

cat(sprintf("    %d loci x %d samples (%d old, %d new)\n",
            nrow(geno_mat), ncol(geno_mat),
            sum(pop_vec == "old"), sum(pop_vec == "new")))

# ---------------------------------------------------------------------------
# 2. OutFLANK genome-wide FST outlier scan
# ---------------------------------------------------------------------------
cat("[2/5] Running OutFLANK FST outlier scan...\n")

snp_mat <- t(geno_mat)
snp_mat[is.na(snp_mat)] <- 9
mode(snp_mat) <- "integer"

fst_df <- MakeDiploidFSTMat(SNPmat = snp_mat,
                             locusNames = locus_id,
                             popNames = pop_vec)

fst_df_clean <- fst_df[!is.na(fst_df$FST) & fst_df$He > 0, ]

of_result <- OutFLANK(fst_df_clean,
                       LeftTrimFraction = 0.05,
                       RightTrimFraction = 0.05,
                       Hmin = 0.1,
                       NumberOfSamples = 2,
                       qthreshold = opt$qthresh)

outflank_table <- of_result$results
outflank_table$OutlierFlag <- outflank_table$OutlierFlag & !is.na(outflank_table$qvalues)
n_outliers <- sum(outflank_table$OutlierFlag, na.rm = TRUE)
cat(sprintf("    %d loci pass OutFLANK trimming; %d outliers at q < %.3f\n",
            nrow(outflank_table), n_outliers, opt$qthresh))

write.csv(outflank_table, paste0(opt$outprefix, "_outflank_results.csv"), row.names = FALSE)

p_of <- OutFLANKResultsPlotter(of_result, withOutliers = TRUE,
                                NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005,
                                Zoom = FALSE, RightZoomFraction = 0.05,
                                titletext = "OutFLANK: old vs new LEPC")
ggsave(paste0(opt$outprefix, "_outflank_fst_fit.png"), plot = p_of,
       width = 7, height = 5, dpi = 300)

# ---------------------------------------------------------------------------
# 3. Temporal Fk and drift-based Ne estimate (Nei & Tajima 1981; Waples 1989)
# ---------------------------------------------------------------------------
cat("[3/5] Estimating genome-wide Fk and effective population size (Ne)...\n")

S0 <- sum(pop_vec == "old")
St <- sum(pop_vec == "new")

allele_freq <- function(geno_mat, pop_vec, pop_label) {
  sub <- geno_mat[, pop_vec == pop_label, drop = FALSE]
  rowSums(sub, na.rm = TRUE) / (2 * rowSums(!is.na(sub)))
}
p_old <- allele_freq(geno_mat, pop_vec, "old")
p_new <- allele_freq(geno_mat, pop_vec, "new")

freq_df <- data.frame(locus = locus_id, p_old = p_old, p_new = p_new,
                       stringsAsFactors = FALSE)
freq_df <- freq_df[complete.cases(freq_df) &
                    !(freq_df$p_old == freq_df$p_new & freq_df$p_old %in% c(0, 1)), ]

freq_df$Fk <- with(freq_df,
  (p_old - p_new)^2 / (((p_old + p_new) / 2) - p_old * p_new)
)
freq_df <- freq_df[is.finite(freq_df$Fk), ]

Fc_mean <- mean(freq_df$Fk, na.rm = TRUE)
S0_genes <- 2 * S0
St_genes <- 2 * St
t_gen <- opt$generations
Ne_est <- t_gen / (2 * (Fc_mean - 1 / S0_genes - 1 / St_genes))

cat(sprintf("    Genome-wide mean Fk = %.5f across %d loci\n", Fc_mean, nrow(freq_df)))
cat(sprintf("    Elapsed generations = %.2f\n", t_gen))
if (!is.finite(Ne_est) || Ne_est <= 0) {
  cat("    NOTE: Ne estimate is negative/non-finite (Fk too low relative to\n")
  cat("          1/S0+1/St). Falling back to Ne=500 for the null simulation;\n")
  cat("          treat the true Ne as effectively large/undefined here.\n")
} else {
  cat(sprintf("    Estimated Ne (temporal method) = %.1f\n", Ne_est))
}

Ne_for_sim <- ifelse(is.finite(Ne_est) && Ne_est > 1, Ne_est, 500)

cat("    Simulating Wright-Fisher drift null (binned by starting frequency)...\n")
freq_df$bin <- cut(freq_df$p_old, breaks = seq(0, 1, length.out = opt$freqbins + 1),
                    include.lowest = TRUE)

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

bin_levels <- levels(freq_df$bin)
null_by_bin <- lapply(bin_levels, function(b) {
  bounds <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", b), ",")[[1]])
  bin_mid <- mean(bounds)
  bin_mid <- min(max(bin_mid, 0.01), 0.99)
  simulate_null_Fk(bin_mid, Ne_for_sim, t_gen, S0_genes, St_genes, opt$nsim)
})
names(null_by_bin) <- bin_levels

freq_df$p_drift_null <- mapply(function(fk, b) {
  null_vec <- null_by_bin[[as.character(b)]]
  if (is.null(null_vec) || length(null_vec) == 0) return(NA_real_)
  mean(null_vec >= fk, na.rm = TRUE)
}, freq_df$Fk, freq_df$bin)

freq_df$q_drift_null <- p.adjust(freq_df$p_drift_null, method = "BH")
n_drift_outliers <- sum(freq_df$q_drift_null < opt$qthresh, na.rm = TRUE)
cat(sprintf("    %d loci deviate from the simulated drift null at q < %.3f\n",
            n_drift_outliers, opt$qthresh))

write.csv(freq_df, paste0(opt$outprefix, "_temporal_Fk_driftnull.csv"), row.names = FALSE)

sink(paste0(opt$outprefix, "_Ne_summary.txt"))
cat("Temporal Ne / Fk summary (Nei & Tajima 1981; Waples 1989)\n")
cat("===========================================================\n")
cat(sprintf("N loci used:          %d\n", nrow(freq_df)))
cat(sprintf("Old sample size:      %d individuals (%d genes)\n", S0, S0_genes))
cat(sprintf("New sample size:      %d individuals (%d genes)\n", St, St_genes))
cat(sprintf("Elapsed generations:  %.2f\n", t_gen))
cat(sprintf("Genome-wide mean Fk:  %.5f\n", Fc_mean))
cat(sprintf("Estimated Ne:         %s\n",
            ifelse(is.finite(Ne_est) && Ne_est > 0, sprintf("%.1f", Ne_est), "undefined/very large")))
cat(sprintf("Ne used for null sim: %.1f (fallback used if estimate was undefined)\n", Ne_for_sim))
sink()

# ---------------------------------------------------------------------------
# 4. Site-class comparison (deleterious vs neutral)
# ---------------------------------------------------------------------------
if (!is.null(opt$siteclass)) {
  cat("[4/5] Comparing allele frequency change between site classes...\n")
  siteclass <- read.table(opt$siteclass, header = FALSE, stringsAsFactors = FALSE,
                           col.names = c("CHROM", "POS", "class"))
  siteclass$locus <- paste(siteclass$CHROM, siteclass$POS, sep = "_")
  siteclass$class <- tolower(siteclass$class)

  merged <- merge(freq_df, siteclass[, c("locus", "class")], by = "locus")
  merged <- merged[merged$class %in% c("deleterious", "neutral"), ]

  if (length(unique(merged$class)) == 2) {
    wt <- wilcox.test(Fk ~ class, data = merged)

    obs_diff <- with(merged, median(Fk[class == "deleterious"]) - median(Fk[class == "neutral"]))
    n_perm <- 10000
    set.seed(1)
    perm_diffs <- replicate(n_perm, {
      shuffled <- sample(merged$class)
      median(merged$Fk[shuffled == "deleterious"]) - median(merged$Fk[shuffled == "neutral"])
    })
    perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))

    cat(sprintf("    Deleterious sites: n=%d, median Fk=%.5f\n",
                sum(merged$class == "deleterious"),
                median(merged$Fk[merged$class == "deleterious"])))
    cat(sprintf("    Neutral sites:     n=%d, median Fk=%.5f\n",
                sum(merged$class == "neutral"),
                median(merged$Fk[merged$class == "neutral"])))
    cat(sprintf("    Wilcoxon rank-sum p = %.4g\n", wt$p.value))
    cat(sprintf("    Permutation p (10,000 perms) = %.4g\n", perm_p))

    merged$outflank_outlier <- merged$locus %in%
      outflank_table$LocusName[outflank_table$OutlierFlag]
    merged$drift_outlier <- merged$q_drift_null < opt$qthresh

    class_summary <- merged %>%
      group_by(class) %>%
      summarise(n = n(),
                median_Fk = median(Fk),
                pct_outflank_outlier = 100 * mean(outflank_outlier, na.rm = TRUE),
                pct_drift_outlier = 100 * mean(drift_outlier, na.rm = TRUE))

    write.csv(class_summary, paste0(opt$outprefix, "_siteclass_summary.csv"), row.names = FALSE)
    write.csv(merged, paste0(opt$outprefix, "_siteclass_merged.csv"), row.names = FALSE)

    p_class <- ggplot(merged, aes(x = class, y = Fk, fill = class)) +
      geom_boxplot(outlier.alpha = 0.3) +
      scale_y_continuous(trans = "sqrt") +
      labs(title = "Temporal allele frequency change (Fk) by site class",
           x = NULL, y = expression(F[k] ~ "(sqrt scale)")) +
      theme_bw() + theme(legend.position = "none")
    ggsave(paste0(opt$outprefix, "_siteclass_Fk_boxplot.png"), plot = p_class,
           width = 5, height = 5, dpi = 300)

    sink(paste0(opt$outprefix, "_siteclass_test.txt"))
    cat("Site-class comparison of temporal Fk (deleterious vs neutral)\n")
    cat("================================================================\n")
    print(wt)
    cat(sprintf("\nPermutation test on median difference (n=%d): p = %.4g\n", n_perm, perm_p))
    sink()
  } else {
    cat("    Skipping: need both 'deleterious' and 'neutral' classes present.\n")
  }
} else {
  cat("[4/5] No --siteclass file provided; skipping site-class comparison.\n")
}

# ---------------------------------------------------------------------------
# 5. Combined outlier summary (OutFLANK outliers x drift-null outliers)
# ---------------------------------------------------------------------------
cat("[5/5] Writing combined outlier summary...\n")

combined <- merge(outflank_table[, c("LocusName", "FST", "He", "qvalues", "OutlierFlag")],
                   freq_df[, c("locus", "p_old", "p_new", "Fk", "p_drift_null", "q_drift_null")],
                   by.x = "LocusName", by.y = "locus", all = TRUE)
combined$outflank_outlier <- ifelse(is.na(combined$OutlierFlag), FALSE, combined$OutlierFlag)
combined$drift_outlier <- combined$q_drift_null < opt$qthresh
combined$both_methods_outlier <- combined$outflank_outlier & combined$drift_outlier

write.csv(combined, paste0(opt$outprefix, "_combined_outlier_summary.csv"), row.names = FALSE)

cat(sprintf("\nDone. %d loci flagged by OutFLANK, %d by drift-null, %d by both.\n",
            sum(combined$outflank_outlier, na.rm = TRUE),
            sum(combined$drift_outlier, na.rm = TRUE),
            sum(combined$both_methods_outlier, na.rm = TRUE)))
cat(sprintf("Outputs written with prefix: %s\n", opt$outprefix))
