## ============================================================================
library(readxl)
library(ggplot2)

COV_FILE   <- "~/final_autosomal_unrel.cov"
ORDER_FILE <- "~/final_cramlist_unrel.txt"
META_FILE  <- "/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity_unrel.xlsx"
QC_FILE    <- "~/vcf_samples.tsv"
META_ID    <- "ID"    # <- name of the sample-ID column in META_FILE (e.g. "F17")
META_DEPTH <- "DOC"   # <- name of the depth column in META_FILE
META_GRP   <- "GRP"   # Past / Present

## Short IDs like "F17" from any naming scheme (paths, "normal_F17", etc.)
short_id <- function(x) regmatches(x, regexpr("F[0-9]+", x))

## --- PCA scores, labelled by sample ----------------------------------------
cov  <- as.matrix(read.table(COV_FILE))
ids  <- short_id(basename(readLines(ORDER_FILE)))
stopifnot(length(ids) == nrow(cov), !anyDuplicated(ids))
eig  <- eigen(cov)
pve  <- eig$values / sum(eig$values) * 100
pcs  <- data.frame(id = ids, PC1 = eig$vectors[, 1], PC2 = eig$vectors[, 2])

## --- Metadata: depth + group ----------------------------------------------
meta <- as.data.frame(read_xlsx(META_FILE))
meta <- data.frame(id = short_id(as.character(meta[[META_ID]])),
                   depth = as.numeric(meta[[META_DEPTH]]),
                   group = meta[[META_GRP]])

## --- VCF-based QC: missingness + heterozygosity ----------------------------
qc <- read.delim(QC_FILE)
qc <- data.frame(id = short_id(qc$sample),
                 missing = qc$missing_frac,
                 het = qc$het_rate_polymorphic_sites)

d <- merge(merge(pcs, meta, by = "id"), qc, by = "id")
if (nrow(d) != nrow(pcs)) {
    warning("Samples lost in merge: ", paste(setdiff(pcs$id, d$id), collapse = ", "))
}
print(d[order(d$depth), ], digits = 3, row.names = FALSE)

## --- Correlations -----------------------------------------------------------
## Spearman is the primary test: with n = 19 and a few extreme points, Pearson is
## driven by the outliers themselves. PC signs are arbitrary, so the direction of
## a correlation means nothing; only its strength matters.
## The leave-outliers-out column asks whether a relationship remains once the
## birds that define the PC (F17, F3 on PC1; F340 on PC2) are removed.
outliers <- c("F17", "F3", "F340")
covars <- c(depth = "Depth of coverage (x)", missing = "Missing genotype rate",
            het = "Heterozygosity (polymorphic sites)")

res <- do.call(rbind, lapply(c("PC1", "PC2"), function(pc) {
    do.call(rbind, lapply(names(covars), function(v) {
        s  <- suppressWarnings(cor.test(d[[pc]], d[[v]], method = "spearman"))
        p  <- cor.test(d[[pc]], d[[v]], method = "pearson")
        dd <- d[!d$id %in% outliers, ]
        s2 <- suppressWarnings(cor.test(dd[[pc]], dd[[v]], method = "spearman"))
        data.frame(PC = pc, covariate = v, n = nrow(d),
                   spearman_rho = round(unname(s$estimate), 3), spearman_p = signif(s$p.value, 3),
                   pearson_r = round(unname(p$estimate), 3), pearson_p = signif(p$p.value, 3),
                   rho_without_outliers = round(unname(s2$estimate), 3),
                   p_without_outliers = signif(s2$p.value, 3), n_without = nrow(dd))
    }))
}))
print(res, row.names = FALSE)
write.csv(res, "pc_quality_correlations.csv", row.names = FALSE)

## --- Plot: each PC against each quality metric -------------------------------
long <- do.call(rbind, lapply(c("PC1", "PC2"), function(pc) {
    do.call(rbind, lapply(names(covars), function(v) {
        data.frame(id = d$id, group = d$group, PC = sprintf("%s (%.1f%%)", pc, pve[as.integer(sub("PC", "", pc))]),
                   covariate = covars[[v]], x = d[[v]], y = d[[pc]])
    }))
}))
ann <- unique(data.frame(
    PC = sprintf("%s (%.1f%%)", res$PC, pve[as.integer(sub("PC", "", res$PC))]),
    covariate = covars[res$covariate],
    label = sprintf("rho = %.2f, P = %.2g", res$spearman_rho, res$spearman_p)))

g <- ggplot(long, aes(x = x, y = y)) +
    geom_smooth(method = "lm", se = FALSE, color = "grey60", linewidth = 0.5) +
    geom_point(aes(fill = group), shape = 21, color = "black", size = 3.5) +
    geom_text(data = subset(long, id %in% outliers), aes(label = id),
              vjust = -1, size = 3) +
    geom_text(data = ann, aes(label = label), x = -Inf, y = Inf,
              hjust = -0.5, vjust = 1.3, size = 3, inherit.aes = FALSE) +
    facet_wrap(PC ~ covariate, scales = "free") +
    scale_fill_manual("Group", values = c(Past = "cadetblue", Present = "black")) +
    labs(x = NULL, y = "PC score") +
    theme_classic() +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
print(g)
ggsave("pc_vs_quality.png", g, width = 11, height = 6, dpi = 300)
