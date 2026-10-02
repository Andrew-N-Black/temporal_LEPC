# =============================================================================
# Genetic load: tables, robustness checks and figures for the USFWS report
# (Tables S2-S3, Figures S1-S3).
#
# Input: genetic_load_per_individual.tsv from analysis/load/load.sh (Step 6),
# which now has three site classes only -- LOF, MISSENSE, NEUTRAL. The
# Grantham radical/conservative missense split was removed (2026-10).
#
# Uses the 19-sample unrelated set (normal_F10 excluded) from
# processing/popmap_unrelated.txt. Run load.sh on the autosomal VCF first
# (see AUDIT.md item 1) so these numbers match the rest of the report.
# =============================================================================
suppressPackageStartupMessages({
    library(tidyverse)
})

LOAD_TSV <- "genetic_load_per_individual.tsv"      # copy from $OUT/load/ on Gautschi
POPMAP   <- "processing/popmap_unrelated.txt"
OUTDIR   <- "load_report"
dir.create(OUTDIR, showWarnings = FALSE)

classes   <- c(LOF = "Loss-of-function", MISSENSE = "Missense", NEUTRAL = "Neutral")
metrics   <- c(total_load = "Total", realized_load = "Realized", masked_load = "Masked")
era_cols  <- c("Past (2019)" = "cadetblue", "Present (2026)" = "black")

pop <- read_tsv(POPMAP, col_names = c("individual", "era"), show_col_types = FALSE)
d <- read_tsv(LOAD_TSV, show_col_types = FALSE) %>%
    inner_join(pop, by = "individual") %>%                 # drops normal_F10
    filter(category %in% names(classes)) %>%
    mutate(era   = factor(ifelse(era == "Past", "Past (2019)", "Present (2026)"),
                          levels = names(era_cols)),
           class = factor(classes[category], levels = classes))
stopifnot(n_distinct(d$individual) == 19)

wt <- function(v, era) {                                   # Present vs Past
    w <- suppressWarnings(wilcox.test(v[era == "Present (2026)"], v[era == "Past (2019)"]))
    c(W = unname(w$statistic), P = w$p.value)
}

# ---- Call rate at classified sites -----------------------------------------
cr <- d %>% group_by(individual, era) %>%
    summarise(call_rate = sum(n_sites) / sum(n_sites_class), .groups = "drop")
d  <- d %>% left_join(cr, by = c("individual", "era"))
cat("Call rate by era:\n"); print(cr %>% group_by(era) %>% summarise(mean = mean(call_rate)))
cat("Call rate Past vs Present: P =", signif(wt(cr$call_rate, cr$era)["P"], 3), "\n")

cr_cor <- d %>% pivot_longer(names(metrics), names_to = "metric", values_to = "load") %>%
    group_by(class, metric) %>%
    summarise(rho = cor(load, call_rate, method = "spearman"), .groups = "drop")
write_tsv(cr_cor, file.path(OUTDIR, "load_vs_callrate_spearman.tsv"))

# ---- Table S2: raw load by class --------------------------------------------
long <- d %>% pivot_longer(names(metrics), names_to = "metric", values_to = "value")
s2 <- long %>% group_by(class, metric) %>%
    summarise(n_sites = max(n_sites_class),
              past = mean(value[era == "Past (2019)"]),
              present = mean(value[era == "Present (2026)"]),
              W = wt(value, era)["W"], P = wt(value, era)["P"], .groups = "drop") %>%
    mutate(metric = metrics[metric])
write_tsv(s2, file.path(OUTDIR, "TableS2_load_by_class.tsv"))

# ---- Table S3: deleterious : neutral ratios (primary comparison) -------------
neut <- long %>% filter(category == "NEUTRAL") %>% select(individual, metric, neutral = value)
ratios <- long %>% filter(category != "NEUTRAL") %>%
    left_join(neut, by = c("individual", "metric")) %>%
    mutate(ratio = value / neutral)

loo_range <- function(df) {                               # leave-one-out P range
    p <- sapply(unique(df$individual), function(i)
        with(filter(df, individual != i), wt(ratio, era)["P"]))
    sprintf("%.3f-%.3f", min(p), max(p))
}
s3 <- ratios %>% group_by(class, metric) %>%
    group_modify(~ tibble(past = mean(.x$ratio[.x$era == "Past (2019)"]),
                          present = mean(.x$ratio[.x$era == "Present (2026)"]),
                          W = wt(.x$ratio, .x$era)["W"], P = wt(.x$ratio, .x$era)["P"],
                          P_loo = loo_range(.x),
                          P_callrate_ge_0.98 = with(filter(.x, call_rate >= 0.98),
                                                    wt(ratio, era)["P"]))) %>%
    ungroup() %>%
    mutate(BH_q = p.adjust(P, "BH"), metric = metrics[metric])   # 6 comparisons
write_tsv(s3, file.path(OUTDIR, "TableS3_load_ratios.tsv"))
print(s3)

# ---- Figures ------------------------------------------------------------------
box <- function(p) p +
    geom_boxplot(aes(colour = era), fill = NA, outlier.shape = NA, width = 0.35,
                 position = position_dodge(0.8)) +
    geom_point(aes(colour = era), size = 2.5,
               position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8)) +
    stat_summary(aes(group = era), fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "white", colour = "black", position = position_dodge(0.8)) +
    scale_colour_manual(NULL, values = era_cols) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

n_lab <- d %>% distinct(class, n_sites_class) %>% group_by(class) %>%
    summarise(lab = paste0("n = ", format(max(n_sites_class), big.mark = ",")))
f1 <- box(ggplot(d, aes(class, total_load))) +
    geom_text(data = n_lab, aes(class, Inf, label = lab), vjust = 1.5, inherit.aes = FALSE) +
    labs(x = NULL, y = "Total load (within-individual derived-allele frequency)")
ggsave(file.path(OUTDIR, "FigS1_load_by_class.png"), f1, width = 8, height = 5.5, dpi = 300)

p_lab <- s3 %>% mutate(lab = sprintf("P = %.3f", P),
                       metric = factor(metric, levels = metrics))
f2 <- box(ggplot(ratios %>% mutate(metric = factor(metrics[metric], levels = metrics)),
                 aes(class, ratio))) +
    geom_text(data = p_lab, aes(class, Inf, label = lab), vjust = 1.5, inherit.aes = FALSE) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = NULL, y = "Deleterious : neutral load ratio")
ggsave(file.path(OUTDIR, "FigS2_load_ratios.png"), f2, width = 10, height = 5, dpi = 300)

rho_lab <- cr_cor %>% filter(metric == "realized_load") %>% mutate(lab = sprintf("rho = %.2f", rho))
f3 <- ggplot(d, aes(call_rate, realized_load)) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dotted", colour = "grey30", linewidth = 0.5) +
    geom_point(aes(colour = era), size = 2.5) +
    geom_text(data = rho_lab, aes(-Inf, Inf, label = lab), hjust = -0.1, vjust = 1.5,
              inherit.aes = FALSE) +
    facet_wrap(~ class, scales = "free_y") +
    scale_colour_manual(NULL, values = era_cols) +
    labs(x = "Genotype call rate at classified sites", y = "Realized load") +
    theme_minimal(base_size = 13) + theme(legend.position = "top")
ggsave(file.path(OUTDIR, "FigS3_realized_vs_callrate.png"), f3, width = 10, height = 4, dpi = 300)

cat("\nWrote tables and figures to", OUTDIR, "\n")
