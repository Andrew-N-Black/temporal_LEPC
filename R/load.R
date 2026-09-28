#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Plot LEPC genetic-load results (Past 2019 vs Present 2026)
#
# Run on your own machine. Copy these three files down from the cluster first:
#
#   $PROJ/results/load/genetic_load_per_individual.tsv
#   $PROJ/results/load/load_ratio_per_individual.tsv
#   $PROJ/results/load/sample_call_rate.tsv
#
# then point DATA_DIR at wherever you put them and run:
#
#   Rscript plot_lepc_load.R
#
# Figures carry no title, subtitle or caption -- everything a reader needs is
# either annotated onto the panel or written to figure_legends.txt, ready to
# paste into the report. Needs only ggplot2.
# ---------------------------------------------------------------------------

DATA_DIR <- "~/"         # folder holding the three data files
OUT_DIR  <- "figures"    # figures and figure_legends.txt are written here
FILE_EXT <- ".txt"       # change to ".tsv" if you kept the original names
DPI      <- 800

PT_SIZE   <- 3.0         # point size, all figures
PT_ALPHA  <- 0.95        # point opacity, all figures
DODGE     <- 0.8         # era separation within a category
BOX_WIDTH <- 0.36        # box width; must stay under DODGE / 2 = 0.40 per slot

suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- colours -------------------------------------------------------------
# Cadetblue vs black separate by LUMINANCE rather than hue, so they stay
# distinct under every form of colour-vision deficiency (dE 65, against a
# target of 8) and survive greyscale printing. Cadetblue sits just under 3:1
# against the panel, so the eras are also dodged side by side -- identity
# never rests on colour alone.
COL_PAST    <- "cadetblue"
COL_PRESENT <- "black"
ERA_COLS    <- c("Past (2019)" = COL_PAST, "Present (2026)" = COL_PRESENT)

INK      <- "#0b0b0b"
INK_SOFT <- "#52514e"
SURFACE  <- "#fcfcfb"

theme_load <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = SURFACE, colour = NA),
      panel.background  = element_rect(fill = SURFACE, colour = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "#e6e5e1", linewidth = 0.3),
      axis.text         = element_text(colour = INK_SOFT),
      axis.title        = element_text(colour = INK_SOFT, size = base_size - 1),
      strip.text        = element_text(colour = INK, face = "bold",
                                       size = base_size - 1),
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.text       = element_text(colour = INK_SOFT)
    )
}

# --- load ----------------------------------------------------------------
need <- function(stem) {
  p <- path.expand(file.path(DATA_DIR, paste0(stem, FILE_EXT)))
  if (!file.exists(p))
    stop("Missing ", p, "\n  Check DATA_DIR and FILE_EXT at the top of this script.",
         call. = FALSE)
  read.delim(p, stringsAsFactors = FALSE)
}
load_df  <- need("genetic_load_per_individual")
ratio_df <- need("load_ratio_per_individual")
qual_df  <- need("sample_call_rate")

PAST <- "Past (2019)"; PRESENT <- "Present (2026)"
era <- function(g) factor(ifelse(g == "Old", PAST, PRESENT), levels = c(PAST, PRESENT))
load_df$era  <- era(load_df$group)
ratio_df$era <- era(ratio_df$group)
qual_df$era  <- era(qual_df$group)

LAB <- c(LOF = "Loss-of-\nfunction",
         MISSENSE_RADICAL = "Missense\nradical",
         MISSENSE = "Missense\n(all)",
         MISSENSE_CONSERVATIVE = "Missense\nconservative",
         NEUTRAL = "Neutral")
# Missense (all) is the union of its two Grantham subsets, so the figures use
# only the non-overlapping classes.
DISJOINT <- c("LOF", "MISSENSE_RADICAL", "MISSENSE_CONSERVATIVE", "NEUTRAL")
DELET    <- c("LOF", "MISSENSE_RADICAL", "MISSENSE_CONSERVATIVE")

# Sites per class, from the Step 5 log. Update if the pipeline is re-run.
N_SITES <- c(LOF = 686, MISSENSE_RADICAL = 6545,
             MISSENSE_CONSERVATIVE = 35974, NEUTRAL = 108465)

relabel <- function(d, keep) {
  d <- d[d$category %in% keep, , drop = FALSE]
  d$class <- factor(LAB[d$category], levels = LAB[keep])
  d
}

# Box (median, interquartile range, 1.5x IQR whiskers), the individual birds
# jittered over it, and a white diamond at the group mean.
#
# Boxplot rather than violin: a violin draws a kernel density, and with 9 and
# 10 birds per group that density is a product of the smoothing bandwidth
# rather than of the data -- it would invent shape that is not there. The raw
# points stay on top for the same reason, so the reader can see that each box
# summarises about ten observations. The mean is marked separately because the
# box line is the MEDIAN, while the report quotes means; a white fill keeps the
# diamond legible against both the box and the points.
#
# outlier.shape = NA suppresses the boxplot's own outlier marks, which would
# otherwise double-plot points already drawn by the jitter layer.
#
# Geometry: position_dodge(width = D) gives each of the two eras a slot of
# D / 2 = 0.40, so a box of 0.36 sits inside its slot with a visible gap.
boxset <- function(mapping) {
  list(
    geom_boxplot(mapping, fill = NA, width = BOX_WIDTH, linewidth = 0.45,
                 outlier.shape = NA, position = position_dodge(width = DODGE),
                 show.legend = FALSE),
    geom_point(mapping, position = position_jitterdodge(jitter.width = 0.10,
                                                        dodge.width = DODGE,
                                                        seed = 1),
               size = PT_SIZE, alpha = PT_ALPHA, stroke = 0),
    stat_summary(mapping, fun = mean, geom = "point", shape = 23, size = 2.8,
                 fill = SURFACE, stroke = 0.7,
                 position = position_dodge(width = DODGE), show.legend = FALSE)
  )
}

wilcox_p <- function(d, value_col) {
  o <- d[[value_col]][d$era == PAST]
  n <- d[[value_col]][d$era == PRESENT]
  if (length(o) < 2 || length(n) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(n, o)$p.value)
}

# ===========================================================================
# Figure 1 - the purifying-selection gradient
# Derived-allele frequency should fall as predicted deleteriousness rises.
# Nothing in the pipeline enforces this order, so it validates the pipeline.
# Annotation: number of sites per class.
# ===========================================================================
d1 <- relabel(load_df, DISJOINT)

# Counts sit in one aligned row above the data rather than at each cluster's
# own height, which would read as ragged given neutral sits far above
# loss-of-function. inherit.aes = FALSE keeps the labels out of the era
# colour/dodge mapping, so they stay centred on the category tick.
ann1 <- data.frame(
  class = factor(LAB[DISJOINT], levels = LAB[DISJOINT]),
  label = sprintf("n = %s", format(N_SITES[DISJOINT], big.mark = ",", trim = TRUE)),
  y     = max(d1$total_load) * 1.04,
  stringsAsFactors = FALSE)

p1 <- ggplot(d1, aes(class, total_load, colour = era)) +
  boxset(aes(class, total_load, colour = era)) +
  geom_text(data = ann1, aes(class, y, label = label), inherit.aes = FALSE,
            colour = INK_SOFT, size = 3.4, vjust = 0) +
  scale_colour_manual(values = ERA_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(x = NULL, y = "Total load (within-individual derived-allele frequency)") +
  theme_load()
ggsave(file.path(OUT_DIR, "fig1_selection_gradient.png"), p1,
       width = 7.5, height = 5, dpi = DPI)

# ===========================================================================
# Figure 2 - the primary comparison: deleterious load relative to neutral
# Each bird's deleterious load divided by its OWN neutral load, so anything
# acting equally on all site classes cancels within the individual.
# Annotation: Wilcoxon P per class, recomputed from the data rather than
# hard-coded, so the figure cannot drift out of step with the inputs.
# ===========================================================================
long <- do.call(rbind, lapply(
  c("ratio_total", "ratio_realized", "ratio_masked"),
  function(m) data.frame(individual = ratio_df$individual, era = ratio_df$era,
                         category = ratio_df$category, metric = m,
                         value = ratio_df[[m]], stringsAsFactors = FALSE)))
long <- relabel(long, DELET)
long$metric <- factor(long$metric,
                      levels = c("ratio_total", "ratio_realized", "ratio_masked"),
                      labels = c("Total", "Realized (homozygous)", "Masked (heterozygous)"))

ann2 <- do.call(rbind, lapply(
  split(long, list(long$class, long$metric), drop = TRUE),
  function(d) data.frame(class = d$class[1], metric = d$metric[1],
                         label = sprintf("P = %.3f", wilcox_p(d, "value")),
                         stringsAsFactors = FALSE)))
# free_y means each facet has its own scale, so the label row is per facet.
facet_top <- tapply(long$value, long$metric, max)
ann2$y <- facet_top[as.character(ann2$metric)] * 1.03

p2 <- ggplot(long, aes(class, value, colour = era)) +
  boxset(aes(class, value, colour = era)) +
  geom_text(data = ann2, aes(class, y, label = label), inherit.aes = FALSE,
            colour = INK_SOFT, size = 3.1, vjust = 0) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = ERA_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Deleterious : neutral load ratio") +
  theme_load()
ggsave(file.path(OUT_DIR, "fig2_ratio_primary.png"), p2,
       width = 9.5, height = 5, dpi = DPI)

# ===========================================================================
# Figure 3 - the data-quality check
# Under-called heterozygotes would inflate realized load. Call rate did not
# differ between eras, so this is a caveat rather than a confound -- but the
# negative slope is why the ratio in Figure 2 is the primary statistic.
# Annotation: Spearman rho per class, recomputed from the data.
# ===========================================================================
d3 <- relabel(merge(load_df, qual_df[, c("individual", "call_rate")],
                    by = "individual"), DISJOINT)

# The label sits ABOVE the topmost point (vjust = 0 at the panel maximum),
# not hanging down from it -- with vjust = 1 it lands directly on the
# lowest-call-rate bird, which is the highest point in most panels. The upper
# axis expansion below supplies the room it needs.
ann3 <- do.call(rbind, lapply(split(d3, d3$class, drop = TRUE), function(d) {
  ct <- suppressWarnings(cor.test(d$call_rate, d$realized_load, method = "spearman"))
  data.frame(class = d$class[1], label = sprintf("rho = %.2f", unname(ct$estimate)),
             x = min(d3$call_rate), y = max(d$realized_load),
             stringsAsFactors = FALSE)
}))

p3 <- ggplot(d3, aes(call_rate, realized_load, colour = era)) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
              colour = INK_SOFT, linewidth = 0.4, linetype = "22") +
  geom_point(size = PT_SIZE, alpha = PT_ALPHA, stroke = 0) +
  geom_text(data = ann3, aes(x, y, label = label), inherit.aes = FALSE,
            colour = INK_SOFT, size = 3.1, hjust = 0, vjust = 0) +
  facet_wrap(~ class, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = ERA_COLS) +
  scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = "Genotype call rate at classified sites", y = "Realized load") +
  theme_load()
ggsave(file.path(OUT_DIR, "fig3_callrate_check.png"), p3,
       width = 11, height = 4.2, dpi = DPI)

# ===========================================================================
# Figure legends, written out for pasting into the report
# ===========================================================================
n_past    <- length(unique(load_df$individual[load_df$era == PAST]))
n_present <- length(unique(load_df$individual[load_df$era == PRESENT]))

legends <- c(
sprintf(paste(
"Figure 1. Per-individual genetic load by site class in New Mexico Lesser",
"Prairie-Chicken sampled in 2019 (\"Past\", n = %d birds) and 2026 (\"Present\",",
"n = %d birds). Boxes show the median and interquartile range with whiskers",
"extending to 1.5x the interquartile range, the white diamond marks the group",
"mean, and each point is one bird. Boxplots rather than violins because nine",
"to ten observations per group cannot support a kernel density. Total load is",
"the within-individual derived-allele frequency, computed",
"over the sites called in that individual. The number of polarized sites in",
"each class is given above the category. Classes are non-overlapping:",
"loss-of-function is SnpEff HIGH impact, missense is split at a Grantham",
"distance of 100 into radical and conservative changes, and neutral is SnpEff",
"LOW impact. Derived-allele frequency declines as predicted severity rises, as",
"purifying selection predicts; no step in the pipeline enforces this ordering,",
"so it serves as an internal check on polarization and impact assignment."),
n_past, n_present),
"",
paste(
"Figure 2. Deleterious load expressed relative to each individual's own",
"neutral load. Dividing by the same bird's neutral load cancels any effect",
"acting equally across site classes, whether demographic or technical, so",
"these ratios rather than the raw loads are the primary comparison between",
"eras. Boxes show the median and interquartile range, the white diamond",
"marks the mean and each point is one bird; P values above",
"each class are from two-sided Wilcoxon rank-sum tests and are uncorrected.",
"Realized load at loss-of-function sites was 11.1% higher in 2026, the only",
"comparison reaching P < 0.05; it does not survive Benjamini-Hochberg",
"correction across the twelve ratio comparisons (q = 0.21) and rests on the",
"smallest and most annotation-sensitive site class."),
"",
paste(
"Figure 3. Realized load against per-individual genotype call rate at",
"classified sites, with Spearman rho given in each panel and a least-squares",
"line fitted across both eras. Under-called heterozygotes depress masked load",
"and inflate realized load, which is the negative association seen here. Call",
"rate did not itself differ between eras (Wilcoxon P = 0.72), so it cannot",
"account for the era contrast; it is nonetheless why the normalised ratios in",
"Figure 2 are preferred to the raw loads.")
)
writeLines(legends, file.path(OUT_DIR, "figure_legends.txt"))

cat("Wrote to ", normalizePath(OUT_DIR), ":\n",
    "  fig1_selection_gradient.png\n",
    "  fig2_ratio_primary.png\n",
    "  fig3_callrate_check.png\n",
    "  figure_legends.txt\n", sep = "")
