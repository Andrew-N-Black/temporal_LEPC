# =============================================================================
# Individual heterozygosity from ANGSD (realSFS) and ROHan, Past vs Present.
# One panel per method; Wilcoxon rank-sum test per method.
# All 20 samples are kept: heterozygosity is an individual-level metric, so the
# related pair (F10, F21) stays in, as in the report.
# =============================================================================
#Load libraries
library(ggplot2)
library(readxl)
library(reshape2)
library(ggpubr)
library(dplyr)

#load metadata
HET_FILT <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity.xlsx")

#Extract relevant information
# heterozygosity20x   = ANGSD/realSFS (heterozygosity4x holds the same values)
# heterozygosity_rohan = ROHan genome-wide estimate (CI in 'low'/'high');
#   a second ROHan column, 'heterozygosity_rohan.1' (CI 'minus'/'plus'), also
#   exists -- swap it in below if that is the estimate you want to report.
het_cols <- c("heterozygosity20x", "heterozygosity_rohan")
sub <- as.data.frame(HET_FILT[, c("ID", "GRP", het_cols)])

#Convert to long format
melt_data <- melt(sub, id.vars = c("ID", "GRP"), measure.vars = het_cols,
                  variable.name = "het_method", value.name = "H")
melt_data$het_method <- factor(melt_data$het_method, levels = het_cols,
                               labels = c("ANGSD", "ROHan"))
melt_data$GRP <- factor(melt_data$GRP, levels = c("Past", "Present"))

#Plot
p <- ggplot(melt_data, aes(x = GRP, y = H, fill = GRP)) +
    # Boxes semi-transparent so same-colored points stay visible on top
    geom_boxplot(outlier.shape = NA) +
    # Points take each group's fill; black outline keeps Present points
    # distinct from the Present box
    geom_jitter(aes(fill = GRP), shape = 21, color = "black", size = 2.5,
                alpha = 0.9, width = 0.12, height = 0) +
    facet_wrap(~ het_method, scales = "free_y") +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x = 1.35, size = 5) +
    scale_fill_manual("", values = c("Past" = "cadetblue", "Present" = "black")) +
    xlab("") + ylab("H") +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", fill = "grey90", linewidth = 1)
    ) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(legend.position = "none") +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 22, face = "italic")) +
    theme(strip.text = element_text(size = 18))
print(p)
ggsave("heterozygosity_angsd_vs_rohan.png", p, width = 8, height = 5, dpi = 300)

#Test for normality, per method (non-normal -> rank-based tests below)
melt_data %>%
    group_by(het_method) %>%
    summarise(shapiro_W = shapiro.test(H)$statistic, shapiro_p = shapiro.test(H)$p.value)

#Wilcoxon rank-sum test of Past vs Present, per method
for (m in levels(melt_data$het_method)) {
    cat("\n==", m, "==\n")
    d <- subset(melt_data, het_method == m)
    print(pairwise.wilcox.test(d$H, d$GRP, p.adjust.method = "BH"))
}

#Summary table (mean +/- SD by method and period)
melt_data %>%
    group_by(het_method, GRP) %>%
    summarise(n = n(), mean_H = mean(H), sd_H = sd(H), median_H = median(H), .groups = "drop")

#Agreement between methods
cor.test(sub$heterozygosity20x, sub$heterozygosity_rohan, method = "pearson")
