# =============================================================================
# fROH per individual from bcftools roh and ROHan, by ROH length class
# Facets: roh_method ~ roh_category; bars = individuals, Past then Present,
# each sorted by ascending fROH within its panel.
# =============================================================================
library(readxl)
library(reshape2)
library(ggplot2)

# Read in metadata
metadata <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity.xlsx")

# Extract fROH columns from both methods
froh_cols <- c("fROH_100kb-1Mb", "fROH_1Mb", "fROH_total",
               "fROH_100kb-1Mb_rohan", "fROH_1Mb_rohan", "fROH_total_rohan")
sub <- as.data.frame(metadata[, c("ID", "GRP", froh_cols)])

# Convert to long format
melt_data <- melt(sub, id.vars = c("ID", "GRP"), measure.vars = froh_cols,
                  variable.name = "variable", value.name = "value")
melt_data$variable <- as.character(melt_data$variable)

# Split each column name into method and length class
melt_data$roh_method <- ifelse(grepl("_rohan$", melt_data$variable), "ROHan", "bcftools roh")
melt_data$roh_category <- sub("_rohan$", "", melt_data$variable)
melt_data$roh_category <- factor(melt_data$roh_category,
                                 levels = c("fROH_100kb-1Mb", "fROH_1Mb", "fROH_total"),
                                 labels = c("100 kb–1 Mb", ">1 Mb", "Total"))
melt_data$roh_method <- factor(melt_data$roh_method, levels = c("bcftools roh", "ROHan"))
melt_data$GRP <- factor(melt_data$GRP, levels = c("Past", "Present"))

# Order bars within each panel: Past first, then Present, each ascending.
# One x level per bar and panel, so free x scales keep each panel's own order.
melt_data$bar <- paste(melt_data$roh_method, melt_data$roh_category, melt_data$ID, sep = "|")
ord <- order(melt_data$GRP, melt_data$value)
melt_data$bar <- factor(melt_data$bar, levels = unique(melt_data$bar[ord]))

# Plot fROH by method and length class
p <- ggplot(melt_data, aes(fill = GRP, y = value, x = bar)) +
    geom_bar(stat = "identity") +
    facet_wrap(roh_method ~ roh_category, scales = "free", nrow = 2) +
    theme_classic() +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(color = "black", fill = "grey90", linewidth = 1),
        axis.line = element_line(color = "black", linewidth = 1)
    ) +
    scale_fill_manual("", values = c("Past" = "cadetblue", "Present" = "black")) +
    ylab("fROH") + xlab(paste0("Sample (N=", length(unique(melt_data$ID)), ")")) +
    theme(legend.position = "top") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title = element_text(size = 22, face = "italic")) +
    theme(strip.text = element_text(size = 14), legend.text = element_text(size = 14)) +
    theme(panel.spacing = unit(0.3, "lines"))
print(p)
ggsave("fROH_bcftools_vs_rohan.png", p, width = 12, height = 7, dpi = 300)
