#Load libraries
library(readxl)
library(ggplot2)

#Read in metadata
HET_FILT <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity_unrel.xlsx")

#Read in covariation matrix (rows follow final_cramlist_unrel.txt, NOT the spreadsheet)
cov <- as.matrix(read.table("~/final_autosomal_unrel.cov"))
stopifnot(nrow(cov) == nrow(HET_FILT))   # 19 and 19

#Check sample order: the .cov rows follow final_cramlist_unrel.txt.
#Copy that file locally, then make sure the metadata rows are in the same order.
#order <- sub("\\..*$", "", basename(readLines("~/final_cramlist_unrel.txt")))  # e.g. "F92"
#HET_FILT <- HET_FILT[match(order, HET_FILT$<ID column>), ]            # set your ID column name
#stopifnot(!anyNA(HET_FILT$GRP))

#Extract and calculate explained variation
axes <- eigen(cov)
pve  <- axes$values / sum(axes$values) * 100
round(head(pve), 2)
#6.95 5.98 5.87 5.72 5.64 5.55

#Bind vectors with metadata and plot
PC1_3 <- as.data.frame(axes$vectors[, 1:3])
x <- cbind(PC1_3, HET_FILT)

ggplot(data = x, aes(x = V1, y = V2)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    # shape 21 is a filled circle: 'fill' sets the interior, 'color' the outline.
    # The default shape (19) has no fill, so color = "black" painted every point black.
    geom_point(aes(fill = GRP), shape = 21, color = "black", size = 6, stroke = 0.5) +
    scale_fill_manual("Group", values = c("cadetblue", "black")) +
    xlab(sprintf("PC1 (%.1f%%)", pve[1])) +
    ylab(sprintf("PC2 (%.1f%%)", pve[2])) +
    theme_classic() +
    theme(legend.position = "right")
#Note: F17 and F38 are the ones to the right along PC1. Not sure why they stand apart from the other ones. Check to see what the DOC is for these two.
