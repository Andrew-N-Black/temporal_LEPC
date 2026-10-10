# =============================================================================
# Reads temporal_metadata.xlsx (all 20 birds, 10 per era). Per-individual
# statistics (heterozygosity, f_ROH, sequencing QC, sampling locations) are
# reported for all 20 birds in the USFWS report (Objective 1, Table 1,
# Figures 2-4, S2); relatedness filtering (F10 removal) applies only to the
# population-level analyses (PCA, fastStructure, F_ST, N_e, load), which use
# the 19-bird set. All values are autosomal (Z scaffolds excluded upstream).
# =============================================================================
#Load libraries
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)

#By DPS and Species
#load metadata
HET_FILT <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files/temporal_metadata.xlsx")
#Plot
ggplot(HET_FILT, aes(x=GRP, y=heterozygosity_angsd, fill=GRP)) +
    geom_boxplot() +
    scale_fill_manual("", values=c("Past"="cadetblue","Present"="black")) +
    xlab("") + ylab("H") +
    theme_classic() +
    theme(axis.text.y = element_text(size=12)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=22, face="italic")) +
    theme(strip.text = element_text(size=18))


#Test for normality
shapiro.test(HET_FILT$heterozygosity_angsd)

#data:  HET_FILT$HET
#W = 0.84473, p-value = 0.004355

#Pairwise test of heterozygosity by ecoregion
pairwise.wilcox.test(HET_FILT$heterozygosity_angsd, HET_FILT$GRP, p.adjust.method = "BH")

#	Pairwise comparisons using Wilcoxon rank sum exact test 

#data:  HET_FILT$heterozygosity_angsd and HET_FILT$GRP 

#        Past
#Present 0.53

#P value adjustment method: BH 
