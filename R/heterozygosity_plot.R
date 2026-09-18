#Load libraries
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)

#By DPS and Species
#load metadata
HET_FILT <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity.xlsx")
#Plot
ggplot(HET_FILT, aes(x=GRP, y=heterozygosity20x, fill=GRP)) +
    geom_boxplot() +
    scale_fill_manual("", values=c("Past"="cadetblue","Present"="black")) +
    xlab("") + ylab("H") +
    theme_classic() +
    theme(axis.text.y = element_text(size=12)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=22, face="italic")) +
    theme(strip.text = element_text(size=18))


#Test for normality
shapiro.test(HET_FILT$HET)

data:  HET_FILT$HET
W = 0.84473, p-value = 0.004355

#Pairwise test of heterozygosity by ecoregion
pairwise.wilcox.test(HET_FILT$heterozygosity20x, HET_FILT$GRP, p.adjust.method = "BH")

	Pairwise comparisons using Wilcoxon rank sum exact test 

data:  HET_FILT$heterozygosity20x and HET_FILT$GRP 

        Past
Present 0.53

P value adjustment method: BH 
