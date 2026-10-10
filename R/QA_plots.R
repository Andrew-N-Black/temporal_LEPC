# =============================================================================
# Reads temporal_metadata.xlsx (all 20 birds, 10 per era). Per-individual
# statistics (heterozygosity, f_ROH, sequencing QC, sampling locations) are
# reported for all 20 birds in the USFWS report (Objective 1, Table 1,
# Figures 2-4, S2); relatedness filtering (F10 removal) applies only to the
# population-level analyses (PCA, fastStructure, F_ST, N_e, load), which use
# the 19-bird set. All values are autosomal (Z scaffolds excluded upstream).
# =============================================================================
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)

#By DPS and Species
#load metadata
QA <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files/temporal_metadata.xlsx")

#depth
ggplot(data=QA, aes(y=DOC, x=reorder(ID,DOC),fill=GRP, label = "mean=18.85, range=9.9-25.4"))+geom_bar(stat="identity")+scale_fill_manual("", values =c("Past"="cadetblue","Present"="black"))+xlab("Sample (N=20)")+ylab("Depth of Coverage")+theme(legend.position="top")+theme_classic()+annotate("text", x = 10, y = 30, label = "mean=18.85, range=9.88-25.4",color = "grey", size = 4, fontface = "italic")+ylim(0,30)
ggsave("~/temporal_DOC.jpeg")
#Mapping
ggplot(data=QA, aes(y=properlyPaired, x=reorder(ID,properlyPaired),fill=GRP, label = "mean=18.85, range=9.9-25.4"))+geom_bar(stat="identity")+scale_fill_manual("", values =c("Past"="cadetblue","Present"="black"))+xlab("Sample (N=20)")+ylab("% Mapped")+theme(legend.position="top")+theme_classic()+annotate("text", x = 10, y = 100, label = "mean=96%, range=84.4%-98.9%",color = "grey", size = 4, fontface = "italic")+ylim(0,100)
ggsave("~/temporal_mapping.jpeg")
