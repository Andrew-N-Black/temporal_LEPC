# =============================================================================
# *** CHECK BEFORE USE: reads old_new_heterozygosity.xlsx, the pre-relatedness-
# filter metadata (all 20 samples, including normal_F10). Three sibling
# scripts in this directory -- plot_pca.R, plot_pca_plink.R, and
# PC_correlations.R -- were already updated to read
# old_new_heterozygosity_unrel.xlsx (19 samples, normal_F10 removed as a
# first-degree relative of normal_F21) plus the autosomal-only PCA/ROH
# inputs. This script was not updated to match, and its metadata file may
# also predate the autosomal (Z-scaffold-excluded) confirmation applied
# elsewhere in this repo. Not changed automatically here: the "_unrel"
# workbook's exact columns were not available to verify this plot's fields
# still exist there. See AUDIT.md before using this script's output.
# =============================================================================
library(ggplot2)
library(readxl)
library(ggpubr)
library(dplyr)

#By DPS and Species
#load metadata
QA <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity.xlsx")

#breadth
ggplot(data=QA, aes(y=DOC, x=reorder(ID,DOC),fill=GRP, label = "mean=18.85, range=9.9-25.4"))+geom_bar(stat="identity")+scale_fill_manual("", values =c("Past"="cadetblue","Present"="black"))+xlab("Sample (N=20)")+ylab("Depth of Coverage")+theme(legend.position="top")+theme_classic()+annotate("text", x = 10, y = 30, label = "mean=18.85, range=9.88-25.4",color = "grey", size = 4, fontface = "italic")+ylim(0,30)

#Mapping
ggplot(data=QA, aes(y=properlyPaired, x=reorder(ID,properlyPaired),fill=GRP, label = "mean=18.85, range=9.9-25.4"))+geom_bar(stat="identity")+scale_fill_manual("", values =c("Past"="cadetblue","Present"="black"))+xlab("Sample (N=20)")+ylab("% Mapped")+theme(legend.position="top")+theme_classic()+annotate("text", x = 10, y = 100, label = "mean=96%, range=84.4%-98.9%",color = "grey", size = 4, fontface = "italic")+ylim(0,100)

