# =============================================================================
# Reads temporal_metadata.xlsx (all 20 birds, 10 per era). Per-individual
# statistics (heterozygosity, f_ROH, sequencing QC, sampling locations) are
# reported for all 20 birds in the USFWS report (Objective 1, Table 1,
# Figures 2-4, S2); relatedness filtering (F10 removal) applies only to the
# population-level analyses (PCA, fastStructure, F_ST, N_e, load), which use
# the 19-bird set. All values are autosomal (Z scaffolds excluded upstream).
# =============================================================================
#Load library
library(reshape2)
library(ggplot2)

#read in metadata
metadata <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/USFWS_REPORTS/files/temporal_metadata.xlsx")

#Extract relevant information
sub<-metadata[,c("ID","GRP","fROH_100kb-1Mb","fROH_1Mb","fROH_total")]

#Convert to long format
melt_data<-melt(sub,id.vars = c("ID","GRP"))
melt_data$GRP <- factor(melt_data$GRP, levels = c("Past","Present"))

#Plot ROH by groups
ggplot(melt_data, aes(fill=GRP, y=value, x=reorder(ID,value))) +
    geom_bar(stat="identity") +
    facet_grid(variable~GRP, scales="free") +
    theme_classic() +
    theme(
        panel.border = element_rect(color="black", fill=NA, linewidth=1),
        strip.background = element_rect(color="black", fill="grey90", linewidth=1),
        axis.line = element_line(color="black", linewidth=1)
    ) +
    scale_fill_manual("", values=c("Past"="cadetblue","Present"="black")) +
    ylab("fROH") + xlab("Sample (N = 20)") +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=22, face="italic")) +
    theme(strip.text = element_text(size=16))+theme(panel.spacing = unit(0, "lines"))

#Test for normality
shapiro.test(sub$`fROH_100kb-1Mb`)

#W = 0.8886, p-value = 0.02535

shapiro.test(sub$`fROH_1Mb`)
#W = 0.23587, p-value = 2.693e-09

shapiro.test(sub$`fROH_total`)
#W = 0.87335, p-value = 0.01346

#Pairwise test
pairwise.wilcox.test(sub$`fROH_100kb-1Mb`, sub$GRP, p.adjust.method = "BH")

#        Past
#Present 0.53

#P value adjustment method: BH 

