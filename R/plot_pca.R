#Load libraries
library(readxl)
library(ggplot2)

#Read in metadata
HET_FILT <- read_xlsx("/Users/andrewblack/Documents/Research/GROUSE/sarek_old_new/old_new_heterozygosity.xlsx")

#Read in covariation matrix
cov<-as.matrix(read.table("~/past_present.cov"))

#Extract and calculate eplained variation
axes<-eigen(cov)
head(axes$values/sum(axes$values)*100)
#6.948021 5.978121 5.870980 5.721814 5.643101 5.545446


#Bind vectors with metadata and plot
PC1_3<-as.data.frame(axes$vectors[,1:3])
x<-cbind(PC1_3,HET_FILT)
 #By species and group
ggplot(data=x, aes(y=V2, x=V1)) +
    geom_point(size=6, color="black", aes(fill=GRP)) +
    theme_classic() +
    xlab("PC1 (6.9%)") + ylab("PC2 (5.9%)") +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    scale_fill_manual("Group", values=c("cadetblue","black"))  +
    theme(legend.position = "right") +
    guides(
        fill = guide_legend(override.aes = list(shape=21, size=5, stroke=0.5)),
        shape = guide_legend(override.aes = list(fill="grey50", size=5))
    )+guides(
        fill = guide_legend(override.aes = list(shape=21, size=5, stroke=0.5)),
        shape = guide_legend(override.aes = list(fill=NA, size=5))
    )

#Note: F10 and F21 are the ones to the right along PC1. Not sure why they stand apart from the other ones. Check to see what the DOC is for these two.
