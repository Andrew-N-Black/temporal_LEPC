# =============================================================================
# plot_faststructure.R -- fastStructure ancestry bar plots (K = 2, 3) for the
# 19 unrelated birds, grouped by era, plus the era-by-cluster Fisher's exact
# test. USFWS report, Objective 1, Figure 6.
# Inputs: output_results_log.<K>.meanQ from analysis/PCA/run_faststructure.sh
# and popmap_unrel.txt (column 1 = era, same sample order as the .fam file).
# =============================================================================
library(readxl)
library(pophelper)
popmap_unrel <- read.csv("~/popmap_unrel.txt", sep="",header = F)
k2<-readQ(files ="~/output_results_log.2.meanQ")
both <- as.data.frame(popmap_unrel)
plotQ(k2,returnplot=T,exportplot=T,clustercol=c("cadetblue","black"),grplab=data.frame(both$V1),ordergrp=T,showlegend=F,height=4,indlabsize=1.2,indlabheight=0.08,indlabspacer=1,barbordercolour="white",divsize = 0.20,grplabsize=2.0,barbordersize=0.1,linesize=0.6,showsp = F,splabsize = 0,outputfilename="past_present_K2",imgtype="pdf",exportpath=getwd(),divtype=2,divcol = "red",splabcol="black",grplabheight=.2)

k3<-readQ(files ="~/output_results_log.3.meanQ")
plotQ(k3,returnplot=T,exportplot=T,clustercol=c("cadetblue","black","grey"),grplab=data.frame(both$V1),ordergrp=T,showlegend=F,height=4,indlabsize=1.2,indlabheight=0.08,indlabspacer=1,barbordercolour="white",divsize = 0.20,grplabsize=2.0,barbordersize=0.1,linesize=0.6,showsp = F,splabsize = 0,outputfilename="past_present_K3",imgtype="pdf",exportpath=getwd(),divtype=2,divcol = "red",splabcol="black",grplabheight=.2)

# Era vs. majority cluster at K = 2 (report: Fisher's exact test, P = 0.65)
clu <- apply(as.matrix(k2[[1]]), 1, which.max)
print(table(era = both$V1, cluster = clu))
print(fisher.test(table(both$V1, clu)))
