#Load libraries
library(readxl)
library(pophelper)
popmap_unrel <- read.csv("~/popmap_unrel.txt", sep="",header = F)
k2<-readQ(files ="~/output_results_log.2.meanQ")
both <- as.data.frame(popmap_unrel)
plotQ(k2,returnplot=T,exportplot=T,clustercol=c("cadetblue","black"),grplab=data.frame(both$V1),ordergrp=T,showlegend=F,height=4,indlabsize=1.2,indlabheight=0.08,indlabspacer=1,barbordercolour="white",divsize = 0.20,grplabsize=2.0,barbordersize=0.1,linesize=0.6,showsp = F,splabsize = 0,outputfilename="past_present_K2",imgtype="pdf",exportpath=getwd(),divtype=2,divcol = "red",splabcol="black",grplabheight=.2)

k3<-readQ(files ="~/output_results_log.3.meanQ")
plotQ(k3,returnplot=T,exportplot=T,clustercol=c("cadetblue","black","grey"),grplab=data.frame(both$V1),ordergrp=T,showlegend=F,height=4,indlabsize=1.2,indlabheight=0.08,indlabspacer=1,barbordercolour="white",divsize = 0.20,grplabsize=2.0,barbordersize=0.1,linesize=0.6,showsp = F,splabsize = 0,outputfilename="past_present_K3",imgtype="pdf",exportpath=getwd(),divtype=2,divcol = "red",splabcol="black",grplabheight=.2)
