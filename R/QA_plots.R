

depth_breadth <- read_excel("metadata.xlsx")
#breadth
ggplot(data=depth_breadth, aes(y=BREADTH_POST, x=reorder(ID,BREADTH_POST),fill=SPECIES))+geom_bar(stat="identity")+scale_fill_manual("", values =c("Tympanuchus pallidicinctus"="brown","Tympanuchus cupido"="goldenrod"))+xlab("Sample (N=481)")+ylab("Genome Breadth")+theme( axis.ticks.x=element_blank(),axis.text.x = element_blank())+ geom_hline(yintercept=80, linetype="dashed", color = "white", linewidth=1)+ylim(0,100)+theme(legend.position="none")
#mapping
ggplot(data=depth_breadth, aes(y=MAPPING_TOTAL, x=reorder(ID,MAPPING_TOTAL),fill=SPECIES))+geom_bar(stat="identity")+scale_fill_manual("", values =c("Tympanuchus pallidicinctus"="brown","Tympanuchus cupido"="goldenrod"))+xlab("Sample (N=481)")+ylab("Reads Aligned (%)")+theme( axis.ticks.x=element_blank(),axis.text.x = element_blank())+ylim(c(0,100))+ geom_hline(yintercept=80, linetype="dashed", color = "white", linewidth=1)+theme(legend.position="top")
#depth
ggplot(data=depth_breadth, aes(y=DEPTH_POST, x=reorder(ID,DEPTH_POST),fill=SPECIES))+geom_bar(stat="identity")+scale_fill_manual("", values =c("Tympanuchus pallidicinctus"="brown","Tympanuchus cupido"="goldenrod"))+xlab("Sample (N=481)")+ylab("Depth Of Coverage")+theme( axis.ticks.x=element_blank(),axis.text.x = element_blank())+ geom_hline(yintercept=3, linetype="dashed", color = "white", linewidth=1)+theme(legend.position="top")
