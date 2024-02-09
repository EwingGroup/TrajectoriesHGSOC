#Alternative figure 2:

#Figure 2 - Recurrent CNV hotspots are associated with HRD and WGD


library(gplots)

## data for CNV figure heatmap
load(file='regionAmpDelSamp_20230814.RData');## array gistic region * sample
## denoting amp/del as 1 / -1;
## GISTIC regions named according to cytoband
load(file='cellularity_20230814.RData');     ## cellularity
load(file="WGD_HRD_20230814.RData");         ## WGD HRD
load(file="subtypes_20230814.RData");        ## subtypes 
load(file="Data_FG.RData");      ## derived from ConsensusDups_CnvkitClimatPurple_withX.bed
load(file="Data_HI.RData");      ## derived from ConsensusDels_CnvkitClimatPurple_withX.bed

## 'deepskyblue' == deletions
## 'darkseagreen' == amplifications
## colour-blind colours for subtypes
Verhaak_colours <- array(c("#E69F00", "#009E73", "#0072B2", "#CC79A7", "white"));
rownames(Verhaak_colours) <- c("IMR", "DIF", "PRO", "MES", "none");

## CNV figure heatmap has 5 parts (has not been assembled in R)
pdf(file='CNV_fig_D_heatmap_20230814.pdf',width=12,height=8);
h1 <- heatmap.2(regionAmpDelSamp,trace='none',col=c("#0072B2",'white',"#009E73"),
                cexRow=0.4,cexCol=0.055,key=F,labRow=F);
dev.off();

rwnms <- array(rev(rownames(regionAmpDelSamp)[h1$rowInd]));

pdf(file='CNV_fig_D_rowcount_20230814.pdf',width=8,height=4);
barplot(apply(rwnms,1,function(x){ c(sum(regionAmpDelSamp[x,]==-1,na.rm=T),
                                     sum(regionAmpDelSamp[x,]==1,na.rm=T))}),border=NA,
        #names=rev(rownames(regionAmpDelSamp)[h1$rowInd]),
        las=2,col=c("#0072B2","#009E73"),cex.names=0.4);
dev.off();

colnms <- array((colnames(regionAmpDelSamp)[h1$colInd]));

pdf(file='CNV_fig_A_colcount_20230814.pdf',width=12,height=4);
barplot(apply(colnms,1,function(x){ c(sum(regionAmpDelSamp[,x]==-1,na.rm=T),sum(regionAmpDelSamp[,x]==1,na.rm=T))}),border=NA,
        las=2,cex.names=0.05,col=c("#0072B2","#009E73"));
dev.off();


#GISTIC hotspots
library(maftools)


all.lesions <- "hgsoc_purple_x.all_lesions.conf_95.txt"
amp.genes <- "hgsoc_purple_x.amp_genes.conf_95.txt"
del.genes <- "hgsoc_purple_x.del_genes.conf_95.txt"
scores.gis <- "hgsoc_purple_x.scores.gistic"

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = FALSE)
pdf("Fig2A_FocalGainsLosses_inclX_recolour.pdf",width=10,height=4)
source("gisticChromPlot2.R")
gisticChromPlot2(gistic = laml.gistic, markBands = "all",color=c(Amp="#009E73",Del="#0072B2"),
                 ref.build="hg38", y.axis.label="GISTIC score of\nhotspot enrichment")
dev.off()

broad_results <- read.table("hgsoc_purple_x.broad_significance_results.txt",sep="\t",header=T)

broad_results$Del.frequency.Inv <- -broad_results$Del.frequency
broad_results_plot <- broad_results[,c("Arm","Amp.frequency","Amp.q.value","Del.frequency.Inv","Del.q.value")]
broad_results_amp <- cbind(dir="Amp",broad_results_plot[,c("Arm","Amp.frequency","Amp.q.value")])
broad_results_del <- cbind(dir="Del",broad_results_plot[,c("Arm","Del.frequency.Inv","Del.q.value")])
colnames(broad_results_amp)<-colnames(broad_results_del)<- c("Dir","Arm","Frequency","Qvalue")
long.br_res <- rbind(broad_results_amp,broad_results_del)

library(ggplot2)
library(ggbreak)
armlevels <- c(paste0(rep(1:22,each=2),rep(c("p","q"),22),sep=""),"Xp","Xq")
long.br_res$Arm <- factor(long.br_res$Arm,levels=armlevels)

pdf("SuppFig9_Broad_armlevel_gainslosses_gistic_barplot_withX.pdf",width=8,height=4)
 library(ggplot2)
 library(ggbreak)
 ggplot(long.br_res,aes(x=Arm,y=Frequency,fill=Dir)) + 
   geom_bar(stat="identity",position="identity",col="black")+
   theme_classic()+
  # scale_y_reverse()+
   scale_y_break(breaks=c(0,0),expand=FALSE,space=0.3)+#ticklabels = c(0.00,0.25,0.5,0.75,1.00)) + 
   ylim(c(-1,1))+theme(legend.position="none",
                       axis.title.y=element_text(size=18),
                       axis.title.x=element_text(size=18),
                       axis.text.x=element_text(size=13),
                       axis.text.y=element_text(size=13))+
   scale_fill_manual(values=c("#009E73","#0072B2"))+
   geom_text( aes(x = Arm, y= Frequency),
              label=ifelse(long.br_res$Qvalue<0.05,"*","") ,
              position = position_nudge(y = ifelse(long.br_res$Dir=="Amp",0.01,-0.04)),
              size=8)+
   annotate("text",x="Xq",y=0.95,label= "q-value <0.05 (*)",vjust=1,hjust=1,fontface=2,cex=8)
dev.off()


load(file="CNV_plot_data.RData");
load(file="SV_plot_data.RData");

library(ggplot2);

## deletions:      'deepskyblue'
## amplifications: 'darkseagreen'

df_amp<- data.frame(Gene=rownames(CNV_plot_data)[which(CNV_plot_data[,3]==1)],Direction="Amp",NumSamps=CNV_plot_data[CNV_plot_data[,3]==1,1],FC=CNV_plot_data[CNV_plot_data[,3]==1,4])
df_del<- data.frame(Gene=rownames(CNV_plot_data)[which(CNV_plot_data[,3]==0)],Direction="Del",NumSamps=CNV_plot_data[CNV_plot_data[,3]==0,2],FC=-CNV_plot_data[CNV_plot_data[,3]==0,4])
df_sv <- data.frame(Gene=SV_plot_data$gene_short_name, Direction=SV_plot_data$SV_event_type,NumSamps=SV_plot_data$percentageSV,FC=SV_plot_data$log2FoldChange)
df<- rbind(df_amp,df_del,df_sv)

lgd1_cols = c("Del" = "#0072B2", "Amp" = "#009E73","SV_del" = "#56B4E9", "SV_dup" = "khaki4", "SV_inv" = "#CC79A7","SV_breakpoint"="#E69F00")
df$Direction <- factor(df$Direction,levels=c("Del","Amp","SV_del","SV_dup","SV_inv","SV_breakpoint"))
pdf("Fig2D_ExpressionFC_athotspots.pdf",width=12,height=6)
ggplot(df, aes(x=NumSamps,y=FC))+geom_point(aes(col=Direction))+theme_classic(base_size=19)+
  xlim(c(0,30))+ylim(c(-1.2,1.2))+
  geom_text(aes(label=Gene,col=Direction),nudge_x=-0.3,size=6,hjust=1,nudge_y=0.05,check_overlap=T,fontface=4,show.legend=F)+
  annotate("text",x=18,y=-0.95,label="WRN",size=6,col="#0072B2",fontface=4)+
  annotate("text",x=25,y=-0.65,label="LEPROTL1",size=6,col="#0072B2",fontface=4)+
  scale_colour_manual(values=lgd1_cols,name="Event type",
                      labels=c("Del"="CNA deletion","Amp"="CNA duplication","SV_del"="SV deletion","SV_dup"="SV duplication",
                               "SV_inv"="Inversion","SV_breakpoint"="Breakpoint disruption"))+
  geom_hline(yintercept = 0,linetype="dashed",col="darkgrey")+
  ylab("Expression fold change with and\nwithout a hotspot event (log2)")+
  xlab("Percentage of samples with hotspot event")+
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=19))+
  guides(color = guide_legend(override.aes = list(size = 5))) 


dev.off()






