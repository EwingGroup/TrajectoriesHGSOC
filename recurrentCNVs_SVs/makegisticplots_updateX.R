library(maftools)

#Data is here: /exports/igmm/eddie/HGS-OvarianCancerA-SGP-WGS/recurrentCNVs/gistic/HGSOC_purple_x

all.lesions <- "~/Documents/HGSOC_marker_paper/ConsensusCNVs/HGSOC_purple_x/hgsoc_purple_x.all_lesions.conf_95.txt"
amp.genes <- "~/Documents/HGSOC_marker_paper/ConsensusCNVs/HGSOC_purple_x/hgsoc_purple_x.amp_genes.conf_95.txt" 
del.genes <- "~/Documents/HGSOC_marker_paper/ConsensusCNVs/HGSOC_purple_x/hgsoc_purple_x.del_genes.conf_95.txt"
scores.gis <- "~/Documents/HGSOC_marker_paper/ConsensusCNVs/HGSOC_purple_x/hgsoc_purple_x.scores.gistic"

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis, isTCGA = FALSE)
#png("FocalGainsLosses_inclX.png",width=13,height=5, unit="in",res=300)
source("ConsensusCNVs/scripts/gisticChromPlot2.R")
gisticChromPlot2(gistic = laml.gistic, markBands = "all",color=c(Amp="tomato1",Del="steelblue3"),
                 ref.build="hg38")
#dev.off()


newrownames <- sapply(strsplit(rownames(laml.gistic@numericMatrix),":"), function(x) x[2])
rownames(laml.gistic@numericMatrix)<-newrownames
rownames(laml.gistic@cnMatrix)<-newrownames

amp_cytobands <- sapply(laml.gistic@cytoband.summary[laml.gistic@cytoband.summary$Variant_Classification=="Amp",'Cytoband'],as.character)
del_cytobands <- sapply(laml.gistic@cytoband.summary[laml.gistic@cytoband.summary$Variant_Classification=="Del",'Cytoband'],as.character)
gisticOncoPlot(gistic = laml.gistic,
               showTumorSampleBarcodes=FALSE,
               bands=amp_cytobands,
               colors=c(Amp="tomato1",Del="steelblue3"))
gisticOncoPlot(gistic = laml.gistic,
               showTumorSampleBarcodes=FALSE,
               bands=del_cytobands,
               colors=c(Amp="tomato1",Del="steelblue3"))

png("FocalGainsLosses_oncoplot.png",width=8,height=6, unit="in",res=300)
gisticOncoPlot(gistic = laml.gistic,
               showTumorSampleBarcodes=FALSE,
               colors=c(Amp="tomato1",Del="steelblue3"))
dev.off()


#Broad events - chromosome arm oncoplot/heatmap - needs rerun to include x if you want it but I think it looks messy.
broad_dat <- read.table("ConsensusCNVs/HGSOC_genes_purple/hgsoc_genes_purple.broad_values_by_arm.txt",sep="\t",header=T)

ordbroad <- hclust( dist(t(broad_dat[,-1]) ))$order
ordbroad

broad_dat<- cbind(broad_dat[,1],broad_dat[,ordbroad])
newcolnames<- gsub("_PrimaryTumour","",colnames(broad_dat))
newcolnames2<- gsub(".purple","",newcolnames)
newcolnames3<- gsub("T","",newcolnames2)
colnames(broad_dat)<- newcolnames3

chr.arms <- paste(rep(1:22,each=2) ,c("p","q"),recycle0 = TRUE,sep="")

library(reshape2)
broad_melt <- melt(broad_dat)
new_value <- rep("",length(broad_melt$value))
new_value[broad_melt$value< -0.1] <- "Arm loss"
new_value[broad_melt$value> 0.1] <- "Arm gain"
broad_melt <- cbind(broad_melt,new_value)
broad_melt$variable <- factor(broad_melt$variable, levels=colnames(broad_dat)[-1])
broad_melt$Chromosome.Arm <- factor(broad_melt$Chromosome.Arm, levels=rev(chr.arms))


library(ggplot2)
png("ChromosomeArmGainsLosses.png",width=10,height=5, unit="in",res=300)
ggp <- ggplot(broad_melt, aes(variable, Chromosome.Arm)) +                     
  geom_tile(aes(fill = as.factor(new_value)))+ 
  scale_fill_manual(values=c("white","tomato1","steelblue3"))+xlab("Sample")+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size=6), legend.title=element_blank())

ggp  
dev.off()
