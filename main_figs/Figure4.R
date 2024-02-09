### HGSOC - Figure 4 ###

# Load data

cSvs_bysample <- read.table("SummaryTableSampleCallOfcSV.tsv",header=T,sep="\t",row.names=1)

#subset for heatmap

cSvs_bysample_subs <- cSvs_bysample[,c("chromoplexy","HRDetectStat","WholeGenomeDuplication","tyfonas","SeismicAmplification","AA_ecDNA",
                                       "bfb","SS_Chromothripsis","rigma","pyrgo")]

cSvs_bysample_subs[cSvs_bysample_subs=="Present"]<- 1
cSvs_bysample_subs[cSvs_bysample_subs=="Absent"]<- 0

cSvs_bysample_subs[cSvs_bysample_subs=="High_Confidence"]<- 1
cSvs_bysample_subs[cSvs_bysample_subs=="Low_Confidence"]<- 1
cSvs_bysample_subs_num <- apply(cSvs_bysample_subs,2,as.numeric)
cSvs_bysample_subs_num <- as.matrix(cSvs_bysample_subs_num)


formatted_cSVnames <- c("Chromoplexy","HRD","WGD","Tyfonas","Seismic Amplification","ecDNA","Breakage Fusion Bridges",
                        "Chromothripsis","Rigma","Pyrgo")
colnames(cSvs_bysample_subs_num)<- formatted_cSVnames
library(ComplexHeatmap)
#Panel A
pdf("Fig4A.pdf",width=15,height=8)
Heatmap(t(cSvs_bysample_subs_num), col = structure(c("grey","#440154"), names = c("0","1")),
        rect_gp = gpar(col = "white", lwd = 1),
        show_column_dend = F,
        row_dend_width = unit(40,"mm"),
        column_title="Samples",column_title_side = "bottom",column_title_gp = gpar(fontsize = 24,colour="black"),
        row_names_side = "left",row_names_max_width = unit(11, "cm"),
        row_names_gp = gpar(fontsize = 26,colour="black"),
        heatmap_legend_param = list( labels=c("Absent","Present"),title="",
                                grid_width=unit(10,"mm"),grid_height=unit(10,"mm"),
                                    labels_gp=gpar(cex=2.2,colour="black"),legend.direction = "horizontal")
)

dev.off()

#Check ecDNA enrichment
chisq.twomat <- matrix(c(46,34,(56+23+106+41),18),nrow=2,byrow=T)
chisq.twomat
#     [,1] [,2]
#[1,]   34   46
#[2,]   18  226
chisq.test(chisq.twomat)
or <- (34/46)/(18/226)

#Load fraction of genome doubled data
mcn_fract <- read.table("mcn_fraction.txt")
colnames(mcn_fract)<- c("Sample","Fraction_genome_doubled")

add_wgd<- merge(cSvs_bysample,mcn_fract,by="Sample")
add_wgd$SS_Chromothripsis[add_wgd$SS_Chromothripsis=="High_Confidence" | add_wgd$SS_Chromothripsis=="Low_Confidence"] <- "Present"

col <- c(
  "AOCS" = "blue",
  "BCCA" = "black",
  "SHGSOC" = "brown",
  "TCGA" = "grey",
  "MDA" = "orange"
)
#Panel C
png("Figure4_C.png",width=7,height=6, unit="in",res=300)

add_wgd$SS_Chromothripsis <- factor(add_wgd$SS_Chromothripsis,levels=c("Present","Absent"))
add_wgd$bfb <- factor(add_wgd$bfb,levels=c("Present","Absent"))
df_plot <- add_wgd[,c("Sample","SS_Chromothripsis","bfb","Fraction_genome_doubled")]

library(reshape2)
df_plot_long<- melt(df_plot,id.var=c("Sample","Fraction_genome_doubled"))
df_plot_long$combo<- paste(df_plot_long$variable,df_plot_long$value,sep="")
df_plot_long$combo <- factor(df_plot_long$combo,levels=c("bfbPresent","bfbAbsent","SS_ChromothripsisPresent",
                                                         "SS_ChromothripsisAbsent"))

pres_cols <- c("Present"="#440154","Absent"="grey")
png("Figure4_EF_combo.png",width=8,height=5, unit="in",res=300)
ggplot(df_plot_long,aes(x=combo,y=Fraction_genome_doubled))+geom_violin()+
  geom_jitter(aes(col=value),width=0.1)+theme_classic()+
  ylab("Fraction of genome duplicated")+
  scale_x_discrete(labels=c("bfbAbsent"="No breakage\nfusion bridges","bfbPresent"="Breakage\nfusion bridges",
                            "SS_ChromothripsisAbsent"="No chromothripsis","SS_ChromothripsisPresent"="Chromothripsis"))+
  scale_color_manual(values=pres_cols)+
  theme(legend.position = "none",
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=16,colour="black"),
        axis.title.x=element_blank())
dev.off()  



#Explore ways of visualising heatmap
corrmat<- cSvs_bysample_subs_num
totals<-rowSums(corrmat)
corrmat<- cbind(corrmat,totals)
rownames(corrmat)<-cSvs_bysample$Sample
corrmat<- corrmat[corrmat[,"totals"]>2,]
corrmat_ord <- corrmat[order(corrmat[,"totals"],decreasing=T),]
corrmat_ord<- data.frame(corrmat_ord[,setdiff(colnames(corrmat_ord),"totals")],Sample=rownames(corrmat_ord))
library(reshape2)

corrmat_long<- melt(corrmat_ord,id.vars="Sample")
corrmat_long$Sample<- factor(corrmat_long$Sample,levels=corrmat_ord$Sample)


#Try Latent class analysis

colnames(corrmat)[5]<-"Seismic_Amp"
colnames(corrmat)[7]<-"Breakage_Fusion_Bridges"

#Try pca again
pca<- prcomp(corrmat[,1:10])
pca$sdev/sum(pca$sdev)

corrmat_nohrdwgd<- corrmat[,c(1,4:10)]
pca2<- prcomp(corrmat_nohrdwgd)
library(reshape2)
library(RColorBrewer)
library(viridis)
reduced_heatmap<- pca2$rotation[,1:5]
heatmap_long<-melt(reduced_heatmap)
heatmap_long$Var2<- factor(heatmap_long$Var2,levels=c("PC1","PC2","PC3","PC4","PC5"))
#heatmap_long$Var1<- factor(heatmap_long$Var1, levels=rev(c("HRD","Chromoplexy","WGD","Chromothripsis","Pyrgo","Rigma","Breakage_Fusion_Bridges",
                                                       "ecDNA","Seismic Amplification","Tyfonas")))

newheatmap<-heatmap_long

hrd_row<- data.frame(Var1=rep("HRD",5),Var2=c("PC1","PC2","PC3","PC4","PC5"),value=0)
wgd_row<- data.frame(Var1=rep("WGD",5),Var2=c("PC1","PC2","PC3","PC4","PC5"),value=0)
newheatmap<- rbind(heatmap_long,hrd_row,wgd_row)
newheatmap$Var1<- factor(newheatmap$Var1, levels=rev(c("HRD","Chromoplexy","WGD","Chromothripsis","Pyrgo","Rigma","Breakage_Fusion_Bridges",
                                                       "ecDNA","Seismic_Amp","Tyfonas")))





png("Reduced_complexSV_profiles.png",height=8,width=10,unit="in",res=300)
ggplot(newheatmap, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white",lwd = 1.5,linetype = 1) +
  scale_fill_viridis(name="Event\ncontribution")+
  xlab("Complex SV signatures\n(principal components)")+
  ylab("")+theme_classic(base_size=24)+
  #scale_y_discrete(labels=rev(c("HRD","Chromoplexy","WGD","Chromothripsis","Pyrgo","Rigma","Breakage\nFusion Bridges",
 # "ecDNA","Seismic\nAmplification","Tyfonas")))+
  geom_rect(xmin =0.52,
            xmax =5.48,
            ymin=7.52,
            ymax=8.48 , fill = "lightgrey")+
  geom_rect(xmin =0.52,
            xmax =5.48,
            ymin=9.52,
            ymax=11.48 , fill = "lightgrey")
dev.off()
corrmat_withpcs<- data.frame(corrmat,PC1=pca2$x[,1],PC2=pca2$x[,2])
wgd_hrd<-paste(corrmat_withpcs$WGD,corrmat_withpcs$HRD)  
corrmat_withpcs<- cbind(corrmat_withpcs, Cluster=wgd_hrd)

wgd_hrd_cols <- c("1 1" = "#0072B2",
          "1 0" = "#009E73",
          "0 1" = "#CC79A7",
          "0 0" = "#E69F00")

png("PCs_complexSVstatus.png",width=7,height=6,unit="in",res=300)
ggplot(corrmat_withpcs,aes(x=PC1, y=PC2,col=Cluster))+geom_point(size=3)+
  #ylim(c(-1,1.3))+xlim(c(-1,1.5))+
  theme_classic(base_size=18)+scale_colour_manual(name="WGD/HRD status",values=wgd_hrd_cols)+
  annotate("text",x=-0.8,y=1.1,label="no WGD,\nHRD",size=6,col="#CC79A7",fontface=4)+
  annotate("text",x=-0.8,y=-1,label="WGD, HRD",size=6,col="#0072B2",fontface=4)+
  annotate("text",x=1.2,y=1.1,label="no WGD,\nno HRD",size=6,col="#E69F00",fontface=4)+
  annotate("text",x=1.2,y=-1,label="WGD, no HRD",size=6,col="#009E73",fontface=4)+
  theme(legend.position = "none")
dev.off()

#Boxplots
forboxes <- data.frame(Sample=rownames(corrmat_withpcs),corrmat_withpcs[,c("PC1","PC2","WGD","HRD","Cluster")])
boxes_long<-melt(forboxes,id.vars=c("Sample","WGD","HRD","Cluster"))
ggplot(boxes_long,aes(x=as.factor(variable),y=value,col=Cluster))+facet_grid(~Cluster)+
  geom_boxplot()+scale_colour_manual(name="WGD/HRD status",values=wgd_hrd_cols)+theme_classic()

g1<-ggplot(corrmat_withpcs,aes(x=as.factor(HRD),y=PC1))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+theme_classic(base_size=18)+
  scale_x_discrete(labels=c("No HRD","HRD"))+theme(axis.title.x=element_blank())
g2<-ggplot(corrmat_withpcs,aes(x=as.factor(WGD),y=PC2))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+theme_classic(base_size=18)+
  scale_x_discrete(labels=c("No WGD","WGD"))+theme(axis.title.x=element_blank())

png("States_bycomplexSVPCs.png",width=8,height=5,unit="in",res=300)
g1+g2
dev.off()
