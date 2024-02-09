#copy number sigs

library(ComplexHeatmap)
library(reshape2)
library(RColorBrewer)
library(viridis)

sig_exposures<-read.table("drews_exposure_per_sample.txt",header=T)
sig_exposures<- reshape(sig_exposures, idvar = c("cohort","sample"), timevar = "CN_sig", direction = "wide")
sig_exposures_ordCX03<- sig_exposures[order(sig_exposures$activity.CX03,decreasing = T),]

sig_mat<- as.matrix(sig_exposures_ordCX03[,3:19])
rownames(sig_mat)<- sig_exposures_ordCX03$sample
colnames(sig_mat) <- gsub("activity.","",colnames(sig_exposures_ordCX03[,3:19]))

png("CNV_sigs_exposures_heatmap.png",width=12,height=6,res=300,unit="in")
  Heatmap(t(sig_mat),col=viridis(100,option = "D"),
          rect_gp = gpar(col = "white", lwd = 0.1),
        show_column_dend = F,
        show_row_dend = F,
        cluster_rows =F,cluster_columns=F,
        column_title="Samples",column_title_side = "bottom",column_title_gp = gpar(fontsize = 13),
        row_names_side = "left",row_names_max_width = unit(11, "cm"),
        row_names_gp = gpar(fontsize = 13),show_column_names = F,
        heatmap_legend_param = list(title="Relative frequency\nof CN signatures"))
dev.off()

library(reshape2)
sig_exposures_ordCX03<- cbind(Sample=sig_exposures_ordCX03$sample,sig_exposures_ordCX03[,3:19])
expo_long <- melt(sig_exposures_ordCX03,id.vars = "Sample")
expo_long$Sample<-factor(expo_long$Sample,levels=sig_exposures_ordCX03$Sample)

library(ggplot2)
sig_cols <- c("activity.CX01"="#009E73",
              "activity.CX03" = "#CC79A7",
              "Other" = "grey")
colour_cat <- as.character(expo_long$variable)
colour_cat[which((colour_cat %in% c("activity.CX01","activity.CX03"))==FALSE)] <- "Other"
expo_long<-cbind(expo_long,colour_cat)

expo_long$colour_cat<- factor(expo_long$colour_cat,levels=c("Other","activity.CX01","activity.CX03"))

png("CNV_sigs_exposures_stackbar.png",width=12,height=4,res=300,unit="in")
ggplot(expo_long, aes(x=Sample,y=value,fill=colour_cat))+
  geom_bar(position="stack",stat="identity")+
  scale_fill_manual(name="CN signature",values=sig_cols,labels=c("CX01 - chromosome missegregation","CX03 - impaired HR","Other"))+
  theme_classic(base_size=13)+theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  ylab("Relative frequency of CN signatures")
dev.off()

#Correlation scatter
png("CNVsigs_scatter.png",width=8,height=8,unit="in",res=300)
ggplot(sig_exposures_ordCX03,aes(x=activity.CX03,y=activity.CX01))+geom_point()+
  theme_classic(base_size = 24)+ylim(c(0,1))+xlim(c(0,1))+ylab("CX01")+xlab("CX03")
dev.off()

#Load HRD/WGD status

hrd <- read.table("data/HRDetect_results_BRCAstatus.txt",sep="\t",header=T)
hrd<- hrd[,c("Sample","HRDetect","BRCA_SNV_status")]
brca_mut <- numeric()
brca_mut[hrd$BRCA_SNV_status=="None"]<- 0
brca_mut[hrd$BRCA_SNV_status!="None"]<- 1
hrd<- cbind(hrd,brca_mut)
hrd$Sample<-gsub("-P$","",hrd$Sample)
cn_withhrd<- merge(sig_exposures_ordCX03,hrd,by="Sample")

#Load WGD
mcn_fract <- read.table("mcn_fraction.txt")
colnames(mcn_fract)<- c("Sample","Fraction_genome_doubled")
wgd <- ifelse(mcn_fract$Fraction_genome_doubled>0.75,1,0)
addwgd<-cbind(mcn_fract,wgd)
cn_with_hrd_wgd <- merge(cn_withhrd,addwgd,by="Sample")
hrd <- ifelse(cn_with_hrd_wgd$HRDetect>0.7,1,0)
cn_with_hrd_wgd_full <- cbind(cn_with_hrd_wgd,hrd) 


cx1<-ggplot(cn_with_hrd_wgd_full,aes(x=as.factor(hrd),y=activity.CX01))+geom_boxplot(outlier.shape=NA)+theme_classic()+geom_jitter(width=0.1)+
  ylab("CX01")+theme(axis.title.x=element_blank())+scale_x_discrete(labels=c("No BRCA1/2 mutation","BRCA1/2 mutation"))

cx3<-ggplot(cn_with_hrd_wgd_full,aes(x=as.factor(hrd),y=activity.CX03))+geom_boxplot(outlier.shape=NA)+theme_classic()+geom_jitter(width=0.1)+
  ylab("CX03")+theme(axis.title.x=element_blank())+scale_x_discrete(labels=c("No BRCA1/2 mutation","BRCA1/2 mutation"))

cx1_wgd<-ggplot(cn_with_hrd_wgd_full,aes(x=as.factor(wgd),y=activity.CX01))+geom_boxplot(outlier.shape=NA)+theme_classic()+geom_jitter(width=0.1)+
  ylab("CX01")+theme(axis.title.x=element_blank())+scale_x_discrete(labels=c("No WGD","WGD"))

cx3_wgd<-ggplot(cn_with_hrd_wgd_full,aes(x=as.factor(wgd),y=activity.CX03))+geom_boxplot(outlier.shape=NA)+theme_classic()+geom_jitter(width=0.1)+
  ylab("CX03")+theme(axis.title.x=element_blank())+scale_x_discrete(labels=c("No WGD","WGD"))

                                                                                                                      
library(patchwork)
png("CNsigs_boxplots.png",width=8,height=4,res=300, unit="in")
(cx1+cx3)/(cx1_wgd+cx3_wgd)
dev.off()
