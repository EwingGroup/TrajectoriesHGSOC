# Supp Figs - complex 

#Load dups at CCNE1
dups <- read.table("CCNE1_dups.bed",sep="\t",stringsAsFactors = F)

dups$V10 <- gsub("_PrimaryTumour","",dups$V10)
dups$V10 <- gsub("T","",dups$V10)
dups[substr(dups$V10,1,2)=="WG","V10"]<- sapply(strsplit(dups[substr(dups$V10,1,2)=="WG","V10"],"-"), function(x) paste(x[1],x[2],x[3],sep="-"))

#Load complex at CCNE1
complex <- read.table("complexSV_CCNE1.bed",sep="\t")
complex[substr(complex$V4,1,2)=="WG","V4"]<- sapply(strsplit(complex[substr(complex$V4,1,2)=="WG","V4"],"-"), function(x) paste(x[1],x[2],x[3],sep="-"))

#Load cn at CCNE1
purple <- read.table("purple_cn_CCNE1.bed")
purple$V5 <- gsub("T","",purple$V5)
purple$V5 <- gsub("_2","",purple$V5)
purple[substr(purple$V5,1,2)=="WG","V5"]<- sapply(strsplit(purple[substr(purple$V5,1,2)=="WG","V5"],"-"), function(x) paste(x[1],x[2],x[3],sep="-"))

#IDs
  usable <-fread("data/Combined_usable_ids_marker.txt", header = F) %>%
  mutate(V1 = gsub("_2$|-P$", "", V1))
usable<-data.frame(usable)
usable<-as.character(usable[,1])
indicator_mat<-mat.or.vec(nr=324,nc=3)
rownames(indicator_mat)<-usable
colnames(indicator_mat)<-c("Simple_dup","Complex_dup","CN")

uni_simple_dups<- unique(dups$V10)
indicator_mat[uni_simple_dups,"Simple_dup"]<-1

for (i in 1:dim(complex)[1]){
  
  if (indicator_mat[complex[i,"V4"],"Complex_dup"] != 0){
      indicator_mat[complex[i,"V4"],"Complex_dup"]<- paste(indicator_mat[complex[i,"V4"],"Complex_dup"],complex[i,"V5"],sep=";")
  }else{
    indicator_mat[complex[i,"V4"],"Complex_dup"]<- complex[i,"V5"]
  }
}
for (i in 1:dim(purple)[1]){
  purple_sample<- purple[i,"V5"]
  purple_cn <- round(mean(purple[purple$V5==purple_sample,"V4"]),0)
  indicator_mat[purple_sample,"CN"]<- purple_cn
}

indicator_df<- data.frame(indicator_mat)
Comb_event<- indicator_df$Simple_dup
Comb_event[indicator_df$Complex_dup!=0]<- indicator_df[indicator_df$Complex_dup!=0,"Complex_dup"]

getUniqCompl<- function(x){
  uniq_complex<-paste(sort(unique(x)),collapse=";")
  return(uniq_complex)
}
Comb_event_uni<-sapply(strsplit(Comb_event,";"),getUniqCompl)

indicator_df<-cbind(indicator_df,Comb_event_uni)
indicator_df$CN<-as.numeric(indicator_df$CN)

indicator_df$Comb_event_uni<- as.character(indicator_df$Comb_event_uni)
indicator_df[indicator_df$Comb_event_uni=="AAecDNA;gGecDNA","Comb_event_uni"]<-"AAecDNA"
indicator_df[indicator_df$Comb_event_uni=="AAecDNA;Tyfonas","Comb_event_uni"]<-"AAecDNA"
indicator_df[indicator_df$Comb_event_uni=="bfb;Chromoplexy","Comb_event_uni"]<-"bfb"
indicator_df[indicator_df$Comb_event_uni=="bfb;HighConfidenceSSChromothripsis","Comb_event_uni"]<-"bfb;Chromothripsis"
indicator_df[indicator_df$Comb_event_uni=="bfb;LowConfidenceSSChromothripsis","Comb_event_uni"]<-"bfb;Chromothripsis"
indicator_df[indicator_df$Comb_event_uni=="HighConfidenceSSChromothripsis","Comb_event_uni"]<-"Chromothripsis"
indicator_df[indicator_df$Comb_event_uni=="HighConfidenceSSChromothripsis;Tyfonas","Comb_event_uni"]<-"Chromothripsis"
indicator_df[indicator_df$Comb_event_uni=="LowConfidenceSSChromothripsis","Comb_event_uni"]<-"Chromothripsis"

indicator_keep <- c("1","AAecDNA","AAecDNA;bfb","bfb","bfb;Chromothripsis","Chromothripsis")
indicator_df2 <- indicator_df[indicator_df$Comb_event_uni %in% indicator_keep,]

indicator_df2$Comb_event_uni<-factor(indicator_df2$Comb_event_uni,levels=c("1","AAecDNA","bfb","Chromothripsis","AAecDNA;bfb","bfb;Chromothripsis"))

pdf("individual_panels/supp_figs/SuppFig7B_CCNE1_CN_byevent.png",width=8,height=4)
ggplot(indicator_df2,aes(x=Comb_event_uni,y=CN))+geom_violin()+
  geom_jitter(width=0.1)+theme(axis.text.x=element_text(angle=60,hjust=1,colour="black",size=12))+theme_classic(base_size=14)+
  scale_x_discrete(labels=c("1"="Simple duplication","AAecDNA"="ecDNA","bfb"="Breakage\nfusion bridges",
                            "Chromothripsis"="Chromothripsis","AAecDNA;bfb"="ecDNA and\nbreakage\nfusion bridges",
                            "bfb;Chromothripsis"="Breakage\nfusion bridges\nand chromothripsis"))+
  ylab("CCNE1 copy number") + xlab("Duplication event at CCNE1")
dev.off()
library(reshape2)
indicator_df2$Sample<- row.names(indicator_df2)
indicator_df2_forplot<- indicator_df2[,c("Comb_event_uni","Sample")]
indicator_df2_forplot$Comb_event_uni<- as.character(indicator_df2_forplot$Comb_event_uni)
indicator_df2_forplot[indicator_df2_forplot[,1]=="1",1]<- "Simple_dup"
write.table(indicator_df2_forplot,file="CCNE1_dup_mechanism.txt",sep="\t",row.names=F,quote=F)

#CN by frequency
simple_or_complex<- as.character(indicator_df2$Comb_event_uni)
simple_or_complex[simple_or_complex!="1" ]<-"Complex"
indicator_df3<- cbind(indicator_df2,simple_or_complex)
table(indicator_df3$CN,indicator_df3$simple_or_complex)                                                      
pdf("SuppFig7A_CCNE_cn_vs_density.pdf",width=10,height=6)
ggplot(indicator_df3,aes(x=CN,fill=simple_or_complex))+geom_density()+theme_classic(base_size = 16)+
  scale_fill_manual(name="Method of duplication",labels=c("Simple","Complex"),values= c( "#56B4E9","#E69F00"))+
  theme(axis.title.y=element_text(size=18,colour="black"),
        axis.title.x=element_text(size=18,colour="black"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"),
        legend.title=element_text(size=18,colour="black"),
        legend.text=element_text(size=16,colour="black"))+
  ylab("Sample density")+xlab("CCNE1 copy number")
dev.off()
pdf("SuppFig7A_extra_CCNE_cn_vs_density_boxplot.pdf",width=4,height=4)
ggplot(indicator_df3,aes(x=simple_or_complex, y=CN,fill=simple_or_complex))+
  geom_boxplot(outlier.shape=NA,width=0.5)+
  geom_jitter(width=0.1)+theme_classic(base_size=16)+
  scale_fill_manual(name="Method of duplication",labels=c("Simple","Complex"),values= c( "#56B4E9","#E69F00"))+
  scale_x_discrete(labels=c("Simple","Complex"))+xlab("Method of duplication")+ylab("CCNE1 copy number")+
  theme(legend.position="none",
        axis.title.y=element_text(size=18,colour="black"),
        axis.title.x=element_text(size=18,colour="black"),
        axis.text.x=element_text(size=16,colour="black"),
        axis.text.y=element_text(size=16,colour="black"))
dev.off()                    
