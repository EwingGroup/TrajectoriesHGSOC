# Stacked SV densities

#Load SV data




#Make x_axis
chrom_lengths <- read.table("data/chrom_lengths.txt",sep="\t")
chrom_lengths<- chrom_lengths[,c(1,2)]
chrom_lengths[,2]<-as.numeric(gsub(",","",chrom_lengths[,2]))

total_lengths <- dim(chrom_lengths)[1]
total_lengths[1] <- chrom_lengths[1,2]

for (i in 2:dim(chrom_lengths)[1]){
  total_lengths[i]<- total_lengths[i-1]+chrom_lengths[i,2]
}
chrom_lengths<- cbind(chrom_lengths,total_lengths)
chrom_lengths[chrom_lengths[,1]=="X",1]<- 23
chrom_lengths[chrom_lengths[,1]=="Y",1]<- 24
chrom_lengths[,1] <- as.numeric(chrom_lengths[,1])

#load svs
svs_all <- read.table("data/svFormatted.bedpe",sep="\t")
svs_all <- svs_all[svs_all$V11!="INS",]
svs1<- svs_all[,c(1:3,7:14)]
svs2 <- svs_all[,c(4:14)]
colnames(svs1)<-colnames(svs2)<- colnames(svs_all)[4:14]
svs<- rbind(svs1,svs2)

svs<- svs[svs[,1] %in% 
            c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"),]
svs[svs[,1]=="X",1]<-"23"
svs[,1]<- as.numeric(svs[,1])

mid_pos<-numeric(dim(svs)[1])
for (i in 1:dim(svs)[1]){
  mid_pos[i] <- mean(c(svs[i,2],svs[i,3]))
}

total_pos<- numeric(dim(svs)[1])
for (i in 1:dim(svs)[1]){
  
  if (svs[i,1]==1){
    total_pos[i]<-mid_pos[i]
  }else{
  total_pos[i] <- mid_pos[i] + chrom_lengths[svs[i,1]-1,"total_lengths"]
  }
}
svs<- cbind(svs,mid_pos,total_pos)

chr_breaks<- c(0,chrom_lengths[,3])
chr_labels<-c(paste0(1:22),"X","")

mid_breaks<-numeric(23)
for (i in 1:23){
  mid_breaks[i] <- (chr_breaks[i]+chr_breaks[i+1])/2
}

all_marks <- sort(c(chr_breaks[1:24],mid_breaks))
chr_labels <-c("","1","","2","","3","","4","","5","","6","","7","","8","","9","","10",
               "","11","","12","","13","","14","","15","","16","","17","","18","","19","","20",
               "","21","","22","","X","")

library(ggplot2)
library(ggridges)
#Plot densities


lgd2_cols =  c("DUP" = "#009E73", "DEL" = "#0072B2", "INV" = "#CC79A7", "CTX" = "#E69F00", "All" = "#242423")

facetlabs <- as_labeller(c(
  "DUP"="Duplication",
  "DEL"="Deletion",
  "INV"="Inversion",
  "CTX"="Translocation",
  "All"="All"
))
ccne1<- c(2654411288+29811991, 2654411288+29824312)
chr19<- c(2654411288,2713028904)
svs2 <- rbind(svs,within(svs,V11 <- "All"))
svs2$V11 <- factor(svs2$V11,levels=c("All","DEL","DUP","INV","CTX"))
png("SV_density_acrossGenome.png",width=12,height=5,unit="in",res=300)
ggplot(svs2,aes(total_pos,fill=V11))+geom_density(adjust=0.1)+theme_classic(base_size =13)+
  facet_grid(V11~.,labeller=facetlabs)+
  scale_x_continuous(breaks=all_marks,labels=chr_labels,limits=c(0,all_marks[47]))+
  xlab("Chromosome")+ylab("SV breakpoint density")+
  theme(legend.position="none",
        strip.text.y = element_text(size=13, face="bold",angle=0,hjust=0),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  scale_fill_manual(labels=c("Duplication","Deletion","Inversion","Translocation","All"),values=lgd2_cols)
  
  dev.off()
  
png("CCNE1_SV.png",width=4,height=5,unit="in",res=300)
  ggplot(svs2,aes(total_pos,fill=V11))+geom_density(adjust=0.2)+theme_classic(base_size=13)+
    facet_grid(V11~.,labeller=facetlabs)+
    scale_x_continuous(breaks=all_marks,labels=chr_labels,limits=chr19)+
    xlab("Chromosome")+ylab("")+
    theme(legend.position="none",
          strip.text.y = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_fill_manual(labels=c("Duplication","Deletion","Inversion","Translocation","All"),values=lgd2_cols)+
    geom_vline(xintercept=ccne1[1])+
    geom_vline(xintercept=ccne1[2])
  #geom_segment(aes(x=0,xend=all_marks[47], y=1/all_marks[47],yend=1/all_marks[47]))
  dev.off()
  
#Complex
complex <- read.table("data/ComplexStructuralVariantRegions.csv",sep=",",header=T,row.names=1)
complex[1,]  
  
complex[complex[,2]=="X",2]<-"23"
complex[,2]<- as.numeric(complex[,2])

mid_pos2<-numeric(dim(complex)[1])
for (i in 1:dim(complex)[1]){
  mid_pos2[i] <- mean(c(complex[i,3],complex[i,4]))
}

total_pos<- numeric(dim(complex)[1])
for (i in 1:dim(complex)[1]){
  
  if (complex[i,2]==1){
    total_pos[i]<-mid_pos2[i]
  }else{
    total_pos[i] <- mid_pos2[i] + chrom_lengths[complex[i,2]-1,"total_lengths"]
  }
}
complex<- cbind(complex,mid_pos2,total_pos)

category <- complex$Feature
category[category=="HighConfidenceSSChromothripsis"]<-"Chromothripsis"
category[category=="LowConfidenceSSChromothripsis"]<-"Chromothripsis"
complex_labels<- c("bfb","AAecDNA","Chromothripsis","Chromoplexy",
  "Rigma","Pyrgo")
complex<-cbind(complex,category)

complex <- complex[complex$category %in% complex_labels,]
complex2 <- rbind(complex,within(complex,category <- "All"))
complex2$category <- factor(complex2$category,levels=c("All",complex_labels))

complex_cols <- c("All"="#242423","bfb"="#CC79A7","AAecDNA"="#009E73","Chromothripsis"="#0072B2","Chromoplexy"="#E69F00",
                     "Rigma"="#56B4E9","Pyrgo"="#8bc34a")

complex_labels_nice<- as_labeller(c("All"="All","bfb"="Breakage fusion bridges",
         "AAecDNA"="ecDNA","Chromothripsis"="Chromothripsis","Chromoplexy"="Chromoplexy",
         "Rigma"="Rigma","Pyrgo"="Pyrgo"))

png("complexSV_density_acrossGenome.png",width=12,height=6.5,unit="in",res=300)
ggplot(complex2,aes(total_pos,fill=category))+geom_density(adjust=0.1)+theme_classic(base_size=16)+
  facet_grid(category~.,labeller=complex_labels_nice)+
  scale_fill_manual(values=complex_cols)+
  scale_x_continuous(breaks=all_marks,labels=chr_labels,limits=c(0,all_marks[47]))+
  xlab("Chromosome")+ylab("Complex SV breakpoint density")+
  theme(legend.position="none",
        strip.text.y = element_text(size=18, face="bold",angle=0,hjust=0),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


#Just chromosome 19
chr19<- c(2654411288,2713028904)
ccne1<- c(2654411288+29811991, 2654411288+29824312)
png("CCNE1_complexSV_density.png",width=3.75,height=6.5,unit="in",res=300)
ggplot(complex2,aes(total_pos,fill=category))+geom_density(adjust=0.2)+theme_classic(base_size=13)+
  facet_grid(category~.)+
  scale_x_continuous(breaks=all_marks,labels=chr_labels,limits=chr19)+
  xlab("Chromosome")+ylab("")+
scale_fill_manual(values=complex_cols)+
  theme(legend.position="none",
        strip.text.y = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_vline(xintercept=ccne1[1])+
  geom_vline(xintercept=ccne1[2])
dev.off()

# Test for chromosomal enrichment

#To get cohort specific results subset svs and complex2 by cohort

svs2[1,]
svs2$Cohort <- substr(svs2$V13,1,2)
aocs_svs2 <- svs2[svs2$Cohort=="AO",]

sv_counts<-table(svs2$V11,svs2$V4)
rownames(sv_counts)[1]<-"All_SV"

complex_counts<-table(complex2[,"category"],complex2[,2])
rownames(complex_counts)[1]<-"All_complex"

all_counts<- rbind(sv_counts,complex_counts)
genome <- rowSums(all_counts)
all_counts<- cbind(all_counts,genome)
calculate_expected_p <- chrom_lengths[,2]/3031042417

class_sv <- character(dim(all_counts)[1]*(dim(all_counts)[2]-1))
chr <- character(dim(all_counts)[1]*(dim(all_counts)[2]-1))
or <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))
lci <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))
uci <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))
p <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))

count<- 1
for (i in 1:dim(all_counts)[1]){
  for (j in 1:(dim(all_counts)[2]-1)){
    
      class_sv[count] <- rownames(all_counts)[i]
      chr[count] <- colnames(all_counts)[j]
      
      b<-binom.test(all_counts[i,j],all_counts[i,"genome"],p=calculate_expected_p[j])
      
      or[count] <- b$estimate/b$null.value
      lci[count] <- b$conf.int[1]/b$null.value
      uci[count] <- b$conf.int[2]/b$null.value
      p[count] <- b$p.value
      
      count<- count+1
  }
}
enrich_table<- data.frame(Class_SV=class_sv,Chromosome=chr,OR=or,LCI=lci,UCI=uci,Pvalue=p)
adjp <- p.adjust(enrich_table$Pvalue,"bonferroni")
enrich_table<- cbind(enrich_table,adjp)
write.table(enrich_table,file="GenomeDist_SVs_enrichment.txt",sep="\t",row.names=F,quote=F)

#Enrichment results by cohort
aocs_svs2 <- svs2[svs2$Cohort=="AO",]
shgsoc_svs2 <- svs2[svs2$Cohort=="SH",]
tcga_svs2 <- svs2[svs2$Cohort=="DO",]
mda_svs2 <- svs2[svs2$Cohort=="WG",]
bc_svs2 <- svs2[svs2$Cohort=="DA"| svs2$Cohort=="DG",]


getChr19SubCohort <- function(cohort_logic_svs,cohort_logic_comp,Cohort){
  svs2 <- svs2[cohort_logic_svs,]
  sv_counts<-table(svs2$V11,svs2$V4)
  rownames(sv_counts)[1]<-"All_SV"
  
  complex2 <- complex2[cohort_logic_comp,]
  complex_counts<-table(complex2[,"category"],complex2[,2])
  rownames(complex_counts)[1]<-"All_complex"
  
  all_counts<- rbind(sv_counts,complex_counts)
  genome <- rowSums(all_counts)
  all_counts<- cbind(all_counts,genome)
  calculate_expected_p <- chrom_lengths[,2]/3031042417
  
  class_sv <- character(dim(all_counts)[1]*(dim(all_counts)[2]-1))
  chr <- character(dim(all_counts)[1]*(dim(all_counts)[2]-1))
  or <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))
  lci <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))
  uci <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))
  p <- numeric(dim(all_counts)[1]*(dim(all_counts)[2]-1))
  
  count<- 1
  for (i in 1:dim(all_counts)[1]){
    for (j in 1:(dim(all_counts)[2]-1)){
      
      class_sv[count] <- rownames(all_counts)[i]
      chr[count] <- colnames(all_counts)[j]
      
      b<-binom.test(all_counts[i,j],all_counts[i,"genome"],p=calculate_expected_p[j])
      
      or[count] <- b$estimate/b$null.value
      lci[count] <- b$conf.int[1]/b$null.value
      uci[count] <- b$conf.int[2]/b$null.value
      p[count] <- b$p.value
      
      count<- count+1
    }
  }
  enrich_table<- data.frame(Class_SV=class_sv,Chromosome=chr,OR=or,LCI=lci,UCI=uci,Pvalue=p)
  adjp <- p.adjust(enrich_table$Pvalue,"bonferroni")
  enrich_table<- cbind(enrich_table,adjp)
  write.table(enrich_table,file=paste(Cohort,"_GenomeDist_SVs_enrichment.txt",sep=""),sep="\t",row.names=F,quote=F)
  
}
getChr19SubCohort(svs2$Cohort=="AO",complex2$Cohort=="AOCS",Cohort="AOCS")
getChr19SubCohort(svs2$Cohort=="SH",complex2$Cohort=="SHGSOC",Cohort="SHGSOC")
getChr19SubCohort(svs2$Cohort=="DO",complex2$Cohort=="TCGA",Cohort="TCGA")
getChr19SubCohort(svs2$Cohort=="WG",complex2$Cohort=="MDA",Cohort="MDA")### - needs checked 
getChr19SubCohort(svs2$Cohort=="DA"| svs2$Cohort=="DG",complex2$Cohort=="BCCA",Cohort="BCCA")

#For MDA no chr21 complex- so add zero column to counts
complex_counts<- cbind(complex_counts,"21"=rep(0,7))
complex_counts<- complex_counts[,c(1:20,23,21:22)]

#Enriched chrs by event
enriched_chrs <- enrich_table[enrich_table$adjp<0.05 & enrich_table$OR >1,]
depleted_chrs <- enrich_table[enrich_table$adjp<0.05 & enrich_table$OR <1,]

enriched_chrs_sv <- enriched_chrs[enriched_chrs[,1] %in% c("All_SV","DEL","DUP","CTX","INV"),]
enriched_chrs_sv[enriched_chrs_sv[,1]=="All_SV",1]<-"All"

xmin <- numeric(dim(enriched_chrs_sv)[1])
xmax <- numeric(dim(enriched_chrs_sv)[1])

for (i in 1:dim(enriched_chrs_sv)[1]){
  if (enriched_chrs_sv[i,2]==1){
    xmin[i]<-0
  }else{
  xmin[i] <- chrom_lengths[as.numeric(enriched_chrs_sv[i,2])-1,3]
  }
  xmax[i] <- chrom_lengths[enriched_chrs_sv[i,2],3]
}

enrich_rect_text <- data.frame(V11 = enriched_chrs_sv[,1],
                        xmin = xmin,
                        xmax   = xmax, 
                        ymin=rep(-Inf,length(xmin)),
                        ymax=rep(Inf, length(xmin)))

depleted_chrs_sv <- depleted_chrs[depleted_chrs[,1] %in% c("All_SV","DEL","DUP","CTX","INV"),]
depleted_chrs_sv[depleted_chrs_sv[,1]=="All_SV",1]<-"All"

xmin <- numeric(dim(depleted_chrs_sv)[1])
xmax <- numeric(dim(depleted_chrs_sv)[1])

for (i in 1:dim(depleted_chrs_sv)[1]){
  if (depleted_chrs_sv[i,2]==1){
    xmin[i]<-0
  }else{
    xmin[i] <- chrom_lengths[as.numeric(depleted_chrs_sv[i,2])-1,3]
  }
  xmax[i] <- chrom_lengths[depleted_chrs_sv[i,2],3]
}

deplete_rect_text <- data.frame(V11 = depleted_chrs_sv[,1],
                               xmin = xmin,
                               xmax   = xmax, 
                               ymin=rep(-Inf,length(xmin)),
                               ymax=rep(Inf, length(xmin)))
pdf("individual_panels/supp_figs/SuppFig6A_SV_density_acrossGenome_withpeaks.png",width=12.2,height=5)
ggplot(svs2,aes(total_pos,fill=V11))+geom_density(adjust=0.1)+theme_classic(base_size =13)+
  facet_grid(V11~.,labeller=facetlabs)+
  scale_x_continuous(breaks=all_marks,labels=chr_labels,limits=c(0,all_marks[47]))+
  xlab("Chromosome")+ylab("SV breakpoint density")+
  theme(legend.position="none",
        strip.text.y = element_text(size=18, face="bold",angle=0,hjust=0,colour="black"),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=16,colour="black"),
        axis.title.x=element_text(size=16,colour="black"),
        axis.text.x=element_text(size=12,colour="black"))+
  scale_fill_manual(labels=c("Duplication","Deletion","Inversion","Translocation","All"),values=lgd2_cols)+
  geom_rect(data=enrich_rect_text, aes(xmin = xmin,xmax = xmax, ymin=ymin,ymax=ymax), fill= "#0072B2",
            alpha=.1,inherit.aes = FALSE )
dev.off()
  #geom_rect(data=deplete_rect_text, aes(xmin = xmin,xmax = xmax, ymin=ymin,ymax=ymax), fill= "#CC79A7",
#          alpha=.1,inherit.aes = FALSE )

#Enriched complex peaks
enriched_chrs_complex <- enriched_chrs[enriched_chrs[,1] %in% c("All_complex","bfb","Chromoplexy",
                                                                "AAecDNA","Chromothripsis","Pyrgo","Rigma"),]
enriched_chrs_complex[enriched_chrs_complex[,1]=="All_complex",1]<-"All"

xmin <- numeric(dim(enriched_chrs_complex)[1])
xmax <- numeric(dim(enriched_chrs_complex)[1])

for (i in 1:dim(enriched_chrs_complex)[1]){
  if (enriched_chrs_complex[i,2]==1){
    xmin[i]<-0
  }else{
    xmin[i] <- chrom_lengths[as.numeric(enriched_chrs_complex[i,2])-1,3]
  }
  xmax[i] <- chrom_lengths[enriched_chrs_complex[i,2],3]
}

enrich_rect_text <- data.frame(category = enriched_chrs_complex[,1],
                               xmin = xmin,
                               xmax   = xmax, 
                               ymin=rep(-Inf,length(xmin)),
                               ymax=rep(Inf, length(xmin)))
enrich_rect_text$category <-as.factor(enrich_rect_text$category)

pdf("individual_panels/supp_figs/SuppFig6B_complexSV_density_acrossGenome_withpeaks.pdf",width=13,height=6.5)

ggplot(complex2,aes(total_pos,fill=category))+geom_density(adjust=0.1)+theme_classic(base_size=13)+
  facet_grid(category~.,labeller=complex_labels_nice)+
  scale_fill_manual(values=complex_cols)+
  scale_x_continuous(breaks=all_marks,labels=chr_labels,limits=c(0,all_marks[47]))+
  xlab("Chromosome")+ylab("Complex SV breakpoint density")+
  theme(legend.position="none",
        strip.text.y = element_text(size=18, face="bold",angle=0,hjust=0,colour="black"),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=16,colour="black"),
        axis.title.x=element_text(size=16,colour="black"),
        axis.text.x=element_text(size=12,colour="black"))+
  geom_rect(data=enrich_rect_text, aes(xmin = xmin,xmax = xmax, ymin=ymin,ymax=ymax), fill= "#0072B2",
            alpha=.1,inherit.aes = FALSE )
dev.off()
#geom_rect(data=deplete_rect_text, aes(xmin = xmin,xmax = xmax, ymin=ymin,ymax=ymax), fill= "#CC79A7",
#          alpha=.1,inherit.aes = FALSE )
