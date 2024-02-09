#Cellularity and coverage plots

#Load data

#Cellularity
cellularity <- read.table("HGSOC_cellularity_estimates.txt",sep="\t")
subco <- substr(cellularity[,1],1,2)
cohort <- sapply(subco, function(x) switch(x,"AO"="AOCS","SH"="SHGSOC","DG"="BCCA","DA"="BCCA","DO"="TCGA","WG"="MDA"))
cellularity<- cbind(cellularity,cohort)

col = c("AOCS" = "blue",
        "BCCA" = "black",
        "SHGSOC" = "brown",
        "TCGA" = "grey",
        "MDA" = "orange")

library(ggplot2)
cell<-  ggplot(cellularity,aes(x=cohort,y=V2))+
  geom_boxplot(fill="#0072B2",width=0.3)+theme_classic()+
  # scale_colour_manual(name="Cohort",values=col)+
  ylim(c(0,1))+ylab("Cellularity")+xlab("Cohort")

#Read depth
library(readxl)
coverage <- read_excel("coverage.xlsx")
coverage <- as.data.frame(coverage)
coverage$cohort_type <- paste(coverage$Cohort,coverage$Type,sep="_")

type_cols = c("Tumour" = "#0072B2", "Normal" = "#009E73")

cov<-ggplot(coverage,aes(x=Cohort,y=Avg_coverage,fill=Type))+
  geom_boxplot()+
  theme_classic()+scale_fill_manual(values=type_cols)+
  ylab("Mean coverage")+xlab("Cohort")

library(cowplot)
plot_grid(cov, cell, labels = c('A', 'B'), label_size = 13)

#Number of SVs/CNVs

#Load SVs
svs_all <- read.table("data/svFormatted.bedpe",sep="\t")
svs_all <- svs_all[svs_all$V11!="INS",]

svs_sample_type <- as.data.frame(table(svs_all$V13,svs_all$V11))
subco <- substr(svs_sample_type$Var1,1,2)
cohort <- sapply(subco, function(x) switch(x,"AO"="AOCS","SH"="SHGSOC","DG"="BCCA","DA"="BCCA","DO"="TCGA","WG"="MDA"))
svs_sample_type <- cbind(svs_sample_type,cohort)

sv_cols =  c("DUP" = "#009E73", "DEL" = "#0072B2", "INV" = "#CC79A7", "CTX" = "#E69F00")

svplot<-ggplot(svs_sample_type,aes(x=cohort,y=log10(Freq), fill=Var2))+
  geom_boxplot()+
  scale_fill_manual(name="Type",values=sv_cols,labels=c(
    "DUP"="Duplication",
    "DEL"="Deletion",
    "INV"="Inversion",
    "CTX"="Translocation"))+
  theme_classic()+xlab("Cohort")+ylab("Number of SVs (log10)")

#Load CNVs
dels <- read.table("data/delFormatted.bed") 
dels$type<- "Del"
dups <- read.table("data/dupFormatted.bed") 
dups$type <- "Dup"

df <- rbind(dels,dups)
d<- as.data.frame(table(df$V5,df$type))
subco <- substr(d$Var1,1,2)
cohort <- sapply(subco, function(x) switch(x,"AO"="AOCS","SH"="SHGSOC","DG"="BCCA","DA"="BCCA","DO"="TCGA","WG"="MDA"))
d <- cbind(d,cohort)

cnv_cols =  c("Dup" = "#009E73", "Del" = "#0072B2")

cnv_plot<- ggplot(d,aes(x=cohort,y=log10(Freq), fill=Var2))+
  geom_boxplot()+
  scale_fill_manual(name="Type",values=cnv_cols,labels=c(
    "Dup"="Duplication",
    "Del"="Deletion"))+
  theme_classic()+xlab("Cohort")+ylab("Number of CNAs (log10)")

pdf("individual_panels/supp_figs/SuppFig2AB_numberSVsCNA.pdf",width=12,height=4)
plot_grid(svplot, cnv_plot, labels = c('A', 'B'), label_size = 16)
dev.off()
