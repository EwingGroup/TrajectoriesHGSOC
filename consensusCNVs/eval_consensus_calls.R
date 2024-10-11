# Compare consensus calls to individual calls

#Load libraries
library(VariantAnnotation)
library(reshape2)
library(ggplot2)

#Set colours
SVtypes_cols<- c("#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
 
#Load SV calls
manta <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/data/Combined_manta_pass.bedpe", sep='\t')
manta$Sample <- gsub("_PrimaryTumour","",manta$V13)
gridss <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/data/Combined_gridss_pass.bedpe", sep='\t')
gridss$Sample <- gsub("_PrimaryTumour","",gridss$V13)

consensus100 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/data/Combined_consensus_100bp.bedpe", sep='\t')
consensus100$Sample <- gsub("_PrimaryTumour","",consensus100$V13)
consensus50 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/data/Combined_consensus_50bp.bedpe", sep='\t')
consensus50$Sample <- gsub("_PrimaryTumour","",consensus50$V13)
consensus150 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/data/Combined_consensus_150bp.bedpe", sep='\t')
consensus150$Sample <- gsub("_PrimaryTumour","",consensus150$V13)

dim(manta)
dim(gridss)
dim(consensus150)
dim(consensus100)
dim(consensus50)

#SV calls by type
manta_by_type <- table(manta$Sample,manta$V11)
manta_by_type<-data.frame(manta_by_type) 
colnames(manta_by_type)<- c("Sample","Type","Manta")

gridss_by_type <- table(gridss$Sample,gridss$V11)
gridss_by_type<-data.frame(gridss_by_type) 
colnames(gridss_by_type)<- c("Sample","Type","Gridss")

consensus100_by_type <- table(consensus100$Sample,consensus100$V11)
consensus100_by_type<-data.frame(consensus100_by_type) 
colnames(consensus100_by_type)<- c("Sample","Type","Consensus100")

consensus50_by_type <- table(consensus50$Sample,consensus50$V11)
consensus50_by_type<-data.frame(consensus50_by_type) 
colnames(consensus50_by_type)<- c("Sample","Type","Consensus50")

consensus150_by_type <- table(consensus150$Sample,consensus150$V11)
consensus150_by_type<-data.frame(consensus150_by_type) 
colnames(consensus150_by_type)<- c("Sample","Type","Consensus150")

two <- merge(manta_by_type,gridss_by_type,by=c("Sample","Type"),all=TRUE)
three <- merge(two,consensus50_by_type,by=c("Sample","Type"),all=TRUE)
four <- merge(three,consensus100_by_type,by=c("Sample","Type"),all=TRUE)
all <- merge(four,consensus150_by_type,by=c("Sample","Type"),all=TRUE)

all[is.na(all)]<-0

#Make comparison plots
all_long <- melt(all,id.vars=c("Sample","Type"))
colnames(all_long) <- c("Sample","Type","Caller","Freq")
logFreq <- log(all_long$Freq)
all_long <- data.frame(all_long,logFreq)

ggplot(all_long, aes(x=Type, y=logFreq, fill=Type)) + facet_wrap(.~Caller) + theme_bw() + geom_boxplot() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab('Counts of SVs (log counts)') + xlab("Type")+scale_fill_manual(values=SVtypes_cols)


# Plot just manta, gridss and consensus 100bp

(144950/139761)*11

all_long_orig <- all_long[all_long$Caller %in% c("Manta","Gridss"),]
#png("Comparison_variants_by_caller.png",res=300, unit="in",height=4, width=8)
ggplot(all_long_orig, aes(x=Type, y=logFreq, fill=Type)) + facet_wrap(.~Caller) + theme_bw() + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab('Counts of SVs (log counts)') + xlab("Type")+scale_fill_manual(values=SVtypes_cols)
#dev.off()

all_long2 <- all_long[all_long$Caller %in% c("Manta","Gridss","Consensus100"),]
#png("Comparison_variants_by_caller_andconsensus.png",res=300, unit="in",height=4, width=8)
ggplot(all_long2, aes(x=Type, y=logFreq, fill=Type)) + facet_wrap(.~Caller) + theme_bw() + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab('Counts of SVs (log counts)') + xlab("Type")+scale_fill_manual(values=SVtypes_cols)
#dev.off()

#Between consensus thresholds - 50bp,100bp,150bp
all_long3 <- all_long[all_long$Caller %in% c("Consensus50","Consensus100","Consensus150"),]
all_long3$Caller <- factor(all_long3$Caller, levels=c("Consensus50","Consensus100","Consensus150"))
#png("Comparison_variants_by_consensusthreshold.png",res=300, unit="in",height=4, width=8)
ggplot(all_long3, aes(x=Type, y=logFreq, fill=Type)) + facet_wrap(.~Caller) + theme_bw() + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab('Counts of SVs (log counts)') + xlab("Type")+scale_fill_manual(values=SVtypes_cols)
#dev.off()

table(all_long3$Caller,all_long3$Type)

##########################################################################################################################
# PCAWG
##########################################################################################################################


#Load PCAWG calls and subset to overlap in AOCS samples
pcawg <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/data/Combined_PCAWG_consensus_svs.bedpe", sep='\t')
colnames(pcawg)[13]<-"pcawgfile"
pcawg_ids <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/data/donor_to_file_prefix.txt",sep="\t")
pcawg_ids[,1] <- gsub("S","S_",pcawg_ids[,1])
colnames(pcawg_ids)<-c("Sample","pcawgfile")

pcawg_withsample <-merge(pcawg,pcawg_ids,by="pcawgfile")
manta_pcawg <- manta[manta$Sample %in% pcawg_withsample$Sample,]
gridss_pcawg <- gridss[gridss$Sample %in% pcawg_withsample$Sample,]
consensus100_pcawg <- consensus100[consensus100$Sample %in% pcawg_withsample$Sample,]

dim(pcawg_withsample)
dim(manta_pcawg)
dim(gridss_pcawg)
dim(consensus100_pcawg)

#SV calls by type
manta_by_type_p <- table(manta_pcawg$Sample,manta_pcawg$V11)
manta_by_type_p<-data.frame(manta_by_type_p) 
colnames(manta_by_type_p)<- c("Sample","Type","Manta")

gridss_by_type_pcawg <- table(gridss_pcawg$Sample,gridss_pcawg$V11)
gridss_by_type_pcawg<-data.frame(gridss_by_type_pcawg) 
colnames(gridss_by_type_pcawg)<- c("Sample","Type","Gridss")

consensus100_by_type_pcawg <- table(consensus100_pcawg$Sample,consensus100_pcawg$V11)
consensus100_by_type_pcawg<-data.frame(consensus100_by_type_pcawg) 
colnames(consensus100_by_type_pcawg)<- c("Sample","Type","Consensus100")

pcawg_withsample_by_type <- table(pcawg_withsample$Sample,pcawg_withsample$V11)
pcawg_withsample_by_type<-data.frame(pcawg_withsample_by_type) 
colnames(pcawg_withsample_by_type)<- c("Sample","Type","PCAWG")

two <- merge(manta_by_type_p,gridss_by_type_pcawg,by=c("Sample","Type"),all=TRUE)
four <- merge(two,consensus100_by_type_pcawg,by=c("Sample","Type"),all=TRUE)
all <- merge(four,pcawg_withsample_by_type,by=c("Sample","Type"),all=TRUE)

all[is.na(all)]<-0

#Make comparison plots
all_long <- melt(all,id.vars=c("Sample","Type"))
colnames(all_long) <- c("Sample","Type","Caller","Freq")
logFreq <- log(all_long$Freq)
all_long <- data.frame(all_long,logFreq)

ggplot(all_long, aes(x=Type, y=logFreq, fill=Type)) + facet_grid(.~Caller) + theme_bw() + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab('Counts of SVs (log counts)') + xlab("Type")+scale_fill_manual(values=SVtypes_cols)

################################################################## CNV consensus calls ######################################################################

#Load CNVkit calls

cnvkit <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_cnvkit_ids.bed",sep="\t")
type<- ifelse(cnvkit$V5>2,"DUP","DEL")
cnvkit <- cbind(cnvkit,type)
size_bin <- ifelse(cnvkit$V7 < 1000,"<1kb",
               ifelse(cnvkit$V7 < 10000,"1kb-10kb",
                      ifelse(cnvkit$V7 < 100000,"10kb-100kb",
                             ifelse(cnvkit$V7 < 1000000, "100kb-1Mb",">1Mb"))))
cnvkit <- cbind(cnvkit,size_bin)
cnvkit.typebysample <-table(cnvkit$V4,cnvkit$type)
cnvkit.typebysample <-data.frame(cnvkit.typebysample,Caller=rep("CNVkit",dim(cnvkit.typebysample)[1]))

summary(cnvkit$V7)

climat <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_CLImAT_ids.bed", sep="\t")
type<- ifelse(climat$V4>2,"DUP","DEL")
climat <- cbind(climat,type)
size <- climat$V3 - climat$V2 +1
climat <- cbind(climat,size)
climat.typebysample <-table(climat$V7,climat$type)
climat.typebysample <-data.frame(climat.typebysample,Caller=rep("CLImAT",dim(climat.typebysample)[1]))

cnvkit.climat <- rbind(cnvkit.typebysample,climat.typebysample)

ggplot(cnvkit.climat,aes(x=Caller,y=Freq,fill=))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+facet_wrap(~Var2)+theme_bw()

cnvkit.climat.0.5 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_consensus_cnvkit_CLImAT_0.5.bed",sep="\t")
cnvkit.climat.0.75 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_consensus_cnvkit_CLImAT_0.75.bed",sep="\t")
cnvkit.climat.0.9 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_consensus_cnvkit_CLImAT_0.9.bed",sep="\t")
cnvkit.climat.0.99 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_consensus_cnvkit_CLImAT_0.99.bed",sep="\t")


cnv_pcawg_0.5 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_PCAWG_consensus_cnvkit_CLImAT_0.5.bed",sep="\t")
cnv_pcawg_0.75 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_PCAWG_consensus_cnvkit_CLImAT_0.75.bed",sep="\t")
cnv_pcawg_0.9 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_PCAWG_consensus_cnvkit_CLImAT_0.9.bed",sep="\t")
cnv_pcawg_0.99 <- read.table("~/Documents/HGSOC_marker_paper_onedrive/SVhotspots/variants/Combined_PCAWG_consensus_cnvkit_CLImAT_0.99.bed",sep="\t")


# Compare CNVkit with CLImAT

#Take a quick look by size 

sizes_all <-c(cnvkit$V7,climat$size)
#,(cnvkit.climat.0.5$V3-cnvkit.climat.0.5$V2 +1),
            #  (cnvkit.climat.0.75$V3-cnvkit.climat.0.75$V2 +1),(cnvkit.climat.0.9$V3-cnvkit.climat.0.9$V2 +1))
              #(cnvkit.climat.0.99$V3-cnvkit.climat.0.99$V2 +1))
callers <- c(rep("CNVkit",length(cnvkit$V7)),rep("CLImAT",length(climat$size)))
             
          #   ,rep("CLImAT + CNVkit, 0.5 overlap",dim(cnvkit.climat.0.5)[1]),
            # rep("CLImAT + CNVkit, 0.75 overlap",dim(cnvkit.climat.0.75)[1]),rep("CLImAT + CNVkit, 0.9 overlap",dim(cnvkit.climat.0.9)[1]))
          #   rep("CLImAT + CNVkit, 0.99 overlap",dim(cnvkit.climat.0.99)[1]))
df <- data.frame(caller=callers,sizes=sizes_all)

ggplot(df,aes(x=log10(sizes),colour=caller,fill=caller))+geom_density(alpha=0.3)+theme_bw()+ 
  ylab("Number of CNVs") + xlab("CNV length (log 10 scale)")+scale_x_continuous(breaks=c(0,2,4,6,8),labels=c("0","100bp","10kb","1Mb","100Mb"))

table(df$caller)
cnvkit[1,]
climat[1,]


