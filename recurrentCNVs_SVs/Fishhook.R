#Look at fishhook results
library(fishHook)
load("~/Documents/HGSOC_marker_paper/RecurrentSVs/fish.supertiles.RData")

fish.supertiles$res[1:10,]
df<-gr2dt(fish.supertiles$res)
df_fdr<- df[df$fdr<0.05,]
df_sort<- df_fdr[order(df_fdr$p),]
df_sort_chr <-df_sort[order(df_sort$seqnames),]

fish.supertiles$qqp(plotly = FALSE)

#Fish at gene level

load("RecurrentSVs/fish.genes.RData")
fish.genes$qqp(plotly=FALSE)

#Fish allowing clusters of brkpts
load("RecurrentSVs/fish.supertiles_allowsclusters.RData")
fish.supertiles$qqp(plotly = FALSE)
fish.supertiles$res[1:10,]
df2<-gr2dt(fish.supertiles$res)
df2_fdr<- df2[df2$fdr<0.05,]
df2_sort<- df2_fdr[order(df2_fdr$p),]
df2_sort_chr <-df2_sort[order(df2_sort$seqnames),]

#Plot results across the genome

hg38_chromsizes <- read.table("~/Documents/HGSOC_marker_paper_onedrive/metrics/hg38.chrom.sizes.txt",sep="\t")
rownames(hg38_chromsizes)<- as.character(hg38_chromsizes[,1])
chrs<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
        "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
hg38_chromsizes <- hg38_chromsizes[c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                     "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"),]
cum_sizes <- numeric(dim(hg38_chromsizes)[1])
cum_sizes[1] <- hg38_chromsizes[1,2]
for (i in 2:length(cum_sizes)){
  cum_sizes[i] <- cum_sizes[i-1]+ hg38_chromsizes[i,2]
}
hg38_chromsizes_cuml <- cbind(hg38_chromsizes,cum_sizes)
edges <- c(1, hg38_chromsizes_cuml[,3])
get_midpoints <- function(x){
  
  midpoints <- numeric(length(x)-1)
  for (i in 1:(length(x)-1)){
    midpoints[i] <- (x[i+1] + x[i])/2
  }
  return(midpoints)
}

mids <- get_midpoints(edges)
df$mid<- df$start + (df$end - df$start)/2 

#newcoords
point_cuml <- df[df$seqnames=="chr1",]
point_cuml$mid_cuml <- point_cuml$mid

for (i in 2:24){
  newchr <- df[df$seqnames==chrs[i],]
  newchr$mid_cuml <- hg38_chromsizes_cuml[i-1,3] + newchr$mid
  point_cuml<- rbind(point_cuml,newchr)
}

point_cuml$sig <- ifelse(point_cuml$fdr<0.05,"sig","notsig")

ggplot(point_cuml)+geom_bar(aes(x = mid_cuml, y = effectsize, col=as.factor(sig),fill=as.factor(sig)),stat="identity") + scale_x_continuous(breaks=mids, labels=chrs) +
  xlab("Chromosome") +theme_bw()

ggplot(point_cuml)+geom_line(aes(x = mid_cuml, y = -log10(p))) + scale_x_continuous(breaks=mids, labels=chrs) +
  xlab("Chromosome") +theme_bw()

#Calculate rearrangement dispersion scores

library(VariantAnnotation)
library(StructuralVariantAnnotation)

vcf_df <- read.table("RecurrentSVs/Consensus100bpSVs_ids_infishHookregions_regions.txt")
vcf_df$region <- paste(vcf_df$V1,vcf_df$V2,vcf_df$V3,sep=":")

vcf <- VariantAnnotation::readVcf("RecurrentSVs/Consensus100bpSVs_ids_infishHookregions_withheader.vcf","hg38")
#add samples
samples <- read.table("RecurrentSVs/samples_infishhookregions.txt",sep="\t")
samples <- as.character(samples[,1])
info(vcf)$SAMPLE <- samples
info(vcf)$region <- vcf_df$region
gr<-breakpointRanges(vcf,info_columns=c("SAMPLE","region"),inferMissingBreakends=TRUE)
sv_recurr_dat <- gr2dt(gr)
sv_recurr_dat <- as.data.frame(sv_recurr_dat)
#Split SVs by recurrent region

sv_recurr_split <- split(sv_recurr_dat,as.factor(sv_recurr_dat$region))
dispersion_score<-numeric(length(sv_recurr_split))
region<-character(length(sv_recurr_split))
for (i in 1:length(sv_recurr_split)){
  #for each region- assign sv length of 10^9 for BNDs,calculate median absolute deviation in abs(svlength), median abs(svlength), ratio of 2.
  sv_recurr_split[[i]][is.na(sv_recurr_split[[i]]$svLen),"svLen"] <- 1e9
  med_svlength <- median(abs(sv_recurr_split[[i]]$svLen))
  mad_svlength <- mad(abs(sv_recurr_split[[i]]$svLen))
  dispersion_score[i] <- mad_svlength/med_svlength 
  region[i]<- sv_recurr_split[[i]][1,"region"]
}

ds <- data.frame(ID=1:length(dispersion_score),DS=sort(dispersion_score),region=region)

#Hartigans' dip test for unimodality / multimodality
#data:  ds$DS, alternative hypothesis: non-unimodal, i.e., at least bimodal
#D = 0.083753, p-value = 0.102

plot(ds$ID,ds$DS,pch=19)
ds <- ds[,c("DS","region")]
#Recurrent bkpts
bkpt_regions <- read.table("RecurrentSVs/Breakpoint_enriched_50kbbins_fishHook_gistic_genes_with_chr.txt",sep="\t")
bkpt_regions$region <- paste(bkpt_regions$V1,bkpt_regions$V2, bkpt_regions$V3,sep="")
bkpt_regions_ds <- merge(bkpt_regions,ds,by="region")
write.table(bkpt_regions_ds,file="RecurrentSVs/Breakpoint_enriched_50kbbins_fishHook_gistic_genes_withchr_ds.txt",sep="\t",row.names=F,col.names=F,quote=F)


#Get SVtype split

for (i in 1:length(sv_recurr_split)){
  print(sv_recurr_split[[i]][1,"region"])
  dist<-table(sv_recurr_split[[i]]$svtype)/dim(sv_recurr_split[[i]])[1]
  print(dist)
}

SVcounts<- mat.or.vec(nc=4,nr=length(sv_recurr_split))
regions<-character(length(sv_recurr_split))
for (i in 1:length(sv_recurr_split)){
  regions[i]<-sv_recurr_split[[i]][1,"region"]
  newrow<-data.frame(table(sv_recurr_split[[i]]$svtype))[,2]
  l<-length(newrow)
  SVcounts[i,(5-l):4]<-data.frame(table(sv_recurr_split[[i]]$svtype))[,2]
}
colnames(SVcounts)<-c("BND","DEL","DUP","INV")
rownames(SVcounts)<-regions
sums<-rowSums(SVcounts)
SVrates<-SVcounts/sums

heatmap(SVcounts)
heatmap(SVrates)

library(reshape2)
library(PNWColors)
pal<-pnw_palette("Bay",4,type="discrete")
long.counts<- melt(SVcounts)
orderregions<-rownames(SVrates)[heatmap(SVrates)$rowInd]
long.counts$Var1<-factor(long.counts$Var1,levels=orderregions)

bar_counts<- ggplot(long.counts, aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+ylab("Number of Breakpoints")+
  xlab("Recurrent SV region")+scale_fill_manual(name="SV type",values=pal)

bar_rates <- ggplot(long.counts, aes(x=Var1,y=value,fill=Var2))+geom_bar(stat="identity",position="fill")+theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))+ylab("Proportion of Breakpoints")+
  xlab("Recurrent SV region")+scale_fill_manual(name="SV type",values=pal)
library(patchwork)

bar_counts+bar_rates

#Samples with breakpoints in recurrent regions
#Load sample data

contrib_samples <- read.table("Samples_contributing_breakpoints_toFishHookpeaks.txt",sep="\t")
contrib_samples$region <- paste(contrib_samples$V12,contrib_samples$V13,contrib_samples$V14,sep=":")
contrib_samples_df <- as.matrix(table(contrib_samples$region,contrib_samples$V11))
contrib_samples_df[contrib_samples_df>1]<-1

write.table(contrib_samples_df,file="RecurrentSVs/Samples_withatleast1Brkpt_inRecurrentSVregions.txt",sep="\t",quote=F)

#Rough summary numbers
rowSums(contrib_samples_df)

#load cn data
dat <- read.table("RecurrentSVs/Consensus_del_calls_atFishHookPeaks.bed",sep="\t")

#remove regions with no cn call overlap
dat_withconsensuscalls <- dat[dat$V5!="-1",]
dat_withconsensuscalls$regions <- paste(dat_withconsensuscalls$V1,dat_withconsensuscalls$V2,dat_withconsensuscalls$V3,sep=":")
dat_withconsensuscalls$region_length <- dat_withconsensuscalls$V3-dat_withconsensuscalls$V2


dat2 <- read.table("RecurrentSVs/Consensus_dup_calls_atFishHookPeaks.bed",sep="\t")


#remove regions with no cn call overlap
dat2_withconsensuscalls <- dat2[dat2$V5!="-1",]
dat2_withconsensuscalls$regions <- paste(dat2_withconsensuscalls$V1,dat2_withconsensuscalls$V2,dat2_withconsensuscalls$V3,sep=":")
dat2_withconsensuscalls$region_length <- dat2_withconsensuscalls$V3-dat2_withconsensuscalls$V2



#For samples mentioned calculate average cn per region and input into matrix
all_samples <- read.table("/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/params/Combined_usable_ids_marker_altdupids.txt",sep="\t")
all_samples <- as.character(all_samples[,1])
all_samples[1:115] <- substr(all_samples[1:115],1,9)
all_samples[299:324] <- sapply(strsplit(all_samples[299:324],"-"), function(x) paste(x[1:3],collapse = "-"))

uni_cn_regions<- unique(c(dat_withconsensuscalls$regions,dat2_withconsensuscalls$regions))
cn_mat <- mat.or.vec(nr=length(uni_cn_regions),nc=length(all_samples))
rownames(cn_mat)<-uni_cn_regions
colnames(cn_mat)<-all_samples
cn_mat[,]<- 2

#Combine dels and dups
gistic_withconsensus <- rbind(dat_withconsensuscalls,dat2_withconsensuscalls)
gistic_withconsensus_ord <- gistic_withconsensus[order(gistic_withconsensus$regions,gistic_withconsensus$V8),]


#Steps
#split by region
consensus_split_region <- split(gistic_withconsensus_ord,as.factor(gistic_withconsensus_ord$regions))

for (i in 1:length(consensus_split_region)){
  #split by sample
  consensus_split_sample <- split(consensus_split_region[[i]],as.factor(consensus_split_region[[i]]$V8))
  
  for (j in 1:length(consensus_split_sample)){
    
    #calculate how much of region is diploid in each sample
    
    nodip_bases <- sum(consensus_split_sample[[j]][consensus_split_sample[[j]]$V7!=2,"V13"])
    dip_bases <- consensus_split_sample[[j]][1,"region_length"]-nodip_bases
    av_cn <- sum(as.numeric(consensus_split_sample[[j]]$V7)*(as.numeric(consensus_split_sample[[j]]$V13)/consensus_split_sample[[j]][1,"region_length"])) + 2*(dip_bases/consensus_split_sample[[j]][1,"region_length"])
    cn_mat[consensus_split_sample[[j]][1,"regions"],consensus_split_sample[[j]][1,"V8"]]<-av_cn
  }
}



#Compare CN between contributing and non-contributing samples
library(MASS)


sv_sample_contributions<- read.table("RecurrentSVs/Samples_withatleast1Brkpt_inRecurrentSVregions.txt",sep="\t",header=T)
rownames(sv_sample_contributions)<- gsub('chr','',rownames(sv_sample_contributions))
#split wgs ids by . replace with - and remove P
colnames(sv_sample_contributions)[287:307]<- sapply(strsplit(colnames(sv_sample_contributions)[287:307],'\\.'), function(x) paste(x[[1]],x[[2]],x[[3]],sep="-"))

absent_samples <- setdiff(colnames(cn_mat),colnames(sv_sample_contributions))

absent_df <- mat.or.vec(nr=dim(sv_sample_contributions)[1],nc=length(absent_samples))
absent_df[,]<- 0
colnames(absent_df)<-absent_samples
sv_sample_contributions<- cbind(sv_sample_contributions,absent_df)

results_df<- data.frame()
for (i in 1:dim(cn_mat)[1]){
  samples_with_nums <- which(sv_sample_contributions[rownames(cn_mat)[i],]==1)
  samples_with_ids<- colnames(sv_sample_contributions)[samples_with_nums]
  samples_without<-setdiff(colnames(cn_mat),samples_with_ids)
  
  cn_with<-cn_mat[i,samples_with_ids]
  df_with<- data.frame(brkpt="With",cn=cn_with)
  cn_without<-cn_mat[i,samples_without]
  df_without<- data.frame(brkpt="Without",cn=cn_without)
  df <- rbind(df_with,df_without)
  df$cn <- round(df$cn,0)
  ggplot(df, aes(x=brkpt,y=cn))+geom_boxplot(outlier.shape=NA)+geom_jitter()+theme_bw()+ylab("Copy number")+xlab("Breakpoint involvement")
  mod<- glm.nb(cn~as.factor(brkpt), data=df)
  or<- signif(exp(summary(mod)$coefficients[2,1]  ),2)
  lci <-signif(exp(summary(mod)$coefficients[2,1] -1.96*summary(mod)$coefficients[2,2] ),2)
  uci <-signif(exp(summary(mod)$coefficients[2,1] +1.96*summary(mod)$coefficients[2,2] ),2)
  pval <- signif(summary(mod)$coefficients[2,4],2)
  newres <- c(Region=rownames(cn_mat)[i],NSamples_withBrkpt=length(samples_with_ids),OR=or,LCL=lci,UCL=uci,Pvalue=pval)
  results_df<-rbind(results_df,newres)
}
colnames(results_df)<-c("Region","NSamples_withBrkpt","OR","LCL","UCL","Pvalue")

write.table(results_df,file="Impact_of_breakpoint_inrecurrentSVregion_onCN_glm.nb_res.txt",sep="\t",quote=F,row.names=F)

