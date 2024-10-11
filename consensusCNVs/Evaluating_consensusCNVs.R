#This script combines the pairwise consensus calls into 3 way consensus calls
# Two way consensus sets are on eddie.
#Also it evaluates the number of calls/concordance with PCAWG

#Load packages
library(reshape2)
library(ggplot2)

# Compare with original CNVkit calls
cnvkit <- read.table("/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/Combined_cnvkit_ids.bed",sep="\t")

#Standardise IDs
cnvkit$cohort <- substr(cnvkit$V4,1,2)
cnvkit$Sample <- rep(NA,dim(cnvkit)[1])

stand_ids <- function(df,sample_var){
  df$cohort <- substr(df[,sample_var],1,2)
  df$Sample <- rep(NA,dim(df)[1])
  df[df$cohort=="SH","Sample"] <-substr(df[df$cohort=="SH",sample_var],1,9) 
  df[df$cohort=="AO","Sample"] <-substr(df[df$cohort=="AO",sample_var],1,8) 
  df[df$cohort=="DO","Sample"] <-gsub('T','',df[df$cohort=="DO",sample_var]) 
  df[df$cohort=="DG","Sample"] <-gsub('T','',df[df$cohort=="DG",sample_var]) 
  df[df$cohort=="DA","Sample"] <-gsub('T','',df[df$cohort=="DA",sample_var]) 
  df[df$cohort=="WG","Sample"] <-sapply(strsplit(df[df$cohort=="WG",sample_var],'-'), function(x) paste(x[1:3],collapse='-')) 
  
  return(df)
}
cnvkit <- stand_ids(cnvkit,"V4")
colnames(cnvkit)[7]<- "Length"
cnvkit_dels <- cnvkit[cnvkit[,"V5"]<2,]
cnvkit_dups <- cnvkit[cnvkit[,"V5"]>2,]

# Compare with original Purple calls
purple <- read.table("/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/Combined_purple_ids_noSVs.bed",sep="\t")
purple[purple[,1]=="X",1]<- 23
purple <- stand_ids(purple,"V17")
purple$Length <- purple$V3 - purple$V2+1
purple <- purple[purple$V4 <1.5 | purple$V4 >2.5,]
purple_dels <- purple[purple[,"V4"]<2,]
purple_dups <- purple[purple[,"V4"]>2,]


#Compare with original CLImAT calls
climat <- read.table("/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/Combined_CLImAT_ids.bed",sep="\t")
climat <- stand_ids(climat,"V7")
climat <- climat[!duplicated(climat),]
climat$Length <- climat$V3 - climat$V2 +1
climat_dels <- climat[climat[,"V4"]<2,]
climat_dups <- climat[climat[,"V4"]>2,]

#Compare with original PCAWG calls
pcawg <- read.table("/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/Combined_pcawg_ids.bed",sep="\t")
pcawg <- stand_ids(pcawg,"V5")
pcawg$Length <- pcawg$V3 - pcawg$V2 +1
pcawg_dels <- pcawg[which(pcawg[,"V4"]<2),]
pcawg_dups <- pcawg[which(pcawg[,"V4"]>2),]


#### Compare Lengths between callers #######

callers_dels <- c(rep("CLImAT",dim(climat_dels)[1]),rep("CNVkit",dim(cnvkit_dels)[1]),rep("Purple",dim(purple_dels)[1]),rep("PCAWG",dim(pcawg_dels)[1]))
lengths_dels <- c(climat_dels$Length, cnvkit_dels$Length, purple_dels$Length,pcawg_dels$Length)
df_dels <- data.frame(Caller=callers_dels, Length = lengths_dels)

del_sizes <-ggplot(df_dels, aes(x=log10(Length),col=Caller)) + geom_density() +theme_minimal() +
  xlab("Length in bp (log10)") + ylab("Density")+ggtitle("Deletions")



callers_dups <- c(rep("CLImAT",dim(climat_dups)[1]),rep("CNVkit",dim(cnvkit_dups)[1]),rep("Purple",dim(purple_dups)[1]),rep("PCAWG",dim(pcawg_dups)[1]))
lengths_dups <- c(climat_dups$Length, cnvkit_dups$Length, purple_dups$Length,pcawg_dups$Length)
df_dups <- data.frame(Caller=callers_dups, Length = lengths_dups)

dup_sizes <- ggplot(df_dups, aes(x=log10(Length),col=Caller)) + geom_density() +theme_minimal() +
  xlab("Length in bp (log10)") + ylab("Density")+ggtitle("Duplications")



library(patchwork)

del_sizes + dup_sizes 

####################################################################################################################################

#Read in variants, remove dups,  add sizes
aggregate_consensus_cnvs <- function(overlap,caller){
  dups <-data.frame()
  dels <-data.frame()
  # Load consensus CNV calls 
  files <- list.files("/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/consensusCNVs_persample",paste(caller,"_",overlap,"_",sep=""),full.names = T)
  
    for  (i in 1:length(files)){
      print(files[i])
      file <- read.table(files[i])
      
      # Remove duplicates if there are any
      file_unique <- file[!duplicated(file[,1:4]),]
  
      # Insert fields for cnvkit, climat sizes and overlap threshold (switch for purple-climat,purple-pcawg, climat-pcawg)
      climat_size <- file_unique$V10 - file_unique$V9
      overlap_threshold <- rep(overlap,dim(file_unique)[1])
      #file_withsizes <- cbind(file_unique,climat_size,overlap_threshold)
      file_withsizes <- cbind(file_unique,overlap_threshold)
      
      # Split into dups/dels - cnvkit
      #dups <- rbind(dups,file_withsizes[which(file_withsizes$V5 > 2),])
      #dels <- rbind(dels,file_withsizes[which(file_withsizes$V5 < 2),])
      
      # Split into dups/dels - purple
      Sample<-gsub('-climat-pcawg_0.5_all_lengths.bed','',sapply(strsplit(files[i],'/'), function(x) x[[7]]))
      file_withsizes$Sample <- rep(Sample, dim(file_withsizes)[1])
      dups <- rbind(dups,file_withsizes[which(file_withsizes$V4 > 2),])
      dels <- rbind(dels,file_withsizes[which(file_withsizes$V4 < 2),])
    }
  output <- list(Dups = dups,Dels =dels)
  return(output)
}

consensusCNVs_0.5 <- aggregate_consensus_cnvs(0.5,"purple")
consensusCNVs_0.5_climat <- aggregate_consensus_cnvs(0.5,"CLImAT")
consensusCNVs_0.5_purple_climat <- aggregate_consensus_cnvs(0.5,"purple-climat")

consensusCNVs_0.5_cnvkit_PCAWG <- aggregate_consensus_cnvs(0.5,"PCAWG")
consensusCNVs_0.5_purple_PCAWG <- aggregate_consensus_cnvs(0.5,"purple-pcawg")
consensusCNVs_0.5_climat_PCAWG <- aggregate_consensus_cnvs(0.5,"climat-pcawg")

#Get PCAWG ids

files <-list.files("/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/consensusCNVs_persample","PCAWG_0.5_",full.names = T)

step1<- sapply(strsplit(files,"/"), function(x) x[[7]])
pcawg_ids <- sapply(strsplit(step1,"-"), function(x) x[[1]])


#Save files
save(cnvkit,file="ConsensusCNVs/cnvkit.RData")
save(purple,file="ConsensusCNVs/purple.RData")
save(consensusCNVs_0.5,file="ConsensusCNVs/cnvkit-purple.RData")
save(consensusCNVs_0.5_climat,file="ConsensusCNVs/cnvkit-climat.RData")
save(consensusCNVs_0.5_purple_climat,file="ConsensusCNVs/purple-climat.RData")

########################################################################################################################

#Create one file with which callers support calls (stick with 0.5 threshold as little difference in dist across samples)

##CNVKIT alone##

#Add ids to cnvkit calls
cnvkit_dels$mergeid <- paste(cnvkit_dels$V1,cnvkit_dels$V2,cnvkit_dels$V3,cnvkit_dels$Sample,sep="&")
cnvkit_dups$mergeid <- paste(cnvkit_dups$V1,cnvkit_dups$V2,cnvkit_dups$V3,cnvkit_dups$Sample,sep="&")

## CNVKIT - CLIMAT

#Split cnvkit-climat calls until dels and dups, add caller field and mergingid
climat_dels_0.5 <- consensusCNVs_0.5_climat[["Dels"]]
climat_dups_0.5 <- consensusCNVs_0.5_climat[["Dups"]]   

climat_dels_0.5 <- stand_ids(climat_dels_0.5,"V4")
climat_dups_0.5 <- stand_ids(climat_dups_0.5,"V4")

climat_dels_0.5$mergeid <- paste(climat_dels_0.5$V1,climat_dels_0.5$V2,climat_dels_0.5$V3,climat_dels_0.5$Sample,sep="&")
climat_dels_0.5$Caller <- "CLImAT"
climat_dups_0.5$mergeid <- paste(climat_dups_0.5$V1,climat_dups_0.5$V2,climat_dups_0.5$V3,climat_dups_0.5$Sample,sep="&")
climat_dups_0.5$Caller <- "CLImAT"

climat_dels_0.5$climat_mergeid <- paste(climat_dels_0.5$V8,climat_dels_0.5$V9,climat_dels_0.5$V10,climat_dels_0.5$Sample,sep="&")
climat_dups_0.5$climat_mergeid <- paste(climat_dups_0.5$V8,climat_dups_0.5$V9,climat_dups_0.5$V10,climat_dups_0.5$Sample,sep="&")

#Throw away fields to annotate orig cnvcalls from cnvkit with which ones were supported by climat
climat_dels_0.5_formerge <- climat_dels_0.5[,c("mergeid","Caller","climat_mergeid")]
climat_dups_0.5_formerge <- climat_dups_0.5[,c("mergeid","Caller","climat_mergeid")]

#Throw away fields to annotate orig climat calls with which ones were supported by cnvkit

climat_dels_0.5_formergeclimat <- climat_dels_0.5[,c("mergeid","Caller","climat_mergeid")]
climat_dels_0.5_formergeclimat$Caller <- "CNVkit"
climat_dups_0.5_formergeclimat <- climat_dups_0.5[,c("mergeid","Caller","climat_mergeid")]
climat_dups_0.5_formergeclimat$Caller <- "CNVkit"

##CNVKIT - PURPLE

#Split into dels and dups and add caller
purple_dels_0.5 <- consensusCNVs_0.5[["Dels"]]
purple_dels_0.5$Caller <- "Purple"
purple_dups_0.5 <- consensusCNVs_0.5[["Dups"]]   
purple_dups_0.5$Caller <- "Purple"

purple_dels_0.5<- stand_ids(purple_dels_0.5,"V4")
purple_dups_0.5<- stand_ids(purple_dups_0.5,"V4")

#Add mergeid
purple_dels_0.5$mergeid <- paste(purple_dels_0.5$V1,purple_dels_0.5$V2,purple_dels_0.5$V3,purple_dels_0.5$Sample,sep="&")
purple_dups_0.5$mergeid <- paste(purple_dups_0.5$V1,purple_dups_0.5$V2,purple_dups_0.5$V3,purple_dups_0.5$Sample,sep="&")

#Simplify sample field - dels then dups
purple_dels_0.5$purple_mergeid <- paste(purple_dels_0.5$V8,purple_dels_0.5$V9,purple_dels_0.5$V10,purple_dels_0.5$Sample,sep="&")
purple_dups_0.5$purple_mergeid <- paste(purple_dups_0.5$V8,purple_dups_0.5$V9,purple_dups_0.5$V10,purple_dups_0.5$Sample,sep="&")


#Throw away fields to annotate orig cnvcalls from cnvkit with which ones were supported by purple
purple_dels_0.5_formerge <- purple_dels_0.5[,c("mergeid","Caller","purple_mergeid")]
purple_dups_0.5_formerge <- purple_dups_0.5[,c("mergeid","Caller","purple_mergeid")]

purple_dels_0.5_withcnvkit<- purple_dels_0.5[,c("purple_mergeid","Caller","mergeid")]
purple_dups_0.5_withcnvkit <- purple_dups_0.5[,c("purple_mergeid","Caller","mergeid")]
purple_dels_0.5_withcnvkit$Caller <- "CNVkit" 
purple_dups_0.5_withcnvkit$Caller <- "CNVkit" 

## PURPLE alone ##

#Add merge ids to original purple calls
purple_dels$purple_mergeid <-paste(purple_dels$V1,purple_dels$V2,purple_dels$V3,purple_dels$Sample,sep="&")
purple_dups$purple_mergeid <-paste(purple_dups$V1,purple_dups$V2,purple_dups$V3,purple_dups$Sample,sep="&")

## PURPLE - CLIMAT##

consensus_DELs_purple_climat <- consensusCNVs_0.5_purple_climat[["Dels"]]
consensus_DUPs_purple_climat <- consensusCNVs_0.5_purple_climat[["Dups"]]
consensus_DELs_purple_climat$origSample <- consensus_DELs_purple_climat$Sample
consensus_DUPs_purple_climat$origSample <- consensus_DUPs_purple_climat$Sample


consensus_DELs_purple_climat <- stand_ids(consensus_DELs_purple_climat,"origSample")
consensus_DUPs_purple_climat <- stand_ids(consensus_DUPs_purple_climat,"origSample")

#Add mergeid to merge with orig purple calls
consensus_DELs_purple_climat$purple_mergeid <- paste(consensus_DELs_purple_climat$V7,consensus_DELs_purple_climat$V8,consensus_DELs_purple_climat$V9,consensus_DELs_purple_climat$Sample,sep="&")
consensus_DUPs_purple_climat$purple_mergeid <- paste(consensus_DUPs_purple_climat$V7,consensus_DUPs_purple_climat$V8,consensus_DUPs_purple_climat$V9,consensus_DUPs_purple_climat$Sample,sep="&")

consensus_DELs_purple_climat$climat_mergeid <- paste(consensus_DELs_purple_climat$V1,consensus_DELs_purple_climat$V2,consensus_DELs_purple_climat$V3,consensus_DELs_purple_climat$Sample,sep="&")
consensus_DUPs_purple_climat$climat_mergeid <- paste(consensus_DUPs_purple_climat$V1,consensus_DUPs_purple_climat$V2,consensus_DUPs_purple_climat$V3,consensus_DUPs_purple_climat$Sample,sep="&")

#Throw away fields to annotate orig cnvcalls from purple with which ones were supported by climat
consensus_DELs_purple_climat$Caller <- "Purple"
purple_climat_dels_formerge <- consensus_DELs_purple_climat[,c("purple_mergeid","Caller","climat_mergeid","Sample")]

consensus_DUPs_purple_climat$Caller <- "Purple"
purple_climat_dups_formerge <- consensus_DUPs_purple_climat[,c("purple_mergeid","Caller","climat_mergeid","Sample")]


#Throw away fields to annotate orig climat calls with which ones were supported by purple

consensus_DELs_purple_climat_formergeclimat <- consensus_DELs_purple_climat[,c("purple_mergeid","Caller","climat_mergeid")]
consensus_DELs_purple_climat_formergeclimat$Caller <- "CLImAT"
consensus_DUPs_purple_climat_formergeclimat <- consensus_DUPs_purple_climat[,c("purple_mergeid","Caller","climat_mergeid")]
consensus_DUPs_purple_climat_formergeclimat$Caller <- "CLImAT"

## CLIMAT alone ##
#Add merge ids to original purple calls
climat_dels$climat_mergeid <-paste(climat_dels$V1,climat_dels$V2,climat_dels$V3,climat_dels$Sample,sep="&")
climat_dups$climat_mergeid <-paste(climat_dups$V1,climat_dups$V2,climat_dups$V3,climat_dups$Sample,sep="&")

## Perform merging - dels

## Annotate CNVkit
# cnvkit alone, cnvkit-climat, cnvkit_purple
# dels
merge_dat1 <- merge(cnvkit_dels, climat_dels_0.5_formerge,by="mergeid",all.x=T)
cnvkit_dels_withconsensus <- merge(merge_dat1, purple_dels_0.5_formerge,by="mergeid",all.x=T)
cnvkit_dels_withconsensus$Callers <- paste("CNVkit",cnvkit_dels_withconsensus$Caller.x, cnvkit_dels_withconsensus$Caller.y,sep=",")
#cnvkit_dels_withconsensus <- cnvkit_dels_withconsensus[,c("mergeid","V1","V2","V3","V4","V5","V6","V7","Callers")]

#find extra climat ids not rep'd by cnvkit calls
climat_toincl<- setdiff(climat_dels$climat_mergeid,unique(cnvkit_dels_withconsensus$climat_mergeid))
rownames(climat_dels) <- climat_dels$climat_mergeid
climat_dels_incl <- climat_dels[climat_toincl,]

#merge extra climat with purple
climat_dels_withconsensus <- merge(climat_dels_incl, purple_climat_dels_formerge,by="climat_mergeid",all.x=T)
climat_dels_withconsensus$Callers <- paste("CLImAT",climat_dels_withconsensus$Caller,sep=",")

#Extras####
all_climat_dels_withconsensus <- merge(climat_dels, purple_climat_dels_formerge,by="climat_mergeid",all.x=T)
all_climat_dels_withconsensus$Callers <- paste("CLImAT",all_climat_dels_withconsensus$Caller,sep=",")

#find extra climat calls not rep'd by purple or cnvkit
combineused_purple <-c(cnvkit_dels_withconsensus$purple_mergeid,climat_dels_withconsensus$purple_mergeid)
purple_toincl <- setdiff(purple_dels$purple_mergeid, combineused_purple)
rownames(purple_dels) <- purple_dels$purple_mergeid
purple_dels_incl <- purple_dels[purple_toincl,]
purple_dels_incl$Callers <- "Purple"

#Stick 3 dfs together

cnvkit_dels_withconsensus$cnvkit_mergeid <- cnvkit_dels_withconsensus$mergeid
cnvkit_dels_withconsensus$Chr <- cnvkit_dels_withconsensus$V1
cnvkit_dels_withconsensus$Start <- cnvkit_dels_withconsensus$V2
cnvkit_dels_withconsensus$End <- cnvkit_dels_withconsensus$V3
cnvkit_dels_withconsensus$CN <- cnvkit_dels_withconsensus$V5
cnvkit_dels_withconsensus<- cnvkit_dels_withconsensus[,c("Sample","Callers","Chr","Start","End","CN","cnvkit_mergeid","climat_mergeid","purple_mergeid")]

climat_dels_withconsensus$cnvkit_mergeid <- rep(NA,dim(climat_dels_withconsensus)[1])
climat_dels_withconsensus$Chr <- climat_dels_withconsensus$V1
climat_dels_withconsensus$Start <- climat_dels_withconsensus$V2
climat_dels_withconsensus$End <- climat_dels_withconsensus$V3
climat_dels_withconsensus$CN <- climat_dels_withconsensus$V4
climat_dels_withconsensus$Sample <- climat_dels_withconsensus$Sample.x
climat_dels_withconsensus<- climat_dels_withconsensus[,c("Sample","Callers","Chr","Start","End","CN","cnvkit_mergeid","climat_mergeid","purple_mergeid")]

climat_dels_withconsensus[1,]

purple_dels_incl$cnvkit_mergeid <- rep(NA,dim(purple_dels_incl)[1])
purple_dels_incl$climat_mergeid <- rep(NA,dim(purple_dels_incl)[1])
purple_dels_incl$Chr <- purple_dels_incl$V1
purple_dels_incl$Start <- purple_dels_incl$V2
purple_dels_incl$End <- purple_dels_incl$V3
purple_dels_incl$CN <- purple_dels_incl$V4
purple_dels_incl<- purple_dels_incl[,c("Sample","Callers","Chr","Start","End","CN","cnvkit_mergeid","climat_mergeid","purple_mergeid")]
purple_dels_incl[1,]


all_dels <- rbind(cnvkit_dels_withconsensus,climat_dels_withconsensus,purple_dels_incl)
all_dels_3way <- all_dels[all_dels$Callers =="CNVkit,CLImAT,Purple",c(3,4,5,6,1,2,7,8,9)]


write.table(all_dels,"/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/ConsensusDels_alllevels.txt",sep="\t",row.names=F,quote=F)
write.table(all_dels_3way,"/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/ConsensusDels_CnvkitClimatPurple_130422.bed",sep="\t",row.names=F,quote=F)

#stats
total_climat <- 22008+7561+43360+6048
all3_climat <- 43360/total_climat
cnvkit_climat <- 22008/total_climat
purple_climat <- 6048/total_climat
justclimat <- 7561/total_climat
climat_stats <- c(all3_climat,cnvkit_climat, purple_climat,justclimat)
climat_stats

total_purple <- 43360+121100+6048+50887
all3_purple <- 43360/total_purple
cnvkit_purple<- 121100/total_purple
purple_climat <- 6048/total_purple
justpurple <- 50887/total_purple
purple_stats <- c(all3_purple,cnvkit_purple, purple_climat,justpurple)
purple_stats

total_cnvkit <- 43360+121100+22008+313552
all3_cnvkit <- 43360/total_cnvkit
cnvkit_purple<- 121100/total_cnvkit
cnvkit_climat <- 22008/total_cnvkit
justcnvkit <- 313552/total_cnvkit
cnvkit_stats <- c(all3_cnvkit,cnvkit_purple, cnvkit_climat,justcnvkit)
cnvkit_stats

## Perform merging - dups
## Annotate CNVkit
# cnvkit alone, cnvkit-climat, cnvkit_purple
# dups
merge_dat1 <- merge(cnvkit_dups, climat_dups_0.5_formerge,by="mergeid",all.x=T)
cnvkit_dups_withconsensus <- merge(merge_dat1, purple_dups_0.5_formerge,by="mergeid",all.x=T)
cnvkit_dups_withconsensus$Callers <- paste("CNVkit",cnvkit_dups_withconsensus$Caller.x, cnvkit_dups_withconsensus$Caller.y,sep=",")
#cnvkit_dels_withconsensus <- cnvkit_dels_withconsensus[,c("mergeid","V1","V2","V3","V4","V5","V6","V7","Callers")]

#find extra climat ids not rep'd by cnvkit calls
climat_toincl_dups<- setdiff(climat_dups$climat_mergeid,unique(cnvkit_dups_withconsensus$climat_mergeid))
rownames(climat_dups) <- climat_dups$climat_mergeid
climat_dups_incl <- climat_dups[climat_toincl_dups,]

#merge extra climat with purple
climat_dups_withconsensus <- merge(climat_dups_incl, purple_climat_dups_formerge,by="climat_mergeid",all.x=T)
climat_dups_withconsensus$Callers <- paste("CLImAT",climat_dups_withconsensus$Caller,sep=",")

#Extras####
all_climat_dups_withconsensus <- merge(climat_dups, purple_climat_dups_formerge,by="climat_mergeid",all.x=T)
all_climat_dups_withconsensus$Callers <- paste("CLImAT",all_climat_dups_withconsensus$Caller,sep=",")

#find extra climat calls not rep'd by purple or cnvkit
combineused_purple_dups <-c(cnvkit_dups_withconsensus$purple_mergeid,climat_dups_withconsensus$purple_mergeid)
purple_toincl_dups <- setdiff(purple_dups$purple_mergeid, combineused_purple_dups)
rownames(purple_dups) <- purple_dups$purple_mergeid
purple_dups_incl <- purple_dups[purple_toincl_dups,]
purple_dups_incl$Callers <- "Purple"

#Stick 3 dfs together

cnvkit_dups_withconsensus$cnvkit_mergeid <- cnvkit_dups_withconsensus$mergeid
cnvkit_dups_withconsensus$Chr <- cnvkit_dups_withconsensus$V1
cnvkit_dups_withconsensus$Start <- cnvkit_dups_withconsensus$V2
cnvkit_dups_withconsensus$End <- cnvkit_dups_withconsensus$V3
cnvkit_dups_withconsensus$CN <- cnvkit_dups_withconsensus$V5
cnvkit_dups_withconsensus<- cnvkit_dups_withconsensus[,c("Sample","Callers","Chr","Start","End","CN","cnvkit_mergeid","climat_mergeid","purple_mergeid")]

climat_dups_withconsensus$cnvkit_mergeid <- rep(NA,dim(climat_dups_withconsensus)[1])
climat_dups_withconsensus$Chr <- climat_dups_withconsensus$V1
climat_dups_withconsensus$Start <- climat_dups_withconsensus$V2
climat_dups_withconsensus$End <- climat_dups_withconsensus$V3
climat_dups_withconsensus$CN <- climat_dups_withconsensus$V4
climat_dups_withconsensus$Sample <- climat_dups_withconsensus$Sample.x
climat_dups_withconsensus<- climat_dups_withconsensus[,c("Sample","Callers","Chr","Start","End","CN","cnvkit_mergeid","climat_mergeid","purple_mergeid")]

climat_dups_withconsensus[1,]

purple_dups_incl$cnvkit_mergeid <- rep(NA,dim(purple_dups_incl)[1])
purple_dups_incl$climat_mergeid <- rep(NA,dim(purple_dups_incl)[1])
purple_dups_incl$Chr <- purple_dups_incl$V1
purple_dups_incl$Start <- purple_dups_incl$V2
purple_dups_incl$End <- purple_dups_incl$V3
purple_dups_incl$CN <- purple_dups_incl$V4
purple_dups_incl<- purple_dups_incl[,c("Sample","Callers","Chr","Start","End","CN","cnvkit_mergeid","climat_mergeid","purple_mergeid")]
purple_dups_incl[1,]


all_dups <- rbind(cnvkit_dups_withconsensus,climat_dups_withconsensus,purple_dups_incl)
table(all_dups$Callers)
table(all_dups$Callers)/sum(table(all_dups$Callers))

all_dups_3way <- all_dups[all_dups$Callers =="CNVkit,CLImAT,Purple",c(3,4,5,6,1,2,7,8,9)]


write.table(all_dups,"/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/ConsensusDups_alllevels.txt",sep="\t",row.names=F,quote=F)
write.table(all_dups_3way,"/Volumes/igmm/HGS-OvarianCancerA-SGP-WGS/consensusSVs/ConsensusDups_CnvkitClimatPurple_130422.bed",sep="\t",row.names=F,quote=F)

## Compare with PCAWG


cnvkit_pcawg_dels_0.5 <- consensusCNVs_0.5_cnvkit_PCAWG[["Dels"]]
cnvkit_pcawg_dups_0.5 <- consensusCNVs_0.5_cnvkit_PCAWG[["Dups"]]   

cnvkit_pcawg_dels_0.5 <- stand_ids(cnvkit_pcawg_dels_0.5,"V4")
cnvkit_pcawg_dups_0.5 <- stand_ids(cnvkit_pcawg_dups_0.5,"V4")

cnvkit_pcawg_dels_0.5$cnvkit_mergeid <- paste(cnvkit_pcawg_dels_0.5$V1,cnvkit_pcawg_dels_0.5$V2,cnvkit_pcawg_dels_0.5$V3,cnvkit_pcawg_dels_0.5$Sample,sep="&")
cnvkit_pcawg_dups_0.5$cnvkit_mergeid <- paste(cnvkit_pcawg_dups_0.5$V1,cnvkit_pcawg_dups_0.5$V2,cnvkit_pcawg_dups_0.5$V3,cnvkit_pcawg_dups_0.5$Sample,sep="&")

cnvkit_pcawg_dels_0.5$PCAWG <- "YES"
cnvkit_pcawg_dups_0.5$PCAWG <- "YES"

cnvkit_pcawg_dels_0.5_formerge <- cnvkit_pcawg_dels_0.5[,c("cnvkit_mergeid","PCAWG")]
cnv_dels_with_pcawg <- merge(cnvkit_dels_withconsensus, cnvkit_pcawg_dels_0.5_formerge, by="cnvkit_mergeid",all.x=T)

cnvkit_pcawg_dups_0.5_formerge <- cnvkit_pcawg_dups_0.5[,c("cnvkit_mergeid","PCAWG")]
cnv_dups_with_pcawg <- merge(cnvkit_dups_withconsensus, cnvkit_pcawg_dups_0.5_formerge, by="cnvkit_mergeid",all.x=T)

climat_pcawg_dels_0.5 <- consensusCNVs_0.5_climat_PCAWG[["Dels"]]
climat_pcawg_dups_0.5 <- consensusCNVs_0.5_climat_PCAWG[["Dups"]]   

climat_pcawg_dels_0.5$origSample <- climat_pcawg_dels_0.5$Sample
climat_pcawg_dups_0.5$origSample <- climat_pcawg_dups_0.5$Sample

climat_pcawg_dels_0.5 <- stand_ids(climat_pcawg_dels_0.5,"origSample")
climat_pcawg_dups_0.5 <- stand_ids(climat_pcawg_dups_0.5,"origSample")

climat_pcawg_dels_0.5$climat_mergeid <- paste(climat_pcawg_dels_0.5$V1,climat_pcawg_dels_0.5$V2,climat_pcawg_dels_0.5$V3,climat_pcawg_dels_0.5$Sample,sep="&")
climat_pcawg_dups_0.5$climat_mergeid <- paste(climat_pcawg_dups_0.5$V1,climat_pcawg_dups_0.5$V2,climat_pcawg_dups_0.5$V3,climat_pcawg_dups_0.5$Sample,sep="&")

climat_pcawg_dels_0.5$PCAWG <- "YES"
climat_pcawg_dups_0.5$PCAWG <- "YES"

climat_pcawg_dels_0.5_formerge <- climat_pcawg_dels_0.5[,c("climat_mergeid","PCAWG")]
climat_dels_with_pcawg <- merge(climat_dels_withconsensus, climat_pcawg_dels_0.5_formerge, by="climat_mergeid",all.x=T)
purple_climat_dels_with_pcawg <- merge(purple_climat_dels_formerge, climat_pcawg_dels_0.5_formerge, by="climat_mergeid",all.x=T)
purple_climat_dels_with_pcawg[purple_climat_dels_with_pcawg$Sample %in% pcawg_ids & is.na(purple_climat_dels_with_pcawg$PCAWG), "PCAWG"] <- "NO"

climat_pcawg_dups_0.5_formerge <- climat_pcawg_dups_0.5[,c("climat_mergeid","PCAWG")]
climat_dups_with_pcawg <- merge(climat_dups_withconsensus, climat_pcawg_dups_0.5_formerge, by="climat_mergeid",all.x=T)
purple_climat_dups_with_pcawg <- merge(purple_climat_dups_formerge, climat_pcawg_dups_0.5_formerge, by="climat_mergeid",all.x=T)
purple_climat_dups_with_pcawg[purple_climat_dups_with_pcawg$Sample %in% pcawg_ids & is.na(purple_climat_dups_with_pcawg$PCAWG), "PCAWG"] <- "NO"


purple_pcawg_dels_0.5 <- consensusCNVs_0.5_purple_PCAWG[["Dels"]]
purple_pcawg_dups_0.5 <- consensusCNVs_0.5_purple_PCAWG[["Dups"]]   

purple_pcawg_dels_0.5$origSample <- purple_pcawg_dels_0.5$Sample
purple_pcawg_dups_0.5$origSample <- purple_pcawg_dups_0.5$Sample

purple_pcawg_dels_0.5 <- stand_ids(purple_pcawg_dels_0.5,"origSample")
purple_pcawg_dups_0.5 <- stand_ids(purple_pcawg_dups_0.5,"origSample")

purple_pcawg_dels_0.5$purple_mergeid <- paste(purple_pcawg_dels_0.5$V1,purple_pcawg_dels_0.5$V2,purple_pcawg_dels_0.5$V3,purple_pcawg_dels_0.5$Sample,sep="&")
purple_pcawg_dups_0.5$purple_mergeid <- paste(purple_pcawg_dups_0.5$V1,purple_pcawg_dups_0.5$V2,purple_pcawg_dups_0.5$V3,purple_pcawg_dups_0.5$Sample,sep="&")

purple_pcawg_dels_0.5$PCAWG <- "YES"
purple_pcawg_dups_0.5$PCAWG <- "YES"

purple_pcawg_dels_0.5_formerge <- purple_pcawg_dels_0.5[,c("purple_mergeid","PCAWG")]
purple_dels_with_pcawg <- merge(purple_dels_incl, purple_pcawg_dels_0.5_formerge, by="purple_mergeid",all.x=T)

purple_pcawg_dups_0.5_formerge <- purple_pcawg_dups_0.5[,c("purple_mergeid","PCAWG")]
purple_dups_with_pcawg <- merge(purple_dups_incl, purple_pcawg_dups_0.5_formerge, by="purple_mergeid",all.x=T)


all_dels_with_pcawg <- rbind(cnv_dels_with_pcawg,purple_dels_with_pcawg,climat_dels_with_pcawg)
all_dels_with_pcawg[all_dels_with_pcawg$Sample %in% pcawg_ids & is.na(all_dels_with_pcawg$PCAWG), "PCAWG"] <- "NO"

all_dups_with_pcawg <- rbind(cnv_dups_with_pcawg,purple_dups_with_pcawg,climat_dups_with_pcawg)
all_dups_with_pcawg[all_dups_with_pcawg$Sample %in% pcawg_ids & is.na(all_dups_with_pcawg$PCAWG), "PCAWG"] <- "NO"

#Individual caller concordance with PCAWG
all_climat_dels_with_pcawg <- merge(climat_dels, climat_pcawg_dels_0.5_formerge, by="climat_mergeid",all.x=T)
all_climat_dels_with_pcawg[all_climat_dels_with_pcawg$Sample %in% pcawg_ids & is.na(all_climat_dels_with_pcawg$PCAWG), "PCAWG"] <- "NO"
table(all_climat_dels_with_pcawg$PCAWG)/sum(table(all_climat_dels_with_pcawg$PCAWG))

all_climat_dups_with_pcawg <- merge(climat_dups, climat_pcawg_dups_0.5_formerge, by="climat_mergeid",all.x=T)
all_climat_dups_with_pcawg[all_climat_dups_with_pcawg$Sample %in% pcawg_ids & is.na(all_climat_dups_with_pcawg$PCAWG), "PCAWG"] <- "NO"
table(all_climat_dups_with_pcawg$PCAWG)

cnvkit_dels$cnvkit_mergeid <- cnvkit_dels$mergeid
cnvkit_dups$cnvkit_mergeid <- cnvkit_dups$mergeid

all_cnvkit_dels_with_pcawg <- merge(cnvkit_dels, cnvkit_pcawg_dels_0.5_formerge, by="cnvkit_mergeid",all.x=T)
all_cnvkit_dels_with_pcawg[all_cnvkit_dels_with_pcawg$Sample %in% pcawg_ids & is.na(all_cnvkit_dels_with_pcawg$PCAWG), "PCAWG"] <- "NO"
table(all_cnvkit_dels_with_pcawg$PCAWG)

all_cnvkit_dups_with_pcawg <- merge(cnvkit_dups, cnvkit_pcawg_dups_0.5_formerge, by="cnvkit_mergeid",all.x=T)
all_cnvkit_dups_with_pcawg[all_cnvkit_dups_with_pcawg$Sample %in% pcawg_ids & is.na(all_cnvkit_dups_with_pcawg$PCAWG), "PCAWG"] <- "NO"
table(all_cnvkit_dups_with_pcawg$PCAWG)

all_purple_dels_with_pcawg <- merge(purple_dels, purple_pcawg_dels_0.5_formerge, by="purple_mergeid",all.x=T)
all_purple_dels_with_pcawg[all_purple_dels_with_pcawg$Sample %in% pcawg_ids & is.na(all_purple_dels_with_pcawg$PCAWG), "PCAWG"] <- "NO"
table(all_purple_dels_with_pcawg$PCAWG)/sum(table(all_purple_dels_with_pcawg$PCAWG))

all_purple_dups_with_pcawg <- merge(purple_dups, purple_pcawg_dups_0.5_formerge, by="purple_mergeid",all.x=T)
all_purple_dups_with_pcawg[all_purple_dups_with_pcawg$Sample %in% pcawg_ids & is.na(all_purple_dups_with_pcawg$PCAWG), "PCAWG"] <- "NO"
table(all_purple_dups_with_pcawg$PCAWG)/sum(table(all_purple_dups_with_pcawg$PCAWG))

cs <- colSums(table(all_dels_with_pcawg$PCAWG, all_dels_with_pcawg$Callers))
df_props <-mat.or.vec(nr=2,nc=7)
for (i in 1:7){
 df_props[,i] <- table(all_dels_with_pcawg$PCAWG, all_dels_with_pcawg$Callers)[,i] /cs[i]
}
colnames(df_props) <- names(cs)
df_props

cs <- colSums(table(all_dups_with_pcawg$PCAWG, all_dups_with_pcawg$Callers))
df_propsdups <-mat.or.vec(nr=2,nc=7)
for (i in 1:7){
  df_propsdups[,i] <- table(all_dups_with_pcawg$PCAWG, all_dups_with_pcawg$Callers)[,i] /cs[i]
}
colnames(df_propsdups) <- names(cs)
df_propsdups

############# Categories by length #############################
all_dels$Length <- all_dels$End - all_dels$Start +1
all_dels$Callers <- factor(all_dels$Callers, levels=c("CNVkit,CLImAT,Purple","CNVkit,NA,Purple","CNVkit,CLImAT,NA","CNVkit,NA,NA","CLImAT,Purple","CLImAT,NA","Purple"))
ggplot(all_dels, aes(x=Callers, y=log10(Length))) + geom_boxplot(outlier.shape = NA) + theme_bw()



all_dups$Length <- all_dups$End - all_dups$Start +1
all_dups$Callers <- factor(all_dups$Callers, levels=c("CNVkit,CLImAT,Purple","CNVkit,NA,Purple","CNVkit,CLImAT,NA","CNVkit,NA,NA","CLImAT,Purple","CLImAT,NA","Purple"))
ggplot(all_dups, aes(x=Callers, y=log10(Length))) + geom_boxplot(outlier.shape = NA) + theme_bw()


# Does size differ between deletions and duplications

all_dels$Type <- "Deletion"
all_dups$Type <- "Duplication"

df_type <- rbind(all_dels,all_dups)
ggplot(df_type, aes(x=Type, y=log10(Length))) + geom_boxplot(outlier.shape = NA) + theme_bw()+geom_jitter(size=1,width=0.05)



# What's the distribution of cnvs across samples by consensus category
library(reshape2)
t <- table(all_dels$Sample, all_dels$Callers)
t2 <- melt(t, id.vars=rownames(t))

ggplot(t2, aes(x=Var2,y=log10(value))) + geom_boxplot() + geom_jitter(width=0.05)

tdups <- table(all_dups$Sample, all_dups$Callers)
t2dups <- melt(tdups, id.vars=rownames(tdups))
ggplot(t2dups, aes(x=Var2,y=log10(value))) + geom_boxplot() + geom_jitter(width=0.05)

# Earlier code-----
#dups
merge_dat1_dups <- merge(cnvkit_dups, climat_dups_0.5_formerge,by="mergeid",all.x=T)
cnvkit_dups_withconsensus <- merge(merge_dat1_dups, purple_dups_0.5_formerge,by="mergeid",all.x=T)
cnvkit_dups_withconsensus$Callers <- paste("CNVkit",cnvkit_dups_withconsensus$Caller.x, cnvkit_dups_withconsensus$Caller.y,sep=",")
cnvkit_dups_withconsensus <- cnvkit_dups_withconsensus[,c("mergeid","V1","V2","V3","V4","V5","V6","V7","Callers")]

## Annotate Purple
#purple alone, purple - climat, purple - cnvkit
merge_dat1_purp_dels <- merge(purple_dels, purple_climat_dels_formerge,by="purple_mergeid",all.x=T)
purple_dels_withconsensus <- merge(merge_dat1_purp_dels, purple_dels_0.5_withcnvkit,by="purple_mergeid",all.x=T)
purple_dels_withconsensus$Callers <- paste("Purple",purple_dels_withconsensus$Caller.x, purple_dels_withconsensus$Caller.y,sep=",")
purple_dels_withconsensus <- purple_dels_withconsensus[,setdiff(colnames(purple_dels_withconsensus),"mergeid")]
purple_dels_withconsensus <- purple_dels_withconsensus[!duplicated(purple_dels_withconsensus),]

merge_dat1_purp_dels <- merge(purple_dels,  purple_dels_0.5_withcnvkit,by="purple_mergeid",all.x=T)
purple_dels_withconsensus <- merge(merge_dat1_purp_dels, purple_climat_dels_formerge,by="mergeid",all.x=T)
purple_dels_withconsensus$Callers <- paste("Purple",purple_dels_withconsensus$Caller.x, purple_dels_withconsensus$Caller.y,sep=",")
purple_dels_withconsensus <- purple_dels_withconsensus[,setdiff(colnames(purple_dels_withconsensus),"mergeid")]
purple_dels_withconsensus <- purple_dels_withconsensus[!duplicated(purple_dels_withconsensus),]

#purple alone, purple - climat, purple - cnvkit
merge_dat1_purp_dups <- merge(purple_dups, purple_climat_dups_formerge,by="purple_mergeid",all.x=T)
purple_dups_withconsensus <- merge(merge_dat1_purp_dups, purple_dups_0.5_withcnvkit,by="purple_mergeid",all.x=T)
purple_dups_withconsensus$Callers <- paste("Purple",purple_dups_withconsensus$Caller.x, purple_dups_withconsensus$Caller.y,sep=",")
purple_dups_withconsensus <- purple_dups_withconsensus[,setdiff(colnames(purple_dups_withconsensus),"mergeid")]
purple_dups_withconsensus <- purple_dups_withconsensus[!duplicated(purple_dups_withconsensus),]


## Annotate CLImAT
#climat alone, cnvkit - climat, purple - climat
merge_dat1_climat_dels <- merge(climat_dels, climat_dels_0.5_formergeclimat ,by="climat_mergeid",all.x=T)
climat_dels_withconsensus <- merge(merge_dat1_climat_dels, consensus_DELs_purple_climat_formergeclimat,by="climat_mergeid",all.x=T)
climat_dels_withconsensus$Callers <- paste("CLImAT",climat_dels_withconsensus$Caller.x, climat_dels_withconsensus$Caller.y,sep=",")
climat_dels_withconsensus <- climat_dels_withconsensus[,setdiff(colnames(climat_dels_withconsensus),"mergeid")]
climat_dels_withconsensus <- climat_dels_withconsensus[!duplicated(climat_dels_withconsensus),]

#climat alone, cnvkit - climat, purple - climat
merge_dat1_climat_dups <- merge(climat_dups, climat_dups_0.5_formergeclimat,by="climat_mergeid",all.x=T)
climat_dups_withconsensus <- merge(merge_dat1_climat_dups, consensus_DUPs_purple_climat_formergeclimat,by="climat_mergeid",all.x=T)
climat_dups_withconsensus$Callers <- paste("Purple",climat_dups_withconsensus$Caller.x, climat_dups_withconsensus$Caller.y,sep=",")
climat_dups_withconsensus <- climat_dups_withconsensus[,setdiff(colnames(climat_dups_withconsensus),"mergeid")]
climat_dups_withconsensus <- climat_dups_withconsensus[!duplicated(climat_dups_withconsensus),]

#Consensus stats
table(cnvkit_dels_withconsensus$Callers)/dim(cnvkit_dels)[1]
table(cnvkit_dups_withconsensus$Callers)/dim(cnvkit_dups)[1]

table(purple_dels_withconsensus$Callers)/dim(purple_dels)[1]
table(purple_dups_withconsensus$Callers)/dim(purple_dups)[1]



#save(cnvkit_dels,file="ConsensusCNVs/cnvkit_dels.RData")
#save(cnvkit_dups,file="ConsensusCNVs/cnvkit_dups.RData")
#save(climat_dels_0.5,file="ConsensusCNVs/consensus_climat_cnvkit_dels_0.5.RData")
#save(climat_dups_0.5,file="ConsensusCNVs/consensus_climat_cnvkit_dups_0.5.RData")
#save(purple_dels_0.5,file="ConsensusCNVs/consensus_purple_cnvkit_dels_0.5.RData")
#save(purple_dups_0.5,file="ConsensusCNVs/consensus_purple_cnvkit_dups_0.5.RData")
#save(cnvkit_dups,file="ConsensusCNVs/cnvkit_dups.RData")

#Checks
all3_climat <- climat_dels_withconsensus[climat_dels_withconsensus$Callers == "CLImAT,CNVkit,Purple",]
all3_cnvkit <- cnvkit_dels_withconsensus[cnvkit_dels_withconsensus$Callers == "CNVkit,CLImAT,Purple",]
all3_purple <- purple_dels_withconsensus[purple_dels_withconsensus$Callers == "Purple,CLImAT,CNVkit",]
cnvkit_dels_withconsensus[1,]

length(unique(all3_cnvkit$mergeid))
length(unique(all3_climat$mergeid))
length(unique(all3_purple$mergeid))

#Just 2

#cnvkit and purple

cnvkit_check <- cnvkit_dels_withconsensus[cnvkit_dels_withconsensus$Callers %in% c("CNVkit,CLImAT,Purple","CNVkit,NA,Purple"),]
purple_check <- purple_dels_withconsensus[purple_dels_withconsensus$Callers %in% c("Purple,CLImAT,CNVkit","Purple,NA,CNVkit"),]

length(unique(cnvkit_check$mergeid))
length(unique(purple_check$mergeid))

#164460
#30652

cnvkit_check <- cnvkit_dels_withconsensus[cnvkit_dels_withconsensus$Callers %in% c("CNVkit,CLImAT,Purple"),]
purple_check <- purple_dels_withconsensus[purple_dels_withconsensus$Callers %in% c("Purple,CLImAT,CNVkit"),]
cnvkit_check <- cnvkit_dels_withconsensus[cnvkit_dels_withconsensus$Callers %in% c("CNVkit,NA,Purple"),]
purple_check <- purple_dels_withconsensus[purple_dels_withconsensus$Callers %in% c("Purple,NA,CNVkit"),]

intersect(cnvkit_check$mergeid,)

length(unique(cnvkit_check$purple_mergeid))
length(unique(purple_check$purple_mergeid))

#164460



