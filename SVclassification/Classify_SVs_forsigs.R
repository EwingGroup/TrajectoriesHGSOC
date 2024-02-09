# SV signatures - classification of SVs

#Load consensus calls with gridss ids
dat_consensus <- read.table("~/Documents/HGSOC_marker_paper/MutSignatures/data/sv_sigs/Combined_consensus_100bp_withgridssid.bedpe",sep="\t")
dat_consensus[,13]<- gsub("_PrimaryTumour","",dat_consensus[,13])
dat_consensus[,13]<- gsub("T","",dat_consensus[,13])
dat_consensus[,1]<- gsub("chr","",dat_consensus[,1])
dat_consensus[,4]<- gsub("chr","",dat_consensus[,4])
dat_consensus[grep('SHGSOC',dat_consensus[,13]),13]<- substr(dat_consensus[grep('SHGSOC',dat_consensus[,13]),13],1,9)
dat_consensus[grep('WGS',dat_consensus[,13]),13]<- sapply( strsplit(dat_consensus[grep('WGS',dat_consensus[,13]),13],"-"), function(x) paste(x[1],x[2],x[3],sep="-"))

mergeid_consensus <- character()

for (i in 1:dim(dat_consensus)[1]){
  mergeid_consensus [i] <- paste(dat_consensus[ i, c(1:6,13) ],collapse = ":" )
}
dat_consensus<- cbind(dat_consensus,mergeid=mergeid_consensus)

#Load data_clusters and complex identity

dat <- read.table("~/Documents/HGSOC_marker_paper/MutSignatures/data/sv_sigs/Combined_consensus.100bp.bedpe.clusters_and_complex_incl.txt",sep="\t",header=T)
dat[grep('WGS',dat[,12]),12]<- sapply( strsplit(dat[grep('WGS',dat[,12]),12],"-"), function(x) paste(x[1],x[2],x[3],sep="-"))

mergeid <- character()

for (i in 1:dim(dat)[1]){
  mergeid[i] <- paste(dat[ i, c(1:6,12) ],collapse = ":" )
}

dat<- cbind(dat,mergeid)
dat_consensus <- dat_consensus[,c("mergeid","V14")]
dat_withgridss <- merge(dat,dat_consensus,by="mergeid")

#Merge foldback data
foldbacks <- read.table("~/Documents/HGSOC_marker_paper/MutSignatures/data/sv_sigs/Foldback_inversions_SHGSOC.txt",sep="\t",header=T)
foldbacks[grep('SHGSOC',foldbacks[,"SAMPLE"]),"SAMPLE"]<- substr(foldbacks[grep('SHGSOC',foldbacks[,'SAMPLE']),"SAMPLE"],1,9)
foldbacks[grep('WGS',foldbacks[,"SAMPLE"]),"SAMPLE"]<- sapply( strsplit(foldbacks[grep('WGS',foldbacks[,"SAMPLE"]),"SAMPLE"],"-"), function(x) paste(x[1],x[2],x[3],sep="-"))

#Get foldback rows
rows <- data.frame()
for (i in 1:dim(foldbacks)[1]){
  newrows<- dat_withgridss[grep(foldbacks[i,"EVENT"],dat_withgridss$V14),]
  newrows <- newrows[newrows$sample==foldbacks[i,"SAMPLE"],]
  if (dim(newrows)[1] < 1){
    print(foldbacks[i,])
  }
  rows <- rbind(rows,newrows)
}

#Annotate SVs with foldback INVs
rownames(dat_withgridss) <- dat_withgridss$mergeid
foldback_inv <- rep(NA,dim(dat_withgridss)[1])
foldback_inv[which(dat_withgridss$sv_type=="INV")]<- 0
dat_withgridss<- cbind(dat_withgridss,foldback_inv)
dat_withgridss[rows$mergeid,"foldback_inv"]<- 1

dat_withgridss<- dat_withgridss[dat_withgridss$sv_length>50,]
consensus_svs_complex_clust_foldbacks <- dat_withgridss[,c(2:23,25)]

clustered <- ifelse(consensus_svs_complex_clust_foldbacks$cluster_size>=3 , 1,0)
anycomplex<- rowSums(consensus_svs_complex_clust_foldbacks[,15:22])
incomplex <- ifelse(anycomplex>=1 , 1,0)
consensus_svs_complex_clust_foldbacks <- cbind(consensus_svs_complex_clust_foldbacks,clustered,incomplex)
#write.table(consensus_svs_complex_clust_foldbacks,"MutSignatures/data/sv_sigs/Combined_consensus.100bp.clusters_foldbackinv_complex.txt",sep="\t",row.names=F,quote=F)

#Classify SVs
whichComplex <- character()
for (i in 1:dim(consensus_svs_complex_clust_foldbacks)[1]){
whichComplex[i] <- paste(ifelse(consensus_svs_complex_clust_foldbacks[i,"Chromoplexy"],"Chromoplexy",""),
                      ifelse(consensus_svs_complex_clust_foldbacks[i,"bfb"],"BFB",""),
                      ifelse(consensus_svs_complex_clust_foldbacks[i,"Rigma"],"Rigma",""),
                      ifelse(consensus_svs_complex_clust_foldbacks[i,"Pyrgo"],"Pyrgo",""),
                      ifelse(consensus_svs_complex_clust_foldbacks[i,"Tyfonas"],"Tyfonas",""),
                      ifelse(consensus_svs_complex_clust_foldbacks[i,"Chromothripsis"],"Chromothripsis",""),
                      ifelse(consensus_svs_complex_clust_foldbacks[i,"ecDNA"],"ecDNA",""),
                      ifelse(consensus_svs_complex_clust_foldbacks[i,"SeismicAmplification"],"SeismicAmplification",""),
                      sep=":")
}
consensus_svs_complex_clust_foldbacks<- cbind(consensus_svs_complex_clust_foldbacks,whichComplex)
complexCat<- character(dim(consensus_svs_complex_clust_foldbacks)[1])
complexCat[which(consensus_svs_complex_clust_foldbacks$clustered==1 & 
                   consensus_svs_complex_clust_foldbacks$incomplex==0)]<-"clustered_not_complex"
complexCat[which(consensus_svs_complex_clust_foldbacks$clustered==0 & 
                   consensus_svs_complex_clust_foldbacks$incomplex==0)]<-"simple"
complexCat[which(consensus_svs_complex_clust_foldbacks$incomplex==1)]<-"uncertain_multiple_complex"
complexCat[which(consensus_svs_complex_clust_foldbacks$whichComplex==":::::Chromothripsis::") ]<-"chromothripsis"
complexCat[which(consensus_svs_complex_clust_foldbacks$whichComplex==":BFB::::::") ]<-"BFB"
complexCat[which(consensus_svs_complex_clust_foldbacks$whichComplex=="Chromoplexy:::::::") ]<-"chromoplexy"
complexCat[which(consensus_svs_complex_clust_foldbacks$whichComplex==":::Pyrgo::::") ]<-"pyrgo"
complexCat[which(consensus_svs_complex_clust_foldbacks$whichComplex=="::Rigma:::::") ]<-"rigma"
complexCat[which(consensus_svs_complex_clust_foldbacks$whichComplex=="::::Tyfonas:::") ]<-"tyfonas"
complexCat[which(consensus_svs_complex_clust_foldbacks$whichComplex=="::::::ecDNA:") ]<-"ecDNA"

consensus_svs_complex_clust_foldbacks<- cbind(consensus_svs_complex_clust_foldbacks,complexCat)

consensus_svs_complex_clust_foldbacks_class <- consensus_svs_complex_clust_foldbacks[consensus_svs_complex_clust_foldbacks$sv_type!="INS",]
#Classify simple
#Classify dels

simpleCat <- character(dim(consensus_svs_complex_clust_foldbacks_class )[1])

#Dels by size
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DEL" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 100]<- "< 100bp del"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DEL" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 1000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 100]<- "< 1kb del"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DEL" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 10000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 1000]<- "< 10kb del"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DEL" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 100000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 10000]<- "< 100kb del"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DEL" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 1000000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 100000]<- "< 1Mb del"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DEL" &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 1000000]<- ">= 1Mb del"

#DUPS
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DUP" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 1000]<- "< 1kb dup"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DUP" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 10000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 1000]<- "< 10kb dup"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DUP" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 100000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 10000]<- "< 100kb dup"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DUP" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 1000000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 100000]<- "< 1Mb dup"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DUP" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 10000000 &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 1000000]<- "< 10Mb dup"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="DUP" &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 10000000]<- ">= 10Mb dup"


#Inv
consensus_svs_complex_clust_foldbacks_class$event_info <- sapply( strsplit(consensus_svs_complex_clust_foldbacks_class$sv_id,":"), 
                          function(x) paste(x[1],x[2],sep=":"))
recipevents_inv <- consensus_svs_complex_clust_foldbacks_class[consensus_svs_complex_clust_foldbacks_class$sv_type=="INV" &
                                                           duplicated(consensus_svs_complex_clust_foldbacks_class$event_info),
                                                           "event_info"]
consensus_svs_complex_clust_foldbacks_class$reciprocal_inv <- "unbalanced"
consensus_svs_complex_clust_foldbacks_class[consensus_svs_complex_clust_foldbacks_class$event_info %in% recipevents_inv,
                                            "reciprocal_inv"] <- "reciprocal"
consensus_svs_complex_clust_foldbacks_class[which(consensus_svs_complex_clust_foldbacks_class$foldback_inv==1),
                                            "reciprocal_inv"] <- "foldback"

simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="INV" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 1000000 & 
            consensus_svs_complex_clust_foldbacks_class$reciprocal_inv=="reciprocal" ]<- "small reciprocal inv"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="INV" &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 1000000 & 
            consensus_svs_complex_clust_foldbacks_class$reciprocal_inv=="reciprocal" ]<- "large reciprocal inv"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="INV" &
            consensus_svs_complex_clust_foldbacks_class$sv_length < 1000000 & 
            consensus_svs_complex_clust_foldbacks_class$reciprocal_inv=="unbalanced" ]<- "small unbalanced inv"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="INV" &
            consensus_svs_complex_clust_foldbacks_class$sv_length >= 1000000 & 
            consensus_svs_complex_clust_foldbacks_class$reciprocal_inv=="unbalanced" ]<- "large unbalanced inv"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="INV" &
            consensus_svs_complex_clust_foldbacks_class$reciprocal_inv=="foldback" ]<- "foldback inv"

#TRA
recipevents_tra <- consensus_svs_complex_clust_foldbacks_class[consensus_svs_complex_clust_foldbacks_class$sv_type=="CTX" &
                                                               duplicated(consensus_svs_complex_clust_foldbacks_class$event_info),
                                                               "event_info"]
consensus_svs_complex_clust_foldbacks_class$reciprocal_tra <- "unbalanced tra"
consensus_svs_complex_clust_foldbacks_class[consensus_svs_complex_clust_foldbacks_class$event_info %in% recipevents_tra,"reciprocal_tra"] <- "reciprocal"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="CTX" &
            consensus_svs_complex_clust_foldbacks_class$reciprocal_tra=="reciprocal" ]<- "reciprocal tra"
simpleCat[consensus_svs_complex_clust_foldbacks_class$sv_type=="CTX" &
            consensus_svs_complex_clust_foldbacks_class$reciprocal_tra=="unbalanced tra" ]<- "unbalanced tra"

simpleCat[consensus_svs_complex_clust_foldbacks_class$complexCat!="simple"]<- "complex_or_clustered"

consensus_svs_complex_clust_foldbacks_simpleCat <- cbind(consensus_svs_complex_clust_foldbacks_class,simpleCat)

#Make signature category
sigCat <- consensus_svs_complex_clust_foldbacks_simpleCat$complexCat
sigCat[sigCat=="simple"] <- consensus_svs_complex_clust_foldbacks_simpleCat[sigCat=="simple","simpleCat"]
consensus_svs_complex_clust_foldbacks_simpleCat <- cbind(consensus_svs_complex_clust_foldbacks_simpleCat,sigCat)

consensus_svs_simple_complex_classes <- consensus_svs_complex_clust_foldbacks_simpleCat[,c(1:26,32)]

write.table(consensus_svs_simple_complex_classes,"MutSignatures/data/sv_sigs/Combined_consensus.100bp.clusters_foldbackinv_complex_withsignatureclasses.txt",sep="\t",row.names=F,quote=F)

#Make sigs counts matrix

consensus_svs_simple_complex_classes$sigCat <- gsub(" ","_",consensus_svs_simple_complex_classes$sigCat)
consensus_svs_simple_complex_classes$sigCat <- gsub("<","lessthan",consensus_svs_simple_complex_classes$sigCat)
consensus_svs_simple_complex_classes$sigCat <- gsub(">=","greaterequal",consensus_svs_simple_complex_classes$sigCat)
sv_sigs_counts_matrix_withcomplex<- table(consensus_svs_simple_complex_classes$sample,consensus_svs_simple_complex_classes$sigCat)

#Add WGS-ER-2 no SVs
extrrow <- rep(0,dim(sv_sigs_counts_matrix_withcomplex)[2])
sv_sigs_counts_matrix_withcomplex<- rbind(sv_sigs_counts_matrix_withcomplex,extrrow)
rownames(sv_sigs_counts_matrix_withcomplex)[dim(sv_sigs_counts_matrix_withcomplex)[1]]<-"WGS-ER-2"

write.table(sv_sigs_counts_matrix_withcomplex,"MutSignatures/data/sv_sigs/SV_signatures_counts_matrix_withcomplex.txt",sep="\t",quote=F)




#####################################################################################################
#Early classification attempts
#Size catagories come from mixture modelling approaches here
#Group all clustered vars together in one catagory. 
#No account of foldback inv or complex SV calls


#Separate clustered SVs
clustdat <- dat[dat$cluster_size>=3,]
simpledat <- dat[dat$cluster_size<3,]
dim(clustdat)[1]/dim(dat)[1]

# 26.6 % of SVs are in clusters

#Explore clusters
summary(clustdat$cluster_size)
hist(clustdat[clustdat$cluster_size<10,"cluster_size"])

#Simple SVs

#Relative numbers by type
table(simpledat$sv_type)

#CTX   DEL   DUP   INS   INV 
#8541 34848 23016    35  5338 

#V small number of insertions.. let's not include in signatures

tabtypes <- table(simpledat$sample,simpledat$sv_type)

library(reshape2)
library(ggplot2)

tabtype_long <- melt(tabtypes,id.vars=rownames(tabtypes))
tabtype_long$Cohort <- substr(tabtype_long[,1],1,2)
colnames(tabtype_long)<- c("Sample","sv_type","Number","Cohort")

tabtype_long[tabtype_long$Number==0,"Number"]<- 1
ggplot(tabtype_long,aes(x=sv_type,y=log(Number)))+ geom_boxplot(outlier.shape=NA)+theme_bw()+xlab("SV type")+
  ylab("Number of SVs per sample (logged)")+geom_jitter(width=0.02,pch=20)
                                                                             
ggplot(tabtype_long,aes(x=sv_type,y=log(Number)))+ geom_violin()+theme_bw()+xlab("SV type")+
  ylab("Number of SVs (logged)")+facet_grid(~Cohort)

#split by sv type
del <- simpledat[simpledat$sv_type=="DEL",]
dup <- simpledat[simpledat$sv_type=="DUP",]
inv <- simpledat[simpledat$sv_type=="INV",]
tra <- simpledat[simpledat$sv_type=="CTX",]

# Evaluate size distributions
summary(del$sv_length)
hist(log10(del$sv_length))


library(ggridges)
notra<-simpledat[simpledat$sv_type!="CTX",]
ggplot(notra, aes(x=log10(sv_length),y=sv_type,group=sv_type,height = ..density..))+
  geom_density_ridges(scale=0.9)+
  #geom_density_ridges(scale=0.9,rel_min_height = 0.01, jittered_points = TRUE,
                    #  position = position_points_jitter(width = 0.5, height = 0),
                     # point_shape = "|", point_size = 2,
                     # alpha = 0.7)+
  theme_bw()+ylab("SV type")+xlab("SV length (log10 bp)")+scale_x_continuous(breaks=c(0,2,4,6,8),labels=c("1bp","100bp","10kb","1Mb","100Mb"))

# Find thresholds in size distributions for sv classification
# Tried using mixture of normals to determine thresholds but no real substructure so arbitrary 
# thresholds on the log10 scale may be the way to go


library(mclust)
mod1 <- densityMclust(log10(dup$sv_length))
summary(mod1)
plot(mod1,what="BIC")
plot(mod1, what = "density", data =log10(dup$sv_length), breaks = 100, xlab="Duplication length")

mod1dr <- MclustDR(mod1)
summary(mod1dr)

par(mfrow=c(2,1))
plot(mod1, what = "density", data =log10(dup$sv_length), breaks = 100, xlab="Deletion length")
plot(mod1dr, what = "pairs")


ggplot(del,aes(x=as.factor(length_class),y=log10(sv_length)))+geom_boxplot(outlier.shape=NA)+
  theme_bw()

#Determine thresholds
#Class 1
max(log10(del[del$length_class==1,"sv_length"]))
min(log10(del[del$length_class==2,"sv_length"])) 

# < 100bp

#Class 2
max(log10(del[del$length_class==2,"sv_length"])) #2.247973
min(log10(del[del$length_class==3,"sv_length"])) #2.25042
10^2.9154
log10(823) # < ~1kb 

#Class 3
max(log10(del[del$length_class==3,"sv_length"])) #3.313656
min(log10(del[del$length_class==4,"sv_length"])) #3.313867
10^5.098197
log10(125371) #[1] < ~ 100kb

#Class 4 > 100kb

#######

#Classify dels

del$class <- character(dim(del)[1])
del[del$sv_length < 100 ,"class"]<- "< 100bp del"
del[del$sv_length < 1000 & del$sv_length >= 100,"class"]<- "< 1kb del"
del[del$sv_length < 10000 & del$sv_length >= 1000,"class"]<- "< 10kb del"
del[del$sv_length < 100000 & del$sv_length >= 10000,"class"]<- "< 100kb del"
del[del$sv_length < 1000000 & del$sv_length >= 100000,"class"]<- "< 1Mb del"
del[del$sv_length >= 1000000,"class"]<- ">= 1Mb del"

#< 100bp < 100kb  < 10kb   < 1kb   < 1Mb  >= 1Mb 
#6114    6098    7984    9299    3337    2016 

#Classify dups
dup$class <- character(dim(dup)[1])
dup[dup$sv_length < 1000 ,"class"]<- "< 1kb dup"
dup[dup$sv_length < 10000 & dup$sv_length >= 1000,"class"]<- "< 10kb dup"
dup[dup$sv_length < 100000 & dup$sv_length>= 10000,"class"]<- "< 100kb dup"
dup[dup$sv_length < 1000000 & dup$sv_length >= 100000,"class"]<- "< 1Mb dup"
dup[dup$sv_length < 10000000 & dup$sv_length >= 1000000,"class"]<- "< 10Mb dup"
dup[dup$sv_length>= 10000000,"class"]<- ">= 10Mb dup"

#< 100kb  < 10kb  < 10Mb   < 1kb   < 1Mb >= 10Mb 
#6670    6402    2122     576    6544     702 

#### Classify inversions ####

#Split into recip and unbalanced, check size thresholds
inv$event_info <- sapply( strsplit(inv$sv_id,":"), 
                          function(x) paste(x[1],x[2],sep=":"))
recipevents <- inv[duplicated(inv$event_info),"event_info"]
inv$reciprocal <- "unbalanced"
inv[inv$event_info %in% recipevents,"reciprocal"] <- "reciprocal"

#classify inversions by size
hist(log10(inv[inv$reciprocal=="reciprocal","sv_length"]))
hist(log10(inv[inv$reciprocal=="unbalanced","sv_length"]))

#recip inversions seem to be bimodal... perhaps two categories... use mixture modelling to check. Threshold approx 1Mb

mod1 <- densityMclust(log10(inv[inv$reciprocal=="reciprocal","sv_length"]),G=2)
summary(mod1)
plot(mod1,what="BIC")
plot(mod1, what = "density", data =log10(inv[inv$reciprocal=="reciprocal","sv_length"]), breaks = 100, xlab="Inversion length")

#what about unbalanced inversions? 1Mb also looks reasonable
mod1 <- densityMclust(log10(inv[inv$reciprocal=="unbalanced","sv_length"]),G=2)
summary(mod1)
plot(mod1,what="BIC")
plot(mod1, what = "density", data =log10(inv[inv$reciprocal=="unbalanced","sv_length"]), breaks = 100, xlab="Inversion length")

inv$class <- character(dim(inv)[1])
inv[inv$sv_length < 1000000 & inv$reciprocal=="reciprocal" ,"class"]<- "small reciprocal inv"
inv[inv$sv_length >= 1000000 & inv$reciprocal=="reciprocal" ,"class"]<- "large reciprocal inv"
inv[inv$sv_length < 1000000 & inv$reciprocal=="unbalanced" ,"class"]<- "small unbalanced inv"
inv[inv$sv_length >= 1000000 & inv$reciprocal=="unbalanced" ,"class"]<- "large unbalanced inv"
inv<- inv[,setdiff(colnames(inv),c("event_info","reciprocal"))]

# Translocations - balanced, unbalanced, intra-chromosomal, inter-chromosomal

tra$event_info <- sapply(strsplit(tra$sv_id, ":"), function(x) paste(x[1],x[2],sep=":"))
recipevents <- tra[duplicated(tra$event_info),"event_info"]
tra$reciprocal <- "unbalanced tra"
tra[tra$event_info %in% recipevents,"reciprocal"] <- "reciprocal tra"
tra$class <- tra$reciprocal
tra<- tra[,setdiff(colnames(tra),c("event_info","reciprocal"))]

#Clustered events - treat at event level not SV level.. v few large clusters, most samples have none so keep them all together

hist(clustdat$cluster_size,breaks=50)
hist(log(clustdat$cluster_size))
mean(clustdat$cluster_size)

mod1 <- densityMclust(clustdat$cluster_size,G=2)
summary(mod1)
plot(mod1,what="BIC")
plot(mod1, what = "density", data =clustdat$cluster_size, breaks = 100, xlab="Cluster size")
abline(h=34, col="blue",lty="dashed")

clustdat[clustdat$cluster_size > mean(clustdat$cluster_size),"size_class"]<- "Large cluster"
clustdat[clustdat$cluster_size <= mean(clustdat$cluster_size),"size_class"]<- "Small cluster"

clustlevel <-clustdat[,c("sample","cluster_id")]
clustlevel <- clustlevel[!duplicated(clustlevel),]
clusttab<-as.data.frame(table(clustlevel$sample))
clusttab$Var2 <- "clustered"

#Combine event classes
combine_simplesig_events <- rbind(del,dup,inv,tra)
combine_simplesig_events[1,]

#Create counts matrix and visualise

counts_matrix <- as.data.frame(table(combine_simplesig_events$sample,combine_simplesig_events$class))

#Add in samples with no SVs

samples<- read.table("~/Documents/HGSOC_marker_paper_onedrive/params/Combined_usable_ids_marker.txt",sep="\t")
samples<- as.character(samples[,1])


setdiff(samples,counts_matrix$Var1)
extradf <- data.frame(Var1=c(rep("SHGSOC107",18),rep("WGS-ER-2-P",18)), Var2=rep(unique(counts_matrix$Var2),2), Freq=rep(0,36))
counts_matrix<- rbind(counts_matrix,extradf)

setdiff(samples,clusttab$Var1)
extradf_clust <- data.frame(Var1=c(rep("SHGSOC084",1),rep("DG1250",1),rep("DG1251",1),rep("WGS-ER-2-P",1)),
                            Var2=rep(unique(clusttab$Var2),4), Freq=rep(0,4))
clusttab<- rbind(clusttab,extradf_clust)


counts_plusclust <- rbind(counts_matrix,clusttab)

svtype<- sapply( strsplit(as.character(counts_plusclust$Var2)," "),function(y) y[length(y)])
counts_plusclust$SVtype <- svtype
cohort_sh <- substr(counts_plusclust$Var1,1,2)
counts_plusclust$Cohort <- character(dim(counts_plusclust)[1])
counts_plusclust[which(cohort_sh=="SH"),"Cohort"] <- "SHGSOC"
counts_plusclust[which(cohort_sh=="AO"),"Cohort"] <- "AOCS"
counts_plusclust[which(cohort_sh=="DO"),"Cohort"] <- "TCGA"
counts_plusclust[which(cohort_sh=="WG"),"Cohort"] <- "MDA"
counts_plusclust[which(cohort_sh=="DG"),"Cohort"] <- "BCCA"
counts_plusclust[which(cohort_sh=="DA"),"Cohort"] <- "BCCA"

#Make wide
widecounts<-reshape(counts_plusclust[,1:3], idvar = "Var1", timevar = "Var2", direction = "wide")
colnames(widecounts)<-c("Sample",unique(as.character(counts_plusclust$Var2)))
#write.table(widecounts, file="MutSignatures/Counts_forSVsigs.txt",sep="\t",quote=F,row.names=F)


#Make some plots
classes <- rev(c("< 100bp del","< 1kb del","< 10kb del","< 100kb del","< 1Mb del", ">= 1Mb del",
             "< 1kb dup","< 10kb dup","< 100kb dup","< 1Mb dup", "< 10Mb dup", ">= 10Mb dup",
             "small reciprocal inv", "large reciprocal inv","small unbalanced inv","large unbalanced inv",
             "reciprocal tra","unbalanced tra","clustered"))
             
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Stratify by cohort          
counts_plusclust$Var2 <- factor(counts_plusclust$Var2,levels=classes)
ggplot(counts_plusclust,aes(x=Var2,y=log(Freq+1),col=SVtype))+geom_boxplot()+scale_color_manual(values=cbPalette[c(7,6,2,8,4)])+
theme_bw()+facet_grid(~Cohort,scales='free')+coord_flip()+xlab("Class of SV")+ylab("Frequency of SV\nlog(freq+1)")

#All cohorts
ggplot(counts_plusclust,aes(x=Var2,y=log(Freq+1),col=SVtype))+geom_boxplot()+scale_color_manual(values=cbPalette[c(7,6,2,8,4)])+
  theme_bw()+coord_flip()+xlab("Class of SV")+ylab("Frequency of SV\nlog(freq+1)")

