#To perform mixture modelling to investigate whether there are multiple recurrence patterns across samples of dels/dups at a gene level. e.g. how recurrent is interesting? 
#Take homes: 3 peaks of del recurrence but only 1 for dups
#Where did this go: Largely exploratory, haven't used it for anything.

recurrently_del_genes <- read.table("ConsensusCNVs/intermediate_tables/genes_recurrently_deleted_gisticpeaks.txt")


par(mfrow=c(2,1))
library(mclust)
mod1 <- densityMclust(recurrently_del_genes$V1,G=3)
summary(mod1)
plot(mod1,what="BIC")
plot(mod1, what = "density", data =recurrently_del_genes$V1, breaks = 100, xlab="Number of samples affected (recurrence)")

mod1dr <- MclustDR(mod1)
summary(mod1dr)
plot(mod1dr, what = "pairs")

most_recurrent<- as.character(recurrently_del_genes[which(mod1dr$classification == 3),2])
mid_recurrent<- as.character(recurrently_del_genes[which(mod1dr$classification == 2),2])
least_recurrent<- as.character(recurrently_del_genes[which(mod1dr$classification == 2),2])

recurrently_dup_genes <- read.table("ConsensusCNVs/intermediate_tables/genes_recurrently_duplicated_gisticpeaks.txt")
mod2 <- densityMclust(recurrently_dup_genes$V1)
summary(mod2)
plot(mod2,what="BIC")
plot(mod2, what = "density", data =recurrently_dup_genes$V1,breaks=100,  xlab="Number of samples affected (recurrence)")


mod2dr <- MclustDR(mod2)
summary(mod2dr)
plot(mod2dr, what = "pairs")
