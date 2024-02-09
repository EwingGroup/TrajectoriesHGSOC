#Create centromere and telomere files

centros  <- read.table("hg38.centromeres",sep="\t")
centros.chrs<-gsub("chr","",centros$V1)
centros.chrs2<-gsub("X","23",centros.chrs)
centros.chrs3<-gsub("Y","24",centros.chrs2)

centros.nopref<-cbind(centros.chrs3,centros[,2:5])
centros.nopref[,1]<-as.numeric(centros.nopref[,1])
centros.sorted <- centros.nopref[order(centros.nopref[,1]),]


chr<-centros.sorted[seq(1,46,2),1]
cen_start<-centros.sorted[seq(1,46,2),2]
cen_end<-centros.sorted[seq(2,46,2),3

centromeres.hg38<-cbind(chr,cen_start,cen_end)


telos  <- read.table("hg38.telomeres",sep="\t")
telos.chrs<-gsub("chr","",telos$V1)
telos.chrs2<-gsub("X","23",telos.chrs)
telos.chrs3<-gsub("Y","24",telos.chrs2)

telos.nopref<-cbind(telos.chrs3,telos[,2:4])
telos.nopref[,1]<-as.numeric(telos.nopref[,1])
telos.sorted <- telos.nopref[order(telos.nopref[,1]),]


chr<-telos.sorted[seq(1,46,2),1]
ptel<-telos.sorted[seq(1,46,2),3]
qtel<-telos.sorted[seq(2,46,2),2]

telomeres.hg38<-cbind(chr,ptel,qtel+1)

hg38_cen_telo<-merge(centromeres.hg38,telomeres.hg38,by="chr")
hg38_cen_telo<-hg38_cen_telo[,c(1,4,2,3,5)]
colnames(hg38_cen_telo)<-c("chr","ptel","cen_start","cen_end","qtel")

hg38_cen_telo$chr<-paste("chr",hg38_cen_telo$chr,sep="")
hg38_cen_telo$chr<-gsub("23","X",hg38_cen_telo$chr)
 write.table(hg38_cen_telo,file="hg38_centromere_and_telomere_coords.txt",sep="\t",row.names=F,quote=F)
