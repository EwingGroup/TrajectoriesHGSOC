args <- commandArgs(trailingOnly=TRUE)

input<-as.character(args[1])
print(input)

output<-as.character(args[2])
print(output)

dat<-read.table(input,sep="\t")
dat <- dat[dat[,4]!=2 & dat[,10] !=2,]
 cnvkit_cndir<-ifelse(dat[,4]>2,"amp","del")
 climat_cndir<-ifelse(dat[,10]>2,"amp","del")
if (length(which(climat_cndir!=cnvkit_cndir))>0){
	dat<-dat[-which(climat_cndir!=cnvkit_cndir),]
}else{
	dat<-dat
}

write.table(dat,file=output,sep="\t",col.names=F,row.names=F,quote=F)
