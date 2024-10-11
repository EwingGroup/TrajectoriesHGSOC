#################################################################
# To Run: 
# Rscript combine datatables

# Read in command line arguments
args <- commandArgs(trailingOnly=T)
datfile=args[1]

dat1 <- read.table(args[1], header=T)
ci=ncol(dat1)
datf <- dat1[,c(2:ci)]

for (f in args[2:(length(args)-1)]){
	tmpdat <- read.table(f, header=T)
	datf <- datf + tmpdat[,c(2:ci)] 
}
datf <- cbind(dat1[,1], datf)
colnames(datf) <- colnames(dat1)
write.table(datf, file=args[length(args)], quote=FALSE, sep="\t", row.names=T, col.names=T)

