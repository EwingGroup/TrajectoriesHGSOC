library(VariantAnnotation)
library(StructuralVariantAnnotation)

args <- commandArgs(trailingOnly=TRUE)

vcffile<-as.character(args[1])
print(vcffile)

outfile<-gsub(".vcf","",vcffile)

samples <- read.table(paste(outfile,"_samples.txt",sep=""))
samples <-as.character(samples[,1])

breakpointgr2bedpeType <- function(gr){

	bedpe <- data.frame(chrom1 = GenomeInfoDb::seqnames(gr), 
        		start1 = start(gr) - 1, 
		end1 = end(gr), 
		chrom2 = GenomeInfoDb::seqnames(partner(gr)), 
        		start2 = start(partner(gr)) - 1, 
		end2 = end(partner(gr)), 
        		name = names(gr), 
		partner.name = names(partner(gr)), 
        		score = gr$QUAL, strand1 = strand(gr), 
		strand2 = strand(partner(gr)),
		svtype = gr$simpleSVtype,
		svlen = gr$svLen,
		sample = gr$SAMPLE)	
    	bedpe <- bedpe[(as.character(bedpe$chrom1) < as.character(bedpe$chrom2)) | 
        (bedpe$chrom1 == bedpe$chrom2 & bedpe$start1 < bedpe$start2), 
        -c(8)]
    return(bedpe)

}




vcf <- VariantAnnotation::readVcf(vcffile,"hg38")
info(vcf)$SAMPLE <- samples
gr<-breakpointRanges(vcf, info_columns = "SAMPLE",inferMissingBreakends=TRUE)
svtype <- simpleEventType(gr)
mcols(gr)$simpleSVtype<-svtype
bedpe<-breakpointgr2bedpeType(gr)

write.table(bedpe,file=paste(outfile,".bedpe",sep=""),sep="\t",row.names=F,quote=F,col.names=F)

