library(VariantAnnotation)
library(StructuralVariantAnnotation)
vcffile="Consensus100bpSVs_ids.vcf"
samples <- read.table("Consensus100bpSVs_ids_samples.txt")
samples <-as.character(samples[,1])

vcf <- VariantAnnotation::readVcf(vcffile,"hg38")
info(vcf)$SAMPLE <- samples
gr<-breakpointRanges(vcf, info_columns = "SAMPLE", inferMissingBreakends=TRUE)

#Tile the genome
tiles <- gr.tile(seqinfo(gr), 49e3)
supertiles <- tiles + 500

#load mask

library(rtracklayer)
eligible.wgs =  import("um75-hs37d5.covered_chr.bed")

#Model breakpoint density using nb model (fishHook)
fish.supertiles = Fish(hypotheses = supertiles, events = gr, eligible = eligible.wgs, idcol = 'SAMPLE', use_local_mut_density = TRUE)
fish.supertiles.snv = fish.supertiles.snv[fish.supertiles.snv$data$frac.eligible > 0.75, ]
fish.supertiles.snv$score()
fish.supertiles$res[1,]
save(fish.supertiles,file="fish.supertiles.RData")
fish.supertiles$qqp(plotly = FALSE)

df<-gr2dt(fish.supertiles$res)
df_fdr<- df[df$fdr<0.25,]
df_sort<- df_fdr[order(df_fdr$p),]
