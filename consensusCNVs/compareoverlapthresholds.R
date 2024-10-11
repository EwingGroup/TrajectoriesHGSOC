### Compare overlap thresholds ######

# Compare number of dels and dups by threshold - no meaningful difference in distribution across samples
dat_consensus_dels <- table(consensus_DELs$V4,consensus_DELs$overlap_threshold)
dat_consensus_dups <- table(consensus_DUPs$V4,consensus_DUPs$overlap_threshold)
dat_cnvkit_dels <- table(cnvkit_dels$V4)
dat_cnvkit_dups <- table(cnvkit_dups$V4)
overlap_threshold <- "CNVkit"
dat_cnvkit_dels <- data.frame(dat_cnvkit_dels,Overlap_threshold=rep(overlap_threshold,dim(dat_cnvkit_dels)[1]))
dat_cnvkit_dups <- data.frame(dat_cnvkit_dups,Overlap_threshold=rep(overlap_threshold,dim(dat_cnvkit_dups)[1]))
colnames(dat_cnvkit_dels)[1] <- "Sample"
colnames(dat_cnvkit_dups)[1] <- "Sample"

wide_dat_consensus_dels<- melt(dat_consensus_dels, id.vars=rownames(dat_consensus_dels))
wide_dat_consensus_dups<- melt(dat_consensus_dups, id.vars=rownames(dat_consensus_dups))
colnames(wide_dat_consensus_dels) <- colnames(wide_dat_consensus_dups) <- c("Sample","Overlap_threshold","Freq")

compare_dels <- rbind(dat_cnvkit_dels,wide_dat_consensus_dels)
compare_dups <- rbind(dat_cnvkit_dups,wide_dat_consensus_dups)
compare_dups$Overlap_threshold <- factor(compare_dups$Overlap_threshold,levels=c("CNVkit","0.5","0.7","0.9"))
compare_dels$Overlap_threshold <- factor(compare_dels$Overlap_threshold,levels=c("CNVkit","0.5","0.7","0.9"))

#png("ConsensusCNVs/Figs/CNVkit_Purple_consensus_dups.png",width=8,height=5,unit="in",res=300)
ggplot(compare_dups, aes(x=Overlap_threshold, y=log10(Freq))) + geom_violin()+geom_jitter(width=0.05)+
  xlab("Proportion of CNVkit call overlapped by Purple")+theme_bw()+ylab("Number of duplications (log10)")
#dev.off()

#png("ConsensusCNVs/Figs/CNVkit_Purple_consensus_dels.png",width=8,height=5,unit="in",res=300)
ggplot(compare_dels, aes(x=Overlap_threshold, y=log10(Freq))) + geom_violin()+geom_jitter(width=0.05)+
  xlab("Proportion of CNVkit call overlapped by Purple")+theme_bw()+ylab("Number of deletions (log10)")
#dev.off()

