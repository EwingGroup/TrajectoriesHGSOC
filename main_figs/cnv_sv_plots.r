
load(file="CNV_plot_data.RData");
load(file="SV_plot_data.RData");
load(file="CNV_gene_sample.RData");
library(gplots);

## deletions:      'deepskyblue'
## amplifications: 'darkseagreen'

## CNV plot
pdf(file='CNV_foldchange_scatterplot.pdf',width=5,height=5);
plot(CNV_plot_data[CNV_plot_data[,3]==1,c(1,4)],pch=19, xlim=c(0,30), ylim=c(0,1.1),
     col='darkseagreen',
     xlab='Percentage of samples with CNV',
     ylab='log2 Fold change')
points(CNV_plot_data[CNV_plot_data[,3]==0,c(2,4)],pch=19,col='deepskyblue')
for(i in 1:dim(CNV_plot_data)[1]) {
    os <- 0.03; colp<-'black';
    if(rownames(CNV_plot_data)[i]=="FGFR1OP") { os <- -os;  }
    if(CNV_plot_data[i,3]==1) { text(rownames(CNV_plot_data)[i],x=CNV_plot_data[i,1],y=CNV_plot_data[i,4]+os,cex=0.6,col=colp); } else {  text(rownames(CNV_plot_data)[i],x=CNV_plot_data[i,2],y=CNV_plot_data[i,4]+os,cex=0.6,col=colp); }
}
legend('topleft',pch=c(19,19),col=c('darkseagreen','deepskyblue'),legend=c('Amplification','Deletion'),cex=0.7);

dev.off();

## SV plot
## colour blind palette for event types
cb8bPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7");
cbp_SV <- cb8bPalette[c(2,4,6,8)];
cols <- array(dim=dim(SV_plot_data)[1],cbp_SV[1]); ## dup
cols[SV_plot_data[,4]=='SV_del'] <- cbp_SV[2];
cols[SV_plot_data[,4]=='SV_inv'] <- cbp_SV[3];
cols[SV_plot_data[,4]=='SV_breakpoint'] <- cbp_SV[4];


pdf(file='SV_foldchange_scatterplot.pdf',width=5,height=5);
plot(SV_plot_data[,c(3,5)],pch=19, xlim=c(0,25), ylim=c(-1.5,1.5),col=cols,
     xlab='Percentage of samples with SV',
     ylab='log2 Fold change')
lines(c(0,100),c(0,0))
for(i in 1:dim(SV_plot_data)[1]) {
    os <- 0.1;
    if(SV_plot_data[i,1]=='CEBPA') { os <- -os; }
    colp <- 'black';
    text(SV_plot_data[i,1],x=SV_plot_data[i,3],y=SV_plot_data[i,5]+os,cex=0.6,col=colp);
}
legend('topleft',col=cbp_SV,pch=19,legend=c('Duplication','Deletion','Inversion','Breakpoint'),cex=0.7)
dev.off();


## CNV heatmap
pdf(file=paste('CNV_heatmap.pdf',sep=''),width=8,height=4)
heatmap.2(CNV_gene_sample,trace='none',col=(c('deepskyblue','white','darkseagreen')),cexCol=0.2,cexRow=0.9,mar=c(8,8),key=F)
dev.off()
