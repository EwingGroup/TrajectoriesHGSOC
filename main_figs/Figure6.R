#Survival_plots


library(glmnet)
library(survival)
library(ggfortify)

load(file="coxph_fit_combined.RData");
load(file="oneFeatureCoefs.RData");
load(file="oneFeaturePvalues.RData");
load(file="nocases.RData");
load(file="cvplotdata.RData");

## part B. Coxph model of all 34 variables
pdf(file='../coxph_combined_variables_forest_plot_20230906_A_ln.pdf',width=5,height=12);
plotd <- cbind((confint(coxph_fit_combined))[,1],(coef(coxph_fit_combined)),(confint(coxph_fit_combined))[,2])
pvalues <- summary(coxph_fit_combined)$coefficients[,5];
nolines <- dim(plotd)[1];
cols <- array(dim=nolines,'black');
cols[pvalues<=0.05] <- 'orange';
cols[pvalues<=0.001] <- 'red';

pchs <- array(dim=nolines,15);
cexs <- array(dim=nolines,1);

par(mar=c(5.1, 8.1, 4.1, 2.1))
par(cex.axis=1, cex.lab=1.2, cex.main=1.2, cex.sub=1)
plot(0,0,xlim=c(min(plotd[,1]),max(plotd[,3])*1.2),ylim=c(1,nolines),
     yaxt='n',xlab='log(HR)',ylab='',col='white');

lines(c(0,0),c(-1,nolines+1),lty=2)
## coxph confidence intervals
for(i in 1:nolines) {
  lines(plotd[i,c(1,3)],c(i,i),lwd=0.5)
}
## coxph symbols
points(cbind(plotd[,2],array(1:nolines)),pch=pchs,cex=cexs,col=cols);
lbl <- sub('x','',names(coef(coxph_fit_combined)));
lbl[lbl=="Old"] <- "Age [old]";
lbl[lbl=="Late"] <- "Stage [late]";
lbl[lbl=="HRD"] <- "HRD [present]";
lbl[lbl=="bfb"] <- "Breakage Fusion Bridge";
lbl[lbl=="eccDNA"] <- "ecDNA";
lbl[lbl=="Chromothripsis"] <- "Severe Chromothripsis";
lbl <- gsub('_',' ',lbl);

#axis(2,at=1:nolines,labels=lbl,las=2,cex=1)
legend('topleft',pch=15,col=c('red','orange','black'),
       legend=c(expression(p<=0.001),expression(p<=0.05),expression(p>0.05)),cex=0.5)
dev.off();


## part A. Coxph models of each of the 34 variables + cohort
pdf(file='../variable_only_forest_plot_20230906_ln.pdf',width=6,height=12);
plotd <- oneFeatureCoefs
pvalues <- oneFeaturePvalues;
nolines <- dim(plotd)[1];
cols <- array(dim=nolines,'black');
cols[pvalues<=0.05] <- 'orange';
cols[pvalues<=0.001] <- 'red';

par(mar=c(5.1, 11.1, 4.1, 2.1))
par(cex.axis=1, cex.lab=1.2, cex.main=1.2, cex.sub=1)
plot(0,0,xlim=c(min(plotd[,1]),max(plotd[,3])*1.2),ylim=c(1,nolines),
     yaxt='n',xlab='log(HR)',ylab='',col='white');

lines(c(0,0),c(-1,nolines+1),lty=2)
## coxph confidence intervals
for(i in 1:nolines) {
  lines(plotd[i,c(1,3)],c(i,i),lwd=0.5)
}
## coxph symbols
points(cbind(plotd[,2],array(1:nolines)),pch=15,col=cols);


axis(2,at=1:nolines,labels=lbl,las=2,cex=1)
legend('topleft',pch=15,col=c('red','orange','black'),
       legend=c(expression(p<=0.001),expression(p<=0.05),expression(p>0.05)),cex=0.5)
dev.off();


## part C. no. cases in A and B
pdf(file='barplot_no_cases_20230906.pdf',width=12,height=5);
par(mar=c(8.1, 4.1, 4.1, 2.1))
barplot(nocases,las=2,ylim=c(0,280),names=rownames(nocases),cex.names=0.9,col='grey80',
        ,yaxt='n',)
lines(c(0,45),c(277,277),lty=1)
lines(c(0,45),c(138,138),lty=2,lwd=0.5)
lines(c(0,45),c(69,69),lty=2,lwd=0.5)
axis(4,at=c(0,50,100,150,200,250),labels=c(0,50,100,150,200,250))
dev.off();

## part D. ten runs of cv.glmnet(), selecting the 8 variables retained in all 10 runs
## plot as box plots horizontally; overlay data from Coxph B.
pdf(file='../cv_coxph_combined_forest_plot_20230906.pdf',width=5,height=4.5);
plotd <- cbind((confint(coxph_fit_combined))[,1],(coef(coxph_fit_combined)),(confint(coxph_fit_combined))[,2])
plotd <- plotd[(xselrns<-paste('x',cvplotdata$sel,sep='')),];
pvalues <- summary(coxph_fit_combined)$coefficients[,5];
pvalues <- pvalues[xselrns]
nolines <- dim(plotd)[1];
cols <- array(dim=nolines,'black');
cols[pvalues<=0.05] <- 'orange';
cols[pvalues<=0.001] <- 'red';

pchs <- array(dim=nolines,15);
cexs <- array(dim=nolines,1);

par(mar=c(5.1, 8.1, 4.1, 2.1))
par(cex.axis=1, cex.lab=1.2, cex.main=1.2, cex.sub=1);
plot(0,0,xlim=c(min(plotd[,1]),max(plotd[,3])*1.2),ylim=c(0.75,nolines+0.25),
     yaxt='n',xlab='log(HR)',ylab='',col='white');
lines(c(0,0),c(-1,nolines+1),lty=2);
## box plots
for(i in 1:length(cvplotdata$sel)) { ## box bounds from boxplot.stats; min, median and max from summary
  confi <- log2(exp(boxplot.stats((cvplotdata$smy[cvplotdata$sel[i],]))$conf));
  smi <- log2(exp(c(summary(cvplotdata$smy[cvplotdata$sel[i],]))[c(1,3,6)]));
  yi1 <- i - 0.1 -0.125;
  yi2 <- i + 0.1 -0.125;
  yi11 <- i - 0.15 -0.125;
  yi12 <- i + 0.15 -0.125;
  lines(c(confi[1],confi[1]),c(yi1,yi2),lwd=0.5); ## lower conf
  lines(c(confi[2],confi[2]),c(yi1,yi2),lwd=0.5); ## upper conf
  lines(c(confi[1],confi[2]),c(yi1,yi1),lwd=0.5);
  lines(c(confi[1],confi[2]),c(yi2,yi2),lwd=0.5);
  lines(c(smi[1],smi[3]),c(i -0.125,i -0.125),lwd=0.5);
  lines(c(smi[1],smi[1]),c(yi1,yi2),lwd=0.5); ## 
  lines(c(smi[3],smi[3]),c(yi1,yi2),lwd=0.5); ## 
  points(smi[2],i -0.125,pch=18,cex=1.5); ## median
}
## coxph confidence intervals
for(i in 1:nolines) {
  lines(plotd[i,c(1,3)],c(i+0.125,i+0.125),lwd=0.5)
}
## coxph symbols
points(cbind(plotd[,2],array(1:nolines)+0.125),pch=pchs,cex=cexs,col=cols);

lbl2 <- sub('x','',names(coef(coxph_fit_combined)[xselrns]));
lbl2[lbl2=="Old"] <- "Age [old]";
lbl2[lbl2=="Late"] <- "Stage [late]";
lbl2[lbl2=="HRD"] <- "HRD [present]";
lbl2[lbl2=="bfb"] <- "Breakage Fusion Bridge";
lbl2[lbl2=="Chromothripsis"] <- "Severe\nChromothripsis";
lbl2[lbl2=="eccDNA"] <- "EC DNA";
lbl2 <- gsub('_',' ',lbl2);

axis(2,at=1:nolines,labels=lbl2,las=2,cex=1)
legend('topleft',pch=15,col=c('red','orange','black'),
       legend=c(expression(p<=0.001),expression(p<=0.05),expression(p>0.05)),cex=0.5)
dev.off();

## part E. no. cases
pdf(file='cv_coxph_combined_barplot_no_cases_20230906.pdf',width=5,height=4);
par(mar=c(8.1, 4.1, 4.1, 2.1));
barplot(rev(nocases[lbl2]),las=2,ylim=c(0,280),names=rev(lbl2),cex.names=0.9,col='grey80',
        ,yaxt='n',)
lines(c(0,22),c(277,277),lty=1)
lines(c(0,22),c(138,138),lty=2,lwd=0.5)
lines(c(0,22),c(69,69),lty=2,lwd=0.5)
axis(4,at=c(0,50,100,150,200,250),labels=c(0,50,100,150,200,250))
dev.off();

