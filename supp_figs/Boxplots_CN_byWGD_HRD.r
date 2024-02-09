


load(file="data/Data_X_20230920.RData");      ## dups/dels overlapping GISTIC regions
load(file="data/Data_Y_20230920.RData");      ## all dups/dels

## 'deepskyblue' == deletions
## 'darkseagreen' == amplifications


## boxplots of ratio GISTIC dups / all dups
## total dups/dels per sample
hrd_wgd_total_dups <- Data_Y[[1]]; ## total dups HRD WGD
hrd_nowgd_total_dups <- Data_Y[[2]];
nohrd_wgd_total_dups <- Data_Y[[3]];
nohrd_nowgd_total_dups <- Data_Y[[4]];
hrd_wgd_total_dels <- Data_Y[[5]]; ## total dels HRD WGD
hrd_nowgd_total_dels <- Data_Y[[6]];
nohrd_wgd_total_dels  <- Data_Y[[7]];
nohrd_nowgd_total_dels <- Data_Y[[8]];

hrd_status <- c(rep(1,length(c(hrd_wgd_total_dups,hrd_nowgd_total_dups))),
                  rep(0,length(c(nohrd_wgd_total_dups,nohrd_nowgd_total_dups))))
wgd_status <- c(rep(1,length(hrd_wgd_total_dups)),rep(0,length(hrd_nowgd_total_dups)),
                rep(1,length(nohrd_wgd_total_dups)),rep(0, length(nohrd_nowgd_total_dups)))

dups <- data.frame(Total_dups=c(hrd_wgd_total_dups, hrd_nowgd_total_dups, nohrd_wgd_total_dups,nohrd_nowgd_total_dups),
                   HRD=hrd_status,WGD=wgd_status)

dels <- data.frame(Total_dels=c(hrd_wgd_total_dels, hrd_nowgd_total_dels, nohrd_wgd_total_dels,nohrd_nowgd_total_dels),
                   HRD=hrd_status,WGD=wgd_status)

#Dups and WGD plot
dups$WGD<- factor(dups$WGD,levels=c("1","0"))
dup_plot<-ggplot(dups,aes(x=as.factor(WGD),y=log10(Total_dups+1)))+geom_boxplot(fill="#009E73")+theme_classic(base_size=16)+
  scale_x_discrete(labels=c("0"="No WGD","1"="WGD"))+ylab("Total CNA duplication (log10)")+xlab("")
wilcox.test(dups$Total_dups~dups$WGD)

#Dels and HRD plot
dels$HRD<- factor(dels$HRD,levels=c("1","0"))
del_plot<-ggplot(dels,aes(x=as.factor(HRD),y=log10(Total_dels+1)))+geom_boxplot(fill="#0072B2")+theme_classic(base_size=16)+
  scale_x_discrete(labels=c("0"="No HRD","1"="HRD"))+ylab("Total CNA deletion (log10)")+xlab("")
wilcox.test(dels$Total_dels~dels$HRD)

hrd_wgd<-dups[,c("HRD","WGD")]
hrd_wgd$Sample<- rownames(hrd_wgd)
library(reshape2)
long_hrdwgd<- melt(hrd_wgd,id.vars="Sample")
library(ggmosaic)
ggplot(long_hrdwgd) + geom_mosaic(aes(x = variable, fill=value))

hrd_wgd$HRD<-factor(hrd_wgd$HRD,levels=c("0","1"))
hrd_wgd$WGD<-factor(hrd_wgd$WGD,levels=c("0","1"))
hrd_wgd %>%
  count(HRD, WGD, name = "freq") %>%
  ggplot() +
  geom_mosaic(aes(weight = freq, 
                  x = product(HRD), 
                  fill = WGD))

library(patchwork)
pdf("individual_panels/supp_figs/SuppFig3BC_WGD_HRD_dups_dels_boxplots.pdf",width=8,height=4)
dup_plot+del_plot
dev.off()

## total amps per sample
boxplot(log10(1+hrd_wgd_total_dups),
        log10(1+hrd_nowgd_total_dups),
        log10(1+nohrd_wgd_total_dups),
        log10(1+xx4),        
        col='darkseagreen',xlab='',ylab='log10 Amplifications per sample',main='',
        names=c('HRD\nWGD','HRD\nno WGD','no HRD\nWGD','no HRD\nno WGD'),las=2);


## total dels per sample 
boxplot(log10(1+yy1),
        log10(1+yy2),
        log10(1+yy3),
        log10(1+yy4),
        col='deepskyblue',xlab='',ylab='log10 Deletions per sample',main='',
        names=c('HRD\nWGD','HRD\nno WGD','no HRD\nWGD','no HRD\nno WGD'),las=2);




