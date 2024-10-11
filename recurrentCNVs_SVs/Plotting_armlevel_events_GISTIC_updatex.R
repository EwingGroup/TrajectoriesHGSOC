#Arm level event plotting

broad_results <- read.table("ConsensusCNVs/HGSOC_purple_x/hgsoc_purple_x.broad_significance_results.txt",sep="\t",header=T)

broad_results$Del.frequency.Inv <- -broad_results$Del.frequency
broad_results_plot <- broad_results[,c("Arm","Amp.frequency","Amp.q.value","Del.frequency.Inv","Del.q.value")]
broad_results_amp <- cbind(dir="Amp",broad_results_plot[,c("Arm","Amp.frequency","Amp.q.value")])
broad_results_del <- cbind(dir="Del",broad_results_plot[,c("Arm","Del.frequency.Inv","Del.q.value")])
colnames(broad_results_amp)<-colnames(broad_results_del)<- c("Dir","Arm","Frequency","Qvalue")
long.br_res <- rbind(broad_results_amp,broad_results_del)

library(ggplot2)
library(ggbreak)
armlevels <- c(paste0(rep(1:22,each=2),rep(c("p","q"),22),sep=""),"Xp","Xq")
long.br_res$Arm <- factor(long.br_res$Arm,levels=armlevels)

#Cite fig2a here: https://www.nature.com/articles/s41467-021-25787-x

#png("Broad_armlevel_gainslosses_gistic_barplot_withX.png",width=10,height=4,unit="in",res=300)
ggplot(long.br_res,aes(x=Arm,y=Frequency,fill=Dir)) + geom_bar(stat="identity",position="identity",col="black")+theme_bw()+
  scale_y_break(breaks=c(0,0),expand=FALSE,space=0.3) + ylim(c(-1,1))+theme(legend.title=element_blank())+scale_fill_manual(values=c("tomato1","steelblue3"))+
  geom_text( aes(x = Arm, y= Frequency),label=ifelse(long.br_res$Qvalue<0.05,"*","") ,position = position_nudge(y = ifelse(long.br_res$Dir=="Amp",0.01,-0.04)),size=8)+
  annotate("text",x="Xq",y=0.95,label= "q-value <0.05 (*)",vjust=1,hjust=1,fontface=2)
#dev.off()

#Sample level data

signif_arms <- as.character(unique(long.br_res[long.br_res$Qvalue<0.05 , "Arm"]))

armleveldat <- read.table("ConsensusCNVs/HGSOC_purple_x/hgsoc_purple_x.broad_values_by_arm.txt",sep="\t",header=T,row.names = 1)
armleveldat[1:4,1:4]
armlevel_cnstatus <- armleveldat
armlevel_cnstatus[armlevel_cnstatus>0] <- "Amp"
armlevel_cnstatus[armlevel_cnstatus<0] <- "Del"
armlevel_cnstatus[armlevel_cnstatus==0] <- "Neutral"
write.table(armlevel_cnstatus,file="ConsensusCNVs/GISTIC_armlevel_cnstatus_update_withX.txt",sep="\t",quote=F)
