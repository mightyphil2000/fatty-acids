load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")

dat<-rbind(dat1,dat2)
# dat<-dat[dat$study != "UCSF_AGS + SFAGS",]

dat<-dat[dat$lnor != "-Inf",]
dat<-dat[dat$lnor != "Inf",]
dat<-dat[dat$lnor != 0,]
dat$bias<-(dat$lnor_sh-dat$lnor )/dat$lnor*100	
IDS<-unique(dat$ID[which(!is.na(dat$HWEp))])
length(IDS)
dat<-dat[dat$ID %in% IDS,]
# IDS<-unique(Dat$ID[which(!is.na(Dat$phet))])
Study<-dat$study

Xlab<-"P value for Hardy Weinberg Equilibrium"
Ylab<-"Absolute bias in predicted log odds ratio"
dat$Y<-abs(dat$bias)
# dat$Y<-dat$bias
dat$X<-dat$HWEp

Title<-""
Title_size<-0
Title_xaxis_size<-10

png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/hwep2.png",width=1000,height=500)
ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Study))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=10))
dev.off()