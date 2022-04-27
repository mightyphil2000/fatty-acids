load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")

dat<-rbind(dat1,dat2)

MAF<-dat$eaf
MAF[MAF>0.5]<-1-MAF[MAF>0.5]
dat$MAC_case<-MAF*dat$ncase*2
dat$MAC_control<-MAF*dat$ncontrol*2
dat<-dat[order(dat$MAC_case),]
mac_case<-dat$MAC_case
mac_control<-dat$MAC_control

dat$Min_mac<-apply(dat[,c("MAC_case","MAC_control")],1, FUN=min)

# head(dat[dat$MAC_case>0,c("MAC_case","MAC_control","Min_mac")])
# dat<-dat[!dat$study %in% c("UCSF_MAYO","UKB","UCSF_AGS + SFAGS"),]
dat<-dat[dat$lnor != "-Inf",]
dat<-dat[dat$lnor != "Inf",]
dat<-dat[dat$lnor != 0,]
dat$bias<-(dat$lnor_sh-dat$lnor )/dat$lnor*100	
IDS<-unique(dat$ID[which(!is.na(dat$eaf))])
dat<-dat[dat$ID %in% IDS,]
dat<-dat[abs(dat$bias)<400,]
dat<-dat[dat$MAC_case<500,]
dat<-dat[dat$ID !=3,]
dat<-dat[dat$Min_mac>0,]
Study<-dat$study
dat$Y<-abs(dat$bias)
# dat$Y<-dat$bias
dat$X<-dat$MAC_case
dat$X<-dat$MAC_control
dat$X<-dat$Min_mac


Model<-summary(lm(Y~X+study,dat))

Xlab<-"Minimum minor allele count in cases or controls"
Ylab<-"Absolute bias in predicted log odds ratio"

png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/mac.png",width=1000,height=500)

ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Study))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=10))
dev.off()
