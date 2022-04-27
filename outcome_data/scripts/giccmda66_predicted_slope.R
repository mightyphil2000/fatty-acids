load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")

Dat<-rbind(dat1,dat2)
Dat<-estimate_bias(dat=Dat)

MAF<-Dat$eaf
MAF[MAF>0.5]<-1-MAF[MAF>0.5]
Dat$MAC_case<-MAF*Dat$ncase*2
Dat$MAC_control<-MAF*Dat$ncontrol*2
Dat<-Dat[Dat$MAC_case>=50 & Dat$MAC_control>=50,]
Dat<-Dat[Dat$lnor < 1 & Dat$lnor> -1,]
Dat<-Dat[Dat$lnor_sh < 1 & Dat$lnor_sh> -1,]

Dat1<-Dat[Dat$ID == 66, ]
Dir<-Dat1$Direction

Dir1<-Dir[1]
Dir2<-gsub("\\?","",Dir)
Dat1$Nstudies<-nchar(Dir2)
dat<-Dat1
dat<-dat[order(dat$Nstudies),]
Nstudies<-dat$Nstudies
Nstudies2<-as.character(dat$Nstudies)
Nstudies2[Nstudies==2]<-"2 studies"
Nstudies2[Nstudies==1]<-"1 study"
table(Nstudies2)
median(Nstudies)
dat$Y<-dat$lnor_sh
dat$X<-dat$lnor

Model1<-summary(lm(Y~X,dat[Nstudies2=="2 studies",]))
Model2<-summary(lm(Y~X,dat))
int1<-Model1$coefficients[1,1]
slope1<-Model1$coefficients[2,1]
int2<-Model2$coefficients[1,1]
slope2<-Model2$coefficients[2,1]

subtitle<-paste0("intercept=",round(int1,3)," | ","slope=",round(slope1,3)," (SNPs present in 2 studies)","\nintercept=",round(int2,3)," | ","slope=",round(slope2,3)," (all SNPs)")

Title<-paste0(unique(dat$study)," | ID: ",unique(dat$ID))
Xlab<-"Reported log odds ratio"
Ylab<-"Predicted log odds ratios"
Title_size<-12
Title_xaxis_size=9

unique(Nstudies2)
png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/giccmda66.png",width=500,height=500)
	ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Nstudies2))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		ggplot2::scale_colour_manual(name = "Nstudies that contributed\nto each SNP",
		labels = c(" 1 study","2 studies"),
		values = c("red","black"))+
		ggplot2::theme(legend.title=ggplot2::element_text(size=9))+
		ggplot2::theme(legend.text=ggplot2::element_text(size=8))
dev.off()

dat$Y<-abs(dat$bias)
dat$X<-dat$Nstudies
Xlab<-"N studies"
Ylab<-"Absolute bias"
png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/giccmda_66_bias.png",width=500,height=500)
	ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Nstudies2))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		ggplot2::scale_colour_manual(name = "Nstudies that contributed\nto each SNP",
		labels = c("<6 studies","6 studies"),
		values = c("red","black"))+
		ggplot2::theme(legend.title=ggplot2::element_text(size=9))+
		ggplot2::theme(legend.text=ggplot2::element_text(size=8))
dev.off()