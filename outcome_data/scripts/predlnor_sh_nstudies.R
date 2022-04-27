load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")

Dat<-rbind(dat1,dat2)
# IDS<-unique(Dat$ID[which(!is.na(Dat$HWEp))])
# IDS<-unique(Dat$ID[which(!is.na(Dat$phet))])
IDS<-unique(Dat$ID[which(!is.na(Dat$Direction))])
Dat<-Dat[Dat$ID %in% IDS,]
Dir<-Dat$Direction
Dir2<-gsub("\\?","",Dir)
Dat$Nstudies<-nchar(Dir2)
Dat<-Dat[Dat$lnor != "-Inf",]
Dat<-Dat[Dat$lnor != "Inf",]
Dat<-Dat[Dat$lnor != 0,]
Dat$bias<-(Dat$lnor_sh-Dat$lnor )/Dat$lnor*100	
Dat$Nstudies_median<-NA
Dat$Nstudies2<-NA
Dat_list<-NULL
for(i in 1:length(IDS)){
	print(i)
	Dat1<-Dat[Dat$ID == IDS[i],]
	# Med<-median(Dat1$Nstudies)
	Sum<-summary(Dat1$Nstudies)
	Min<-Sum[1]
	p25<-Sum[2]
	Med<-Sum[3]
	p75<-Sum[5]
	Max<-Sum[6]
	Dat1$Nstudies2[which(Dat1$Nstudies == Min )]<- 0
	Dat1$Nstudies2[which(Dat1$Nstudies == p25 )]<- 0.25
	Dat1$Nstudies2[which(Dat1$Nstudies == Med )]<- 0.5
	Dat1$Nstudies2[which(Dat1$Nstudies == p75 )]<- 0.75
	Dat1$Nstudies2[which(Dat1$Nstudies == Max )]<- 1
	# Dat1$Nstudies_median[which(Dat1$Nstudies>=Med)]<-1
	# Dat1$Nstudies_median[which(Dat1$Nstudies<Med)]<-0

	Dat_list[[i]]<-Dat1
}

Dat2<-do.call(rbind,Dat_list)
table(Dat2$Nstudies2)

dat<-Dat2
# dat<-dat[dat$study != "UCSF_AGS + SFAGS",]
Study<-dat$study

Xlab<-"SNP present above or below median number of studies"
Ylab<-"Absolute bias in predicted log odds ratio"
dat$Y<-abs(dat$bias)
dat$Y<-dat$bias
dat$X<-dat$Nstudies


Title<-""
Title_size<-0
Title_xaxis_size<-5

ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Study))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))


