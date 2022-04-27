load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")

Dat<-rbind(dat1,dat2)
Dat<-estimate_bias(dat=Dat)

# 0
Dat<-Dat[Dat$lnor_sh <= 1.999 & Dat$lnor_sh>= -1.999,] #lnor_sh ==1.999 is an artifiact. 

#1 
MAF<-Dat$eaf
MAF[MAF>0.5]<-1-MAF[MAF>0.5]
Dat$MAC_case<-MAF*Dat$ncase*2
Dat$MAC_control<-MAF*Dat$ncontrol*2
Dat<-Dat[Dat$MAC_case>=50 & Dat$MAC_control>=50,]
dat<-Dat1

# 2
Dat<-Dat[Dat$study != "UKB",]
ID_outliers<-c(3,59,63,64,65,66,70,71,#datasets with slope <0.8 25 explained by N studies 
	25,104,132)#imputation info/r2

Dat<-Dat[!Dat$ID %in% ID_outliers,]


# 3
Dat<-Dat[Dat$lnor < 1 & Dat$lnor> -1,]

dat<-Dat[Dat$ID == 8, ]

Bias<-abs(dat$bias)
Bias2<-as.character(Bias)
summary(Bias)
Bias2[Bias>=50]<-">50"
Bias2[Bias<50 & Bias>=25]<-"25-50"
Bias2[Bias<25]<-"<25"

table(Bias2)

dat$Y<-dat$lnor_sh
dat$X<-dat$lnor


Model1<-summary(lm(Y~X,dat[Bias < 50 ,]))
Model2<-summary(lm(Y~X,dat))
int1<-Model1$coefficients[1,1]
slope1<-Model1$coefficients[2,1]
int2<-Model2$coefficients[1,1]
slope2<-Model2$coefficients[2,1]

subtitle<-paste0("intercept=",round(int1,3)," | ","slope=",round(slope1,3)," (Bias<50%)","\nintercept=",round(int2,3)," | ","slope=",round(slope2,3)," (All SNPs)")
Title<-paste0(unique(dat$study)," | ID: ",unique(dat$ID))
Xlab<-"Reported log odds ratio"
Ylab<-"Predicted log odds ratios"
Title_size<-12
Title_xaxis_size=9

labels_colour<-unique(Bias2)
values_colour<-unique(Bias2)
values_colour[values_colour==">50"]<-"red"
values_colour[values_colour=="25-50"]<-"orange"
values_colour[values_colour=="<25"]<-"black"
Pos<-order(values_colour)
Pos<-order(labels_colour)
values_colour<-values_colour[Pos]
labels_colour<-labels_colour[Pos]


png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/chc23.png",width=500,height=500)
	ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Bias2))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		ggplot2::scale_colour_manual(name = "Bias",
		labels = labels_colour,
		values = values_colour)+
		ggplot2::theme(legend.title=ggplot2::element_text(size=9))+
		ggplot2::theme(legend.text=ggplot2::element_text(size=8))
dev.off()
# ACCC
# other studies with DIrection 
# Dat1<-Dat[Dat$ID == 130, ]
# Dat1<-Dat[Dat$ID == 71, ]
# Dat1<-Dat[Dat$ID == 70, ]
# Dat1<-Dat[Dat$ID == 66, ]
# Dat1<-Dat[Dat$ID == 3, ]
# Dir<-Dat1$Direction

# Dir1<-Dir[1]
# Dir2<-gsub("\\?","",Dir)
# Dat1$Nstudies<-nchar(Dir2)
# summary(Dat1$Nstudies)
# # Dat1<-Dat1[Dat1$Nstudies==2,]
# Dat1<-Dat1[Dat1$Nstudies>13,]
# # Dat1<-Dat1[Dat1$Nstudies==3,]
# # Dat1<-Dat1[Dat1$Nstudies==3,]
# Dat<-Dat1