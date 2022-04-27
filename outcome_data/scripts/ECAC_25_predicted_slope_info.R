source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

# load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
# load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")

Dat<-load_data_info()

MAF<-Dat$eaf
MAF[MAF>0.5]<-1-MAF[MAF>0.5]
Dat$MAC_case<-MAF*Dat$ncase*2
Dat$MAC_control<-MAF*Dat$ncontrol*2
Dat<-Dat[Dat$MAC_case>=50 & Dat$MAC_control>=50,]
Dat<-Dat[Dat$lnor < 1 & Dat$lnor> -1,]
Dat<-Dat[Dat$lnor_sh < 1 & Dat$lnor_sh> -1,]

# 23,97
Dat<-Dat[Dat$ID == 25, ]

# min(Dat1$HWEp,na.rm=TRUE)

dat<-Dat
dat<-dat[order(dat$mean_info),]
info<-dat$mean_info
info2<-as.character(info)

info2[info>=0.90]<-"≥0.90"
info2[info<0.90 & info>=0.80 ]<-"0.80-0.90"
info2[info<0.80 & info>=0.70 ]<-"0.70-0.80"
# info2[info<0.70 ]<-"<0.70"
info2[info<0.70 & info>=0.60 ]<-"0.60-0.70"
info2[info<0.60 & info>=0.50 ]<-"0.50-0.60"
info2[info<0.50 ]<-"<0.50"
# info2[info<0.50 & info>=0.40 ]<-"0.40-0.50"
table(info2)
unique(info2)

dat$Y<-dat$lnor_sh
dat$X<-dat$lnor

# dat$Y<-dat$bias
# dat$X<-dat$mean_info

Model1<-summary(lm(Y~X,dat[info>0.90,]))
Model2<-summary(lm(Y~X,dat))
int1<-Model1$coefficients[1,1]
slope1<-Model1$coefficients[2,1]
int2<-Model2$coefficients[1,1]
slope2<-Model2$coefficients[2,1]

subtitle<-paste0("intercept=",round(int1,3)," | ","slope=",round(slope1,3)," (SNPs with info>0.90)","\nintercept=",round(int2,3)," | ","slope=",round(slope2,3)," (All SNPs)")
Title<-paste0(unique(dat$study)," | ID: ",unique(dat$ID))
Xlab<-"Reported log odds ratio"
Ylab<-"Predicted log odds ratios"
Title_size<-12
Title_xaxis_size=9

Colour<-info2
Colour[Colour == "≥0.90" ]<-"black"
Colour[Colour ==  "0.80-0.90" ]<-"blue"
Colour[Colour ==  "0.70-0.80" ]<-"purple"
Colour[Colour ==  "0.60-0.70" ]<-"green"
Colour[Colour == "0.50-0.60" ]<-"orange"
Colour[Colour == "<0.50" ]<-"red"
labels_colour<-unique(info2)
values_colour<-unique(info2)
values_colour[values_colour == "≥0.90" ]<-"black"
values_colour[values_colour ==  "0.80-0.90" ]<-"blue"
values_colour[values_colour ==  "0.70-0.80" ]<-"purple"
values_colour[values_colour ==  "0.60-0.70" ]<-"green"
values_colour[values_colour == "0.50-0.60" ]<-"orange"
values_colour[values_colour == "<0.50" ]<-"red"
Pos<-order(values_colour)
values_colour<-values_colour[Pos]
labels_colour<-labels_colour[Pos]


png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/ECAC_25.png",width=500,height=500)
	ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Colour))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size)) +
					ggplot2::scale_colour_manual(name = "Imputation info/r2",
		labels = labels_colour,
		values = values_colour)+
		ggplot2::theme(legend.title=ggplot2::element_text(size=9))+
		ggplot2::theme(legend.text=ggplot2::element_text(size=8))
dev.off()


dat$Y<-abs(dat$bias)
dat$X<-dat$mean_info
Xlab<-"Imputation info/r2"
Ylab<-"Absolute bias"

png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/ecac_25_bias.png",width=500,height=500)
		ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Colour))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size)) +
					ggplot2::scale_colour_manual(name = "Imputation info/r2",
		labels = labels_colour,
		values = values_colour)+
		ggplot2::theme(legend.title=ggplot2::element_text(size=9))+
		ggplot2::theme(legend.text=ggplot2::element_text(size=8))
dev.off()

		# 		+
		# ggplot2::scale_colour_manual(name = "HWEp",
		# labels = c("<5 studies","≥14 studies","10-14 studies","5-10 studies"),
		# values = c("red","black","blue","orange"))+
		# ggplot2::theme(legend.title=ggplot2::element_text(size=9))+
		# ggplot2::theme(legend.text=ggplot2::element_text(size=8))
# dev.off()
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