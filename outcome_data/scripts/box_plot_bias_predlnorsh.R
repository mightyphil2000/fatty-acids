
library(ggplot2)


#########################################
# Predicted log odds ratios versus reported effect sizes #
########################################
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")
Dat<-rbind(dat1,dat2)
Dat<-estimate_bias(dat=Dat)
# Dat<-fix_info(dat=Dat)# include the SNPs with imputation scores< 0.8  

Dat$ID2<-as.factor(Dat$ID)
Dat<-Dat[!Dat$study %in% "UKB",]
Dat<-Dat[abs(Dat$bias)<100,]
Dat$bias<-abs(Dat$bias)
png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/box_plot_infofiltered_bias100.png",width=2000,height=500)
ggplot(Dat, aes(x=ID2, y=bias,fill=study)) + 
    geom_boxplot()+
    theme(axis.text.x = element_text(size = 8,angle=90)) 
dev.off()

#flags: lnor_sh<=1.999, mac >50, info>0.8, UKB, sample size bias, bias<100
load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")
Dat<-rbind(dat1,dat2)
Dat<-estimate_bias(dat=Dat)
# 0
Dat<-Dat[Dat$lnor_sh <= 1.999 & Dat$lnor_sh>= -1.999,] #lnor_sh ==1.999 is an artifiact. 
# 1
MAF<-Dat$eaf
MAF[MAF>0.5]<-1-MAF[MAF>0.5]
Dat$MAC_case<-MAF*Dat$ncase*2
Dat$MAC_control<-MAF*Dat$ncontrol*2
# Dat<-Dat[order(Dat$MAC_case),]
Dat<-Dat[Dat$MAC_case>=50 | Dat$MAC_control >= 50, ]

# 2
Dat<-Dat[Dat$study != "UKB",]
ID_outliers<-c(3,59,63,64,65,66,70,71,#datasets with slope <0.8 25 explained by N studies 
	25,104,132)#imputation info/r2
Dat<-Dat[!Dat$ID %in% ID_outliers,]

# 3
Dat<-Dat[Dat$lnor < 1 & Dat$lnor> -1,]

Dat<-Dat[abs(Dat$bias)<100,]
Dat$ID2<-as.factor(Dat$ID)
Dat<-Dat[!Dat$study %in% "UKB",]
Dat$bias<-abs(Dat$bias)

png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/box_plot_flags123_filtered_bias100.png",width=2000,height=500)
ggplot(Dat, aes(x=ID2, y=bias,fill=study)) + 
    geom_boxplot()+
    theme(axis.text.x = element_text(size = 8,angle=90)) 
dev.off()

