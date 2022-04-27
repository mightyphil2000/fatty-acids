source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")
IDS<-unique(Dat$ID)
# i<-which(IDS ==74)
# Plot1[i]
Plot1<-NULL 
for(i in 1:length(IDS)){
	Dat1<-Dat[Dat$ID == IDS[i],]
	# Dat1$z.p<-qnorm(Dat1$p/2,lower.tail=F)
	# Dat1$z.lnor<-abs(Dat1$lnor/Dat1$se)
	# Dat1<-Dat1[Dat1$z.lnor<20,]
	print(i)
	Dat1<-Dat[Dat$ID == 74,]
	Title<-unique(paste0(Dat1$study," | ID: ",Dat1$ID))
	Plot1[[i]]<-predz_vs_obsz(dat=Dat1,Title=Title,Xlab="",Ylab="",Title_xaxis_size=0,Title_size=9)
}

make_cow_plot2(width=2000,height=1500,Title_axis_size=15,bycols=FALSE,
		Plot_list=Plot1,
		Title="",
		Title_size=0,
		Ylab="Z score from log odds ratio and standard error",
		Xlab="Z score from P value",
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/zzplots/cowplot_zzplots_preqc.png"))

Dat1<-Dat[Dat$ID == 74,]
Dat1$z.p<-qnorm(Dat1$p/2,lower.tail=F)
Dat1$z.lnor<-abs(Dat1$lnor/Dat1$se)
Dat1[Dat1$z.lnor>10 & Dat1$z.p < 15,]
cor(Dat1$z.p,Dat1$z.lnor)

Dat1<-Dat[Dat$ID == 133,]
Dat1$z.p<-qnorm(Dat1$p/2,lower.tail=F)
Dat1$z.lnor<-abs(Dat1$lnor/Dat1$se)
Dat1<-Dat1[Dat1$z.lnor<20,]
cor(Dat1$z.p,Dat1$z.lnor)