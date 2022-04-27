source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")

################################
# Maf plots for schs and CHARGE#
################################

Dat<-Dat[Dat$population == "European",]
# Dat<-Dat[Dat$population == "East Asian",]

IDS<-unique(Dat$ID)
# i<-1
Plot1<-NULL
for(i in 1:length(IDS)){
	print(i)
	Dat1<-Dat[Dat$ID == IDS[i],]	
	if(any(!is.na(Dat1$eaf)))
	{
		ref_studies<-c("CHARGE","SCHS")
		ref_pop<-c("European","East Asian")
		Pop<-unique(Dat1$population)
		ref_study<-ref_studies[which(ref_pop %in% Pop)]
		Plot1[[i]]<-make_plot_maf(refstudy=ref_study,target_dat=Dat1,Title_xaxis_size=0,Title_size=9,Title=paste0(unique(Dat1$study)," | ID: ",unique(Dat1$ID)),Ylab="",Xlab="")
	}
}


make_cow_plot2(width=1300,height=1000,Title_axis_size=0,bycols=FALSE,
		Plot_list=Plot1,
		Title="",
		Title_size=0,
		Ylab="",
		Xlab="",
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/cowplot_charge_preqc2.png"))


#########################################
# Maf plots for 1000 G super populations#
#########################################

load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")
# Dat<-Dat[Dat$population == "East Asian",]
Dat<-Dat[Dat$population == "European",]
Dat<-Dat[Dat$all_summary_stats,]
Dat<-Dat[!is.na(Dat$eaf),]
# Dat<-Dat[Dat$population == "European",]
IDS<-unique(Dat$ID)
IDS1<-IDS[1:34]
IDS2<-IDS[35:68]
IDS3<-IDS[69:103]
# length(IDS)
Plot1<-NULL
for(i in 1:length(IDS2))
{
	print(i)
	Dat1<-Dat[Dat$ID == IDS2[i],]
	Plot1[[i]]<-make_plot_maf(refstudy="ALL_1000G",target_dat=Dat1,Title_xaxis_size=8,Title_size=12,cowplot_title=paste0(unique(Dat1$study)," | ID: ",unique(Dat1$ID)),Ylab="",Xlab="")
}

make_cow_plot2(width=1700,height=1100,Title_axis_size=0,bycols=FALSE,
		Plot_list=Plot1,
		Title="",
		Title_size=0,
		Ylab="",
		Xlab="",
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/cowplot_eur_all_1000G_preqc2_2.png"))



