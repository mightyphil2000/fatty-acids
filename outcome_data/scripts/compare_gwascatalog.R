source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")

Dat<-Dat[Dat$all_summary_stats, ]
IDS<-unique(Dat$ID)
# i<-which(IDS==133)
Plot1<-NULL 
# Plot1[i]	
for(i in 1:length(IDS)){
# for(i in 1:10){
	Dat1<-Dat[Dat$ID == IDS[i],]	
	print(i)
	Title<-unique(paste0(Dat1$study," | ID: ",Dat1$ID))	
	Plot1[[i]]<-make_plot_gwas_catalog_zscores(dat=Dat1,efo=unique(Dat1$efo),Title_size_subplot=12,Title=Title,Ylab="",Xlab="",Title_xaxis_size=0,legend=FALSE)
	#gwas_catalog_ancestral_group=unique(Dat1$population) 
}


make_cow_plot2(width=2000,height=1000,Title_axis_size=15,bycols=FALSE,
		Plot_list=Plot1[1:64],
		Title="",
		Title_size=0,
		Ylab="Z score in target study",
		Xlab="Z score in GWAS catalog",
		Legend_z=TRUE,
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/gwas_catalog_hits/cowplot_gwascatalog_zscores1_preqc.png"))


make_cow_plot2(width=2000,height=1000,Title_axis_size=15,bycols=FALSE,
		Plot_list=Plot1[65:length(Plot1)],
		Title="",
		Title_size=0,
		Ylab="Z score in target study",
		Xlab="Z score in GWAS catalog",
		Legend_z=TRUE,
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/gwas_catalog_hits/cowplot_gwascatalog_zscores2_preqc.png"))


source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")
Dat<-Dat[Dat$all_summary_stats & !is.na(Dat$eaf), ]
IDS<-unique(Dat$ID)
# i<-which(IDS==118)
Plot1<-NULL 
# Plot1[i]
for(i in 1:length(IDS)){
	# for(i in 1:3){
	print(i)
	Dat1<-Dat[Dat$ID == IDS[i],]
	Title<-unique(paste0(Dat1$study," | ID: ",Dat1$ID))
	Plot1[[i]]<-make_plot_gwas_catalog_eaf(dat=Dat1,efo=unique(Dat1$efo),Title_size_subplot=12,Title=Title,Ylab="",Xlab="",Title_xaxis_size=0,legend=FALSE)
}

# length(Plot1)/2

make_cow_plot2(width=2000,height=1000,Title_axis_size=0,bycols=FALSE,
		Plot_list=Plot1[1:61],
		Title="",
		Title_size=0,
		Ylab="",
		Xlab="",
		Legend_eaf=TRUE,
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/gwas_catalog_hits/cowplot_gwascatalog_eaf1_preqc.png"))


make_cow_plot2(width=2000,height=1000,Title_axis_size=0,bycols=FALSE,
		Plot_list=Plot1[62:length(Plot1)],
		Title="",
		Title_size=0,
		Ylab="",
		Xlab="",
		Legend_eaf=TRUE,
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/gwas_catalog_hits/cowplot_gwascatalog_eaf2_preqc.png"))



