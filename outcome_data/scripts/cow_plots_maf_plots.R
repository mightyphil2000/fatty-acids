source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")

Dat1<-Dat[!is.na(Dat$Direction),]
length(unique(Dat1$ID))
IDS1<-unique(Dat1$ID)
IDS<-unique(Dat$ID)
length(IDS[!IDS %in% IDS1])
length(IDS[IDS %in% IDS1])


# Dat[Dat$ID == 58,]

# length(unique(Dat$ID))
Dat1$Index<-paste(Dat1$ID,Dat1$ncase)
Dat1<-Dat1[!duplicated(Dat1$Index),]
Dat1[,c("ID","ncase")]
# Dat[Dat$ID == 95,]
Res<-table(Dat1$ID)
length(which(Res > 1))
Res1<-data.frame(Res)
length(which(Res1$Freq>1))
length(which(Res1$Freq==1))

# Dat<-Dat[Dat$population=="European",]
Dat<-Dat[Dat$population=="East Asian",]
IDS<-unique(Dat$ID)
Plot1<-NULL
Plot2<-NULL
for(i in 1:length(IDS))
{

	print(IDS[i])
	Dat1<-Dat[Dat$ID == IDS[i],]	
	head(Dat1)
	
	ID<-unique(Dat1$ID)
	
	# if(any(!is.na(Dat1$eaf)))
	# {
	# 	ref_studies<-c("CHARGE","SCHS")
	# 	ref_pop<-c("European","East Asian")
	# 	Pop<-unique(Dat1$population)
	# 	ref_study<-ref_studies[which(ref_pop %in% Pop)]
	# 	Plot1[[ID]]<-make_plot_maf(refstudy=ref_study,target_dat=Dat1,Title_xaxis_size=12,Title_size=12,Title="",Ylab="",Xlab="")
	# }

	if(any(Dat1$all_summary_stats))
	{
		if(any(!is.na(Dat1$eaf)))
		{
			Plot2[[i]]<-make_plot_maf(refstudy="ALL_1000G",target_dat=Dat1,Title_xaxis_size=8,Title_size=12,cowplot_title=unique(paste0(Dat1$study," | ID:", Dat1$ID)),Ylab="",Xlab="")
		}
	}
}

Pos<-which(unlist(lapply(1:length(Plot2),FUN=function(x) !is.null(unlist(Plot2[x])))))

# for European studies
length(Pos)/4
length(Plot2)/4

make_cow_plot(Plot_list=Plot2[Pos[1:25]],out_file="~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/cowplot_eur_all_1000G_1.png",Xlab="Minor allele frequency in 1000 genomes",Ylab="Allele frequency in cancer dataset",width=1500)
make_cow_plot(Plot_list=Plot2[Pos[26:50]],out_file="~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/cowplot_eur_all_1000G_2.png",Xlab="Minor allele frequency in 1000 genomes",Ylab="Allele frequency in cancer dataset",width=1500)
make_cow_plot(Plot_list=Plot2[Pos[51:75]],out_file="~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/cowplot_eur_all_1000G_3.png",Xlab="Minor allele frequency in 1000 genomes",Ylab="Allele frequency in cancer dataset",width=1500)
make_cow_plot(Plot_list=Plot2[Pos[76:103]],out_file="~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/cowplot_eur_all_1000G_4.png",Xlab="Minor allele frequency in 1000 genomes",Ylab="Allele frequency in cancer dataset",width=1500)

# for East Asian studies
make_cow_plot(Plot_list=Plot2[Pos],out_file="~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/cowplot_eas_all_1000G.png",Xlab="Minor allele frequency in 1000 genomes",Ylab="Allele frequency in cancer dataset",width=1500)