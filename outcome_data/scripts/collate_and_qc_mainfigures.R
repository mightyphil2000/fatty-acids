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

# IDS<-74
# IDS<-c(161,3,6,5,132,967,67,133,106,164,134,64,162)
IDS<-c(967,67,106,5,133,68,132,74)
i<-which(IDS==967)
IDS<-unique(Dat$ID)
# IDS<-132
for(i in 1:length(IDS))
{

	print(IDS[i])
	Dat1<-Dat[Dat$ID == IDS[i],]	
	head(Dat1)
	
	ID<-unique(Dat1$ID)
	
	if(any(!is.na(Dat1$eaf)))
	{
		ref_studies<-c("CHARGE","SCHS")
		ref_pop<-c("European","East Asian")
		Pop<-unique(Dat1$population)
		ref_study<-ref_studies[which(ref_pop %in% Pop)]
		Plot1<-make_plot_maf(refstudy=ref_study,target_dat=Dat1,Title_xaxis_size=12,Title_size=12,Title=paste0("Comparison of MAF between target study and ",ref_study),Ylab="MAF in target study",Xlab=paste0("MAF in ",ref_study))
	}

	if(any(Dat1$all_summary_stats))
	{
		if(any(!is.na(Dat1$eaf)))
		{
			Plot2<-make_plot_maf(refstudy="ALL_1000G",target_dat=Dat1,Title_xaxis_size=11,Title_size=12,cowplot_title="Comparison of MAF between target study and 1000 genomes project",Ylab="MAF in target study",Xlab="")
		}
	}

	if(any(Dat1$all_summary_stats))
	{
		Plot3<-make_plot_gwas_catalog_zscores(dat=Dat1,efo=unique(Dat1$efo),Title_size_subplot=12,Title="Comparison of Z scores between target study and GWAS catalog",Ylab="Z score in target study",Xlab="Z score in GWAS catalog",Title_xaxis_size=12)
	}
	
	if(any(Dat1$all_summary_stats))
	{
		if(any(!is.na(Dat1$eaf)))
		{
			Plot4<-make_plot_gwas_catalog_eaf(dat=Dat1,efo=unique(Dat1$efo),Title_size_subplot=12,Title="Comparison of EAF between target study and GWAS catalog",Ylab="EAF in target study",Xlab="EAF in GWAS catalog",Title_xaxis_size=12,	gwas_catalog_ancestral_group=unique(Dat1$population))
		}
	}

	if(any(!is.na(Dat1$eaf)))
	{
		# Dat1<-Dat1[which(Dat1$lnor > -99),]
		# sort(Dat1$lnor)
		# Dat1[Dat1$lnor < -0.3,]
		# Dat2<-Dat1[Dat1$rsid != "rs4927355",]
		Dat2<-prep_data_plot_predlnor_sh(ID=unique(Dat1$ID))
		if(nrow(Dat2)>0){
			Plot5<-make_plot_predlnor2(dat=Dat2,Xlab="Reported log odds ratio",Ylab="Predicted log odds ratio",Title_size=12,Title="Predicted versus reported log odds ratios",Title_xaxis_size=12)

		# MAF<-Dat1$eaf
		# MAF[MAF>0.5]<-1-MAF[MAF>0.5]
		# Dat1$maf<-MAF
		# Res<-pred_lnor_sh(dat2=Dat1)

		}
	}
	# Plot6<-make_plot_predlnor(dat=Dat1,threshold=0.1,maf_filter=TRUE,bias=TRUE,Xlab="Reported log odds ratio",Ylab="% deviation of predicted log odds ratio from reported log odds ratio",Title_size=12,Title="Percentage deviation of predicted log odds ratio from reported log odds ratio",Title_xaxis_size=12)
	Plot7<-predz_vs_obsz(dat=Dat1,Title="ZZ plot",Xlab="Z score inferred from p value",Ylab="Z score inferred from log odds ratio and standard error",Title_xaxis_size=12,Title_size=12)
	# Dat1$z<-abs(Dat1$lnor/Dat1$se)
	# Dat1$zp<-qnorm(Dat1$p/2,lower.tail=F)
	# plot(zp,z)
	# Dat1[Dat1$z>30,]
	# # rm(Plot6)
	Plot_list2<-ls()[grep("Plot[0-9]",ls())] 

	# Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
	# 	!is.null(nrow(eval(parse(text=trait_list[x]))))))
	# rm(list=trait_list[!Test]) #remove the objects with no data
	# # ,envir=.GlobalEnv
	Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
	out_file<-unique(paste0(Dat1$study,"_",Dat1$ID))
	out_file<-gsub("/","_",out_file)
	make_cow_plot2(width=800,height=1000,Title_axis_size=0,
		Plot_list=Plot_list,
		Title="",
		Title_size=0,
		Ylab="",
		Xlab="",
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/cowplot_",out_file,"v2.png"))
	rm(list=Plot_list2)
	rm(Plot_list)
	rm(Plot_list2)
}

Res[Res$ID == 134,]
Res[Res$ID == 161,]

log(0.20)
# Slope of predicted versus reported log odds ratio
IDS<-unique(Dat$ID)
Metrics<-NULL
for(i in 1:length(IDS))
{
	Dat1<-Dat[Dat$ID == IDS[i],]
	if(any(!is.na(Dat1$eaf)))
	{
		Dat1<-Dat1[!is.na(Dat1$eaf),]
		maf<-Dat1$eaf
		maf[maf>0.5]<-1-maf[maf>0.5]
		Dat1<-Dat1[maf>threshold,]
		threshold<-0.1
		Res<-predlnor_model(dat=Dat1)
		
		Res<-c(Res,IDS[i])
		Metrics[[i]]<-Res
	}
}

Res<-do.call(rbind,Metrics)
Res<-data.frame(Res)
names(Res)<-c("intercept","slope","ID")
Res<-merge(Res,unique(Dat[,c("ID","study")]),by="ID")
hist(Res$slope,main="Slope for predicted regressed on reported log odds ratios",xlab="Slope for predicted log odds ratio regressed on reported log odds ratio")

Res1<-Res[Res$ID<500,]
Labels<-Res1$ID
Labels[Res1$slope<1.5 & Res1$slope >0.5 ]<-""
Res1[Res1$slope < 0.8,]
Res1$lnslope<-log(Res1$slope)
# which(Res1$slope<0.8)

png("~/fatty-acids/outcome_data/results/plots/plot_predor.png", width = 500, height = 500)
	plot(lnslope~ID,data=Res1,col="lightblue", pch=19, cex=2,main="Unusual log odds ratios by study ID",ylab="Predicted lnOR regressed on reported lnOR (slope)",yaxt="n")
	abline(h=log(1), col="red")
	axis(2, at=c(-2.30,0,2.3,4.6,6.9),labels=c(0.1,1,10,100,1000))
	text(lnslope~ID, labels=Labels,data=Res1,cex=1)
dev.off()	


