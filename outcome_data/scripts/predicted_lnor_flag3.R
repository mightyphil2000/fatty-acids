#########################################
# Predicted log odds ratios versus reported effect sizes #
########################################
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")
Dat<-rbind(dat1,dat2)
Dat<-estimate_bias(dat=Dat)
Dat<-fix_info(dat=Dat) # include the SNPs with imputation scores< 0.8  

# 0
Dat<-Dat[Dat$lnor_sh <= 1.999 & Dat$lnor_sh>= -1.999,] #lnor_sh ==1.999 is an artifiact. 

#1 
MAF<-Dat$eaf
MAF[MAF>0.5]<-1-MAF[MAF>0.5]
Dat$MAC_case<-MAF*Dat$ncase*2
Dat$MAC_control<-MAF*Dat$ncontrol*2
Dat<-Dat[Dat$MAC_case>=50 & Dat$MAC_control>=50,]

# 2
Dat<-Dat[Dat$study != "UKB",]
ID_outliers<-c(3,59,63,64,65,66,70,71,#datasets with slope <0.8 25 explained by N studies 
	25,104,132)#imputation info/r2

Dat<-Dat[!Dat$ID %in% ID_outliers,]


# 3
Dat<-Dat[Dat$lnor < 1 & Dat$lnor> -1,]

# Dat<-Dat[Dat$lnor_sh < 1 & Dat$lnor_sh> -1,]
# Dat<-Dat[abs(Dat$bias)<10,]
# head(Dat[Dat$ID == 8,])
# predict_lnor(lnor=0.5524,se=0.2903,n=2810+40941,p=0.1704,cases=2810,controls=40941)
# MAF<-Dat$eaf
# MAF[MAF>0.5]<-1-MAF[MAF>0.5]
# MAC<-MAF*Dat$ncase*2
# Dat<-Dat[MAC>=10,]
# Dat<-Dat[MAF>=0.05,]
IDS<-unique(Dat$ID)
# i<-which(IDS == 23)
Int_list<-NULL
Slope_list<-NULL
ID_list<-NULL
for(i in 1:length(IDS)){
	print(i)	
	# i<-which(IDS==109)
	Dat1<-Dat[Dat$ID == IDS[i],]
	# Dat1<-Dat[Dat$ID == 3,]
	# Dat1<-Dat1[which(Dat1$lnor<0.5),]
	Model<-summary(lm(lnor_sh~lnor,Dat1))
	Int_list[[i]]<-Model$coefficients[1,1]
	Slope_list[[i]]<-Model$coefficients[2,1]
	ID_list[[i]]<-IDS[i]
}
Int<-unlist(Int_list) 
ID<-unlist(ID_list)	
Slope<-unlist(Slope_list)

Pos<-which(Slope > 1.20 | Slope < 0.8)
length(Pos)

IDS2<-ID[Pos]
Dat1<-Dat[Dat$ID %in% IDS2,]

# Dat<-Dat[Dat$summary_set ==  "full summary stats",]
# i<-1
IDS<-unique(Dat1$ID)
IDS<-IDS[order(IDS)]
Plot1<-NULL
for(i in 1:length(IDS)){
	print(i)
	Dat2<-Dat1[Dat1$ID == IDS[i],]
	# Dat2<-Dat1
	Plot1[[i]]<-make_plot_predlnor2(dat=Dat2,Xlab="",Ylab="",Title_size=12,Title=paste0(unique(Dat2$study)," | ID: ",unique(Dat2$ID)),Title_xaxis_size=0,maf_filter=FALSE,legend=FALSE,standard_errors=FALSE)
}

make_cow_plot2(width=2000,height=1000,Title_axis_size=0,bycols=FALSE,
		Plot_list=Plot1,
		Title="",
		Title_size=0,
		Ylab="",
		Xlab="",
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/predicted_lnor/cowplot_shmethod_flagged3.png"))
