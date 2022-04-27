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
length(which(abs(Dat$bias)>100))
summary(abs(Dat$bias))
Dat$ID2<-as.factor(Dat$ID)

# 0
Dat<-Dat[Dat$lnor_sh <= 1.999 & Dat$lnor_sh>= -1.999,] #lnor_sh ==1.999 is an artifiact. 

# Dat<-Dat[!Dat$study %in% "UKB",]

# # 1
# # MAF<-Dat$eaf
# # MAF[MAF>0.5]<-1-MAF[MAF>0.5]
# # Dat<-Dat[MAF>=0.01,]

# MAF<-Dat$eaf
# MAF[MAF>0.5]<-1-MAF[MAF>0.5]
# Dat$MAC_case<-MAF*Dat$ncase*2
# Dat$MAC_control<-MAF*Dat$ncontrol*2
# Dat<-Dat[Dat$MAC_case>=50 & Dat$MAC_control>=50,]

IDS<-unique(Dat$ID)
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

Pos<-which(Slope > 1.25 | Slope < 0.8)
length(Pos)

IDS2<-ID[Pos]
Dat1<-Dat[Dat$ID %in% IDS2,]

IDS<-unique(Dat1$ID)
Plot1<-NULL
for(i in 1:length(IDS)){
	print(i)
	Dat2<-Dat1[Dat1$ID == IDS[i],]
	Dat2<-Dat2[Dat2$lnor> -99,]
	Dat2[Dat2$lnor< -0.30,"rsid"]
	# Dat2<-Dat1
	Plot1[[i]]<-make_plot_predlnor2(dat=Dat2,Xlab="",Ylab="",Title_size=12,Title=paste0(unique(Dat2$study)," | ID: ",unique(Dat2$ID)),Title_xaxis_size=0,maf_filter=FALSE,legend=FALSE,standard_errors=FALSE)
}

make_cow_plot2(width=2000,height=1000,Title_axis_size=0,bycols=FALSE,
		Plot_list=Plot1,
		Title="",
		Title_size=0,
		Ylab="",
		Xlab="",
		out_file=paste0("~/fatty-acids/outcome_data/results/plots/predicted_lnor/cowplot_shmethod_flagged0.png"))




# # plot of slopes by ID
# load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
# load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")
# Dat<-rbind(dat1,dat2)
# Dat<-Dat[which(Dat$lnor < 1 & Dat$lnor > -1),]
# MAF<-Dat$eaf
# MAF[MAF>0.5]<-1-MAF[MAF>0.5]
# Dat<-Dat[MAF>0.1,]
# IDS<-unique(Dat$ID)
# Slope_list<-NULL
# Int_list<-NULL
# ID_list<-NULL
# max(Dat$eaf)
# for(i in 1:length(IDS)){
# 	print(i)	
# 	Dat1<-Dat[Dat$ID == IDS[i],]
# 	Model<-summary(lm(lnor_sh~lnor,Dat1))
# 	Int_list[[i]]<-Model$coefficients[1,1]
# 	Slope_list[[i]]<-Model$coefficients[2,1]
# 	ID_list[[i]]<-IDS[i]
# }
# Int<-unlist(Int_list) 
# ID<-unlist(ID_list)	
# Slope<-unlist(Slope_list)

# # Dat1<-Dat[Dat$ID<500,]

# Pos<-which(ID<500)
# ID<-ID[Pos]
# Slope<-Slope[Pos]
# Int<-Int[Pos]
# Labels<-ID
# Labels[Slope<1.25 & Slope >0.80 ]<-""
# # Dat1[Dat1$slope < 0.7,]
# lnslope<-log(Slope)
# # png("~/fatty-acids/outcome_data/Datults/plots/plot_predor.png", width = 500, height = 500)
# 	plot(lnslope~ID,col="lightblue", pch=19, cex=2,main="Unusual log odds ratios by study ID",ylab="Predicted lnOR regDatsed on reported lnOR (slope)",yaxt="n")
# 	abline(h=log(1), col="red")
# 	axis(2, at=c(-2.30,0,2.3,4.6,6.9),labels=c(0.1,1,10,100,1000))
# 	text(lnslope~ID, labels=Labels,cex=1)
# # dev.off()	


# # load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")
# # # Dat<-Dat[Dat$population == "East Asian",]
# # # Dat<-Dat[Dat$population == "European",]
# # Dat<-Dat[!is.na(Dat$eaf),]
# # IDS<-unique(Dat$ID)
# # Plot1<-NULL
# # for(i in 1:length(IDS)){
# # 	print(i)
# # 	Dat1<-Dat[Dat$ID == IDS[i],]
# # 	Plot1[[i]]<-make_plot_predlnor(dat=Dat1,thDathold=0.1,maf_filter=TRUE,Xlab="Reported log odds ratio",Ylab="Predicted log odds ratio",Title_size=12,Title="Predicted versus reported log odds ratios",Title_xaxis_size=12)
# # }

# # Plot1[1]
