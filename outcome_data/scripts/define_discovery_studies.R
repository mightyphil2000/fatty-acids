# #####################
# Combine all datasets#
#######################
library(plyr)
# library(TwoSampleMR)
# download google sheet
# install.packages('gsheet')
library(gsheet)
# ieugwasr::get_access_token()
# ao<-available_outcomes()
# ao$ncase[ao$trait == "Lung cancer"]

url<-"https://docs.google.com/spreadsheets/d/1sTliqObb08y3rdfB-vsYecjveud1pPk4FTWSYcpx3Xk/edit?usp=sharing"
meta.tab<-data.frame(gsheet2tbl(url))
meta.tab1<-format_dat(dat=meta.tab) #meta data after some basic formatting and excluding lowest quality studies
meta.tab2<-format_dat2(dat=meta.tab1)
meta.tab3<-combine_cancers(dat=meta.tab2)
meta.tab4<-meta.tab3[meta.tab3$cases>1000,]
meta.tab5<-combine_cancer_groups2(dat=meta.tab2,cancers=meta.tab4$cancer,IDS=meta.tab4$ID)
meta.tab6<-rbind(meta.tab4,meta.tab5)
meta.tab7<-combine_nervous2(dat=meta.tab2,IDS=meta.tab6$ID)
meta.tab8<-rbind(meta.tab6,meta.tab7)
lymphoid_lymphomas<-combine_cancers(dat=meta.tab1,IDS=c(44,150)) 
nasop_carcinomas<-combine_cancers(dat=meta.tab1,IDS=c(97,132))
extra_dat<-rbind(lymphoid_lymphomas,nasop_carcinomas) #these cancers have less than 1000 cases
meta.tab8<-rbind(meta.tab8,extra_dat)
disc.tab9<-fads_power2(dat=meta.tab8)
meta.tab9<-format_dat3(dat=meta.tab1) #further formatting to remove studies with include set to 0 but keeping CRC in BJ and UKB

save(list=c("disc.tab9","meta.tab9"),file="~/MR_FattyAcids/data/summary_data/meta_data.Rdata")


# alternative strategy defining discovery studies in terms of expected effect size
# meta.tab<-data.frame(gsheet2tbl(url))
# meta.tab1<-format_dat(dat=meta.tab)
# meta.tab2<-combine_cancers(dat=meta.tab1)
# disc.tab1<-fads_power(dat=meta.tab2,OR=1.10)
# dim(disc.tab1)
# meta.tab3<-combine_cancer_groups(dat=meta.tab1,cancer=disc.tab1$cancer,IDS=disc.tab1$ID)
# dim(meta.tab3)
# disc.tab2<-fads_power(dat=meta.tab3,OR=1.10,drop80=FALSE)
# dim(disc.tab2)
# disc.tab3<-rbind(disc.tab1,disc.tab2)
# disc.tab3<-fads_power2(dat=disc.tab3)
# disc.tab3[,c("cancer","cases","power","power","power10","power15","power20")]
# dim(disc.tab3)
# combine_cancer_groups(dat=meta.tab1,cancers=disc.tab3$cancer,IDS=disc.tab3$ID)

combine_cancers<-function(dat=NULL,IDS=NULL){
	# c(44,150) Lymphoid leukaemia
	# c(97,132) Nasopharyngeal carcinoma 
	dat1<-dat[dat$ID %in% IDS ,]
	dat1<-dat1[order(dat1$cases,decreasing=T),]
	dat1<-dat1[!duplicated(dat1$study.abbreviation),]		
	dat1$cases<-sum(as.numeric(dat1$cases))
	dat1$controls<-sum(as.numeric(dat1$controls))
	dat1$total<-sum(as.numeric(dat1$total))
	dat1$cancer<-unique(dat1$cancer)
	for(j in 1:ncol(dat1)){
		dat1[,j]<-paste(unique(dat1[,j]),collapse="; ")
	}
	dat1<-dat1[!duplicated(dat1),]
	return(dat1)
}

power_dist<-function(dat=disc.tab){
	Temp<-dat[order(dat$power,decreasing=T),c("cancer","power","cases","controls","total","ID")]
	Temp1<-Temp[Temp$power == 1,]
	Temp1<-Temp1[order(Temp1$cases,decreasing=T),]
	Temp2<-Temp[Temp$power != 1,]
	Temp<-rbind(Temp1,Temp2)
	max(Temp$cases)
	min(Temp$cases)
	median(Temp$cases)
	median(Temp$controls)
	Q<-nrow(Temp)/4
	Index<-c(1,Q,Q*2,Q*2+1,Q*3,nrow(Temp))
	Temp2<-Temp[Index,]

	b1<-log(1.05)
	b2<-log(1.10)
	b3<-log(1.15)
	b4<-log(1.20)
	# b2<-log(1.20)
	sig<-0.05 #alpha
	rsq<-c(0.34)
	Temp2$total<-as.numeric(Temp2$total)
	Temp2$cases<-as.numeric(Temp2$cases)
	Temp2$controls<-as.numeric(Temp2$controls)
	n<-Temp2$total
	ratio<-Temp2$cases/Temp2$controls

	Temp2$power05<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))
	Temp2$power10<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b2-qnorm(1-sig/2))
	Temp2$power15<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b3-qnorm(1-sig/2))
	Temp2$power20<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b4-qnorm(1-sig/2))
	return(Temp2)	
}

# the only remaining underpowered cancers (OR=1.15) it makes sense to combine on a system level is nervous system cancers
combine_nervous<-function(){
# cancers<-disc.tab1$cancer
	IDS<-unlist(strsplit(disc.tab$ID,split=";"))
	IDS<-as.numeric(trimws(IDS))
	dat<-meta.tab1[!meta.tab1$ID %in% IDS,]
	# dim(dat)
	dat1<-dat[dat$system=="Nervous",]
	dat1<-dat1[order(dat1$cases,decreasing=T),]
	dat1<-dat1[!duplicated(dat1$study.abbreviation),]		
	dat1$cases<-sum(as.numeric(dat1$cases))
	dat1$controls<-sum(as.numeric(dat1$controls))
	dat1$total<-sum(as.numeric(dat1$total))
	dat1$cancer<-"Central nervous system and eye cancer"
	for(j in 1:ncol(dat1)){
		dat1[,j]<-paste(unique(dat1[,j]),collapse="; ")
	}
	dat1<-dat1[!duplicated(dat1),]
	return(dat1)
}

combine_nervous2<-function(dat=NULL,IDS=NULL){
	IDS<-unlist(strsplit(IDS,split=";"))
	IDS<-as.numeric(trimws(IDS))
	dat1<-dat[!dat$ID %in% IDS,]
	dat1<-dat1[dat1$system=="Nervous",]	
	dat1<-dat1[order(dat1$cases,decreasing=T),]
	dat1<-dat1[!duplicated(dat1$study.abbreviation),]		
	dat1$cases<-sum(as.numeric(dat1$cases))
	dat1$controls<-sum(as.numeric(dat1$controls))
	dat1$total<-sum(as.numeric(dat1$total))
	dat1$cancer<-"Central nervous system and eye cancer"
	for(j in 1:ncol(dat1)){
		dat1[,j]<-paste(unique(dat1[,j]),collapse="; ")
	}
	dat1<-dat1[!duplicated(dat1),]
	dat1$cases<-as.numeric(dat1$cases)
	dat1$controls<-as.numeric(dat1$controls)
	dat1$total<-as.numeric(dat1$total)	
	dat1<-dat1[dat1$cases>1000,]		
	return(dat1)
}



# combine cancers with fewer than 1000 cases together within the same cancer group 
combine_cancer_groups2<-function(dat=NULL,cancers=NULL,IDS=NULL){
	IDS<-unlist(strsplit(IDS,split=";"))
	IDS<-as.numeric(trimws(IDS))
	dat1<-dat[!dat$ID %in% IDS,]
	dat2<-dat1[!dat1$Cancer.Group %in% cancers,]
	dat1<-dat2
	Dups<-unique(dat1$Cancer.Group[duplicated(dat1$Cancer.Group)])
		# disc.tab9$cancer[disc.tab9$cancer %in% Dups]
	dat2<-dat1[dat1$Cancer.Group %in% Dups,]
	# dat2[order(dat2$Cancer.Group),c("study.abbreviation","cancer","Cancer.Group","cases")]
	groups<-unique(dat2$Cancer.Group)
	Dat_list<-NULL
	for(i in 1:length(groups)){
		print(groups[i])
		Dat<-dat2[dat2$Cancer.Group == groups[i],]	
		Dat<-Dat[order(Dat$cases,decreasing=T),]
		Dat<-Dat[!duplicated(Dat$study.abbreviation),]		
		Dat$cases<-sum(as.numeric(Dat$cases))
		Dat$controls<-sum(as.numeric(Dat$controls))
		Dat$total<-sum(as.numeric(Dat$total))
		Dat$cancer<-groups[i]
		for(j in 1:ncol(Dat)){
			Dat[,j]<-paste(unique(Dat[,j]),collapse="; ")
		}
		Dat_list[[i]]<-Dat[!duplicated(Dat),]		
	}
	Dat1<-do.call(rbind,Dat_list)
	Dat1$cases<-as.numeric(Dat1$cases)
	Dat1$controls<-as.numeric(Dat1$controls)
	Dat1$total<-as.numeric(Dat1$total)	
	Dat1<-Dat1[Dat1$cases>1000,]		
	return(Dat1)
}

combine_cancer_groups<-function(dat=NULL,cancers=NULL,IDS=NULL){	
	IDS<-unlist(strsplit(IDS,split=";"))
	IDS<-as.numeric(trimws(IDS))
	dat<-dat[!dat$ID %in% IDS,]
	dat1<-dat[!dat$Cancer.Group %in% cancers,]
	oth.dat<-dat[dat$Cancer.Group %in% cancers,]

	Dups<-unique(dat1$Cancer.Group[duplicated(dat1$Cancer.Group)])
	dat2<-dat1[dat1$Cancer.Group %in% Dups,]
	oth.dat<-dat1[!dat1$Cancer.Group %in% Dups,]
	dat2[order(dat2$Cancer.Group),c("study.abbreviation","cancer","Cancer.Group","cases")]
	groups<-unique(dat2$Cancer.Group)
	Dat_list<-NULL
	for(i in 1:length(groups)){
		print(groups[i])
		Dat<-dat2[dat2$Cancer.Group == groups[i],]	
		Dat<-Dat[order(Dat$cases,decreasing=T),]
		Dat<-Dat[!duplicated(Dat$study.abbreviation),]		
		Dat$cases<-sum(as.numeric(Dat$cases))
		Dat$controls<-sum(as.numeric(Dat$controls))
		Dat$total<-sum(as.numeric(Dat$total))
		Dat$cancer<-groups[i]
		for(j in 1:ncol(Dat)){
			Dat[,j]<-paste(unique(Dat[,j]),collapse="; ")
		}
		Dat_list[[i]]<-Dat[!duplicated(Dat),]		
	}
	Dat1<-do.call(rbind,Dat_list)
	return(Dat1)
}

fads_power<-function(dat=NULL,OR=NULL,drop80=TRUE){
	b1<-log(OR)
	# b2<-log(1.20)
	sig<-0.05 #alpha
	rsq<-c(0.34)
	dat$total<-as.numeric(dat$total)
	dat$cases<-as.numeric(dat$cases)
	dat$controls<-as.numeric(dat$controls)
	n<-dat$total
	ratio<-dat$cases/dat$controls

	dat$power<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))
	if(drop80){
		dat<-dat[dat$power>=0.8,]
	}
	return(dat)
		
}

fads_power2<-function(dat=NULL){
	# dat<-dat[as.numeric(dat$cases)>ncases,]
	b1<-log(1.25)
	b2<-log(1.20)
	b3<-log(1.15)
	b4<-log(1.10)
	b5<-log(1.05)
	b6<-log(1.01)
	# b2<-log(1.20)
	sig<-0.05 #alpha
	rsq<-c(0.34)
	dat$total<-as.numeric(dat$total)
	dat$cases<-as.numeric(dat$cases)
	dat$controls<-as.numeric(dat$controls)
	n<-dat$total
	ratio<-dat$cases/dat$controls

	dat$power25<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))
	dat$power20<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b2-qnorm(1-sig/2))
	dat$power15<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b3-qnorm(1-sig/2))
	dat$power10<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b4-qnorm(1-sig/2))
	dat$power05<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b5-qnorm(1-sig/2))
	dat$power.order<-pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b6-qnorm(1-sig/2))
	dat<-dat[order(dat$power.order,decreasing=T),]
	dat$order<-1:nrow(dat)	
	# dat1<-dat[dat$power>=0.8,]
	return(dat)
}

First.upper<-function(name=NULL){
		paste(toupper(substr(name, 1, 1)), substr(name, 2, nchar(name)), sep="")
	}

format_dat<-function(dat=NULL){
	remove_columns<-c("include_site_discovery", "meta_analyse_cancer",    "site.ind",  "power80","power15","power20","power_include","power","INITIAL.SAMPLE.SIZE"  ,  "MAPPED_TRAIT"   ,        "MAPPED_TRAIT_URI" )	
	dat<-dat[,!names(dat) %in% remove_columns]
	dat$include[dat$cancer=="Male genital cancer"  ]<-0 #controls include males and females
	dat$include[dat$ID == 43]<-0 #FInGen Malignant neoplasm of bronchus and lung 43 ID 42 slightly larger 
	dat$include[dat$cancer=="Male genital cancer"  ]<-0 #controls include males and females
	dat$include[dat$cancer=="Malignant neoplasm of breast" & dat$study.abbreviation == "FinnGen"]<-0 #includes males and females

	dat$cancer<-tolower(dat$cancer)
	dat$cancer<-First.upper(dat$cancer)	
	dat$Cancer.Group<-tolower(dat$Cancer.Group)
	dat$Cancer.Group<-First.upper(dat$Cancer.Group)	
	dat$site<-tolower(dat$site)
	dat$site<-First.upper(dat$site)	
	dat$system<-tolower(dat$system)
	dat$system<-First.upper(dat$system)	
	dat$cell<-tolower(dat$cell)
	dat$cell<-First.upper(dat$cell)			
	dat<-dat[dat$MAC100rs174546 == "PASS",]
	dat$Cancer.Group[dat$Cancer.Group=="Stomach cancer"]<-"Gastric adenocarcinoma"
	dat$study[dat$study.abbreviation == "PanScan I+II"]<-"Pancreatic Cancer Cohort Consortium I+II"
	dat$study[dat$study.abbreviation == "PanScan III"]<-"Pancreatic Cancer Cohort Consortium III"
	names(dat)[names(dat) == "PMID"]<-"pmid"
	dat$cancer[dat$cancer == "Er- breast cancer"]<-"ER- breast cancer"
	dat$cancer[dat$cancer == "Er+ breast cancer"]<-"ER+ breast cancer"
	 

	# cancer<-c("Glioblastoma","Non-glioblastoma glioma")
	# study<-c("GICC/MDA","GICC/MDA")
	# pmid<-c(28346443,28346443)

	# outcomes<-c("Prostate cancer Gleason score","Prostate cancer (high aggressive vs low/intermediate aggressive)","Prostate cancer (high aggressive vs low aggressive)" ,"Prostate cancer (advanced vs non-advanced)","Breast cancer (Survival)","ER+ Breast cancer (Survival)","ER- Breast cancer (Survival)")


	return(dat)
}

# remove duplicates and overlapping studies prior to combining cancers and estimate power 
# only want tumours or tumour subtypes not stratified by sex or smoking status
format_dat2<-function(dat=NULL){
	dat<-dat[dat$include ==1,]	
	dat<-dat[order(dat$cases,decreasing=T),]
	dat1<-dat[dat$cancer == "Female genital cancer",]
	dat2<-dat[dat$cancer != "Female genital cancer",]
	dat2<-dat2[grep("male",dat2$cancer,ignore.case=T,invert=T),]
	dat2<-dat2[grep("smokers",dat2$cancer,invert=T),]
	dat3<-rbind(dat2,dat1)
	# paste(dat$cancer,dat$study)[duplicated(paste(dat$cancer,dat$study))]
	dat3<-dat3[!duplicated(paste(dat3$cancer,dat3$study)),]		
	return(dat3)
}

combine_cancers<-function(dat=NULL){
	cancers<-unique(dat$cancer)
	Dat_list<-NULL
	for(i in 1:length(cancers)){
		print(cancers[i])
		Dat<-dat[dat$cancer == cancers[i],]
		Dat[,c("cancer","study.abbreviation","cases")]
		Dat$cases<-sum(Dat$cases)
		Dat$controls<-sum(Dat$controls)
		Dat$total<-sum(Dat$total)
		for(j in 1:ncol(Dat)){
			Dat[,j]<-paste(unique(Dat[,j]),collapse="; ")
		}
		Dat_list[[i]]<-Dat[!duplicated(Dat),]		
	}
	Dat1<-do.call(rbind,Dat_list)
	Dat1$cases<-as.numeric(Dat1$cases)
	Dat1$controls<-as.numeric(Dat1$controls)
	Dat1$total<-as.numeric(Dat1$total)
	return(Dat1)
}


format_dat3<-function(dat=NULL){
	# IDS<-as.numeric(IDS[grep(";",IDS,invert=T)])
	# dat<-dat[!dat$ID %in% IDS, ]
	dat1<-dat[dat$ID %in% c(12,143),] #CRC in Biobank Japan [12] very likely overlaps with ACCC. UK Biobank cRC [143] overlaps with GECCO. excluded from discovery studies but include in followup analyses
	dat<-dat[order(dat$cases,decreasing=T),]
	dat<-dat[dat$include ==1,]	
	Dups<-paste(dat$cancer,dat$study)[duplicated(paste(dat$cancer,dat$study))] #FinGen DUp for overall cancer. include for sensitivty analyses because huge discrepancy in number of cases 
	# dat2<-dat1[which(!paste(dat1$cancer,dat1$study) %in% Dups),]	
	dat2<-rbind(dat,dat1)	
	return(dat2)
}
