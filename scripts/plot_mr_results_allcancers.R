rm(list=ls())
# source("~/fatty-acids-mr/mr/plot_mr_results_functions.R")
library(ggforestplot)
library(ggplot2)
# install.packages("meta")
library(meta)

##########################################################
#exclude subtypes and use GWAS sig plus replicated instruments####
############################################



# lapply(1:length(unique(Res$exposure)),FUN=function(x)
# 	summary(Res$power[which(as.numeric(as.factor(Res$exposure))==x)]))

# unique(Res$exposure)
# Dat<-Res

# exposure<-exposures[1]
# plot_function(Dat=Res,exposure=exposures[1],prune_ncases=100,exclude=c("COR:49","COR:31"),prune_duplicate_consortia=TRUE,prune_duplicates_outcome=TRUE,Order="b")

# plot_function<-function(Dat=NULL,exposure=NULL,prune_ncases=NULL,power=NULL,exclude=NULL,Name=NULL,prune_duplicates_outcome=FALSE,force_include=NULL,Order=NULL,prune_duplicate_consortia=FALSE){

load("~/fatty-acids/results/mr/mr_results_Cancer_dat_cleaned.Rdata")
# r2_tab$snps

exposures<-c("AA_to_DGLA","POA_to_PA","GLA_to_LA","DHA_to_DPA_n3")

Exclude<-c("Non-glioblastoma glioma COR:49", "Glioblastoma COR:31")

Res1<-format_plot_dat2(Dat=Res,exposure=exposures[1:4], exclude=c("Lung cancer || id:966","Lung cancer adjusted for chip COR:41"),Order="ncase",cancer_site= "Upper gastrointestinal cancer")
Cancer_sites<-unique(Res1$cancer_site)
# Pos<-Res1$outcome=="Acute lymphoblastic leukemia COR:1" & Res1$exposure=="DHA_to_DPA_n3"
# Res1<-Res1[!Pos,]
# Res1[Res1$original_outcome=="Breast cancer",]
Plot<-make_plot(Dat=Res1,Colour="exposure",meta_analysis=TRUE)
pdf(unlist(Plot[2]))
	Plot[1]
dev.off()

Meta<-Plot[3]
Meta

Tab<-table(Res1$original_outcome)
which(Tab>4)

Res1$outcome_study

# Res1$b[Res1$id.outcome == "COR:47"]<-Res1$b[Res1$id.outcome == "COR:47"]*-1
# Res1<-Res1[Res1$id.outcome %in% c("COR:34","COR:3","COR:13"),]
Exposures<-unique(Res1$exposure)
Meta<-lapply(1:4,FUN=function(x)
	metagen(TE=Res1$b[Res1$exposure==Exposures[x]],seTE=Res1$se[Res1$exposure==Exposures[x]],comb.fixed=T,sm="MD"))

# Haemotological cancers


# all cancers single largest per site
Res1<-format_plot_dat2(Dat=Res,exposure=exposures[1:4], exclude=c("Non-glioblastoma glioma COR:49", "Glioblastoma COR:31","Lung cancer || id:966","Lung cancer adjusted for chip COR:41"),prune_duplicates_outcome=TRUE,Order="b",prune_duplicate_consortia=TRUE)

########################
# lung cancer results###
###########################

# one exposure
Res1<-format_plot_dat2(Dat=Res,exposure=exposures[1], exclude=c("Non-glioblastoma glioma COR:49", "Glioblastoma COR:31",       "Lung cancer || id:966","Lung cancer adjusted for chip COR:41"),prune_duplicates_outcome=FALSE,Order="b",prune_duplicate_consortia=FALSE,outcomes="Lung cancer",Name="outcome_study")

# all exposures all studies
Res1<-format_plot_dat2(Dat=Res,exclude=c("Non-glioblastoma glioma COR:49", "Glioblastoma COR:31",       "Lung cancer || id:966","Lung cancer adjusted for chip COR:41"),prune_duplicates_outcome=FALSE,Order="b",prune_duplicate_consortia=FALSE,outcomes="Lung cancer",Name="outcome_study")
Plot<-make_plot(Dat=Res1,Colour="exposure",meta_analysis=TRUE)
pdf(unlist(Plot[2]))
	Plot[1]
dev.off()

# all exposures all subtypes
Res1<-format_plot_dat2(Dat=Res,exclude=c("Non-glioblastoma glioma COR:49", "Glioblastoma COR:31","Lung cancer || id:966","Lung cancer adjusted for chip COR:41"),prune_duplicates_outcome=TRUE,Order="ncase",prune_duplicate_consortia=FALSE,cancer_site="Lung cancer",Name="original_outcome")

Plot<-make_plot(Dat=Res1,Colour="exposure",meta_analysis=TRUE)
pdf(unlist(Plot[2]))
	Plot[1]
dev.off()



make_plot<-function(Dat=NULL,Title.plot=NULL,Colour=NULL,Shape=NULL,text.title=1,text.names=7,meta_analysis=FALSE){
	
	if(length(unique(Dat$exposure))==4){
		Exposure<-"allexp"
	}else{
		Exposure<-paste(unique(Dat$exposure),collapse="_")
	}

	if(!is.null(unique(Dat$analysisinfo))){
		Info<-unlist(strsplit(unique(Dat$analysisinfo),split=" "))
		Info<-paste(unique(Info),length(which(Info==Info[!duplicated(Info)])),sep="_v")
		Exposure<-paste(Exposure,Info,sep="_")
	}

	Exposure<-gsub(" ","_",Exposure)
	Exposure<-gsub("/","_",Exposure)
	Exposure<-gsub("__","_",Exposure)

	Title.file<-paste("~/fatty-acids/results/mr/plots/",Exposure,"_",unique(Dat$power_prune),".pdf",sep="")
	Title.file<-gsub("_.pdf",".pdf",Title.file)
	Title.file<-gsub("__","_",Title.file)

	if(meta_analysis){
		Meta<-metagen(TE=Dat$b,seTE=Dat$se,comb.fixed=T,sm="MD")
	}
	# text.title<-""
	# pdf(Title.file)
	if(!is.null(Colour)){
		Plot<-forestplot(df = Dat,
			 logodds = TRUE,
			 # name=original_outcome,
						  estimate=b,
						  se=se,
						  shape=Shape,
						  colour = eval(parse(text=Colour)),
						   xlab = "")+
				 # labs(title=Title.plot,size=1)+
				 theme(plot.title = element_text(size = text.title))+
				 theme(text = element_text(size=text.names))
	}else{
		Plot<-forestplot(df = Dat,
		logodds = TRUE,
		# name=original_outcome,
			  estimate=b,
			  se=se,
			  shape=Shape,
			  colour = NULL,
			   xlab = "")+
		# labs(title=Title.plot,size=1)+
		theme(plot.title = element_text(size = text.title))+
		theme(text = element_text(size=text.names))
	}
	if(!meta_analysis){
		Meta<-"no meta analysis"
	}
	 return(list(Plot,Title.file,Meta))
	}
	
}

# # prune_ncases=1000,power=0.5

# Res2<-Res1[Res1$exposure == "AA_to_DGLA", ]
# Res3<-Res1[Res1$exposure ==  "GLA_to_LA", ]
# Res4<-Res1[Res1$exposure == "DHA_to_DPA_n3", ]
# Res5<-Res1[Res1$exposure == "POA_to_PA" , ]

# Res1[Res1$original_outcome == "Pancreatic cancer",]
# unique(Res1$outcome[which(Res1$ncase<1000)])
# sort(Res2$power)
# sort(Res3$power)
# sort(Res4$power)
# sort(Res5$power)

# Traits.include<-Res$exposure[2]
# Force_include = c("ER+ Breast cancer" ,"ER- Breast cancer")
# # unique(Res$original_outcome[Res$outcome_study == "BCAC"])
# Res<-Res[!Res$id.outcome %in% c("COR:49","COR:31"),]
# Res<-format_plot_data(Dat=Res,prune_duplicates=TRUE,Name="original_outcome",include_cases_in_name=TRUE,traits_include=Traits.include,prune_duplicates_outcome=TRUE,prune_duplicate_consortia=TRUE,force_include=Force_include)
# Plot<-power_prune(mr_res=Res,threshold=c("min_power",0.5))
# Plot[order(Plot$pval2),c("original_outcome","pval2","exposure","nsnps")]

# length(unique(Plot$outcome))



power_prune<-function(mr_res=NULL,threshold=NULL){
	Outcomes<-unique(mr_res$id.outcome)
	Dat<-NULL

	for(i in  1:length(Outcomes)){
		print(Outcomes[i])
		mr_res2<-mr_res[mr_res$id.outcome == Outcomes[i],]
		Med<-round(median(mr_res2$power),2)
		Min<-round(min(mr_res2$power),2)
		Max<-round(max(mr_res2$power),2)
		Dat[[i]]<-data.frame(matrix(c(Min,Med,Max,Outcomes[i]),nrow=1,ncol=4,byrow=F),stringsAsFactors=F)
	}
	Dat<-do.call(rbind,Dat)
	names(Dat)<-c("min_power","med_power","max_power","id.outcome")
	Dat$min_power<-as.numeric(Dat$min_power)
	Dat$max_power<-as.numeric(Dat$max_power)
	Dat$med_power<-as.numeric(Dat$med_power)
	mr_res3<-merge(mr_res,Dat,by="id.outcome")
	mr_res3<-mr_res3[!is.na(mr_res3[,threshold[1]]),] #some studies are missing numnber of cases and controls. CRC in females and males. And oral cancer and oropharyngeal

	mr_res3<-mr_res3[mr_res3[,threshold[1]] > as.numeric(threshold[2]),]
	mr_res3<-mr_res3[order(mr_res3$ncase,decreasing=T),]
	return(mr_res3)
}



format_plot_dat2<-function(Dat=NULL,exclude=NULL,exposure=NULL,Name=NULL,prune_duplicates_outcome=FALSE,force_include=NULL,Order=NULL,prune_duplicate_consortia=FALSE,prune_ncases=NULL,prune_power=NULL,prune_ocac=TRUE,min_power=NULL,organ.system=NULL,outcomes=NULL,cancer_site=NULL,cancer_site2=NULL){

	Pos<-which(is.na(Dat$outcome_study))
	Dat$outcome_study[is.na(Dat$outcome_study)]<-paste("study",999+1:length(Pos),sep="")
	Dat$ncase[is.na(Dat$ncase)]<-9.9
	if(prune_ocac){
		Dat1<-Dat[which(Dat$outcome_study == "OCAC"),  ]
		Dat1<-Dat1[Dat1$original_outcome %in% c("Ovarian cancer","Endometrioid ovarian cancer", "High grade and low grade serous ovarian cancer", " Mucinous ovarian cancer: invasive and low malignant potential","Clear cell ovarian cancer","Low malignant potential ovarian cancer"),]
		Dat2<-Dat[which(Dat$outcome_study != "OCAC"),  ]
		Dat<-rbind(Dat1,Dat2)
	}
	
	if(!is.null(exclude)){
		Dat<-Dat[!Dat$outcome %in% exclude,]
	}

	Dat$exposure<-as.character(Dat$exposure)
	
	Dat$snp<-as.character(Dat$snp)
	if(!is.null(prune_ncases)){
		Dat<-Dat[which(Dat$ncase>prune_ncases),]
	}

	if(!is.null(exposure)){
		Dat<-Dat[Dat$exposure %in% exposure,]
	}

	Dat$original_outcome[which(Dat$outcome_study == "GAMEON")]<-"Cancer (5 sites)"

	Outcomes.exclude<-c("Illnesses of siblings: Breast cancer || id:UKB-b:12227","Illnesses of mother: Breast cancer || id:UKB-b:13584","Illnesses of father: Lung cancer || id:UKB-b:14521","Illnesses of siblings: Lung cancer || id:UKB-b:15826","Illnesses of siblings: Prostate cancer || id:UKB-b:16522","Illnesses of mother: Lung cancer || id:UKB-b:20176","Illnesses of father: Prostate cancer || id:UKB-b:7773","Prostate cancer Gleason score || id:1242","Breast cancer (Survival) || id:1165","ER- Breast cancer (Survival) || id:1163","ER+ Breast cancer (Survival) || id:1164","Interstitial lung disease COR:39")
	Dat<-Dat[!Dat$outcome %in% Outcomes.exclude,]

	Dat$name<-paste(Dat$original_outcome,"\nN cases=",Dat$ncase,sep="")
	if(!is.null(Name)){
		Dat$name<-paste(Dat[,Name],"\nN cases=",Dat$ncase,sep="")	
	}

	if(prune_duplicates_outcome){
		Dat<-Dat[order(Dat$ncase,decreasing=T),]
		Dat<-Dat[!duplicated(paste(Dat$exposure,Dat$original_outcome)),]
		# Dat<-Dat[!Dat$original_outcome %in% c("Lung cancer in father","Breast cancer in mother","Prostate cancer in father","Lung cancer in mother","Breast cancer in siblings","Lung cancer in siblings","Prostate cancer in siblings"),]
# unique(Cancer_mrbase$outcome)
	}	

	if(prune_duplicate_consortia){ #main purpose is to prune out subtypes (works if they come from the same consortium
		Dat<-Dat[order(Dat$ncase,decreasing=T),]
		Dat1<-Dat[Dat$outcome_study %in% c("OCAC","BCAC","PRACTICAL","TRICL","GECCO","GICC/MDA","INHANCE"),]
		# "UK Biobank", "MRC-IEU"
		Dat1<-Dat1[!duplicated(paste(Dat1$outcome_study,Dat1$exposure)),]
		Dat2<-Dat[!Dat$outcome_study %in% c("OCAC","BCAC","PRACTICAL","TRICL","GECCO","GICC/MDA","INHANCE"),]
		if(!is.null(force_include)){
			Dat3<-Dat[Dat$original_outcome %in% force_include,]
			Dat<-rbind(do.call(rbind,list(Dat1,Dat2,Dat3)))			
		}else{
			Dat<-rbind(do.call(rbind,list(Dat1,Dat2)))			
		}			
	}	

	if(!is.null(Order)){
		Dat<-Dat[order(Dat[,Order],decreasing=T),]
	}

    ########################
	# power pruning section#
	########################
	
	if(!is.null(prune_ncases)){
		Dat<-Dat[which(Dat$ncase>prune_ncases),]
	}

	if(is.null(min_power)){
		if(!is.null(prune_power)){
			Dat<-Dat[which(Dat$power>prune_power),]
		}
	}else{
		Dat<-Dat[which(Dat$power>min_power),]
		Tab<-table(Dat$outcome)
		N<-max(Tab)
		Names_keep<-names(Tab)[Tab==N]
		Dat<-Dat[Dat$outcome %in% Names_keep,]
	}
		
	if(!is.null(prune_power)){
		Dat$power_prune=paste("power",prune_power,sep="_")	
	}

	if(!is.null(prune_ncases)){
		Dat$power_prune=paste("ncases",prune_ncases,sep="_")		
	}

	if(!is.null(prune_power) & !is.null(prune_ncases)){
		Dat$power_prune=paste("power_",prune_power,"_ncase_",prune_ncases,sep="")	
	}

	if(!is.null(min_power) & is.null(prune_ncases)){
		Dat$power_prune=paste("minpower_",min_power,sep="")
	}

		if(!is.null(min_power) & !is.null(prune_ncases)){
		Dat$power_prune=paste("minpower_",min_power,"_ncase_",prune_ncases,sep="")
	}


	Dat$gene<-NA
	Dat$gene[Dat$exposure == "POA_to_PA"]<-"SCD"
	Dat$gene[Dat$exposure == "GLA_to_LA"]<-"FADS2"
	Dat$gene[Dat$exposure == "AA_to_DGLA"]<-"FADS1"
	Dat$gene[Dat$exposure == "DHA_to_DPA_n3"]<-"ELOVL2"

	# systems

	Dat$system<-NA
	Dat$system[Dat$original_outcome %in% c("Acute lymphoblastic leukemia","B cell non-Hodgkin lymphoma" , "Chronic myeloid leukemia","Hematological malignancy" ,"Hodgkin lymphoma" , "Multiple myeloma")]<-"Blood"

	Dat$system[Dat$original_outcome %in% c("Basal cell skin cancer" ,"Melanoma","Squamous cell skin cancer" ,"Uveal melanoma")]<-"Skin"

	Dat$system[Dat$original_outcome %in% c("Biliary tract cancer", "Bladder cancer"  , 
		"Colorectal cancer","Colorectal cancer (distal)" , "Colorectal cancer in males" , "Colorectal cancer in females" , "Colon cancer"    ,"Rectal cancer" ,"Colorectal cancer (proximal)"  ,
		"Esophageal adenocarcinoma",    "Esophageal cancer","Gallbladder cancer","Gastric cancer","Head and neck cancer" ,"Oropharyngeal cancer","Oral cancer"   ,
		"Hepatocellular carcinoma","Kidney cancer", "Pancreatic cancer")]<-"Digestive system"

	Dat$system[Dat$original_outcome %in% c("Lung cancer","Lung cancer in ever smokers" ,"Lung adenocarcinoma"  ,"Squamous cell lung cancer","Small cell lung carcinoma"    ,"Lung cancer in never smokers"            )]<-"Lung"

	Dat$system[Dat$original_outcome %in% c("Breast cancer","ER- Breast cancer", "ER+ Breast cancer" ,
		"Cervical cancer" ,"Endometrial cancer",
		"Ovarian cancer","Clear cell ovarian cancer" , "Endometrioid ovarian cancer"  , "Low malignant potential ovarian cancer"  ,"High grade and low grade serous ovarian cancer",
		"Prostate cancer","Prostate cancer (early-onset)"  ,"Prostate cancer (high aggressive vs low/intermediate aggressive)", "Prostate cancer (advanced)"     ,"Prostate cancer (advanced vs non-advanced)" ,"Prostate cancer (high aggressive vs low aggressive)"   ,
		 "Thyroid cancer") ]<-"Reproductive/endocrine system" 

	Dat$system[Dat$original_outcome %in% c("Glioma" ,"Meningioma" ,"Neuroblastoma")]<-"Nervous system"

	Dat$system[Dat$original_outcome %in% c("Cancer","Cancer (5 sites)")]<-"Other"
	
	Dat$cancer_site<-Dat$original_outcome
	Dat$cancer_site[grep("prostate",Dat$cancer_site,ignore.case=T)]<-"Prostate cancer"
	Dat$cancer_site[grep("ovarian",Dat$cancer_site,ignore.case=T)]<-"Ovarian cancer"
	Dat$cancer_site[grep("lung",Dat$cancer_site,ignore.case=T)]<-"Lung cancer"
	Dat$cancer_site[grep("colorectal",Dat$cancer_site,ignore.case=T)]<-"Colorectal cancer"
	Dat$cancer_site[grep("colon",Dat$cancer_site,ignore.case=T)]<-"Colorectal cancer"
	Dat$cancer_site[grep("rectal",Dat$cancer_site,ignore.case=T)]<-"Colorectal cancer"
	Dat$cancer_site[grep("breast",Dat$cancer_site,ignore.case=T)]<-"Breast cancer"
	Dat$cancer_site[which(Dat$system == "Blood")]<-Dat$system[which(Dat$system == "Blood")]
	Dat$cancer_site[which(Dat$system == "Skin")]<-Dat$system[which(Dat$system == "Skin")]
	Dat$cancer_site[grep("Oropharyngeal",Dat$cancer_site,ignore.case=T)]<-"Head & neck"
	Dat$cancer_site[grep("oral",Dat$cancer_site,ignore.case=T)]<-"Head & neck"
	Dat$cancer_site[grep("glioblastoma",Dat$cancer_site,ignore.case=T)]<-"Brain"
	Dat$cancer_site[grep("glioma",Dat$cancer_site,ignore.case=T)]<-"Brain"
	Dat$cancer_site[grep("meningioma",Dat$cancer_site,ignore.case=T)]<-"Brain"
	Dat$cancer_site[unlist(lapply(c("biliary","gallbladder","liver","hepatocellular"),FUN=function(x)
		grep(x,Dat$cancer_site,ignore.case=T)))]<-"Liver cancer"
	Dat$cancer_site[unlist(lapply(c("Esophageal adenocarcinoma","Esophageal cancer","Gastric cancer"),FUN=function(x)
		grep(x,Dat$cancer_site,ignore.case=T)))]<-"Upper gastrointestinal cancer"
	Dat$cancer_site2<-Dat$cancer_site
	Dat$cancer_site2[Dat$cancer_site2 %in% c("Upper gastrointestinal cancer","Head & neck cancer")]<-"Upper gastrointestinal, head & neck"
	Dat$cancer_site2[Dat$cancer_site2 %in% c("Ovarian cancer","Cervical cancer","Endometrial cancer")]<-"Female reproductive cancer"

	if(!is.null(cancer_site)){
		Dat<-Dat[Dat$cancer_site == cancer_site,]
	}

	if(!is.null(cancer_site2)){
		Dat<-Dat[Dat$cancer_site2 == cancer_site2,]
	}	
	
	if(!is.null(organ.system)){
		Dat<-Dat[which(Dat$system %in% organ.system),]
	
	}
	
	# if(!is.null(outcomes)){
	# 	Dat$analysisinfo<-paste(Dat$analysisinfo,"_",unique(outcomes),sep="")
	# }

	if(!is.null(outcomes)){
		Dat<-Dat[Dat$original_outcome %in% outcomes,]
	}
	Info<-c(organ.system,cancer_site,cancer_site2,outcomes)
	Dat$analysisinfo<-paste(gsub(" ","_",Info),collapse=" ")
	return(Dat)
	# if(!include_cases_in_name & !is.null(Name)){
	# 	Dat$name<-Dat[,Name]
	# 	Dat$outcome_study<-paste(Dat$outcome_study," \nNcases=",Dat$ncase,sep="")
	# }

	# if(Name == "exposure"){
	# 	Dat$name<-Dat$exposure
	# }

	# if(include_pval_name){
	# 	Dat$name<-paste(Dat$name,"\nP=",Dat$pval2,sep="")
	# }
}


# b<-temp$beta.sd
# 	se<-temp$se.sd
# 	# p<-temp$p
# 	w<-1/se^2
# 	temp$beta.sd<-sum(b*w)/(sum(w))
# 	temp$se.sd<-sqrt(sum(w)^-1)
# 	z<-abs(temp$beta.sd/temp$se.sd)
# 	temp$p<-pnorm(z,lower.tail=F)*2
# 	temp$pval<-pnorm(z,lower.tail=F)*2
# 	b<-temp$beta
# 	se<-temp$se
# 	w<-1/se^2
# 	temp$beta<-sum(b*w)/(sum(w))
# 	temp$se<-sqrt(sum(w)^-1)
