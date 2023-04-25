
clean_data<-function(){
		dat<-dat[dat$FAMRC.identifier !=96.2,]
		dat<-dat[!dat$FAMRC.identifier %in% c(74,106,5),] #QC pipeline
		dat<-dat[dat$FAMRC.identifier!= 58,] #gallbladder cancer study with 41 cases with minor count in cases <50 across all SNPs. this is the only Dataset where all SNPs had mac<50
		# rename lung cancer unadjusted for chip as lung cancer, so that it's not counted as a unique cancer type, as per definition 1 of unique cancer types. see cancer_type()
		dat$cancer[dat$cancer=="Lung cancer unadjusted for chip" ]<-"Lung cancer"
		return(dat)
	}

cancer_system<-function(){
	systems<-unique(dat$system)
	dat2<-NULL
	for(i in 1:length(systems)){
		# i<-1
		dat1<-dat[dat$system == systems[i],]
		dat1$N<-nrow(dat1)
		dat1<-unique(dat1[,c("N","system")])
		dat2[[i]]<-dat1
	}
	dat3<-do.call(rbind,dat2)
	dat3$system[dat3$system == "integumentary"]<-"skin"
	dat3$N<-as.numeric(dat3$N)
	dat3<-dat3[order(dat3$N,decreasing=TRUE),]
	return(dat3)
}

cancer_group<-function(){
	group<-unique(dat$Cancer.group)
	dat2<-NULL
	for(i in 1:length(group)){
		# i<-1
		dat1<-dat[dat$Cancer.group == group[i],]
		dat1$N<-nrow(dat1)
		dat1<-unique(dat1[,c("N","Cancer.group")])
		dat2[[i]]<-dat1
	}
	dat3<-do.call(rbind,dat2)	
	return(dat3)
}

cancer_source<-function(add_nature_editor=FALSE,miss_effect_allele=FALSE){
	# miss_effect_allele include CGEMS and N-UGC datasets that were missing effect allele, other allele and eaf
	dat[dat$study.abbreviation == "UKB" & is.na(dat$Open.GWAS.identifier), ]

	dat$source<-NA
	dat$source[dat$study.abbreviation=="UKB"]<-"Open GWAS (UK Biobank)"
	# dat$source[dat$FAMRC.identifier== 163]<-"GWAS catalog"
	# dat$Downloaded.from.GWAS.catalog.FTP.site[dat$FAMRC.identifier == 163]
	dat$source[dat$Obtained.by.correspondence ]<-"Application/correspondence"
	dat$source[!is.na(dat$Open.GWAS.identifier) & is.na(dat$source)]<-"Open GWAS (consortia/other)"
	dat$source[grep("finn",dat$Open.GWAS.identifier)]<-"OpenGWAS (FinnGen)"
	dat$source[grep("bbj",dat$Open.GWAS.identifier)]<-"OpenGWAS (Biobank Japan)"
	dat$source[dat$Downloaded.from.GWAS.catalog.FTP.site]<-"GWAS catalog"
	dat$source[dat$FAMRC.identifier==128]<-"Open GWAS (consortia/other)"
	source<-unique(dat$source)
	dat2<-NULL
	for(i in 1:length(source)){
		# i<-1
		dat1<-dat[dat$source == source[i],]
		dat1$N<-nrow(dat1)
		dat1<-unique(dat1[,c("N","source")])
		dat2[[i]]<-dat1
	}
	dat3<-do.call(rbind,dat2)	
	if(add_nature_editor){
		dat3<-rbind(dat3,c(5,"Nature editor"))
	}
	dat3$N<-as.numeric(dat3$N)
	dat3<-dat3[order(dat3$N,decreasing=TRUE),]
	if(miss_effect_allele){
		dat3$N[dat3$source == "Application/correspondence"]<-dat3$N[dat3$source == "Application/correspondence"]+2
	}
	return(dat3)
}


cancer_scope<-function(miss_effect_allele=FALSE){
	a<-data.frame(matrix(c(c(130,34),c("All","Subset")),nrow=2))
	names(a)<-c("N","Summary_data_shared")
	a$N<-as.numeric(a$N)
	if(miss_effect_allele){
		a$N[a$Summary_data_shared=="All"]<-a$N[a$Summary_data_shared=="All"]+2 #CGEMS and N_UGC datasets with missing effect allele, other allele and eaf
	}
	return(a)
}

cancer_ancestry<-function(){
	dat<-dat[order(dat$cases,decreasing=TRUE),]
	dat<-dat[!duplicated(dat$study.abbreviation),]
	population<-unique(dat$population)
	dat2<-NULL
	for(i in 1:length(population)){
		# i<-1
		dat1<-dat[dat$population == population[i],]
		dat1$N<-nrow(dat1)
		dat1<-unique(dat1[,c("N","population")])
		dat2[[i]]<-dat1
	}
	dat3<-do.call(rbind,dat2)	
	return(dat3)
}

cancer_metadata<-function(){
	a<-data.frame(matrix(c(c("beta/se/effect_allele","eaf","imputation quality","P_HWE","P_het"),c(164,155,53,6,18)),nrow=5))
	names(a)<-c("data","N")
	a$N<-as.numeric(a$N)
	return(a)
}

cancer_issue<-function(){
	a<-data.frame(matrix(c(c(124,5,10,25),
	c("major issue not identified","metadata error","beta not log odds ratio","unusual log odds ratios")),nrow=4))
	names(a)<-c("N","QC_issue")
	a$N<-as.numeric(a$N)
	return(a)
}

cancer_issue2<-function(){
	a<-data.frame(matrix(c(c(32,5,1,15),
	c("major issue not identified","metadata error","beta not log odds ratio","unusual log odds ratios")),nrow=4))
	names(a)<-c("N","QC_issue")
	a$N<-as.numeric(a$N)
	a<-a[order(a$N,decreasing=TRUE),]
	return(a)
}


cancer_source2<-function(add_nature_editor=FALSE,miss_effect_allele=FALSE){
	if(qc_pipeline_fail){
		dat<-dat[!dat$FAMRC.identifier %in% c(967,106,5),] #QC pipeline
		dat<-dat[dat$FAMRC.identifier!= 58,] #gallbladder cancer study with 41 cases with minor count in cases <50 across all SNPs
	}
	dat<-dat[order(dat$cases,decreasing=TRUE),]
	dat<-dat[!duplicated(dat$study.abbreviation),]
	dat$source<-NA
	dat$source[dat$study.abbreviation=="UKB"]<-"Open GWAS (UK Biobank)"
	dat$source[dat$Obtained.by.correspondence ]<-"Application/correspondence"
	dat$source[!is.na(dat$Open.GWAS.identifier) & is.na(dat$source)]<-"Open GWAS (consortia/other)"
	dat$source[grep("finn",dat$Open.GWAS.identifier)]<-"OpenGWAS (FinnGen)"
	dat$source[grep("bbj",dat$Open.GWAS.identifier)]<-"OpenGWAS (Biobank Japan)"
	dat$source[dat$Downloaded.from.GWAS.catalog.FTP.site]<-"GWAS catalog"
	s2<-dat$study.abbreviation[grep("Open",dat$source)]
	# dat$Open.GWAS.identifier[dat$FAMRC.identifier==106]

	# s1[!s1 %in% s2]
	# s2[!s2 %in% s1]
	source<-unique(dat$source)
	dat2<-NULL

	for(i in 1:length(source)){
		# i<-1
		dat1<-dat[dat$source == source[i],]
		dat1$N<-nrow(dat1)
		dat1<-unique(dat1[,c("N","source")])
		dat2[[i]]<-dat1
	}
	dat3<-do.call(rbind,dat2)	

	if(add_nature_editor){
		dat3<-rbind(dat3,c(5,"Nature editor"))
		dat3$N[which(dat3$source == "Application/correspondence")]<-
		dat3$N[which(dat3$source == "Application/correspondence")]-5
	}

	dat3$N<-as.numeric(dat3$N)
	dat3<-dat3[order(dat3$N,decreasing=TRUE),]
	
	if(miss_effect_allele){
		dat3$N[dat3$source == "Application/correspondence"]<-dat3$N[dat3$source == "Application/correspondence"]+1 #add the CGEMS study 
	}
	return(dat3)
}
# dat2<-dat
cancer_consortia2<-function(qc_pipeline_fail=FALSE){
	# dat<-dat[order(dat$cases,decreasing=TRUE),]
	# dat1<-dat[!duplicated(dat$study.abbreviation),]
	# dat1<-dat[,c("cancer","cases","study.abbreviation")]
	# dat1<-dat1[dat1$cases>500,]	
	length(unique(dat$study.abbreviation))
	if(qc_pipeline_fail){
		dat<-dat[!dat$FAMRC.identifier %in% c(74,106,5),] #QC pipeline
		dat<-dat[dat$FAMRC.identifier!= 58,] #gallbladder cancer study with 41 cases with minor count in cases <50 across all SNPs
	}


	study<-unique(dat$study.abbreviation)
	dim(dat)
	length(study)

	dat3<-NULL
	for(i in 1:length(study)){
		# i<-1
		dat2<-dat[dat$study.abbreviation == study[i],]
		dat2$N<-nrow(dat2)
		dat2<-unique(dat2[,c("study.abbreviation","N")])
		# dat2<-dat2[!duplicated(dat2$Cancer.group),]
		dat3[[i]]<-dat2
	}
	dat4<-do.call(rbind,dat3)	
	dat4<-dat4[,names(dat4) !="Cancer.group"]
	dat4$study.abbreviation[dat4$study.abbreviation == "UKB"]<-"UK Biobank"
	dat4$study.abbreviation[dat4$study.abbreviation == "BJ"]<-"Biobank Japan"
	dat4$N<-as.numeric(dat4$N)
	dat4<-dat4[order(dat4$N,decreasing=TRUE),]
	names(dat4)[names(dat4) == "study.abbreviation"]<-"Study"

	return(dat4)
}

cancer_consortia<-function(qc_pipeline_fail=FALSE){
	if(qc_pipeline_fail){
		dat<-dat[!dat$FAMRC.identifier %in% c(74,106,5),] #QC pipeline
		dat<-dat[dat$FAMRC.identifier!= 58,] #gallbladder cancer study with 41 cases with minor count in cases <50 across all SNPs. this is the only dataset where all SNPs had mac<50
	}

	dat<-dat[order(dat$cases,decreasing=TRUE),]
	# We assume that cases from IL and BJ are all independent within study. but assume overlap in controls and only take the largest instance of controls within study
	dat2<-dat[dat$study.abbreviation %in% c("InterLymph","BJ"),]
	Pos<-dat2$study.abbreviation == "InterLymph"
	cases<-sum(dat2$cases[Pos])
	dat2$cases[Pos]<-cases
	dat2$cancer[Pos]<-"Blood cancers"
	Pos<-dat2$study.abbreviation == "BJ"
	cases<-sum(dat2$cases[Pos])
	dat2$cases[Pos]<-cases
	dat2$cancer[Pos]<-"Multiple cancers"
	dat2$site[Pos]<-"Multiple"
	dat2<-dat2[order(dat2$controls,decreasing=TRUE),]
	dat2<-dat2[!duplicated(dat2$study.abbreviation),]

	# we assume that cases are nested within the largest case group within study, e.g. UK biobank all site spceific cancer cases are nested within overall cancer. 
	dat1<-dat[!duplicated(dat$study.abbreviation),]
	dat1<-dat1[!dat1$study.abbreviation %in% c("InterLymph","BJ"),]


	# dat1<-dat[,c("cancer","cases","study.abbreviation")]
	# dat1<-dat1[dat1$cases>500,]	
	dat1<-rbind(dat1,dat2)
	study<-unique(dat1$study.abbreviation)
	dat3<-NULL
	for(i in 1:length(study)){
		# i<-1
		dat2<-dat1[dat1$study.abbreviation == study[i],]
		dat2$N<-nrow(dat2)
		dat2<-dat2[,c("cancer","site","study.abbreviation","N","cases","controls","Cancer.group")]
		# dat2<-dat2[!duplicated(dat2$Cancer.group),]
		dat3[[i]]<-dat2
	}
	dat4<-do.call(rbind,dat3)	
	dat4<-dat4[,names(dat4) !="Cancer.group"]
	dat4<-dat4[order(dat4$cases,decreasing=TRUE),]
	dat4$study.abbreviation[dat4$study.abbreviation == "UKB"]<-"UK Biobank"
	dat4$study.abbreviation[dat4$study.abbreviation == "BJ"]<-"Biobank Japan"
	names(dat4)<-c("Cancer","site","study","N","cases","controls")
	dat4$site<-gsub("(^[[:alpha:]])", "\\U\\1", dat4$site, perl=TRUE)
	return(dat4)
}

cancer_consortia3<-function(qc_pipeline_fail=FALSE){
	if(qc_pipeline_fail){
		dat<-dat[!dat$FAMRC.identifier %in% c(74,106,5),] #QC pipeline
		dat<-dat[dat$FAMRC.identifier!= 58,] #gallbladder cancer study with 41 cases with minor count in cases <50 across all SNPs. this is the only dataset where all SNPs had mac<50
	}

	dat<-dat[order(dat$cases,decreasing=TRUE),]
	# We assume that cases from IL and BJ are all independent
	dat2<-dat[dat$study.abbreviation %in% c("InterLymph","BJ"),]
	
	# we assume that cases are nested within the largest case group within study, e.g. UK biobank all site spceific cancer cases are nested within overall cancer
	dat1<-dat[!duplicated(dat$study.abbreviation),]
	# dat1<-dat[,c("cancer","cases","study.abbreviation")]
	# dat1<-dat1[dat1$cases>500,]	
	dat1<-rbind(dat1,dat2)
	study<-unique(dat1$study.abbreviation)
	dat3<-NULL
	for(i in 1:length(study)){
		# i<-1
		dat2<-dat1[dat1$study.abbreviation == study[i],]
		dat2$N<-nrow(dat2)
		dat2<-dat2[,c("cancer","site","study.abbreviation","N","cases","controls","Cancer.group")]
		# dat2<-dat2[!duplicated(dat2$Cancer.group),]
		dat3[[i]]<-dat2
	}
	dat4<-do.call(rbind,dat3)	
	dat4<-dat4[,names(dat4) !="Cancer.group"]
	dat4<-dat4[order(dat4$cases,decreasing=TRUE),]
	dat4$study.abbreviation[dat4$study.abbreviation == "UKB"]<-"UK Biobank"
	dat4$study.abbreviation[dat4$study.abbreviation == "BJ"]<-"Biobank Japan"
	names(dat4)<-c("Cancer","site","study","N","cases","controls")
	dat4$site<-gsub("(^[[:alpha:]])", "\\U\\1", dat4$site, perl=TRUE)
	return(dat4)
}

# Dat=dat
cancer_site<-function(Dat=dat){
	Dat<-Dat[order(Dat$cases,decreasing=TRUE),]
	# Dat1<-Dat[!duplicated(Dat$study.abbreviation),]
	# Dat1<-Dat[,c("cancer","cases","study.abbreviation")]
	# Dat1<-Dat1[Dat1$cases>500,]	
	
	sites<-unique(Dat$site)
	Dat3<-NULL
	for(i in 1:length(sites)){
		# i<-1
		Dat2<-Dat[Dat$site == sites[i],]
		Dat2$N<-length(unique(Dat2$cancer))
		Dat2<-unique(Dat2[,c("site","N")])
		# Dat2<-Dat2[!duplicated(Dat2$Cancer.group),]
		Dat3[[i]]<-Dat2
	}
	Dat4<-do.call(rbind,Dat3)	
	Dat4<-Dat4[order(Dat4$N,decreasing=TRUE),]
	Dat4$site<-gsub("(^[[:alpha:]])", "\\U\\1", Dat4$site, perl=TRUE)
	Dat4$site[grep("Respiratory",Dat4$site) ]<-"chest"
	return(Dat4)
}

c("BCAC","PRACTICAL","GECCO/CORECT/CCFR","ILCCO","OCAC","ACCC","ECAC","MMAC","KidRISK","PanScan I+II+PanC4","INHANCE","InterLymph", "EPITHYR")

# Dat=dat
cancer_types<-function(Dat=NULL,qc_pipeline_fail=FALSE){
	if(qc_pipeline_fail){
		Dat<-Dat[!Dat$FAMRC.identifier %in% c(74,106,5),] #QC pipeline
		Dat<-Dat[Dat$FAMRC.identifier!= 58,] #gallbladder cancer study with 41 cases with minor count in cases <50 across all SNPs. this is the only Dataset where all SNPs had mac<50
	}
	
	# how many unique cancer types are there? 
	# "Lung cancer unadjusted for chip" is not a different type to Lung cancer
	excl1<-"Lung cancer unadjusted for chip"  

	#Cancer 5 sites (combined analysis of colorectal cancer, lung cancer, ovarian cancer, breast cancer and prostate cancer). Is this different to overall cancer or any of the include site specific cancers? Not identical  
	excl2<-c("Cancer (5 sites)"  ,"Lung cancer unadjusted for chip")    

	# Add cancers assessed in population subgroups: males, females, never smokers, ever smokers
	# lung cancer in never smokers could be conisdered a different disease to lung cancer in ever smokers
	# colorectal cancer in males different to colorectal cancer females?
	excl3<-c("Cancer (5 sites)"  ,"Lung cancer unadjusted for chip","Lung cancer in ever smokers","Kidney cancer in females", "Colorectal cancer in females","Colorectal cancer in males","Kidney cancer in males","Lung cancer in never smokers")    

	Pos1<-!Dat$cancer %in% excl1
	Pos2<-!Dat$cancer %in% excl2
	Pos3<-!Dat$cancer %in% excl3
	
	N1<-length(unique(Dat$cancer[Pos1]))
	N2<-length(unique(Dat$cancer[Pos2]))	
	N3<-length(unique(Dat$cancer[Pos3]))

	Dat<-Dat[Pos1,]

	cancer_replicates<-as.numeric(table(Dat$cancer))
	length(cancer_replicates[cancer_replicates==2])
	length(cancer_replicates[cancer_replicates==1])
	length(cancer_replicates[cancer_replicates>2])

	Dat<-Dat[order(Dat$cases,decreasing=TRUE),]
	Dat<-Dat[!duplicated(Dat$cancer),]

	return(list(dfn1=N1,dfn2=N2,dfn3=N3,med_dfn1=median(Dat$cases),min_dfn1=min(Dat$cases),	max_dfn1=max(Dat$cases)))
}

total_cases_independent<-function(){
	dat3<-cancer_consortia(qc_pipeline_fail=TRUE) #
	# exclude overlapping studies. GAMEON overlaps with ILCCO, PRACTICAL, BCAC, OCAC, and GECCO
	# PanScan/PanC4 overlap with PanScan I+II+II+PanC4
	a<-dat3[!dat3$study %in% c("PanScan I","PanScan I+II","PanScan III","PanC4","GAME-ON"),]
	dim(a)
	# a<-dat3[!dat3$study %in% c("PanScan I+II+PanC4","GAME-ON"),]
	ca<-sum(a$cases)
	med_ca<-median(a$cases)
	min_ca<-min(a$cases)
	max_ca<-max(a$cases)
	co<-sum(a$controls)
	med_co<-median(a$controls)
	min_co<-min(a$controls)
	max_co<-max(a$controls)
	return(list(total_cases=ca,med_ca=med_ca,min_ca=min_ca,max_ca=max_ca,total_controls=co,med_co=med_co,min_co=min_co,max_co=max_co))
}