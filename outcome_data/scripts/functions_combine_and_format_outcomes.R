
format_final<-function(dat=NULL,dat2=NULL){
	N<-length(dat2$ID[is.na(dat2$ID )])-1
	dat2$ID[is.na(dat2$ID )]<-(999-N):999 #cancer survival studies 
	IDS<-as.numeric(trimws(unlist(strsplit(dat$ID,split=";"))))
	# length(dat2$ID)
	dat3<-dat2[!dat2$ID %in% IDS,]  
	# dat<-rbind.fill(dat,dat3)
	dat<-rbind.fill(dat,dat2)
	dat$proxy[which(dat$proxy == "NA")]<-NA
	for(i in 1:nrow(dat)){
		Temp<-dat[i,grep("proxy",names(dat),ignore.case=T)]
		for(j in 1:length(names(Temp))){
			print(names(Temp[j]))
			Temp[,j][Temp[,j] == "NA"]<-NA
		}
		if(any(!is.na(Temp))) dat$proxy[i] <- TRUE
	}
	dat$proxy[is.na(dat$proxy)]<-FALSE
	dat$proxy<-as.logical(dat$proxy)
	dat$cancer[is.na(dat$cancer)]<-dat$outcome[is.na(dat$cancer)]
	# dat[which(dat$proxy),grep("proxy",names(dat),ignore.case=T)]
	return(dat)
}

meta_analysis<-function(dat=NULL,IDS=IDS){
	IDS<-IDS[grep(";",IDS)]
	Res_list<-NULL	
	for(i in 1:length(IDS)){
		print(IDS[i])
		ID<-trimws(unlist(strsplit(IDS[i],split=";")))
		if(length(ID)==1) stop("ID=1")
		dat1<-dat[dat$ID %in% ID,]
		dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
		# if(sum(grep("FinnGen",dat1$study.abbreviation))!=0) stop("fingen")
		b<-dat1$lnor
		se<-dat1$se
		# p<-temp$p
		w<-1/se^2
		b.fixed<-sum(b*w)/(sum(w))
		se.fixed<-sqrt(sum(w)^-1)
		z<-abs(b.fixed/se.fixed)
		p.fixed<-pnorm(z,lower.tail=F)*2
		nstudies.fixed<-length(b)
		cancer<-unique(dat1$cancer)
		if(length(cancer)!=1) cancer<-unique(dat1$Cancer.Group)
		if(length(cancer)!=1) 	cancer<-unique(paste(dat1$system,"system cancers"))
		if(length(cancer)!=1) stop("length of cancer not 1")
		ids.fixed<-IDS[i]
		cases<-sum(dat1$cases)
		controls<-sum(dat1$controls)
		study<-"Overall fixed effect"
		Q<-sum((b.fixed-b)^2*w)
		df.Q<-length(b)-1		
		Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
		# names(Meta)
		# Meta$Q
		# Meta$pval.Q
		# Meta<-metagen(TE=b,seTE=se,comb.fixed=T,sm="MD")
		# # Meta$TE.fixed
		# Meta$seTE.fixed
		# Meta$pval.fixed
		# Q.p
		# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
		EA<-unique(dat1$Effect.Allele)
		OA<-unique(dat1$Other.Allele)
		EAF<-round(sum((dat1$eaf*w))/sum(w),3)
		if(length(EA)>1) stop("effect allele not consistent across studies")

		# studies<-paste(dat1$study.abbreviation,collapse="; ")
		Cols<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
		for(j in 1:length(Cols)){
			Cols[j]
			Cols[j]<-paste(unique(dat1[,Cols[j]]),collapse="; ")
		}
		dat.matrix<-c(cancer,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,cases,controls,study,Q.p,EA,OA,EAF)
		Res<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
		names(Res)<-c("cancer","lnor","se","p","nstudies","ID","cases","controls","study","Q.p","Effect.Allele","Other.Allele","eaf")
		Col.dat<-data.frame(matrix(Cols,ncol=length(Cols),nrow=1),stringsAsFactors=F)
		names(Col.dat)<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
		Res.c<-cbind(Res,Col.dat)		
		Res_list[[i]]<-Res.c
	}
	Res1<-do.call(rbind,Res_list)
	Res1$p<-as.numeric(Res1$p)
	return(Res1)
	# i4<-which(as.numeric(Res1$Q.p)<0.10 & as.numeric(Res1$Q.p)>0.05)
}




format_meta_dat4<-function(){ #merge meta data with additional uk biobank studies from automated phenotype pipeline (phesant)
	load("~/MR_FattyAcids/data/summary_data/mrbase_extra/ukb_mrbase.Rdata")
	ukb_dat1<-format_columns(dat=ukb_dat,study="UK Biobank") #from bolt lmm without correction of the betas 
	ukb_dat1$outcome[ukb_dat1$outcome ==  "Cancer code  self-reported: small intestine/small bowel cancer"]<- "Cancer code self-reported: small intestine/small bowel cancer"
	ukb_dat1$outcome[ukb_dat1$outcome ==  "Cancer code  self-reported: squamous cell carcinoma"]<-"Cancer code self-reported: squamous cell carcinoma"
	meta.ukb<-meta.tab9[meta.tab9$study.abbreviation == "UKB",]
	ukb.m<-merge(ukb_dat1,meta.ukb,by.x="outcome",by.y="original_name",all.x=T)


	col.keep2<-c("proxy.outcome","target_snp.outcome","proxy_snp.outcome","target_a1.outcome","target_a2.outcome","proxy_a1.outcome","proxy_a2.outcome","pmid.y","id.outcome")

	col.keep1<-c("proxy","target_snp","proxy_snp","target_a1","target_a2","proxy_a1","proxy_a2" ,"population.y")

	col.keep<-c("outcome","UKbiobank","rsid","Effect.Allele","Other.Allele","lnor" ,"se","eaf","p","effect_allele_confirmed","cancer","study.abbreviation","site","include","cases","system" ,"cell","Cancer.Group","controls","original_name","study.y","ID","subpopulation","total","Note","overlap","MAC_rs603424","MAC_rs3734398","MAC_rs174546","MAC100rs174546")

	ukb.m<-ukb.m[,names(ukb.m) %in% c(col.keep2,col.keep1,col.keep)]

	names(ukb.m)[names(ukb.m) == "study.y"]<-"study" 
	names(ukb.m)[names(ukb.m) == "pmid.y"]<-"pmid" 
	names(ukb.m)[names(ukb.m) == "population.y"]<-"population" 
	names(ukb.m)[names(ukb.m) == "id.outcome"]<-"id.mrbase"

	ukb.m$set<-"ukb2"
	ukb.m$Effect.Allele<-toupper(ukb.m$Effect.Allele)
	ukb.m$Other.Allele<-toupper(ukb.m$Other.Allele)

	# transform_betas
	# formula: log OR = beta / (u(1-u)); where u=ncases/(ncases + ncontrol) REPEAT with SE 	
	beta<-ukb.m$lnor
	se<-ukb.m$se
	u<-ukb.m$cases/(ukb.m$cases+ukb.m$controls)
	ukb.m$lnor <- beta / (u * (1 - u))
	ukb.m$se<-se / (u * (1 - u)) 	
	return(ukb.m)
}


format_meta_dat3<-function(){ #merge meta data with fingen studies
	load("~/MR_FattyAcids/data/summary_data/mrbase_extra/fin_mrbase.Rdata")
	
	fin_dat1<-format_columns(dat=fin_dat,study="FinnGen") #estimates for finGen from SAIGE, which works for binary traits
	meta.fin<-meta.tab9[meta.tab9$study.abbreviation == "FinnGen",]
	fin.m<-merge(fin_dat1,meta.fin,by.x="outcome",by.y="original_name")

	col.keep2<-c("proxy.outcome","target_snp.outcome","proxy_snp.outcome","target_a1.outcome","target_a2.outcome","proxy_a1.outcome","proxy_a2.outcome","pmid.y","id.outcome")

	col.keep1<-c("proxy","target_snp","proxy_snp","target_a1","target_a2","proxy_a1","proxy_a2" ,"population.y")

	col.keep<-c("outcome","UKbiobank","rsid","Effect.Allele","Other.Allele","lnor" ,"se","eaf","p","effect_allele_confirmed","cancer","study.abbreviation","site","include","cases","system" ,"cell","Cancer.Group","controls","original_name","study.y","ID","subpopulation","total","Note","overlap","MAC_rs603424","MAC_rs3734398","MAC_rs174546","MAC100rs174546")
	
	fin.m<-fin.m[,names(fin.m) %in% c(col.keep2,col.keep1,col.keep)] 
	names(fin.m)[names(fin.m) == "pmid.y"]<-"pmid"
	names(fin.m)[names(fin.m) == "population.y"]<-"population"
	names(fin.m)[names(fin.m) == "study.y"]<-"study"
	names(fin.m)[names(fin.m) == "id.outcome"]<-"id.mrbase"
	fin.m$set<-"fingen"
	fin.m$Effect.Allele<-toupper(fin.m$Effect.Allele)
	fin.m$Other.Allele<-toupper(fin.m$Other.Allele)
	
	return(fin.m)
}


format_meta_dat2<-function(){ #merge meta data with studies obtained from mrabse
	
	mrbase<-format_mrbase_data()
	mrbase$outcome[mrbase$outcome == "ER- Breast cancer"]<-"ER- breast cancer"
	mrbase$outcome[mrbase$outcome == "ER+ Breast cancer"]<-"ER+ breast cancer"

	mrbase.m1<-merge(mrbase,meta.tab9,by.x=c("outcome","pmid"),by.y=c("cancer","pmid"))
	mrbase.m1$cancer<-mrbase.m1$outcome
	mrbase.m2<-merge(mrbase,meta.tab9,by.x=c("outcome","pmid"),by.y=c("original_name","pmid"))
	mrbase.m2$original_name <- mrbase.m2$outcome

	mrbase.m<-rbind(mrbase.m1,mrbase.m2)	
	mrbase.m<-mrbase.m[!duplicated(mrbase.m[,c("outcome","pmid","study.abbreviation","rsid")]),]		
	mrbase2<-mrbase[!mrbase$outcome %in% mrbase.m$outcome,]
	mrbase2<-mrbase2[!mrbase2$outcome %in% c("Neuroblastoma","Gallbladder cancer"),]
	mrbase.m<-rbind.fill(mrbase.m,mrbase2)
	names(mrbase.m)[grep("id",names(mrbase.m))]
	
	col.keep1<-c("proxy","target_snp","proxy_snp","target_a1","target_a2","proxy_a1","proxy_a2" ,"population.y","id_outcome")

	col.keep<-c("pmid","outcome","UKbiobank","rsid","Effect.Allele","Other.Allele","lnor" ,"se","eaf","p","effect_allele_confirmed","cancer","study.abbreviation","site","include","cases","system" ,"cell","Cancer.Group","controls","original_name","study.y","ID","subpopulation","total","Note","overlap","MAC_rs603424","MAC_rs3734398","MAC_rs174546","MAC100rs174546")
	mrbase.m<-mrbase.m[,names(mrbase.m) %in% c(col.keep,col.keep1),] 
	# names(mrbase.m)[!names(mrbase.m) %in% col.keep]

	names(mrbase.m)[names(mrbase.m) == "population.y"]<-"population"
	names(mrbase.m)[names(mrbase.m) == "study.y"]<-"study"
	names(mrbase.m)[names(mrbase.m) == "id_outcome"]<-"id.mrbase"

	mrbase.m$set<-"mrbase"
	mrbase.m$Effect.Allele<-toupper(mrbase.m$Effect.Allele)
	mrbase.m$Other.Allele<-toupper(mrbase.m$Other.Allele)
	return(mrbase.m)
}

format_meta_dat1<-function(){ #merge meta data with studies obtained by correspondence
	load("~/MR_FattyAcids/data/summary_data/outcomes_obtainedby_correspondence_V2.Rdata")
	Can_corr$outcome[Can_corr$outcome == "Esophageal cancer" & Can_corr$study == "Japanese Biobank"]<-"Esophageal squamous cell carcinoma"
	# unique(Can_corr[Can_corr$outcome == "Oesoph_cancer",c("outcome","study")])
	
	meta.tab9<-meta.tab9[order(meta.tab9$cancer),]
	dups<-unique(meta.tab9$pmid[duplicated(meta.tab9$pmid)])
	pmids<-meta.tab9$pmid[!meta.tab9$pmid %in% dups]
	pmids2<-pmids[pmids %in% Can_corr$pmid ]
	meta1 <- meta.tab9[meta.tab9$pmid %in% pmids2,]
	corr1 <- Can_corr[Can_corr$pmid %in% pmids2,]
	corr2 <- Can_corr[!Can_corr$pmid %in% pmids2,]
	corr1.m<-merge(corr1,meta1,by="pmid")	
	meta2<-meta.tab9[meta.tab9$pmid %in% unique(corr2$pmid) ,]
	Dups<-paste(corr2$outcome,corr2$pmid,corr2$rsid)[duplicated(paste(corr2$outcome,corr2$pmid,corr2$rsid))]
	corr3<-corr2[!paste(corr2$outcome,corr2$pmid,corr2$rsid) %in% Dups,]
	corr4<-corr2[paste(corr2$outcome,corr2$pmid,corr2$rsid) %in% Dups,]
	meta2$pmid[meta2$study=="UK Biobank"]<-"unpublished_ukb"
	corr5<-corr3[!paste(corr3$outcome,corr3$pmid) %in% paste(meta2$cancer,meta2$pmid),]
	corr4<-rbind(corr4,corr5)
	corr3.m<-merge(corr3,meta2,by.x=c("outcome","pmid"),by.y=c("cancer","pmid"))	
	corr3.m$cancer<-corr3.m$outcome
	corr2.m<-rbind(corr1.m,corr3.m)	
	corr5<-corr4[corr4$study %in% c("UK Biobank","GECCO"),]
	corr6<-corr4[!corr4$study %in% c("UK Biobank","GECCO"),]

	corr5$ID<-NA
	corr5$ID[corr5$outcome == "Colorectal cancer (distal)"]<-63
	corr5$ID[corr5$outcome == "Colorectal cancer (proximal)"]<-64
	corr5$ID[corr5$outcome == "Brain_cancer"]<-138
	corr5$ID[corr5$outcome == "Breast_cancer"]<-139
	corr5$ID[corr5$outcome == "Colorectal_cancer"]<-143
	corr5$ID[corr5$outcome == "Haem_cancer"]<-137
	corr5$ID[corr5$outcome == "Headneck_cancer"]<-157
	corr5$ID[corr5$outcome == "Leuk_cancer"]<-146
	corr5$ID[corr5$outcome == "Liver_bile_cancer"]<-147
	corr5$ID[corr5$outcome == "Liver_cell_cancer"]<-148
	corr5$ID[corr5$outcome == "Lung_cancer"]<-149
	corr5$ID[corr5$outcome == "Lung_cancer_unadj"]<-149
	corr5$ID[corr5$outcome == "Lymph_leuk_cancer"]<-150
	corr5$ID[corr5$outcome == "Mult_myel_cancer"]<-154
	corr5$ID[corr5$outcome == "Myel_leuk_cancer"]<-155
	corr5$ID[corr5$outcome == "Nm_skin_cancer"]<-156
	corr5$ID[corr5$outcome == "Oesoph_cancer"]<-144
	corr5$ID[corr5$outcome == "Ovarian_cancer"]<-158
	corr5$ID[corr5$outcome == "Pan_cancer"]<-141
	corr5$ID[corr5$outcome == "Pan_inclc44_cancer"]<-140
	corr5$ID[corr5$outcome == "Prostate_cancer"]<-159
	corr5$ID[corr5$outcome == "Skin_cancer"]<-153

	unique(corr5[,c("outcome","study","pmid","ID")])
	corr5.m<-merge(corr5,meta.tab9,by="ID")
	names(corr5.m)[names(corr5.m) =="pmid.x" ]<-"pmid"

	corr5.m$cancer[corr5.m$outcome == "Lung_cancer_unadj"]<-"Lung cancer adjusted for array"
	corr5.m$ID[corr5.m$outcome == "Lung_cancer_unadj"]<-1499

	# meta.tab9[meta.tab9$study == "UK Biobank",c("cancer","ID")]
	corr6.m<-merge(corr6,meta.tab9,by.x=c("outcome","study","pmid"),by.y=c("cancer","study.abbreviation","pmid"))	
	corr6.m$cancer<-corr6.m$outcome
	corr6.m$study.abbreviation<-corr6.m$study
	corr7.m<-rbind.fill(corr6.m,corr5.m)
	corr3.m<-rbind.fill(corr2.m,corr7.m)
	names(corr3.m)[names(corr3.m) == "population.y"]<-"population"

	col.keep<-c("pmid","outcome","UKbiobank","rsid","Effect.Allele","Other.Allele","lnor" ,"se","eaf","p","effect_allele_confirmed","cancer","study.abbreviation","site","include","cases","system" ,"cell","Cancer.Group","controls","original_name","study.y","ID","population","subpopulation","total","Note","overlap","MAC_rs603424","MAC_rs3734398","MAC_rs174546","MAC100rs174546")
	
	corr3.m<-corr3.m[,names(corr3.m) %in% col.keep]	
	names(corr3.m)[names(corr3.m) == "study.y"]<-"study"
	corr3.m$set<-"corr"

	corr3.m$Effect.Allele<-toupper(corr3.m$Effect.Allele)
	corr3.m$Other.Allele<-toupper(corr3.m$Other.Allele)
	return(corr3.m)
}




format_columns<-function(dat=NULL,study=NULL){
	names(dat)[names(dat)=="outcome"]<-"outcome.mrbaseid"
	names(dat)[names(dat)=="originalname.outcome"]<-"outcome"
	names(dat)[names(dat)=="SNP"]<-"rsid"
	names(dat)[names(dat)=="effect_allele.outcome"]<-"Effect.Allele"
	names(dat)[names(dat)=="other_allele.outcome"]<-"Other.Allele"
	names(dat)[names(dat)=="eaf.outcome"]<-"eaf"
	names(dat)[names(dat)=="beta.outcome"]<-"lnor"
	names(dat)[names(dat)=="se.outcome"]<-"se"
	names(dat)[names(dat)=="pval.outcome"]<-"p"
	dat$effect_allele_confirmed<-TRUE	
	dat$population<-"European"
	dat$pmid<-"unpublished"
	dat$study<-study
	dat$UKbiobank<-TRUE
	return(dat)
}

# combine and format mrbase cancer datasets and datasets obtained by correspondence 
format_all<-function(){
	Cancer_dat<-rbind.fill(Dat,Can_corr)
}


# Cancer_mrbase<-format_mrbase_data()

format_mrbase_data<-function(){
	load("~/MR_FattyAcids/data/summary_data/cancers_mrbase.Rdata")
	Dat1<-Cancer_mrbase[Cancer_mrbase$outcome == "Prostate cancer (overall) || id:1174",]
	Dups<-unique(Dat1$rsid[duplicated(Dat1$rsid)])
	Dat3<-Dat1[!Dat1$rsid %in% Dups,]
	Dat2<-Cancer_mrbase[Cancer_mrbase$outcome != "Prostate cancer (overall) || id:1174",]
	Dat<-rbind(Dat3,Dat2)

	Pos<-is.na(Dat$study)
	Tempid<-paste(Dat$pmid[Pos],Dat$outcome[Pos])
	Dat$study[Pos]<-paste("unknown",as.numeric(as.factor(Tempid)),sep="")
	Dat$id_outstudy<-paste(Dat$outcome,Dat$study)
	Dat$id_outstudysnp<-paste(Dat$outcome,Dat$study,Dat$rsid)
	Dat<-Dat[order(Dat$ncase,decreasing=T),]
	Dups<-unique(Dat$id_outstudysnp[duplicated(Dat$id_outstudysnp)])

	Temp<-Dat[which(Dat$id_outstudysnp %in% Dups),c("outcome","study","ncase","pmid","outcome.deprecated","rsid","lnor","Effect.Allele","Other.Allele","eaf")]

	if(nrow(Temp)>0) stop("duplicate SNPS within study present") 

		
	names(Dat)[names(Dat) == "outcome"]<-c("id_outcome","outcome")
	
	if(is.null(nrow(grep("ieu",Dat$id)))){
		Dat$id<-paste0("ieu-a-",Dat$id)
	}
	# tricl.id<-paste0("ieu-a-",c(984,985,986,987,988,989))
	# data.frame(ao[ao$id %in% tricl.id,])

	Dat$outcome[which(Dat$id == "ieu-a-985")]<-"Lung cancer in ever smokers"
	Dat$outcome[which(Dat$id == "ieu-a-986")]<-"Lung cancer in never smokers"

	# ao1<-ao[ao$id %in% Dat$id,]	
	# data.frame(ao[ao$id %in% tricl.id,])
	# unique(ao$id[grep("984",ao$id)])
	# data.frame(ao[which(ao$trait == "Lung cancer"),])
	# unique(Dat$id)
	# Dat2<-merge(Dat,ao1[,c("id","note")])
	
	unique(Dat3$outcome[grep("survival",Dat3$outcome,ignore.case=T)])
	Dat2<-Dat
	Dat2$id_outidsnp<-paste(Dat2$id_outcome,Dat2$rsid)
	Dups<-Dat2$id_outidsnp[duplicated(Dat2$id_outidsnp)]
	Dat3<-Dat2[!duplicated(Dat2$id_outidsnp),] #2 duplicate SNPs in PRACTICAL. same rsid but different strands
	# Dat3<-Dat3[!Dat3$outcome %in% c("Prostate cancer Gleason score","Prostate cancer (high aggressive vs low/intermediate aggressive)","Prostate cancer (high aggressive vs low aggressive)" ,"Prostate cancer (advanced vs non-advanced)" ),]

	Dat3$id_outsnp<-paste(Dat3$outcome,Dat3$rsid)

	Dat3$pmid[Dat3$study == "TRICL"]<-28604730
	Dat3<-Dat3[which(Dat3$pmid!=25751625),] #This study is the older version of the BCAC study. 29059683 is more recent with 
	Dat3<-Dat3[!Dat3$id %in% c("UKB-a:205", "UKB-a:204" , "UKB-a:213"),] #these Ben Neale ids are same outcomes as UKB-b:14521, UKB-b:7773 and UKB-b:13584, which have more cases.
	Dat3<-Dat3[order(Dat3$ncase,decreasing=T),]
	Dat3<-Dat3[which(Dat3$study != "ILCCO"),]

	Dat3$study[grep("UKB",Dat3$id)]<-"UK Biobank"
	Dat3$id_outstudysnp<-paste(Dat3$outcome,Dat3$study,Dat3$rsid)
		
	# Temp2<-Dat3[Dat3$id_outstudysnp %in% Dups,c("outcome","id","study","ncase","rsid","id_outsnp","Effect.Allele","Other.Allele","eaf","lnor")]
	Dat3<-Dat3[order(Dat3$ncase,decreasing=T),]
	Dat4<-Dat3[!duplicated(Dat3$id_outstudysnp),]
	Dups<-unique(Dat4$id_outstudysnp[which(duplicated(Dat4$id_outstudysnp))])
	Temp<-Dat4[Dat4$id_outstudysnp %in% Dups,c("outcome","id","study","ncase")]
	if(nrow(Temp) > 0){
	    stop("duplicates SNPs present")
	}
	
	
	# Dat4<-Dat3[which(!duplicated(Dat3$id_outsnp)),] # This step not needed because there are no duplicates 
	Dat<-Dat4
	Dat$outcome[Dat$outcome == "Prostate cancer (overall)"]<-"Prostate cancer"
	Dat$outcome[Dat$outcome == "Illnesses of father: Lung cancer"]<-"Lung cancer in father"   
	Dat$outcome[Dat$outcome == "Illnesses of mother: Breast cancer"  ]<-"Breast cancer in mother"
	Dat$outcome[Dat$outcome == "Illnesses of father: Prostate cancer"    ]<-"Prostate cancer in father"
	Dat$outcome[Dat$outcome == "Illnesses of mother: Lung cancer" ]<-"Lung cancer in mother"
	Dat$outcome[Dat$outcome == "Illnesses of siblings: Breast cancer"  ]<-"Breast cancer in siblings"
	Dat$outcome[Dat$outcome == "Illnesses of siblings: Lung cancer"      ]<-"Lung cancer in siblings"
	Dat$outcome[Dat$outcome =="Illnesses of siblings: Prostate cancer"  ]<-"Prostate cancer in siblings"
	# unique(Dat$outcome)
	                                               
	Dat$outcome[Dat$outcome == "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"Breast cancer"
	Dat$outcome[Dat$outcome == "ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"ER+ Breast cancer" 
	Dat$outcome[Dat$outcome == "Breast cancer (Oncoarray)"]<-"Breast cancer"
	Dat$outcome[Dat$outcome =="Breast cancer (iCOGS)"]<-"Breast cancer"
	Dat$outcome[Dat$outcome == "ER+ Breast cancer (Oncoarray)"]<-"ER+ Breast cancer"
	Dat$outcome[Dat$outcome == "ER+ Breast cancer (iCOGS)"]<-"ER+ Breast cancer"
	Dat$outcome[Dat$outcome == "ER+ Breast cancer (GWAS)"]<-"ER+ Breast cancer"
	Dat$outcome[Dat$outcome == "ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"ER- Breast cancer" 
	Dat$outcome[Dat$outcome == "Breast cancer (GWAS)"]<-"Breast cancer" 
	Dat$outcome[Dat$outcome == "Breast cancer (GWAS) "]<-"Breast cancer" 
	Dat$outcome[Dat$outcome == "ER- Breast cancer (Oncoarray)"]<-"ER- Breast cancer" 
	Dat$outcome[Dat$outcome == "ER- Breast cancer (iCOGS)"]<-"ER- Breast cancer" 
	Dat$outcome[Dat$outcome == "ER- Breast cancer (GWAS)"]<-"ER- Breast cancer" 
	Dat$outcome[Dat$outcome == "Breast cancer (Survival)"]<-"Breast cancer survival" 
	Dat$outcome[Dat$outcome == "ER+ Breast cancer (Survival)"]<-"ER+ Breast cancer survival" 
	Dat$outcome[Dat$outcome == "ER- Breast cancer (Survival)"]<-"ER- Breast cancer survival" 

	Dat$population<-"European"
	# unique(Dat$pmid[Dat$id %in% ID.other])
	Dat$population[which(Dat$pmid == 23504502)]<-"East Asian"
	Dat$population[which(Dat$pmid == 22318345)]<-"East Asian"
	Dat$UKbiobank<-FALSE #confirmed no overlap with PRACTICAL, TRICL, BCAC, OCAC. Of the large 
	# names(Dat)[which(names(Dat)=="outcome")]<-c("id_outcome","outcome")
	# names(Dat)[names(Dat)=="id"]<-"id_mrbase"
	Dat<-Dat[,!names(Dat) %in% c("category","subcategory","units","year","outcome.deprecated", "mr_keep","data_source")]

	Dat$Effect.Allele<-toupper(Dat$Effect.Allele)
	Dat$Other.Allele<-toupper(Dat$Other.Allele)

	unique(Dat$pmid[which(Dat$Effect.Allele == "")]) #PMID 21372204 prostate cancer study of  1175     cases  and 1100 controls missing effect allele informatino. 23504502 UGI cancers dbGAP. outcome name wrong anyway. Drop
	Dat<-Dat[Dat$Effect.Allele != "",] #PMID 21372204 

	Dat$id_snp_pmid_outcome<-paste(Dat$rsid,Dat$pmid,Dat$outcome)
	mrc_ids<-unique(Dat$id[grep("UKB-b",Dat$id)]) #ids from MRC IEU automated pipeline. cancer in relative
	ben_ids<-unique(Dat$id[grep("UKB-a",Dat$id)]) #ids from ben neale automated pipeline. cancer in relative
	# unique(Dat[Dat$study == "OCAC",c("outcome","ncase","ncontrol") ])
	# unique(Dat[Dat$study == "OCAC",c("outcome") ])
	# unique(Dat[Dat$outcome %in% c("Prostate cancer (early-onset)","Prostate cancer (advanced)"),c("outcome","ncase","ncontrol") ])

	#drop gliomoscan. the version in mrbase is wrong
	Dat<-Dat[Dat$study != "GliomaScan",]
	
	# ID.other<-c(ao$id[ao$trait=="Glioma" & ao$filename=="for_Philip_data_delivery.txt.new.tab"],ao$id[ao$trait=="Neuroblastoma"],ao$id[ao$trait=="Upper gastrointestinal cancers"],ao$id[ao$trait=="Thyroid cancer"],ao$id[ao$trait=="Gallbladder cancer"])

	#drop MRC IEU cancers. 

	Dat<-Dat[!Dat$id %in% c(mrc_ids,ben_ids),]
	
	# drop duplicates
	Dat<-Dat[order(Dat$ncase,decreasing=T),]
	# Dups<-unique(Dat$id_snp_pmid_outcome[duplicated(Dat$id_snp_pmid_outcome)])
	# Temp<-Dat[Dat$id_snp_pmid_outcome %in% Dups,]
	# head(Temp[order(Temp$id_snp_pmid_outcome),c("rsid","pmid","outcome","original_outcome")])
	Dat<-Dat[!duplicated(Dat$id_snp_pmid_outcome),] #duplicates all seem to correspond to BCAC breast cancer from GWAS, icogs and oncoarray
	Dat<-Dat[which(Dat$id_outcome != "Prostate cancer (overall) || id:1174"),] #drop overall prostate cancer which is missing SNPs that shouldn't be missing. Something went wrong in the original upload to MR-Base. I extracted the SNPs directly from the raw files and added to the Corr set
	return(Dat)	
}



format_proxies<-function(Dat=NULL,snp=NULL){
    load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata") #this file contains the proxy information
    Pops<-c("EUR","EAS")
    Populations<-c("European","East Asian")
    ids<-unique(Dat$ID)
    ids<-ids[!is.na(ids)]
    # ids2<-Dat$ID[which(Dat$rsid == snp)]
    id.miss<-unique(Dat$ID)
    print(paste("SNP=",snp))
# which(id.miss %in% 99:103)
    L.proxy<-NULL
    for(i in 1:length(id.miss)){
        print(i)
        print(paste("ID=",unique(id.miss[i])))
        print(unique(Dat$cancer[which(Dat$ID == id.miss[i])]))


        miss.dat<-Dat[which(Dat$ID %in% id.miss[i]),]
        Pop<-Pops[which( Populations == unique(miss.dat$population))]

        if(sum(which(fa.reg2$SNP == snp))==0){ 
            stop("snp not present in fa.reg2") 
        }
        alias_snps<-unique(fa.reg2$alias_rs[which(fa.reg2$target.alias == snp )])
        if(any(miss.dat$rsid %in% alias_snps)){
            stop("alias rsid found")
            miss.dat2<-miss.dat[miss.dat$rsid %in% alias_snps,]
        }
        
        snp.proxies1<-unique(unlist(strsplit(fa.reg2$SNP[which(fa.reg2$ProxySNP == snp & fa.reg2$pop == Pop )],split=" ")))
        snp.proxies2<-unique(unlist(strsplit(fa.reg2$ProxySNP[which(fa.reg2$SNP == snp & fa.reg2$pop == Pop )],split=" ")))


        if(all(!miss.dat$rsid %in% snp.proxies2)){
            if(any(miss.dat$rsid %in% snp.proxies1)){
                message("proxies present in snp.proxies1 but not snp.proxies2")
            }

            if(all(!miss.dat$rsid %in% snp.proxies1)){
                message("no proxies found in dataset")                              
            }
        }
        if(any(miss.dat$rsid %in% snp.proxies2)){
            miss.dat2<-miss.dat[miss.dat$rsid %in% snp.proxies2,]
            Datp<-unique(fa.reg2[which(fa.reg2$ProxySNP %in% miss.dat2$rsid & fa.reg2$SNP == snp & fa.reg2$pop == Pop ),c("SNP","ProxySNP","source","Coord.proxy","Alleles","MAF","Distance","Dprime","R2","Correlated_Alleles","Ndistance0","pop","SNP.GRCh37","SNP.GRCh38.p12","Chr","Proxy.GRCh38.p12","Proxy.GRCh37","chr_other.ens", "SNP.GRCh38.p12.other.ens","chr_other.proxy.ens", "Proxy.GRCh38.p12.other.ens")])
            Datp2<-unique(Datp[,!names(Datp) %in% c("chr_other.ens", "SNP.GRCh38.p12.other.ens","chr_other.proxy.ens", "Proxy.GRCh38.p12.other.ens")])
        

            Dat.m<-merge(miss.dat2,Datp2,by.x="rsid",by.y="ProxySNP") 
            if(nrow(Dat.m)>1){ 
                Dat.m[order(Dat.m$R2,decreasing=T),]
                Dat.m<-Dat.m[!duplicated(Dat.m$SNP),]
            }#only want one proxy SNP . choose proxy SNP with strongest r2 with target SNP

            names(Dat.m)[names(Dat.m) == "rsid"]<-"ProxySNP"
            names(Dat.m)[names(Dat.m) == "SNP"]<-"rsid"
            Correlated_Alleles<-Dat.m$Correlated_Alleles
            alleles<-unlist(strsplit(Correlated_Alleles,split=","))
            allele1<-alleles[1]
            allele1_1<-unlist(strsplit(allele1,split="="))[1]
            allele1_2<-unlist(strsplit(allele1,split="="))[2]
            allele2<-alleles[2]
            allele2_1<-unlist(strsplit(allele2,split="="))[1]
            allele2_2<-unlist(strsplit(allele2,split="="))[2]
            alleleA<-c(allele1_1,allele2_1)
            alleleB<-c(allele1_2,allele2_2)
            ea<-Dat.m$Effect.Allele
            oa<-Dat.m$Other.Allele

            if(is.na(ea)){
                ea<-alleleB[alleleB != oa]

            }
            if(is.na(oa)){
                oa<-alleleB[alleleB != ea]
            }
            print(Correlated_Alleles)
    
            # Assuming alleleB are the proxy alleles. Allele B always seems to be proxy Allele, i.e. in Correlated_Alleles the second allele seems to always be the proxy allele
            # sometimes alleles are same for proxy and target SNP, hence can't use above script for alleleA
            if(any(alleleB == ea) & any(alleleB == oa)){
            # if((ea %in% alleleB &  oa %in% alleleB)) {#ea and oa are in alleleB 
                # print(ea %in% alleleB &  oa %in% alleleB)
                message("ea or oa in AlleleB set")
                Dat.m$ea_new<-alleleA[match(ea,alleleB)]
                Dat.m$oa_new<-alleleA[match(oa,alleleB)]
                # L.proxy[[i]]<-Dat.m
            }   
            # Dat.m[,c("rsid","ProxySNP","Effect.Allele","Other.Allele","ea_new","oa_new","Correlated_Alleles","id")]
            Dat.m2<-Dat.m[,!names(Dat.m) %in% c("Effect.Allele","Other.Allele")]
            # Dat.m[,c("ea_new","oa_new","Effect.Allele","Other.Allele")]
            names(Dat.m2)[names(Dat.m2) %in% c("ea_new","oa_new")]<-c("Effect.Allele","Other.Allele")
            Dat.m2$proxy_outcome_snp<-TRUE          
        }else{next}         
        L.proxy[[i]]<-Dat.m2                    
        # rm(Dat.m2)
    }
    if(!is.null(L.proxy)){
        L.dat<-do.call(rbind,L.proxy)
        return(L.dat)
    }
}



get_missing_allele<-function(Dat=NULL,snp=NULL){
    load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata") #this file contains the proxy information
    Traits<-unique(Dat$cancer)
    L<-NULL
    for(i in 1:length(Traits)){
        print(i)
        print(Traits[i])
        # Dat2[Dat2$outcome== Traits[i] & Dat2$rsid == snp,]
        Temp<-Dat[Dat$cancer== Traits[i] & Dat$rsid == snp,]
        ref.snps<-unique(fa.reg2[fa.reg2$SNP == snp,c("SNP","effect_allele","other_allele")])
        ref.snps<-ref.snps[!is.na(ref.snps$other_allele),]
        ref.snps$effect_allele<-toupper(ref.snps$effect_allele)
        ref.snps$other_allele<-toupper(ref.snps$other_allele)
        # Temp2<-Temp[Temp$rsid == target.snps[1],]
        alleles<-paste(ref.snps[,c("effect_allele","other_allele")])
        Temp$Other.Allele<-alleles[alleles != Temp$Effect.Allele ]
        L[[i]]<-Temp
    }
    L.dat<-do.call(rbind,L)
    return(L.dat)
}





harmonise_ea<-function(dat=NULL,ea=NULL,oa=NULL){
	dat$eaf<-round(dat$eaf,2)
	strand1<-c("C","G","A","T")
	strand2<-c("G","C","T","A")
	ea_altstrand<-strand2[which(strand1 %in% ea)]
	dat1<-dat[dat$Effect.Allele == ea,]
	dat2<-dat[dat$Effect.Allele != ea,]
	dat3<-dat2[dat2$Effect.Allele == ea_altstrand,]
	dat3$Effect.Allele <- ea
	dat3$Other.Allele <- oa
	
	dat4<-dat2[dat2$Effect.Allele != ea_altstrand,]
	dat4$Effect.Allele <- ea 
	dat4$Other.Allele <- oa
	dat4$lnor<-dat4$lnor*-1
	dat4$eaf<-1-dat4$eaf

	dat5<-do.call(rbind,list(dat1,dat3,dat4))
	return(dat5)
}



make_refdat2<-function(){
	ref.dat<-data.frame(matrix(c(c("rs174546","European","C","T", "0.3592"), #consistent amongst different european subpops in ensmbl (comes close to 0.5 in fins)
		c("rs174546","East Asian","C","T","0.354"), #curiously the minor alleles are different amongst different east asian populatiouns, including between regions within China) but is consistent between han Chiniese and Japanese
        c("rs174546","North East Asian","C","T","0.354"),
        c("rs174546","South East Asian","C","T","0.586"),
		c("rs2236212","European","G" ,"C","0.4168"), #consistent amongst eur subpops in ensembl but comes to close 0.5 in some
		c("rs2236212","East Asian","G","C","0.709"), #consistent amongst different east asian sub pops in ensembl
		c("rs603424","European","G","A","0.1632"), #consistent amongst sub eur pops
		c("rs603424","East Asian","G","A","0.088"), #consistent amongst esat asians sub pops
		c("rs174576","European","C","A","0.3638"),
		c("rs174576","East Asian","C","A","0.354")), #freq varies between east asi sub pops (incl between chinese reginos). but consistent between japanese and han chinese 
		nrow=10,ncol=5,byrow=TRUE),stringsAsFactors=F)


	ref.dat$X5<-as.numeric(ref.dat$X5)
	Pos<-ref.dat$X5<0.5
	ref.dat$minor.allele<-ref.dat$X3
	ref.dat$minor.allele[Pos]<-ref.dat$X4[Pos]	
	Pos<-ref.dat$X5>0.5
	ref.dat$major.allele<-ref.dat$X3
	ref.dat$major.allele[Pos]<-ref.dat$X4[Pos]
	maf<-ref.dat$X5
	Pos<-which(maf>0.5)
	maf[Pos]<-1-maf[Pos]
	ref.dat$maf<-maf
	strand1<-c("C","G","A","T")
	strand2<-c("G","C","T","A")
	ref.dat$minor.allele2<-strand2[match(ref.dat$minor.allele,strand1)] 
	ref.dat$major.allele2<-strand2[match(ref.dat$major.allele,strand1)] 	
	names(ref.dat)<-c("rsid","population","allele1","allele2","a2freq","minor.allele","major.allele","maf","minor.allele2","major.allele2")
	return(ref.dat)
}



 

find_allele_error2<-function(Dat1=NULL){	
	Dat1<-make_ea_ma(dat=Dat1,eaf=target_dat_eaf_col)
	ref.dat<-make_refdat()
	Dat1$subpopulation2<-Dat1$subpopulation
	Dat1$subpopulation2[which(Dat1$subpopulation=="Central and South (68%) East Asian") ]<-"South East Asian"
	Dat1$subpopulation2[which(Dat1$subpopulation=="South & North (69%) East Asian" ) ]<-"North East Asian"
	Dat1$subpopulation2[which(Dat1$subpopulation=="South (23%), Central & North East Asian" ) ]<-"North East Asian"
	Dat.m<-merge(Dat1,ref.dat,by.x=c("rsid","subpopulation2"),by.y=c("rsid","population"),all.x=T)
	Temp<-Dat.m[which(Dat.m$Effect.Allele != Dat.m$minor.allele),c("Effect.Allele","Other.Allele","eaf","minor.allele","minor.allele2","subpopulation")]
	Temp<-Temp[which(Temp$Effect.Allele != Temp$minor.allele2),]
	nrow(Temp)
	return(Temp)	
}
 

make_ea_ma<-function(dat=NULL,eaf=NULL,effect=NULL,effect_allele=NULL,other_allele=NULL){
	# dat1<-dat[dat$consortium == study,]
	dat[,eaf]<-as.numeric(dat[,eaf])
	dat1<-dat[!is.na(dat[,eaf]),]
	dat2<-dat[is.na(dat[,eaf]),]
	maf<-dat1[,eaf]
	Pos.change<-which(maf>0.5)
	Beta<-dat1[,effect][Pos.change]*-1
	dat1[,effect][Pos.change]<-Beta
	# Beta.sd<-dat1$beta.sd[Pos.change]*-1
	# dat1$beta.sd[Pos.change]<-Beta.sd
	Eaf<-1-dat1[,eaf][Pos.change]
	dat1[,eaf][Pos.change]<-Eaf
	EA<-dat1[,effect_allele][Pos.change]
	OA<-dat1[,other_allele][Pos.change]
	dat1[,effect_allele][Pos.change]<-OA
	dat1[,other_allele][Pos.change]<-EA
	dat1$target_dat_maf<-dat1[,eaf]
	return(dat1)
	# dat2<-dat[dat$consortium != study,]
	# unique(dat2$consortium)
	# dat3 <- rbind(dat1,dat2) 
}


