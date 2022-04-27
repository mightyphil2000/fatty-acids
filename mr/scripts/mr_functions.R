

harmonise_dat<-function(Dat=NULL,ref_dat=NULL,marker="SNP",effect_allele="effect_allele",other_allele="other_allele",effect="beta",EAF="eaf",assume_same_strand=TRUE,assume_same_strand_palindromic_and_non_palindromic_SNPs=FALSE,drop_palindromic_SNPs=FALSE){

	if(assume_same_strand_palindromic_and_non_palindromic_SNPs) assume_same_strand<-FALSE
	ref_dat[,effect_allele]<-toupper(ref_dat[,effect_allele])
	ref_dat[,other_allele]<-toupper(ref_dat[,other_allele])
	names(ref_dat)[names(ref_dat) == effect_allele]<-"ref_ea"
	names(ref_dat)[names(ref_dat) == other_allele]<-"ref_oa"
	names(ref_dat)[names(ref_dat) == EAF]<-"ref_eaf"
	Dat[,effect_allele]<-toupper(Dat[,effect_allele])
	Dat[,other_allele]<-toupper(Dat[,other_allele])	
	ref_dat<-ref_dat[,c(marker,"ref_ea","ref_oa","ref_eaf")]
	Dat<-merge(Dat,ref_dat,by=marker)	

	Alleles<-paste(Dat[,effect_allele],Dat[,other_allele],sep="")	
	Pos<-Alleles %in% c("GC","CG","TA","AT")
	Dat2<-Dat[Pos,] #palindromic SNPs
	Dat3<-Dat[!Pos,] #non palindromic SNPs
	
	# deal with palindromic SNPs
	# check that palindromic SNPs "appear to be on same strand". 
	ea<-Dat2[,effect_allele]
	# if(any(ea != Dat2$ref_ea & ea != Dat2$ref_oa)) 
	Pos<-ea != Dat2$ref_ea & ea != Dat2$ref_oa
	Dat2_1<-Dat2[!Pos,]# drop palindromic SNPs that are not palindromic in both the test and reference datasets. E.g. one rs11231053 is C/G in CHARGE GLA:LA but is A/G UK biobank bim file and dbSNP.  
	# ukb[which(ukb$V2 == "rs11231053"),]

	Dat2_2<-Dat2_1[is.na(Dat2_1[,EAF]),] #palindromic SNPS with missing eaf (gtex)
	Dat2_1<-Dat2_1[!is.na(Dat2_1[,EAF]),] #palindromic SNPs not missing eaf
	
	# exclude SNPs with MAF >=0.43
	Pos<-Dat2_1[,EAF]<0.43 | Dat2_1[,EAF]>0.57
	Dat2_1<-Dat2_1[Pos,]
	Dat2_1$direction<-NA
	Pos<-Dat2_1[,EAF]<0.5
	Dat2_1$direction[Pos]<--1
	Dat2_1$direction[!Pos]<-1
	Pos<-Dat2_1$ref_eaf<0.5
	Dat2_1$ref_direction<-NA
	Dat2_1$ref_direction[Pos]<--1
	Dat2_1$ref_direction[!Pos]<-1
	# if(any(is.na(Dat2_1$direction))) stop("palindromic SNPs with MAF = 0.5 present")
	ea<-Dat2_1[,effect_allele]
	oa<-Dat2_1[,other_allele]
	eaf<-Dat2_1[,EAF]
	beta<-Dat2_1[,effect]
	Pos<-Dat2_1$direction != Dat2_1$ref_direction
	# Dat2_1[Dat2_1$marker=="rs12787928",]
	Dat2_1[,EAF][Pos]<-1-eaf[Pos]
	Dat2_1[,effect][Pos]<-beta[Pos]*-1
	Dat2_1[,effect_allele][Pos]<-oa[Pos]
	Dat2_1[,other_allele][Pos]<-ea[Pos]
	# Dat2_1[Dat2_1[,effect_allele] != Dat2_1$ref_ea,c("effect_allele","other_allele","ref_ea","ref_oa","eaf","ref_eaf")]
	if(!assume_same_strand){
		ea<-Dat3[,effect_allele]
		oa<-Dat3[,other_allele]
		eaf<-Dat3[,EAF]
		beta<-Dat3[,effect]
		Pos1<-ea != Dat3$ref_ea & ea != Dat3$ref_oa 
		Pos2<-oa != Dat3$ref_ea & oa != Dat3$ref_oa 	
		if(any(Pos1 & Pos2)) stop("non palindromic SNPs on different strands")
		Pos<-Pos1 | Pos2
		Dat3_1<-Dat3[!Pos,]
		ea<-Dat3_1[,effect_allele]
		oa<-Dat3_1[,other_allele]
		eaf<-Dat3_1[,EAF]
		beta<-Dat3_1[,effect]
		Pos<-which(ea != Dat3_1$ref_ea) #positions where effect allele is different from effect allele in reference set
		Dat3_1[,effect_allele][Pos]<-oa[Pos]
		Dat3_1[,other_allele][Pos]<-ea[Pos]		
		Dat3_1[,EAF][Pos]<-1-eaf[Pos]
		Dat3_1[,effect][Pos]<-beta[Pos]*-1

		if(any(Dat3_1[,effect_allele] != Dat3_1$ref_ea)) stop("effect alleles not harmonised")
		if(any(Dat3_1[,other_allele] != Dat3_1$ref_oa)) stop("effect alleles not harmonised")

		# strand1<-c("A","T","G","C")
		# strand2<-c("T","A","C","G")
		# ea<-Dat3_1[,effect_allele]
		# oa<-Dat3_1[,other_allele]
		# Pos<-match(ea,strand1)
		# Dat3_1[,effect_allele]<-strand2[Pos]
		# Pos<-match(oa,strand1)
		# Dat3_1[,other_allele]<-strand2[Pos]
		# ea<-Dat3_1[,effect_allele]
		# oa<-Dat3_1[,other_allele]	

		Dat<-rbind.fill(Dat3_1,Dat2_1)
	}
	
	if(assume_same_strand){
		if(!all(Dat2_1[,effect_allele] == Dat2_1$ref_ea)) stop("palindromic SNPs are on different strands. You should harmonise using maf, or drop palindromic SNPs. Check that non-palindomic SNPs are also on different strands") #this steps using Dat2_1 which are palindromic SNPs harmonised just on allele frequency. May not be a reliable way to harmonise datasets from North and South East Asia, where allele frequency is sometimes in conflict directions. Perhaps suggests should not be running MR of exposures and outcomes derived from different East Asian regions? eQTLs and cancer outcomes derived from Japan (north). Fatty acid data derived from Singapore (south)
		if(all(Dat2_1[,effect_allele] == Dat2_1$ref_ea)) message("palindromic SNPs are on the same strand")
		# harmonise non palindromic and palindromic SNPs assuming same strand
			ea<-Dat[,effect_allele]
			oa<-Dat[,other_allele]
			eaf<-Dat[,EAF]
			beta<-Dat[,effect]
			if(any(ea != Dat$ref_ea & ea != Dat$ref_oa)) stop("some SNPs are on different strands to the reference")
			if(all(ea == Dat$ref_ea | ea == Dat$ref_oa)) message("all non palindromic SNPs are on the same strand")

			Pos<-which(ea != Dat$ref_ea) #positions where effect allele is different from effect allele in reference set
			Dat[,effect_allele][Pos]<-oa[Pos]
			Dat[,other_allele][Pos]<-ea[Pos]
			Dat[,EAF][Pos]<-1-eaf[Pos]
			Dat[,effect][Pos]<-beta[Pos]*-1
			all(Dat[,effect_allele] == Dat$ref_ea)	
	}
	
	# This is necessary to harmonise the East Asian datasets on fatty acids, eQTLs and cancer. The eQTL data in BBJ does not have allele frequency. Allele frequency in FA and Cancer datasets in East Asians sometimes in conflicting directinos but it is known that allele frequency is in opposing directions for main FADS variant between South East and North East Asians. Japan is in north east asia. Singapore is in South East Asian
	# appears that palindromic SNPS are already removed from BBJ eQTL datasets (but are present in other,non-eQTL BBJ GWAS datasets). BBJ also known as BJ
	if(assume_same_strand_palindromic_and_non_palindromic_SNPs){
		ea<-Dat[,effect_allele]
		oa<-Dat[,other_allele]
		eaf<-Dat[,EAF]
		beta<-Dat[,effect]
		Pos1<-ea != Dat$ref_ea & ea != Dat$ref_oa 
		Pos2<-oa != Dat$ref_ea & oa != Dat$ref_oa 	
		
		if(any(Pos1) | any(Pos2)) stop("non palindromic SNPs may be on different strands")

		if(all(ea == Dat$ref_ea | ea == Dat$ref_oa)) message("non palindromic SNPs are on the same strand as ref dat")

		Pos<-which(ea != Dat$ref_ea) #positions where effect allele is different from effect allele in reference set
		Dat[,effect_allele][Pos]<-oa[Pos]
		Dat[,other_allele][Pos]<-ea[Pos]
		Dat[,EAF][Pos]<-1-eaf[Pos]
		Dat[,effect][Pos]<-beta[Pos]*-1
		all(Dat[,effect_allele] == Dat$ref_ea)
	}

	# Dat2_1[,c("eaf","ref_eaf")]
	# Pos<-Dat2_1[,EAF]<0.5
	# Dat2_1$direction[Pos]<--1
	# Dat2_1$direction[!Pos]<-1
	# Pos<-Dat2_1$ref_eaf<0.5
	# Dat2_1$ref_direction<-NA
	# Dat2_1$ref_direction[Pos]<--1
	# Dat2_1$ref_direction[!Pos]<-1
	# Pos<-Dat2_1$direction != Dat2_1$ref_direction
	# all(!Pos)

	# join non palindromic and palindromic SNPs together
	# Dat<-do.call(rbind.fill,list(Dat2_1,Dat2_2,Dat3))	 	
	if(drop_palindromic_SNPs)
	{
		Alleles<-paste(Dat[,effect_allele],Dat[,other_allele],sep="")	
		Dat<-Dat[which(!Alleles %in% c("AT","TA","GC","CG")),]
	}
	return(Dat)
}

mr_meta_analysis<-function(Dat=NULL,beta.col="b",se.col="se",outcome="cancer",ncase="cases",ncontrol="controls"){
	id.ma<-disc.tab9$ID[grep(";",disc.tab9$ID)] #identify IDs for outcomes that should be meta analysed
	IDS<-unique(id.ma)

	# i<-which(IDS == "86; 5; 49")
	# i<-which(IDS=="118; 120; 158; 18; 51")
	L<-NULL
	for(i in 1:length(IDS)){
		ID<-trimws(unlist(strsplit(IDS[i],split=";")))
		Dat2<-Dat[Dat$id.outcome %in% ID,]
		# Dat2[order(Dat2$exposure),c("exposure","id.exposure","id.outcome","population")]
		# Meta_analysis<-meta_analysis(dat=Dat2,beta.col=beta.col,se.col=se.col,outcome=outcome,ncase=ncase,ncontrol=ncontrol)
		Meta_analysis<-meta_analysis_v2(dat=Dat2,beta.col=beta.col,se.col=se.col,outcome=outcome,ncase=ncase,ncontrol=ncontrol) #updated so that only GLA:LA is used in East Asian studies and only AA:DGLA used for European populations. This reflects the definition of PUFA desaturase activity as AA:DGLA in European studies and as GLA:LA in East Asian studies 
		# Dat3<-rbind.fill(Dat2,Meta_analysis)
		L[[i]]<-Meta_analysis
	}
	Dat4<-do.call(rbind,L)
	return(Dat4)
}


format_exposure2<-function(dat=NULL,standardise_beta=FALSE,beta="beta",se="se",pval="pval",samplesize="samplesize",effect_allele="effect_allele",other_allele="other_allele",eaf="eaf",rsid="rsid",ID="population",snps=NULL,exposure="exposure"){
	names(dat)[names(dat) ==beta]<-"beta.exposure"
	names(dat)[names(dat) ==se]<-"se.exposure"
	names(dat)[names(dat) ==pval]<-"pval.exposure"
	names(dat)[names(dat) == samplesize]<-"samplesize.exposure"
	names(dat)[names(dat) ==effect_allele]<-"effect_allele.exposure"
	names(dat)[names(dat) ==other_allele]<-"other_allele.exposure"
	names(dat)[names(dat) ==eaf]<-"eaf.exposure"
	names(dat)[names(dat) ==rsid]<-"SNP"		
	names(dat)[names(dat) == exposure]<-"exposure"		
	if(any(dat$effect_allele.exposure == TRUE)) {
		warning("effect alleles set to TRUE. Changing to T")
		dat$effect_allele.exposure[dat$effect_allele.exposure == TRUE]<-"T"
	}
	if(!is.null(snps)){
		dat<-dat[dat$SNP %in% snps,]
	}

	
	if(!is.null(ID)){
		dat$id.exposure<-paste0("exp",as.numeric(as.factor(paste(dat$exposure,dat$SNP,dat[,ID],sep="_"))))
		# dat$id.exposure<-gsub(" ","_",dat$id.exposure)
	}else{
		dat$id.exposure<-paste(dat$exposure,dat$SNP,sep="_")
	}

	if(standardise_beta){
		b<-dat$beta.exposure
		se<-dat$se.exposure
		z<-b/se
		n<-dat$samplesize.exposure
		maf<-dat$eaf.exposure
		maf[maf>0.5]<-1-maf[maf>0.5]
		b_sd2<-b_sd(z=z,maf=maf,n=n)
		b
		b_sd2
		dat$se.exposure<-b_sd2/z
		dat$beta.exposure<-b_sd2

	}
	return(dat)
}


format_exposure3<-function(dat=NULL,standardise_beta=FALSE,beta="beta",se="se",pval="pval",samplesize="samplesize",effect_allele="effect_allele",other_allele="other_allele",eaf="eaf",rsid="rsid",ID="population",snps=NULL,exposure="exposure"){
	names(dat)[names(dat) ==beta]<-"beta.exposure"
	names(dat)[names(dat) ==se]<-"se.exposure"
	names(dat)[names(dat) ==pval]<-"pval.exposure"
	names(dat)[names(dat) == samplesize]<-"samplesize.exposure"
	names(dat)[names(dat) ==effect_allele]<-"effect_allele.exposure"
	names(dat)[names(dat) ==other_allele]<-"other_allele.exposure"
	names(dat)[names(dat) ==eaf]<-"eaf.exposure"
	names(dat)[names(dat) ==rsid]<-"SNP"		
	names(dat)[names(dat) == exposure]<-"exposure"		
	if(any(dat$effect_allele.exposure == TRUE)) {
		warning("effect alleles set to TRUE. Changing to T")
		dat$effect_allele.exposure[dat$effect_allele.exposure == TRUE]<-"T"
	}
	if(!is.null(snps)){
		dat<-dat[dat$SNP %in% snps,]
	}

	
	if(!is.null(ID)){
		dat$id.exposure<-paste0("exp",as.numeric(as.factor(dat[,ID])))
	}else{
		dat$id.exposure<-paste0("exp",as.numeric(as.factor(dat$exposure)))
	}

	if(standardise_beta){
		b<-dat$beta.exposure
		se<-dat$se.exposure
		z<-b/se
		n<-dat$samplesize.exposure
		maf<-dat$eaf.exposure
		maf[maf>0.5]<-1-maf[maf>0.5]
		b_sd2<-b_sd(z=z,maf=maf,n=n)
		b
		b_sd2
		dat$se.exposure<-b_sd2/z
		dat$beta.exposure<-b_sd2

	}
	return(dat)
}
format_outcomes2<-function(dat=NULL,all_cols.keep=FALSE){

	names(dat)[names(dat) =="lnor"]<-"beta.outcome"
	names(dat)[names(dat) =="se"]<-"se.outcome"
	names(dat)[names(dat) =="p"]<-"pval.outcome"
	names(dat)[names(dat) %in% c("cases","ncase")]<-"ncase.outcome"
	names(dat)[names(dat) %in% c("controls","ncontrol")]<-"ncontrol.outcome"
	names(dat)[names(dat) =="Effect.Allele"]<-"effect_allele.outcome"
	names(dat)[names(dat) =="Other.Allele"]<-"other_allele.outcome"
	names(dat)[names(dat) =="eaf"]<-"eaf.outcome"
	names(dat)[names(dat) =="rsid"]<-"SNP"
	names(dat)[names(dat) =="ID"]<-"id.outcome"

	# dat[dat$population == "European; East Asian",c("id.outcome")]
	# exclude meta-analysed outcome data, to avoid MR of mixed European and East Asian populations. Need to MR East Asian and European studies separately prior to meta analysis
	dat<-dat[grep(";",dat$id.outcome,invert=T),]	
	dat$population[dat$id.outcome %in% 993:999]<-"European"
	dat$study.abbreviation[is.na(dat$study.abbreviation)]<-"survival"
	dat$outcome2<-paste(dat$outcome2,dat$id.outcome)
	if(any(names(dat) == "cancer")){
		dat$outcome2<-paste(dat$cancer,dat$id.outcome)
	}

	Dups<-unique(dat$outcome2[duplicated(dat$outcome2)])
	Pos<-dat$outcome2 %in% Dups
	if(any(Pos)) stop("duplicate outcomes when none expected")
	# dat$outcome2[Pos]<-paste(dat$outcome2[Pos],dat$id.outcome[Pos])
	outcome<-dat$outcome2
	dat$outcome2<-dat$cancer
	dat$outcome<-outcome

	# if(all(c("cancer","outcome") %in% names(dat))) warning("two column names for cancer outcome when only one expected")
	# names(dat)[names(dat) %in% c("cancer","outcome")]<-"outcome"
	dat$se.outcome<-as.numeric(dat$se.outcome)
	dat$beta.outcome<-as.numeric(dat$beta.outcome)
	dat$eaf.outcome<-as.numeric(dat$eaf.outcome)
	dat$pval.outcome<-as.numeric(dat$pval.outcome)
	# names(dat)[names(dat) =="cancer"]<-"outcome"
	Cols.keep<-c("SNP","outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","ncase.outcome","ncontrol.outcome","effect_allele.outcome","other_allele.outcome","id.outcome","population")

	dat<-dat[, Cols.keep]
	
	return(dat)
}

format_outcomes<-function(dat=NULL,IDS=NULL){
	names(dat)[names(dat) =="lnor"]<-"beta.outcome"
	names(dat)[names(dat) =="se"]<-"se.outcome"
	names(dat)[names(dat) =="p"]<-"pval.outcome"
	names(dat)[names(dat) =="cases"]<-"ncase.outcome"
	names(dat)[names(dat) =="controls"]<-"ncontrol.outcome"
	names(dat)[names(dat) =="Effect.Allele"]<-"effect_allele.outcome"
	names(dat)[names(dat) =="Other.Allele"]<-"other_allele.outcome"
	names(dat)[names(dat) =="eaf"]<-"eaf.outcome"
	names(dat)[names(dat) =="rsid"]<-"SNP"
	names(dat)[names(dat) =="cancer"]<-"outcome"
	dat$se.outcome<-as.numeric(dat$se.outcome)
	dat$beta.outcome<-as.numeric(dat$beta.outcome)
	dat$eaf.outcome<-as.numeric(dat$eaf.outcome)
	dat$pval.outcome<-as.numeric(dat$pval.outcome)

	# names(dat)[names(dat) =="cancer"]<-"outcome"
	dat$id.outcome<-dat$ID

	dat1<-dat[dat$ID %in% IDS,]
	dat2<-dat1[as.numeric(dat1$ncase.outcome)>1000,]
	dat3<-dat[!dat$ID %in%	IDS,]
	dat4<-dat1[as.numeric(dat1$ncase.outcome)<=1000,]
	dat5<-rbind(dat3,dat4)
	dat5$study.abbreviation[is.na(dat5$study.abbreviation)]<-"survival"
	dat5$outcome2<-paste(dat5$outcome,dat5$study.abbreviation)
	# Temp<-dat5[dat5$cancer2 %in% Dups,c("cancer","rsid","ID","lnor","se","Effect.Allele","Other.Allele","eaf","cases","study.abbreviation")]
	# Temp[order(Temp$cancer),]
	Dups<-unique(dat5$outcome2[duplicated(dat5$outcome2)])
	Pos<-dat5$outcome2 %in% Dups
	dat5$outcome2[Pos]<-paste(dat5$outcome2[Pos],dat5$ID[Pos])
	outcome<-dat5$outcome2
	dat5$outcome2<-dat5$outcome
	dat5$outcome<-outcome
	Cols.keep<-c("SNP","outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","ncase.outcome","ncontrol.outcome","effect_allele.outcome","other_allele.outcome")
	head(dat2[, Cols.keep])
	write.table(dat2[, Cols.keep],"~/fatty-acids/mr/data/outcome_discovery.txt",sep="\t",col.names=T,row.names=F,quote=F)
	write.table(dat5[,Cols.keep],"~/fatty-acids/mr/data/outcome_replication.txt",sep="\t",col.names=T,row.names=F,quote=F)

	return(list(dat2,dat5))
}


format_exposure<-function(dat=NULL,exposure=NULL,snp=NULL){
	# snps<-c("rs2236212", "rs174576",  "rs603424",  "rs174546" )
	# Exposure<-c("rs174546 AA_to_DGLA","rs2236212 DHA_to_DPA_n3","rs603424 POA_to_PA","rs174576 GLA_to_LA")
	# exposure_dat<-dat[paste(dat$SNP,dat$exposure) %in% Exposure,]
	exposure_dat<-dat[dat$exposure == "AA_to_DGLA" ,]
	dat<-dat[dat$SNP == snp & dat$exposure == exposure,]
	return(exposure_dat)
}

tsmr2<-function(outcome_dat=NULL,exposure_dat=NULL,action=2){
	outcome_dat$population[outcome_dat$id.outcome %in% 993:999]<-"European"
	outcome_rep$population[outcome_rep$id.outcome %in% 993:999]<-"European"
	dat <- harmonise_data(exposure_dat = exposure_dat,outcome_dat =outcome_dat,action=action)
	mr_res<-mr(dat,method_list=c("mr_wald_ratio"))		
	# meta.dat1<-exposure_dat[,c("id.exposure","SNP","population","beta.exposure","se.exposure")]
	# meta.dat2<-outcome_dat[,names(outcome_dat) != "outcome"]
	mr_res1<-merge(mr_res,dat[,!names(dat) %in% c("exposure","outcome")],by=c("id.outcome","id.exposure"))
	mr_res1$OR<-exp(mr_res1$b) 
	mr_res1$LCI<-exp(mr_res1$b-1.96*mr_res1$se)
	mr_res1$UCI<-exp(mr_res1$b+1.96*mr_res1$se)
	return(mr_res1)
}


tsmr<-function(outcome_dat=NULL,exposure_dat=NULL,meta.dat=NULL){	
	
	outcome_dat1<-outcome_dat[outcome_dat$population == "European",]
	outcome_dat2<-outcome_dat[outcome_dat$population == "East Asian",]

	exposure_dat1<-exposure_dat[exposure_dat$population == "European", ]
	exposure_dat2<-exposure_dat[exposure_dat$population == "East Asian",]

	dat1 <- harmonise_data(exposure_dat = exposure_dat1,outcome_dat =outcome_dat1,action=2)
	dat2 <- harmonise_data(exposure_dat = exposure_dat2,outcome_dat =outcome_dat2,action=2)
	dat<-rbind(dat1,dat2)
	# dat$exposure<-dat$exposure
	# dat[dat$exposure == exposures[1],c("exposure","beta.exposure","se.exposure","beta.outcome","se.outcome","effect_allele.exposure","effect_allele.outcome","other_allele.exposure","other_allele.outcome")]
	mr_res<-mr(dat,method_list=c("mr_wald_ratio"))		
	mr_res1<-merge(mr_res,meta.dat ,by.x="id.outcome",by.y="ID",all.x=T)		
	mr_res1$population[mr_res1$id.outcome %in% c("1499" ,"993",  "994",  "995",  "996",  "997",  "998",  "999" )]<-"European"
	dat_q<-dat_outcomes_final[,c("ID","Q.p"  ) ]
	mr_res2<-merge(mr_res1,dat_q,by.x="id.outcome",by.y="ID")	
	mr_res2$OR<-exp(mr_res2$b)
	mr_res2$LCI<-exp(mr_res2$b-1.96*mr_res2$se)
	mr_res2$UCI<-exp(mr_res2$b+1.96*mr_res2$se)	
	mr_res3<-merge(mr_res2,exposure_dat[,c("SNP","id.exposure","samplesize.exposure")],by="id.exposure")	
	mr_res4<-merge(mr_res3,outcome_dat[,c("id.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome")],by="id.outcome")	
	return(mr_res4)
}


b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}


format_mr_res<-function(dat=NULL,snp_keep=NULL,keep_adjrs174546=TRUE){
	dat<-dat[dat$SNP == snp_keep,]
	dat<-dat[dat$population.x == dat$population.y,]
	dat<-dat[order(dat$exposure),]
	if(!keep_adjrs174546) dat<-dat[dat$exposure != "GLA:LAadj_rs174546",]
	return(dat)
}

meta_analysis_v2<-function(dat=NULL,beta.col=NULL,se.col=NULL,outcome="outcome",ncase="ncase",ncontrol="ncontrol"){
	# not necessary to meta analyse by exposure because there is only one exposure
	# exposures<-unique(dat$exposure)
	# meta_results<-NULL
	# for(i in 1:length(exposures)){
		# print(i)
		# dat1<-dat[dat$exposure==exposures[i],]
	if(any(duplicated(dat$id.outcome))) stop("duplicate ids present")
	b<-dat[,beta.col]
	se<-dat[,se.col]	
	# p<-temp$p
	w<-1/se^2
	b.fixed<-sum(b*w)/(sum(w))
	se.fixed<-sqrt(sum(w)^-1)
	z<-abs(b.fixed/se.fixed)
	p.fixed<-pnorm(z,lower.tail=F)*2
	nstudies.fixed<-length(b)
	if(all(dat$id.outcome %in% c("49","86"))){
		outcome<-"Cancer.Group"
	}
	cancer<-unique(dat[,outcome])
	
	
	if(all(dat$id.outcome %in% c("138", "165", "35"))){
		cancer<-"Central nervous system and eye cancer"	
	}	
	
	if(all(dat$id.outcome %in% c("49", "86", "5" ))){
		cancer<-unique(dat$Cancer.Group)
		system<-unique(dat$system)
		site<-unique(dat$site)
		cancer.group<-unique(dat$Cancer.Group)
		cell<-"lymphocytes"
	}	
	
	cases<-sum(dat[,ncase])
	controls<-sum(dat[,ncontrol])
	# samplesize.exposure<-sum(dat1$samplesize.exposure)
	study<-"Overall fixed effect"
	Q<-sum((b.fixed-b)^2*w)
	df.Q<-length(b)-1		
	Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
	EA<-unique(dat$effect_allele.outcome)
	OA<-unique(dat$other_allele.outcome)
	EAF<-round(sum((dat$eaf.outcome*w))/sum(w),3)
	exposure<-unique(dat$exposure)
	exposure<-exposure[order(exposure,decreasing=TRUE)] #ensure is in same order as population. East Asian comes before European, GLA:LA always corresponds to East Asian and AA:DGLA always corresponds to Eurpean
	exposure<-paste(exposure,collapse="/")
	population<- paste(sort(unique(dat$population)),collapse="; ")
	id.outcome<-paste(sort(as.numeric(dat$id.outcome)),collapse="; ")
	# exposure<-unique(dat1$exposure)
	dat.matrix<-c(population,cancer,id.outcome,b.fixed,se.fixed,p.fixed,nstudies.fixed,cases,controls,study,Q.p,EA,OA,EAF,exposure)
	meta_results<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
	Res<-meta_results
	# }


	# Res<-do.call(rbind,meta_results)	

	# if(sum(grep("exposure",se.col))==1) pval<-"pval.exposure"
	# if(sum(grep("outcome",se.col))==1) pval<-"pval.outcome"
	# if(se.col=="se") pval<-"pval"

	names(Res)<-c("population","outcome","id.outcome",beta.col,se.col,"pval","nstudies",ncase,ncontrol,"study","Q.p","effect_allele.outcome","other_allele.outcome","eaf.outcome","exposure")	
	meta.tab<-unique(meta.tab9[,c("cancer","system","site","cell","Cancer.Group")])
	meta.tab1<-meta.tab[meta.tab$cancer != "Acute lymphoblastic leukaemia",]
	meta.tab2<-meta.tab[meta.tab$cancer == "Acute lymphoblastic leukaemia",]
	meta.tab2$cell<-"B lymphocytes"
	meta.tab<-rbind(meta.tab1,meta.tab2)	

	if(all(dat$id.outcome %in% c("49", "86", "5" ))){
		meta.tab$system<-system
		meta.tab$site<-site
		meta.tab$Cancer.Group<-cancer.group
		meta.tab$cell<-cell 
		meta.tab$cancer<-cancer
		meta.tab<-unique(meta.tab)
		if(nrow(meta.tab)>1) stop("nrows greater than 1")
	}

	which(meta.tab$cancer == Res$outcome)
	Res<-merge(Res,meta.tab,by.x="outcome",by.y="cancer")
	if(nrow(Res)==0) stop("outcome missing from meta table")
	return(Res)
}


meta_analysis<-function(dat=NULL,beta.col=NULL,se.col=NULL,outcome="outcome",ncase="ncase",ncontrol="ncontrol"){
	exposures<-unique(dat$exposure)
	meta_results<-NULL
	for(i in 1:length(exposures)){
		# print(i)
		dat1<-dat[dat$exposure==exposures[i],]
		if(any(duplicated(dat$id.outcome))) stop("duplicate ids present")
		b<-dat1[,beta.col]
		se<-dat1[,se.col]	
		# p<-temp$p
		w<-1/se^2
		b.fixed<-sum(b*w)/(sum(w))
		se.fixed<-sqrt(sum(w)^-1)
		z<-abs(b.fixed/se.fixed)
		p.fixed<-pnorm(z,lower.tail=F)*2
		nstudies.fixed<-length(b)
		if(all(dat1$id.outcome %in% c("49","86"))){
			outcome<-"Cancer.Group"
		}

		cancer<-unique(dat1[,outcome])
		if(all(dat1$id.outcome %in% c("138", "165", "35"))){
			cancer<-"Central nervous system and eye cancer"	
		}	
		cases<-sum(dat1[,ncase])
		controls<-sum(dat1[,ncontrol])
		# samplesize.exposure<-sum(dat1$samplesize.exposure)
		study<-"Overall fixed effect"
		Q<-sum((b.fixed-b)^2*w)
		df.Q<-length(b)-1		
		Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
		EA<-unique(dat1$effect_allele.outcome)
		OA<-unique(dat1$other_allele.outcome)
		EAF<-round(sum((dat1$eaf.outcome*w))/sum(w),3)
		exposure<-exposures[i]
		population<- paste(unique(dat1$population),collapse="; ")
		id.outcome<-paste(dat1$id.outcome,collapse="; ")
		# exposure<-unique(dat1$exposure)
		dat.matrix<-c(exposure,population,cancer,id.outcome,b.fixed,se.fixed,p.fixed,nstudies.fixed,cases,controls,study,Q.p,EA,OA,EAF)
		meta_results[[i]]<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
	}

	Res<-do.call(rbind,meta_results)	

	if(sum(grep("exposure",se.col))==1) pval<-"pval.exposure"
	if(sum(grep("outcome",se.col))==1) pval<-"pval.outcome"
	if(se.col=="se") pval<-"pval"

	names(Res)<-c("exposure","population","outcome","id.outcome",beta.col,se.col,pval,"nstudies",ncase,ncontrol,"study","Q.p","effect_allele.outcome","other_allele.outcome","eaf.outcome")	
	meta.tab<-unique(meta.tab9[,c("cancer","system","site","cell","Cancer.Group")])
	meta.tab1<-meta.tab[meta.tab$cancer != "Acute lymphoblastic leukaemia",]
	meta.tab2<-meta.tab[meta.tab$cancer == "Acute lymphoblastic leukaemia",]
	meta.tab2$cell<-"B lymphocytes"
	meta.tab<-rbind(meta.tab1,meta.tab2)		
	Res<-merge(Res,meta.tab,by.x="outcome",by.y="cancer")

	return(Res)
}




format_plot_data<-function(mr_res=mr_res1,beta=NULL,se=NULL,modify_weight=NULL,beta_reverse=TRUE,exp=FALSE){
	Dat<-meta_analysis(dat=mr_res,beta.col=beta,se.col=se)
	if(sum(grep("exposure",se.col))==1) pval<-"pval.exposure"
	if(sum(grep("outcome",se.col))==1) pval<-"pval.outcome"
	if(se.col=="se") pval<-"pval"
	Dat<-rbind.fill(Dat,mr_res[,c("exposure","outcome",beta,se,pval,"ncase.outcome","ncontrol.outcome","study")])
	Dat[,beta]<-as.numeric(Dat[,beta])
	Dat[,se]<-as.numeric(Dat[,se])
	if(beta_reverse) Dat[,beta]<-Dat[,beta]*-1
	Dat[,se]<-as.numeric(Dat[,se])
	Dat$weight<-1/Dat[,se]/modify_weight
	Dat$population<-"European"
	Dat$population[Dat$study == "ACCC"]<-"East Asian"
	Dat$population[Dat$study == "Overall fixed effect"]<-"Combined"
	Dat1<-Dat[Dat$population=="Combined",]
	Dat2<-Dat[Dat$population!="Combined",]	
	Dat<-rbind(Dat2,Dat1)
	Dat$lci<-Dat[,beta]-Dat[,se]*1.96
	Dat$uci<-Dat[,beta]+Dat[,se]*1.96
	if(exp) {
		Dat$OR<-exp(Dat[,beta]) 
		Dat$lci<-exp(Dat$lci)
		Dat$uci<-exp(Dat$uci)
	}
	return(Dat)	
}

effect_allele_gtex_data<-function(dat=NULL){
	variant_id<-unlist(strsplit(dat$variant_id,split="_"))
	dat$effect_allele<-variant_id[seq(4,length(variant_id),by=5)]
	dat$other_allele<-variant_id[seq(3,length(variant_id),by=5)]
	return(dat)
}

effect_allele_bbj_data<-function(dat=NULL){
	dat$effect_allele<-dat$ALT
	dat$other_allele<-dat$REF
	return(dat)
}


format_ens<-function(dat=NULL,ens=c("ENSG00000149485","ENSG00000134824")){
	dat<-dat[dat$Ens %in% ens,]	 
	dat$gene[dat$Ens == "ENSG00000149485"] <- "FADS1"
	dat$gene[dat$Ens == "ENSG00000134824"] <- "FADS2"	
	return(dat)	
	# 2  FADS1  rs174528 ENSG00000197977
	# 3  FADS2  rs174528 ENSG00000134824
}



# format_metareg2<-function(){
# 	load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
# 	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
# 	"Mouth & throat cancer"
# 	"Esophageal cancer",
# }

format_metareg2<-function(){
	load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")

	mr_res1<-fix_ids(Dat=mr_res1,id_col="id.outcome")

	disc.tab9<-fix_ids(Dat=disc.tab9,id_col="ID")
	meta1<-disc.tab9[,c("ID2","study.abbreviation")]
	meta2<-	meta.tab9[c("ID","study.abbreviation")]
	names(meta1)[names(meta1) == "ID2"]<-"ID"
	meta_dat<-unique(rbind(meta1,meta2))
	
	mr_res1<-merge(mr_res1,meta_dat,by.x="ID2",by.y="ID")

	# mr_res$ID2[is.na(mr_res1$study.abbreviation.y)]
	mr_res2<-mr_res1[mr_res1$population == "European",]
	mr_res2<-mr_res2[mr_res2$exposure == "AA:DGLA",]
	mr_res3<-mr_res1[mr_res1$population != "European",]
	mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	mr_res1<-rbind(mr_res2,mr_res3)
	mr_res1$cases<-as.numeric(mr_res1$cases)	
	mr_res1$cancer[is.na(mr_res1$cancer)]<-mr_res1$outcome[is.na(mr_res1$cancer)]
	mr_res1<-mr_res1[!mr_res1$id.outcome %in% c("993","994","995","996","997","998","999","1499"),]
	mr_res1$b<-as.numeric(mr_res1$b)
	mr_res1$se<-as.numeric(mr_res1$se)
	mr_res1<-mr_res1[order(mr_res1$cases,decreasing=T),]
	Infl_cancers<-c("Uveal melanoma","Diffuse large b cell lymphoma","Melanoma","Gastric adenocarcinoma","Pancreatic cancer", "Oral cavity and pharyngeal cancer" , "Lung cancer","Esophageal adenocarcinoma","Esophageal squamous cell carcinoma","Colorectal cancer",	"Bladder cancer","Pleural mesothelioma", "Marginal zone lymphoma","Malignant skin cancer")
	sort(unique(mr_res1$cancer))
	
	mr_res2<-mr_res1[mr_res1$cancer %in% Infl_cancers,]
	mr_res2$infl<-1
	mr_res2<-mr_res2[!duplicated(mr_res2$cancer),]
	mr_res3<-mr_res1[!mr_res1$cancer %in% Infl_cancers,]
	mr_res3$infl<-0
	mr_res3<-mr_res3[!duplicated(mr_res3$Cancer.Group),]
	mr_res3<-mr_res3[!mr_res3$cancer %in% c("Cancer (all cause)","Cancer of digestive organs", "Kidney cancer in males","Respiratory and intrathoracic cancer","Kidney cancer in females","Colorectal cancer in males"  , "Lung cancer in ever smokers" ,"Oral cancer" ,"Noncardia gastric adenocarcinoma" ,"Liver & bile duct cancer"    ) ,]
	mr_res1<-rbind(mr_res2,mr_res3)
	# unique(mr_res3$cancer)
	# cancer.group2<-c(
	# 	"Bladder cancer",
	# 	"Colorectal cancer",
	# 	"Esophageal cancer",
	# 	"Lung cancer",
	# 	"Mouth & throat cancer",
	# 	"Pancreatic cancer",
	# 	"Gastric adenocarcinoma")
		
		
	# unique(mr_res1$Cancer.Group)

	
	# mr_res1$infl<-0
	# mr_res1$infl[mr_res1$Cancer.Group %in% cancer.group2]<-1
	# mr_res1$infl[mr_res1$cancer %in% c("Melanoma","Uveal melanoma")]<-1
	
	# mr_res2<-mr_res1[mr_res1$cancer %in% c("Melanoma","Uveal melanoma") ,]
	# mr_res2<-mr_res2[!duplicated(mr_res2$cancer),]
	# mr_res3<-mr_res1[mr_res1$Cancer.Group == "Esophageal cancer",]
	# mr_res3<-mr_res3[!duplicated(mr_res3$cancer),]
	# mr_res4<-mr_res1[mr_res1$Cancer.Group != "Esophageal cancer" & !mr_res1$cancer %in% c("Melanoma","Uveal melanoma"),]
	# mr_res4<-mr_res4[!duplicated(mr_res4$Cancer.Group),]
	# mr_res4<-mr_res4[!mr_res4$cancer %in% c("Cancer (all cause)","Cancer of digestive organs", "Kidney cancer in males","Respiratory and intrathoracic cancer","Kidney cancer in females"  ) ,]
	# mr_res2<-rbind.fill(mr_res2,mr_res4)	
	# mr_res2<-rbind.fill(mr_res2,mr_res3)	
	
	# mr_res2$cancer[mr_res2$infl==0]
	# id(prune_blood_cancers){
	# 	mr_res2[mr_res2$cancer]


	# }
	
	return(mr_res1)
}

format_metareg<-function(prune_blood_cancers=FALSE,drop_skin_cancers=FALSE){
	load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")

	mr_res1<-fix_ids(Dat=mr_res1,id_col="id.outcome")

	disc.tab9<-fix_ids(Dat=disc.tab9,id_col="ID")
	meta1<-disc.tab9[,c("ID2","study.abbreviation")]
	meta2<-	meta.tab9[c("ID","study.abbreviation")]
	names(meta1)[names(meta1) == "ID2"]<-"ID"
	meta_dat<-unique(rbind(meta1,meta2))
	
	mr_res1<-merge(mr_res1,meta_dat,by.x="ID2",by.y="ID")
	table(mr_res1$Cancer.Group)
	# mr_res$ID2[is.na(mr_res1$study.abbreviation.y)]
	mr_res2<-mr_res1[mr_res1$population == "European",]
	mr_res2<-mr_res2[mr_res2$exposure == "AA:DGLA",]
	mr_res3<-mr_res1[mr_res1$population != "European",]
	mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	mr_res1<-rbind(mr_res2,mr_res3)
	mr_res1$cases<-as.numeric(mr_res1$cases)	
	mr_res1$cancer[is.na(mr_res1$cancer)]<-mr_res1$outcome[is.na(mr_res1$cancer)]
	mr_res1<-mr_res1[!mr_res1$id.outcome %in% c("993","994","995","996","997","998","999","1499"),]

	cancer.group<-c(
		"Bladder cancer",
		"Colorectal cancer",
		"Esophageal cancer",
		"Kidney cancer",
		"Liver & biliary tract cancer",
		"Lung cancer",
		"Mouth & throat cancer",
		"Pancreatic cancer",
		"Gastric adenocarcinoma")

	cancer.group2<-c(
		"Bladder cancer",
		"Colorectal cancer",
		"Esophageal cancer",
		"Liver & biliary tract cancer",
		"Lung cancer",
		"Mouth & throat cancer",
		"Pancreatic cancer",
		"Gastric adenocarcinoma","Skin cancer",
		"Non-hodgkin's lymphoma" ,"Ovarian cancer"   )

	# unique(mr_res1$Cancer.Group)


	mr_res1$infl<-0
	mr_res1$infl[mr_res1$Cancer.Group %in% cancer.group2]<-1
	mr_res1$smoking1<-0 #all cancers with no evidence or unclear evidence
	mr_res1$smoking1[mr_res1$Cancer.Group %in% cancer.group]<-1
	mr_res1$smoking1[mr_res1$cancer %in% c("Cervical cancer") ]<-1
	# mr_res1$smoking1[mr_res1$cancer %in% c("Cervical cancer","Myeloid leukaemia") ]<-1
	# mr_res1$smoking1[mr_res1$cancer=="Endometrial cancer"]<-1
	# mr_res1$smoking1[mr_res1$cancer %in% c("Cervical cancer","Myeloid leukaemia") ]<-1

	mr_res1$smoking2<-0 #all cancers with unclear evidence
	mr_res1$smoking2[mr_res1$Cancer.Group %in% cancer.group]<-1
	mr_res1$smoking2[mr_res1$cancer == "Cervical cancer"]<-1
	mr_res1$smoking2[which(mr_res1$site == "Brain")]<-2
	mr_res1$smoking2[which(mr_res1$Cancer.Group %in% c("Breast cancer","Prostate cancer"))]<-2 #no evidence

	mr_res1$smoking3<-0 #all cancers with unclear evidence
	mr_res1$smoking3[mr_res1$Cancer.Group %in% cancer.group]<-1
	mr_res1$smoking3[mr_res1$cancer == "Endometrial cancer"]<-1 #1 = cancers with causal association but endometrial cancer shows protective effect 
	mr_res1$smoking3[which(mr_res1$site == "Brain")]<-2
	mr_res1$smoking3[which(mr_res1$Cancer.Group %in% c("Breast cancer","Prostate cancer"))]<-2

	Cancer.group2<-c(c("Breast cancer","Prostate cancer"),cancer.group)

	mr_res2<-mr_res1[mr_res1$Cancer.Group %in% Cancer.group2,]
	mr_res2<-mr_res2[order(mr_res2$cases,decreasing=T),]
	mr_res4<-mr_res2[!mr_res2$Cancer.Group %in% c("Esophageal cancer"),]
	mr_res3<-mr_res2[mr_res2$Cancer.Group %in% c("Esophageal cancer"),]
	# mr_res3$cancer[duplicated(mr_res3$cancer)]
	mr_res3<-mr_res3[!duplicated(mr_res3$cancer),]
	# Dups<-unique(mr_res3$cancer[duplicated(mr_res3$cancer)])
	# mr_res3$cancer[mr_res3$cancer %in% Dups]
	# unique(mr_res4$cancer[mr_res4$Cancer.Group %in% Dups])
	mr_res4<-mr_res4[!duplicated(mr_res4$Cancer.Group),]
	mr_res2<-rbind(mr_res3,mr_res4)
	
	mr_res3<-mr_res1[!mr_res1$Cancer.Group %in% Cancer.group2,]	

	mr_res3<-mr_res3[order(mr_res3$cases,decreasing=T),]
	mr_res3<-mr_res3[!duplicated(mr_res3$cancer),]
	
	# mr_res3<-mr_res3[!duplicated(mr_res3$Cancer.Group),]

	mr_res<-rbind(mr_res3,mr_res2)
	
	# table(mr_res3$site)
	# mr_res3[mr_res3$site == "Ovary",c("cancer","study","id.outcome","Cancer.Group","cases","study.abbreviation")]
	# unique(meta.tab9[,c("ID","study.abbreviation")])
	# unique(disc.tab9[,c("ID","study.abbreviation")])
	# # mr_res3$cancer %in% c(""Cancer (all cause)","Cancer (excluding non-melanoma skin cancer)" )

	# unique(mr_res3[,c("Cancer.Group","study")])
	mr_res$b<-as.numeric(mr_res$b)
	mr_res$se<-as.numeric(mr_res$se)

	if(prune_blood_cancers){
		mr_res_blood<-mr_res[mr_res$site == "Blood",]
		mr_res<-mr_res[mr_res$site != "Blood",]
	 	mr_res_blood<-mr_res_blood[!duplicated(mr_res_blood$study.abbreviation.y),]
	 	mr_res_blood<-mr_res_blood[!mr_res_blood$cancer %in% c("Follicular lymphoma","Multiple myeloma","Non-follicular lymphoma"),]
		mr_res<-rbind(mr_res,mr_res_blood)	
	}			
	
	mr_res<-mr_res[!mr_res$Cancer.Group %in% c("Digestive system cancer","All cause"),]
	
	mr_res$Cancer.Group[mr_res$cancer == "Kidney cancer in males"]
	mr_res<-mr_res[!mr_res$cancer %in% c("Kidney cancer in males","Kidney cancer in females" ,"Urinary tract cancer", "Respiratory and intrathoracic cancer","Lymphoid leukaemia","Leukaemia","Malignant non-melanoma skin cancer","Malignant skin cancer","Blood cancer","Lymphoma","Brain cancer"),]	
	# "Endometrial cancer"
	
	mr_res1<-mr_res[mr_res$Cancer.Group == "Ovarian cancer",]
	mr_res1<-mr_res1[!duplicated(mr_res1$Cancer.Group),]
	mr_res2<-mr_res[mr_res$Cancer.Group != "Ovarian cancer",]
	mr_res<-rbind(mr_res1,mr_res2)	
	# table(mr_res$site)
	# mr_res<-fix_ids(Dat=mr_res,id_col="id.outcome")
	# disc.tab9<-fix_ids(Dat=disc.tab9,id_col="ID")
	
	# Temp<-mr_res[mr_res$site == "Blood",c("study","outcome","id.outcome","study.abbreviation","ID2")]
	# Temp2<-merge(Temp,disc.tab9,by="ID2",all.x=T)
	# Temp2[order(Temp2$cases,decreasing=T),c("cancer","outcome","study.y","study.abbreviation.y",c("cases","id.outcome"))]
	if(drop_skin_cancers){
		mr_res<-mr_res[mr_res$site!="Skin",]
	}

	mr_res$Colour<-mr_res$smoking1
	mr_res$Colour[mr_res$cancer %in% c("Meningioma","Glioma","Central nervous system and eye cancer","Breast cancer","Prostate cancer")]<-2
	mr_res$Colour[mr_res$cancer == "Endometrial cancer"]<-3
	mr_res$Colour[mr_res$Colour == 0]<-4
	mr_res$weight<-1/mr_res$se/5
	mr_res$Colour1[mr_res$Colour == 1]<-"Increases risk"
	mr_res$Colour1[mr_res$Colour == 2]<-"No relationship"
	mr_res$Colour1[mr_res$Colour == 3]<-"Decreases risk"
	mr_res$Colour1[mr_res$Colour == 4]<-"Unclear relationship"

	mr_res<-mr_res[order(mr_res$cases,decreasing=T),]
	mr_res1<-mr_res[mr_res$Colour == 1,]
	mr_res2<-mr_res[mr_res$Colour == 2,]
	mr_res3<-mr_res[mr_res$Colour == 3,]
	mr_res4<-mr_res[mr_res$Colour == 4,]
	mr_res<-do.call(rbind,list(mr_res1,mr_res2,mr_res3,mr_res4))
	# mr_res[,c("cases","Colour1")]
	mr_res$weight<-1/mr_res$se^2
	mr_res$plot_name<-paste0(mr_res$cancer,"\n",mr_res$cases)
	mr_res$smoking <- mr_res$Colour1
	mr_res$Shape<-"square"
	return(mr_res)
}

# https://www.ncbi.nlm.nih.gov/books/NBK294317/table/ch4.t1/
#evidence sufficient to infer a caual relationship
	# bladder cancer
	# cervical cancer
	# colorectal cancer
	# Esophageal cancer
	# kidney cancer
	# laryngeal cancer 
	# acute myeloid leukemia
	# liver cancer (hepatocellular carcinoma)
	# lung cancer 
	# oral cavity and pharyngeal cancer
	# pancreatic cancer
	# gastric cancer


# evidence is suggestive of no causal effect 
	# brain cancer
	# breast cancer
	# prostate cancer

# The evidence is sufficient to infer that current smoking reduces the risk
	# endometrial cancer in postmenopausal women

# unclear evidence
	# ovarian cancer 
	# everything else 



# smoking_cancers<-c("Colorectum","Esophagus","Kidney and renal pelvis","Larynx","Liver and intrahepatic bile duct","Lung, bronchus, and trachea","Myeloid leukemia","Oral cavity and pharynx","Pancreas","Stomach","Urinary bladder","Uterine cervix")





fix_ids<-function(Dat=NULL,id_col=NULL){
	Dat$ID2<-Dat[,id_col]
	for(i in 1:length(Dat[,id_col])){
		IDS<-as.numeric(unlist(strsplit(Dat[,id_col][i],split=";")))
		IDS<-IDS[order(IDS)]
		IDS<-paste(IDS,collapse="; ")
		Dat$ID2[i]<-IDS
	}
	return(Dat)
}



format_mrresults<-function(){	
	load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")	
 	mr_res1<-mr_res1[!mr_res1$id.outcome %in% c("993","994","995","996","997","998","999","1499"),]
	mr_res1<-mr_res1[mr_res1$population == "European",]
	mr_res1<-mr_res1[mr_res1$exposure == "AA:DGLA",]
	mr_res3<-mr_res1[mr_res1$population != "European",]
	mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	mr_res1<-rbind(mr_res1,mr_res3)
	mr_res1$cases<-as.numeric(mr_res1$cases)	
	mr_res1$cancer[is.na(mr_res1$cancer)]<-mr_res1$outcome[is.na(mr_res1$cancer)]
	mr_res1$OR<-round(exp(as.numeric(mr_res1$b)),2)
	mr_res1$LCI<-round(exp(as.numeric(mr_res1$b)-1.96*as.numeric(mr_res1$se)),2)
	mr_res1$UCI<-round(exp(as.numeric(mr_res1$b)+1.96*as.numeric(mr_res1$se)),2)
	mr_res1<-mr_res1[order(as.numeric(mr_res1$cases),decreasing=T),]
	mr_res1$b<-as.numeric(mr_res1$b)	
	mr_res1$se<-as.numeric(mr_res1$se)	
	mr_res1$plot_name<-paste0(mr_res1$cancer,"\n",mr_res1$cases)
	return(mr_res1)
}


format_metareg3<-function(){
	load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")

	mr_res1<-fix_ids(Dat=mr_res1,id_col="id.outcome")

	disc.tab9<-fix_ids(Dat=disc.tab9,id_col="ID")
	meta1<-disc.tab9[,c("ID2","study.abbreviation")]
	meta2<-	meta.tab9[c("ID","study.abbreviation")]
	names(meta1)[names(meta1) == "ID2"]<-"ID"
	meta_dat<-unique(rbind(meta1,meta2))
	
	mr_res1<-merge(mr_res1,meta_dat,by.x="ID2",by.y="ID")

	# mr_res$ID2[is.na(mr_res1$study.abbreviation.y)]
	mr_res2<-mr_res1[mr_res1$population == "European",]
	mr_res2<-mr_res2[mr_res2$exposure == "AA:DGLA",]
	mr_res3<-mr_res1[mr_res1$population != "European",]
	mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	mr_res1<-rbind(mr_res2,mr_res3)
	mr_res1$cases<-as.numeric(mr_res1$cases)	
	mr_res1$cancer[is.na(mr_res1$cancer)]<-mr_res1$outcome[is.na(mr_res1$cancer)]
	mr_res1<-mr_res1[!mr_res1$id.outcome %in% c("993","994","995","996","997","998","999","1499"),]
	mr_res1$b<-as.numeric(mr_res1$b)
	mr_res1$se<-as.numeric(mr_res1$se)
	mr_res1<-mr_res1[order(mr_res1$cases,decreasing=T),]
	


	# cancers associated with infectious agents Nature. 2002 December 19; 420(6917): 860–867. doi:10.1038/nature01322.Coussens
		# Cholangiosarcoma, 
		# colon carcnoma, 
		# Gall bladder cancer
		# Gastric adenocarcinoma, 
		# MALT (Low-grade primary MALT (mucosa-associated lymphoid tissue) lymphoma of the stomach )
		# Hepatocellular carcinoma
		# B-cell non-Hodgkin’s lymphoma, 
		# Burkitts lymphoma,
		# Non-Hodgkin’s lymphoma, 
		# squamous cell carcinomas, 
		# Kaposi’s sarcoma
		# Skin carcinoma in draining sinuses
		# Ovarian carcinoma, 
		# cervical/anal carcinoma 
		# Bladder, 
		# liver, 
		# rectal carcinoma, 
		# follicular lymphoma of the spleen

	# how infectious agents cancers matched to cancers in FAMRC
	agents<-c("Biliary tract cancer",# Cholangiosarcoma, 
		"Colon cancer",# colon carcnoma, 
		"Gastric adenocarcinoma",# Gastric adenocarcinoma, 
		 "Liver cancer",# Hepatocellular carcinoma
		 "B cell non-hodgkin lymphoma",   # B-cell non-Hodgkin’s lymphoma,     
		 "Non-hodgkin lymphoma unspecified",# Non-Hodgkin’s lymphoma, 
		 "Ovarian cancer",     # Ovarian carcinoma, 
		 "Cervical cancer",   # cervical/anal carcinoma   
		 "Bladder cancer",# Bladder, 
		  "Rectal cancer") # rectal carcinoma, 

	Infl_cancers<-c("Uveal melanoma","Diffuse large b cell lymphoma","Melanoma","Gastric adenocarcinoma","Pancreatic cancer", "Oral cavity and pharyngeal cancer" , "Lung cancer","Esophageal adenocarcinoma","Esophageal squamous cell carcinoma","Colorectal cancer",	"Bladder cancer","Pleural mesothelioma", "Marginal zone lymphoma","Malignant skin cancer")

	sort(unique(mr_res1$cancer))
	mr_res1<-mr_res1[!duplicated(mr_res1$cancer),]
	mr_res2<-mr_res1[mr_res1$cancer %in% unique(c(agents,Infl_cancers)),]
	mr_res2$infl<-1
	
	mr_res3<-mr_res1[!mr_res1$cancer %in% unique(c(agents,Infl_cancers)),]
	mr_res3$infl<-0
	mr_res3<-mr_res3[!duplicated(mr_res3$Cancer.Group),]
	mr_res3<-mr_res3[!mr_res3$cancer %in% c("Cancer (all cause)","Cancer of digestive organs", "Kidney cancer in males","Respiratory and intrathoracic cancer","Kidney cancer in females","Colorectal cancer in males"  , "Lung cancer in ever smokers" ,"Oral cancer" ,"Noncardia gastric adenocarcinoma" ,"Liver & bile duct cancer"    ) ,]
	mr_res1<-rbind(mr_res2,mr_res3)
	# unique(mr_res3$cancer)
	# cancer.group2<-c(
	# 	"Bladder cancer",
	# 	"Colorectal cancer",
	# 	"Esophageal cancer",
	# 	"Lung cancer",
	# 	"Mouth & throat cancer",
	# 	"Pancreatic cancer",
	# 	"Gastric adenocarcinoma")
		
		
	# unique(mr_res1$Cancer.Group)

	
	# mr_res1$infl<-0
	# mr_res1$infl[mr_res1$Cancer.Group %in% cancer.group2]<-1
	# mr_res1$infl[mr_res1$cancer %in% c("Melanoma","Uveal melanoma")]<-1
	
	# mr_res2<-mr_res1[mr_res1$cancer %in% c("Melanoma","Uveal melanoma") ,]
	# mr_res2<-mr_res2[!duplicated(mr_res2$cancer),]
	# mr_res3<-mr_res1[mr_res1$Cancer.Group == "Esophageal cancer",]
	# mr_res3<-mr_res3[!duplicated(mr_res3$cancer),]
	# mr_res4<-mr_res1[mr_res1$Cancer.Group != "Esophageal cancer" & !mr_res1$cancer %in% c("Melanoma","Uveal melanoma"),]
	# mr_res4<-mr_res4[!duplicated(mr_res4$Cancer.Group),]
	# mr_res4<-mr_res4[!mr_res4$cancer %in% c("Cancer (all cause)","Cancer of digestive organs", "Kidney cancer in males","Respiratory and intrathoracic cancer","Kidney cancer in females"  ) ,]
	# mr_res2<-rbind.fill(mr_res2,mr_res4)	
	# mr_res2<-rbind.fill(mr_res2,mr_res3)	
	
	# mr_res2$cancer[mr_res2$infl==0]
	# id(prune_blood_cancers){
	# 	mr_res2[mr_res2$cancer]


	# }
	
	return(mr_res1)
}

format_results_mvmr<-function(res=NULL)
{
		load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
		out<-unique(outcome_dat[,c("id.outcome","population","ncase.outcome","ncontrol.outcome","outcome2")])		
		out$ncase.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
			median(out$ncase.outcome[out$id.outcome==x])))
		out$ncontrol.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
			median(out$ncontrol.outcome[out$id.outcome==x])))
		out<-unique(out)

		res.m<-merge(res,out,by.x="result.id.outcome",by.y="id.outcome")
		res.m<-merge(res.m,meta.tab9[,c("ID","Cancer.Group","system","study","study.abbreviation")],by.x="result.id.outcome",by.y="ID")

		# if(exclude_east_asians)
		# {
		# 	res.m<-res.m[res.m$population!="East Asian",]
		# }
		return(res.m)
} 



format_results2<-function(){
	ids<-res2$id.outcome	
	# ID<-dat.meta$ID
	# ID<-trimws(unlist(strsplit(IDS,split=";")))
	ID<-dat.meta$id.outcome
	ID<-trimws(unlist(strsplit(ID,split=";")))
	ids_keep<-ids[!ids %in% ID]
	res3<-res2[res2$id.outcome %in% ids_keep,]
	res3$nstudies <- 1
	res3$Q.p<-NA 
	# names(res3)[!names(res3) %in% names(dat.meta)]
	res4<-plyr::rbind.fill(res3,dat.meta)
	res4<-res4[order(as.numeric(res4$ncase.outcome),decreasing=TRUE),]
	res4<-res4[!duplicated(res4$outcome2),]
	return(res4)
}


format_results4<-function(dat=NULL,dat.meta=NULL,exclude_finngen=FALSE){
	ids<-dat$id.outcome	
	# ID<-dat.meta$ID
	# ID<-trimws(unlist(strsplit(IDS,split=";")))
	ID<-dat.meta$id.outcome
	if(exclude_finngen){
		if(sum(grep("FinnGen",dat.meta$studies))!=0) stop("FinnGen pdatent in dat.meta")
		dat<-dat[dat$study.abbreviation!="FinnGen",]
	}
	ID<-unique(as.numeric(trimws(unlist(strsplit(ID,split=";")))))
	ids_keep<-unique(ids[!ids %in% ID])
	dat1<-dat[dat$study.abbreviation ==  "GECCO/CORECT/CCFR",]
	dat<-dat[dat$outcome2 !=  "Colorectal cancer",]
	dat.meta<-dat.meta[dat.meta$outcome2 != "Colorectal cancer",] #the only independent CRC analysis comes from GECCO/CORECT/CCFR (which already includes UKB)
	
	dat_bcc<-dat[dat$outcome2 == "Basal cell carcinoma" & dat$exposure == "Alpha-linolenic acid (18:3n3)",] #this result is only available in the 23andMe study and is missing from UKB. ALA therefore gets excluded entirely. A downside is that this result is 23andME only while other PUFas results are from meta analysis of 23andMe and UKB (HNMSC excluded entirely because too many missing SNPs) 
	dat3<-dat[dat$id.outcome %in% ids_keep,]
	dat3<-rbind(dat3,dat_bcc)
	dat3$nstudies <- 1
	dat3$Q.p<-NA 

	# names(dat3)[!names(dat3) %in% names(dat.meta)]
	dat4<-plyr::rbind.fill(dat3,dat.meta)
	dat4<-plyr::rbind.fill(dat1,dat4)
	dat4<-dat4[order(as.numeric(dat4$ncase.outcome),decreasing=TRUE),]
	dat4$study.abbreviation[is.na(dat4$study.abbreviation)]<-dat4$studies[is.na(dat4$study.abbreviation)]
	id<-paste0(dat4$outcome2,dat4$exposure)
	dups<-unique(id[duplicated(id)])	
	dat4[id == dups[1],]
	if(any(duplicated(id))) stop("duplicate analyses pdatent")
	return(dat4)
}


format_results5<-function(dat=NULL,dat.meta=NULL,exclude_finngen=FALSE){
	
	if(exclude_finngen){
		if(sum(grep("FinnGen",dat.meta$studies))!=0) stop("FinnGen present in dat.meta")
		dat<-dat[dat$study.abbreviation!="FinnGen",]
	}
	
	# res3<-dat[dat$id.outcome %in% ids_keep,]
	res3<-dat
	res3$nstudies <- 1
	res3$Q.p<-NA 
	# names(res3)[!names(res3) %in% names(dat.meta)]
	res4<-plyr::rbind.fill(res3,dat.meta)
	res4<-res4[order(as.numeric(res4$ncase.outcome),decreasing=TRUE),]
	res4$study.abbreviation[is.na(res4$study.abbreviation)]<-res4$studies[is.na(res4$study.abbreviation)]

	id<-paste0(res4$outcome2,res4$exposure)
	dups<-unique(id[duplicated(id)])	
	return(res4)
}


format_results3<-function(res=res2)
{
		
		out<-unique(outcome_dat[,c("id.outcome","population","ncase.outcome","ncontrol.outcome","outcome2")])		
		out$ncase.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
			median(out$ncase.outcome[out$id.outcome==x])))
		out$ncontrol.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
			median(out$ncontrol.outcome[out$id.outcome==x])))
		out<-unique(out)

		head(res)
		res.m<-merge(res,out,by="id.outcome")
		res.m<-merge(res.m,meta.tab9[,c("ID","Cancer.Group","system","study","study.abbreviation")],by.x="id.outcome",by.y="ID")

		# if(exclude_east_asians)
		# {
		# 	res.m<-res.m[res.m$population!="East Asian",]
		# }
		return(res.m)
} 

format_results<-function(res=NULL)
{
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")		
	out<-unique(outcome_dat[,c("id.outcome","population","ncase.outcome","ncontrol.outcome","outcome2")])		
	out$ncase.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
		median(out$ncase.outcome[out$id.outcome==x])))
	out$ncontrol.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
		median(out$ncontrol.outcome[out$id.outcome==x])))
	out<-unique(out)

	res.m<-merge(res,out,by="id.outcome")
	res.m<-merge(res.m,meta.tab9[,c("ID","Cancer.Group","system","study","study.abbreviation")],by.x="id.outcome",by.y="ID")

	# if(exclude_east_asians)
	# {
	# 	res.m<-res.m[res.m$population!="East Asian",]
	# }
	return(res.m)
} 

format_outcomes3<-function(dat=NULL,all_cols.keep=FALSE){

	names(dat)[names(dat) =="lnor"]<-"beta.outcome"
	names(dat)[names(dat) =="lnor_se"]<-"se.outcome"
	names(dat)[names(dat) =="p"]<-"pval.outcome"
	names(dat)[names(dat) %in% c("cases","ncase")]<-"ncase.outcome"
	names(dat)[names(dat) %in% c("controls","ncontrol")]<-"ncontrol.outcome"
	names(dat)[names(dat) =="effect_allele"]<-"effect_allele.outcome"
	names(dat)[names(dat) =="other_allele"]<-"other_allele.outcome"
	names(dat)[names(dat) =="eaf"]<-"eaf.outcome"
	names(dat)[names(dat) =="rsid"]<-"SNP"
	names(dat)[names(dat) =="ID"]<-"id.outcome"

	# dat[dat$population == "European; East Asian",c("id.outcome")]
	# exclude meta-analysed outcome data, to avoid MR of mixed European and East Asian populations. Need to MR East Asian and European studies separately prior to meta analysis
	dat<-dat[grep(";",dat$id.outcome,invert=T),]	
	dat$population[dat$id.outcome %in% 993:999]<-"European"
	# dat$study.abbreviation[is.na(dat$study.abbreviation)]<-"survival"
	dat$outcome2<-paste(dat$outcome,dat$id.outcome)
	if(any(names(dat) == "cancer")){
		dat$outcome2<-paste(dat$cancer,dat$id.outcome)
	}

	# Dups<-unique(dat$outcome2[duplicated(dat$outcome2)])
	# Pos<-dat$outcome2 %in% Dups
	# if(any(Pos)) stop("duplicate outcomes when none expected")
	# dat$outcome2[Pos]<-paste(dat$outcome2[Pos],dat$id.outcome[Pos])
	# outcome<-dat$outcome2
	# dat$outcome2<-dat$cancer
	# dat$outcome<-outcome

	# if(all(c("cancer","outcome") %in% names(dat))) warning("two column names for cancer outcome when only one expected")
	# names(dat)[names(dat) %in% c("cancer","outcome")]<-"outcome"
	dat$se.outcome<-as.numeric(dat$se.outcome)
	dat$beta.outcome<-as.numeric(dat$beta.outcome)
	dat$eaf.outcome<-as.numeric(dat$eaf.outcome)
	dat$pval.outcome<-as.numeric(dat$pval.outcome)
	# names(dat)[names(dat) =="cancer"]<-"outcome"
	Cols.keep<-c("SNP","outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","ncase.outcome","ncontrol.outcome","effect_allele.outcome","other_allele.outcome","id.outcome","population","study")
	# "proxy"

	dat<-dat[, Cols.keep]
	# Cols.keep[!Cols.keep %in% names(dat) ]
	# head(dat)
	dat$outcome2<-dat$outcome
	dat$outcome<-paste0(dat$outcome," | ",dat$id.outcome)	
	
	dat$index<-paste0(dat$outcome,dat$SNP)
	Dups<-unique(dat$id.outcome[which(duplicated(dat$index))])
	if(any(duplicated(dat$index))) warning(paste0("duplicate SNPs present, IDs=: ",paste(Dups,collapse=" | ")))

	return(dat)
}


format_outcomes4<-function(dat=NULL,all_cols.keep=FALSE){

	names(dat)[names(dat) =="lnor"]<-"beta.outcome"
	names(dat)[names(dat) =="lnor_se"]<-"se.outcome"
	names(dat)[names(dat) =="p"]<-"pval.outcome"
	names(dat)[names(dat) %in% c("cases","ncase")]<-"ncase.outcome"
	names(dat)[names(dat) %in% c("controls","ncontrol")]<-"ncontrol.outcome"
	names(dat)[names(dat) =="effect_allele"]<-"effect_allele.outcome"
	names(dat)[names(dat) =="other_allele"]<-"other_allele.outcome"
	names(dat)[names(dat) =="eaf"]<-"eaf.outcome"
	names(dat)[names(dat) =="rsid"]<-"SNP"
	names(dat)[names(dat) =="id"]<-"id.outcome"

	# dat[dat$population == "European; East Asian",c("id.outcome")]
	# exclude meta-analysed outcome data, to avoid MR of mixed European and East Asian populations. Need to MR East Asian and European studies separately prior to meta analysis
	dat<-dat[grep(";",dat$id.outcome,invert=T),]	
	dat$population[dat$id.outcome %in% 993:999]<-"European"
	# dat$study.abbreviation[is.na(dat$study.abbreviation)]<-"survival"
	dat$outcome2<-paste(dat$outcome,dat$id.outcome)
	if(any(names(dat) == "cancer")){
		dat$outcome2<-paste(dat$cancer,dat$id.outcome)
	}

	# Dups<-unique(dat$outcome2[duplicated(dat$outcome2)])
	# Pos<-dat$outcome2 %in% Dups
	# if(any(Pos)) stop("duplicate outcomes when none expected")
	# dat$outcome2[Pos]<-paste(dat$outcome2[Pos],dat$id.outcome[Pos])
	# outcome<-dat$outcome2
	# dat$outcome2<-dat$cancer
	# dat$outcome<-outcome

	# if(all(c("cancer","outcome") %in% names(dat))) warning("two column names for cancer outcome when only one expected")
	# names(dat)[names(dat) %in% c("cancer","outcome")]<-"outcome"
	dat$se.outcome<-as.numeric(dat$se.outcome)
	dat$beta.outcome<-as.numeric(dat$beta.outcome)
	dat$eaf.outcome<-as.numeric(dat$eaf.outcome)
	dat$pval.outcome<-as.numeric(dat$pval.outcome)
	# names(dat)[names(dat) =="cancer"]<-"outcome"
	Cols.keep<-c("SNP","outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","ncase.outcome","ncontrol.outcome","effect_allele.outcome","other_allele.outcome","id.outcome","population","proxy","study")

	dat<-dat[, Cols.keep]
	dat$outcome2<-dat$outcome
	dat$outcome<-paste0(dat$outcome," | ",dat$id.outcome)	
	
	dat$index<-paste0(dat$outcome,dat$SNP)
	Dups<-unique(dat$id.outcome[which(duplicated(dat$index))])
	if(any(duplicated(dat$index))) warning(paste0("duplicate SNPs present, IDs=: ",paste(Dups,collapse=" | ")))

	return(dat)
}


meta_analysis2<-function(dat=NULL,IDS=IDS){
	IDS<-IDS[grep(";",IDS)]
	Res_list<-NULL	
	for(i in 1:length(IDS)){
		print(IDS[i])
		ID<-trimws(unlist(strsplit(IDS[i],split=";")))
		if(length(ID)==1) stop("ID=1")
		dat1<-dat[dat$id.outcome %in% ID,]
		dat1<-dat1[dat1$population=="European",]
		if(nrow(dat1)>0)
		{

			
			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
			# if(sum(grep("FinnGen",dat1$study.abbreviation))!=0) stop("fingen")
			b<-dat1$b
			se<-dat1$se
			# p<-temp$p
			w<-1/se^2
			b.fixed<-sum(b*w)/(sum(w))
			se.fixed<-sqrt(sum(w)^-1)
			z<-abs(b.fixed/se.fixed)
			p.fixed<-pnorm(z,lower.tail=F)*2
			nstudies.fixed<-length(b)
			cancer<-unique(dat1$outcome2)
			if(length(cancer)!=1) cancer<-unique(dat1$Cancer.Group)
			if(length(cancer)!=1) 	cancer<-unique(paste(dat1$system,"system cancers"))
			if(length(cancer)!=1) stop("length of cancer not 1")
			if(length(unique(dat1$population))!=1) stop("population not same across studies")
			ids.fixed<-paste(dat1$id.outcome,collapse="; ")
			cases<-sum(dat1$ncase.outcome)
			controls<-sum(dat1$ncontrol.outcome)
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
			# EA<-unique(dat1$Effect.Allele)
			# OA<-unique(dat1$Other.Allele)
			# EAF<-round(sum((dat1$eaf*w))/sum(w),3)
			# if(length(EA)>1) stop("effect allele not consistent across studies")

			# studies<-paste(dat1$study.abbreviation,collapse="; ")
			# Cols<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
			# for(j in 1:length(Cols)){
			# 	Cols[j]
			# 	Cols[j]<-paste(unique(dat1[,Cols[j]]),collapse="; ")
			# }
			dat.matrix<-c(cancer,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,cases,controls,study,Q.p)
			Res<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
			names(Res)<-c("outcome2","b","se","pval","nstudies","id.outcome","ncase.outcome","ncontrol.outcome","study","Q.p")
			# Col.dat<-data.frame(matrix(Cols,ncol=length(Cols),nrow=1),stringsAsFactors=F)
			# names(Col.dat)<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
			# Res.c<-cbind(Res,Col.dat)		
			round(as.numeric(Res$lnor),3)
			round(as.numeric(Res$se),3)
			# round(as.numeric(Res$p),5)
			# options("scipen"=2, "digits"=3)
			# format(as.numeric(Res$p), digits=3)		
			Res_list[[i]]<-Res
		}
	}
	Res1<-do.call(rbind,Res_list)
	Res1$pval<-as.numeric(Res1$pval)
	Res1$population<-"European"
	return(Res1)
	# i4<-which(as.numeric(Res1$Q.p)<0.10 & as.numeric(Res1$Q.p)>0.05)
}

meta_analysis3<-function(dat=NULL,exclude_finngen=FALSE){
	if(exclude_finngen) dat<-dat[dat$study.abbreviation!="FinnGen",] #this script is designed to metaanalyse results across studies for secondary PUFAs. A problem with FinnGen is that most of the SNPs are missing, e.g. 7 for LA (compared to 26-34 for other studies), which could introduce some bias. Would not introduce bias if all SNPs are valid instruments
	Res_list<-NULL	
	outcomes<-unique(dat$outcome2)
	# i<-7
	for(i in 1:length(outcomes)){
		print(outcomes[i])
		dat1<-dat[dat$outcome2 %in% outcomes[i],]
		dat1<-dat1[dat1$population=="European",]
		if(nrow(dat1)>1)
		{

			exposure<-unique(dat$exposure)
			cancer<-outcomes[i]

			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
			# if(sum(grep("FinnGen",dat1$study.abbreviation))!=0) stop("fingen")
			b<-dat1$b
			se<-dat1$se
			# p<-temp$p
			w<-1/se^2
			b.fixed<-sum(b*w)/(sum(w))
			se.fixed<-sqrt(sum(w)^-1)
			z<-abs(b.fixed/se.fixed)
			p.fixed<-pnorm(z,lower.tail=F)*2
			nstudies.fixed<-length(b)
			if(length(unique(dat1$population))!=1) stop("population not same across studies")
			ids.fixed<-paste(dat1$id.outcome,collapse="; ")
			studies<-paste(dat1$study.abbreviation,collapse="; ")
			cases<-sum(dat1$ncase.outcome)
			controls<-sum(dat1$ncontrol.outcome)
			study<-"Overall fixed effect"
			Q<-sum((b.fixed-b)^2*w)
			df.Q<-length(b)-1		
			Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
			nsnp_avg<-mean(dat1$nsnp)
			dat1<-dat1[order(dat1$ncase.outcome,decreasing=TRUE),]
			nsnp_maxn<-dat1$nsnp[1]

			# names(Meta)
			# Meta$Q
			# Meta$pval.Q
			# Meta<-metagen(TE=b,seTE=se,comb.fixed=T,sm="MD")
			# # Meta$TE.fixed
			# Meta$seTE.fixed
			# Meta$pval.fixed
			# Q.p
			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
			# EA<-unique(dat1$Effect.Allele)
			# OA<-unique(dat1$Other.Allele)
			# EAF<-round(sum((dat1$eaf*w))/sum(w),3)
			# if(length(EA)>1) stop("effect allele not consistent across studies")

			# studies<-paste(dat1$study.abbreviation,collapse="; ")
			# Cols<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
			# for(j in 1:length(Cols)){
			# 	Cols[j]
			# 	Cols[j]<-paste(unique(dat1[,Cols[j]]),collapse="; ")
			# }
			dat.matrix<-c(cancer,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,studies,cases,controls,study,Q.p,nsnp_maxn,nsnp_avg)
			Res<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
			names(Res)<-c("outcome2","b","se","pval","nstudies","id.outcome","studies","ncase.outcome","ncontrol.outcome","study","Q_between_study.p","nsnp_maxn","nsnp_avg")
			Res$exposure<-exposure

			# Col.dat<-data.frame(matrix(Cols,ncol=length(Cols),nrow=1),stringsAsFactors=F)
			# names(Col.dat)<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
			# Res.c<-cbind(Res,Col.dat)		
			# Res$lnor<-round(as.numeric(Res$b),3)
			# Res$se<-round(as.numeric(Res$se),3)
			# round(as.numeric(Res$p),5)
			# options("scipen"=2, "digits"=3)
			# format(as.numeric(Res$p), digits=3)		
			Res_list[[i]]<-Res
		}
	}
	Res1<-do.call(rbind,Res_list)
	# Res1$pval<-as.numeric(Res1$pval)
	Res1$population<-"European"
	return(Res1)
	# i4<-which(as.numeric(Res1$Q.p)<0.10 & as.numeric(Res1$Q.p)>0.05)
}


harmonise_dat2<-function(Dat=NULL,ref_dat=NULL,marker="SNP",effect_allele="effect_allele.outcome",other_allele="other_allele.outcome",effect="beta.outcome",EAF="eaf.outcome"){

	Alleles<-unique(paste(Dat[,effect_allele],Dat[,other_allele],sep="")	)
	if(any(Alleles %in% c("GC","CG","TA","AT"))) warning("palindromic SNPs present. This function assumes SNPs are not palindromic")
	
	ref_dat[,effect_allele]<-toupper(ref_dat[,effect_allele])
	ref_dat[,other_allele]<-toupper(ref_dat[,other_allele])
	names(ref_dat)[names(ref_dat) == effect_allele]<-"ref_ea"
	names(ref_dat)[names(ref_dat) == other_allele]<-"ref_oa"
	names(ref_dat)[names(ref_dat) == EAF]<-"ref_eaf"
	Dat[,effect_allele]<-toupper(Dat[,effect_allele])
	Dat[,other_allele]<-toupper(Dat[,other_allele])	
	ref_dat<-ref_dat[,c(marker,"ref_ea","ref_oa","ref_eaf")]
	Dat<-merge(Dat,ref_dat,by=marker)	
	
	ea<-Dat[,effect_allele]
	if(any(ea != Dat$ref_ea & ea != Dat$ref_oa)) stop("strand conflict")
	oa<-Dat[,other_allele]
	eaf<-Dat[,EAF]
	beta<-Dat[,effect]
	if(any(ea != Dat$ref_ea & ea != Dat$ref_oa)) stop("some SNPs are on different strands to the reference")
	if(all(ea == Dat$ref_ea | ea == Dat$ref_oa)) message("all SNPs are on the same strand")

	Pos<-which(ea != Dat$ref_ea) #positions where effect allele is different from effect allele in reference set
	Dat[,effect_allele][Pos]<-oa[Pos]
	Dat[,other_allele][Pos]<-ea[Pos]
	Dat[,EAF][Pos]<-1-eaf[Pos]
	Dat[,effect][Pos]<-beta[Pos]*-1
	Dat<-Dat[,!names(Dat) %in% c("ref_ea",  "ref_oa",  "ref_eaf")]
	# dat1[2:nrow(dat1),]
	all(Dat[,effect_allele] == Dat$ref_ea)	
	return(Dat)
}


meta_analysis_snp<-function(dat=NULL,include_east_asians=TRUE){
	if(!include_east_asians)
	{
		dat<-dat[dat$population == "European",]
	}
	
	Res_list<-NULL	
	snps<-unique(dat$SNP)
	for(i in 1:length(snps)){
		print(snps[i])
		# unique(crc[,c("effect_allele.outcome","other_allele.outcome")])
		# assumes outcome unique
		dat1<-dat[dat$SNP %in% snps[i],]
		# dat1<-dat1[dat1$population=="European",]
		if(length(unique(dat$outcome2))>1) stop("more than one outcome")
		if(nrow(dat1)>1)
		{
			dat2<-harmonise_dat2(Dat=dat1[2:nrow(dat1),],ref_dat=dat1[1,])
			dat1<-rbind(dat1[1,],dat2)
			# exposure<-unique(dat$exposure)
			cancer<-unique(dat1$outcome2)

			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
			# if(sum(grep("FinnGen",dat1$study.abbreviation))!=0) stop("fingen")
			if(length(unique(dat1$effect_allele.outcome))!=1) stop("effect alleles not harmonised")
			if(length(unique(dat1$other_allele.outcome))!=1) stop("effect alleles not harmonised")
			allele<-unique(paste0(dat1$effect_allele.outcome,dat1$other_allele.outcome))
			if(allele %in% c("AT","TA","CG","GC")) stop("palindromic SNP")
			print(dat1[,c("study","SNP","effect_allele.outcome","other_allele.outcome")])
			b<-dat1$b
			se<-dat1$se
			# p<-temp$p
			w<-1/se^2
			b.fixed<-sum(b*w)/(sum(w))
			se.fixed<-sqrt(sum(w)^-1)
			z<-abs(b.fixed/se.fixed)
			p.fixed<-pnorm(z,lower.tail=F)*2
			nstudies.fixed<-length(b)
			ids.fixed<-paste(dat1$id.outcome,collapse="; ")
			studies<-paste(dat1$study,collapse="; ")
			cases<-sum(dat1$ncase.outcome)
			controls<-sum(dat1$ncontrol.outcome)
			study<-"Overall fixed effect"
			population<-paste(dat1$population,collapse="; ")
			Q<-sum((b.fixed-b)^2*w)
			df.Q<-length(b)-1		
			Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
			effect_allele<-unique(dat1$effect_allele.outcome)
			other_allele<-unique(dat1$other_allele.outcome)
			eaf<-sum(dat1$eaf.outcome*dat1$ncase.outcome)/(sum(dat1$ncase.outcome))
			# names(Meta)
			# Meta$Q
			# Meta$pval.Q
			# Meta<-metagen(TE=b,seTE=se,comb.fixed=T,sm="MD")
			# # Meta$TE.fixed
			# Meta$seTE.fixed
			# Meta$pval.fixed
			# Q.p
			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
			# EA<-unique(dat1$Effect.Allele)
			# OA<-unique(dat1$Other.Allele)
			# EAF<-round(sum((dat1$eaf*w))/sum(w),3)
			# if(length(EA)>1) stop("effect allele not consistent across studies")

			# studies<-paste(dat1$study.abbreviation,collapse="; ")
			# Cols<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
			# for(j in 1:length(Cols)){
			# 	Cols[j]
			# 	Cols[j]<-paste(unique(dat1[,Cols[j]]),collapse="; ")
			# }
			dat.matrix<-c(snps[i],effect_allele,other_allele,eaf,cancer,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,studies,cases,controls,study,Q.p,population)
			Res<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
			names(Res)<-c("SNP","effect_allele","other_allele","eaf","outcome2","b","se","pval","nstudies","id.outcome","studies","ncase.outcome","ncontrol.outcome","study","Q.p","population")	
			Res_list[[i]]<-Res
		}
	}
	Res1<-do.call(rbind,Res_list)
	# Res1$pval<-as.numeric(Res1$pval)
	if(include_east_asians){
		Res1$population<-"European & East Asian"
	}
	if(!include_east_asians){
		Res1$population<-"European"
	}
	return(Res1)
	# i4<-which(as.numeric(Res1$Q.p)<0.10 & as.numeric(Res1$Q.p)>0.05)
}


meta_analysis_wr<-function(dat=NULL){
	Res_list<-NULL	
	snps<-unique(dat$SNP)
	for(i in 1:length(snps)){
		print(snps[i])
		# unique(crc[,c("effect_allele.outcome","other_allele.outcome")])
		# assumes outcome unique
		dat1<-dat[dat$SNP %in% snps[i],]
		# dat1<-dat1[dat1$population=="European",]
		if(length(unique(dat$outcome2))>1) stop("more than one outcome")
		if(nrow(dat1)>1)
		{
			# dat2<-harmonise_dat2(Dat=dat1[2:nrow(dat1),],ref_dat=dat1[1,])
			# dat1<-rbind(dat1[1,],dat2)
			# exposure<-unique(dat$exposure)
			cancer<-unique(dat1$outcome2)

			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
			# if(sum(grep("FinnGen",dat1$study.abbreviation))!=0) stop("fingen")
			
			b<-dat1$b
			se<-dat1$se
			# p<-temp$p
			w<-1/se^2
			b.fixed<-sum(b*w)/(sum(w))
			se.fixed<-sqrt(sum(w)^-1)
			z<-abs(b.fixed/se.fixed)
			p.fixed<-pnorm(z,lower.tail=F)*2
			nstudies.fixed<-length(b)
			ids.fixed<-paste(dat1$id.outcome,collapse="; ")
			studies<-paste(dat1$study.abbreviation,collapse="; ")
			cases<-sum(dat1$ncase.outcome)
			controls<-sum(dat1$ncontrol.outcome)
			study<-"Overall fixed effect"
			population<-paste(unique(dat1$population),collapse="; ")
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
			# EA<-unique(dat1$Effect.Allele)
			# OA<-unique(dat1$Other.Allele)
			# EAF<-round(sum((dat1$eaf*w))/sum(w),3)
			# if(length(EA)>1) stop("effect allele not consistent across studies")

			# studies<-paste(dat1$study.abbreviation,collapse="; ")
			# Cols<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
			# for(j in 1:length(Cols)){
			# 	Cols[j]
			# 	Cols[j]<-paste(unique(dat1[,Cols[j]]),collapse="; ")
			# }
			dat.matrix<-c(snps[i],cancer,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,studies,cases,controls,study,Q.p,population)
			Res<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
			names(Res)<-c("SNP","outcome2","b","se","pval","nstudies","id.outcome","studies","ncase.outcome","ncontrol.outcome","study","Q.p","population")	
			Res_list[[i]]<-Res
		}
	}
	Res1<-do.call(rbind,Res_list)
	return(Res1)
	# i4<-which(as.numeric(Res1$Q.p)<0.10 & as.numeric(Res1$Q.p)>0.05)
}


meta_analysis_1snp_moutcomes<-function(dat=NULL,wr=FALSE){
	if(length(unique(dat$effect_allele))!=1) stop("effect alleles not harmominsed") 
	if(length(unique(dat$other_allele)) !=1) stop("effect alleles not harmominsed") 
	# dat2<-harmonise_dat2(Dat=dat[2:nrow(dat1),],ref_dat=dat1[1,])
	

	# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
	# if(sum(grep("FinnGen",dat1$study.abbreviation))!=0) stop("fingen")
	allele<-unique(paste0(dat$effect_allele,dat$other_allele))
	if(allele %in% c("AT","TA","CG","GC")) stop("palindromic SNP")
	print(dat[,c("study","SNP","effect_allele","other_allele")])
	b<-as.numeric(dat$b)
	se<-as.numeric(dat$se)
	if(wr==TRUE)
	{
		b<-as.numeric(dat$b.elovl2)
		se<-as.numeric(dat$se.elovl2)
	}
	# p<-temp$p
	w<-1/se^2
	b.fixed<-sum(b*w)/(sum(w))
	se.fixed<-sqrt(sum(w)^-1)
	z<-abs(b.fixed/se.fixed)
	p.fixed<-pnorm(z,lower.tail=F)*2
	nstudies.fixed<-length(b)
	IDS<-as.numeric(unlist(strsplit(dat$id.outcome,split=";")))
	Pos<-order(IDS)
	IDS<-IDS[Pos]
	ids.fixed<-paste(IDS,collapse="; ")
	studies<-unlist(strsplit(dat$studies,split=";"))
	studies<-gsub(" ","",studies)
	studies<-studies[Pos]
	studies<-paste(dat$studies,collapse="; ")
	cancers<-paste(dat$outcome2,collapse="; ")
	cases<-sum(as.numeric(dat$ncase.outcome))
	controls<-sum(as.numeric(dat$ncontrol.outcome))
	study<-"Overall effect"
	population<-paste(dat$population,collapse="; ")
	Q<-sum((b.fixed-b)^2*w)
	df.Q<-length(b)-1		
	Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
	effect_allele<-unique(dat$effect_allele)
	other_allele<-unique(dat$other_allele)
	eaf<-sum(as.numeric(dat$eaf)*as.numeric(dat$ncase.outcome))/(sum(as.numeric(dat$ncase.outcome)))
	snp<-unique(dat$SNP)
	dat.matrix<-c(snp,effect_allele,other_allele,eaf,cancers,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,studies,cases,controls,study,Q.p,population,study)
	Res<-data.frame(matrix(dat.matrix,nrow=length(snp),ncol=length(dat.matrix)),stringsAsFactors=F)
	names(Res)<-c("SNP","effect_allele","other_allele","eaf","cancers","b","se","pval","nstudies","id.outcome","studies","ncase.outcome","ncontrol.outcome","study","Q.p","population","study")	
	return(Res)
}

meta_analysis_1wr_moutcomes<-function(dat=NULL){
	b<-as.numeric(dat$b)
	se<-as.numeric(dat$se)	
	w<-1/se^2
	b.fixed<-sum(b*w)/(sum(w))
	se.fixed<-sqrt(sum(w)^-1)
	z<-abs(b.fixed/se.fixed)
	p.fixed<-pnorm(z,lower.tail=F)*2
	nstudies.fixed<-length(b)
	IDS<-as.numeric(unlist(strsplit(dat$id.outcome,split=";")))
	Pos<-order(IDS)
	IDS<-IDS[Pos]
	ids.fixed<-paste(IDS,collapse="; ")
	studies<-unlist(strsplit(dat$studies,split=";"))
	studies<-gsub(" ","",studies)
	studies<-studies[Pos]
	studies<-paste(dat$studies,collapse="; ")
	cancers<-paste(unique(dat$outcome2),collapse="; ")
	cases<-sum(as.numeric(dat$ncase.outcome))
	controls<-sum(as.numeric(dat$ncontrol.outcome))
	study<-"Overall effect"
	population<-paste(unique(dat$population),collapse="; ")
	Q<-sum((b.fixed-b)^2*w)
	df.Q<-length(b)-1		
	Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
	snp<-unique(dat$SNP)
	if(length(snp)!=1) stop("more than one SNP")
	dat.matrix<-c(snp,cancers,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,studies,cases,controls,study,Q.p,population,study)
	Res<-data.frame(matrix(dat.matrix,nrow=length(snp),ncol=length(dat.matrix)),stringsAsFactors=F)
	names(Res)<-c("SNP","cancers","b","se","pval","nstudies","id.outcome","studies","ncase.outcome","ncontrol.outcome","study","Q.p","population","study")	
	return(Res)
}


format_meta_analysis_rs3734398<-function(){
	meta_dat<-meta_dat2[meta_dat2$SNP=="rs3734398",] #seems like rs3734398 not present in FinnGen

	# meta_dat2[meta_dat2$outcome2 == "Esophageal squamous cell carcinoma",]
	meta_dat1<-meta_dat[meta_dat$study=="UKB" & meta_dat$outcome2 == "Overall cancer" ,]
	# meta_dat3<-meta_dat2[meta_dat2$study!="UKB",]
	meta_dat<-meta_dat[meta_dat$outcome2 != "Basal cell carcinoma" ,] #weak evidence that ara/D5D causes bcc 	
	meta_dat<-meta_dat[meta_dat$study == "Overall fixed effect",]
	meta_dat<-rbind(meta_dat1,meta_dat)
	# meta_dat<-ara_increasing_effect_allele(Plot_dat=meta_dat)		
	meta_dat<-elovl2_function(Plot_dat=meta_dat)
	meta_dat$study2<-meta_dat$study
	meta_dat$study<-meta_dat$outcome2
	meta_dat9<-meta_analysis_1snp_moutcomes(dat=meta_dat)
	# names(meta_dat9)[!names(meta_dat9) %in% names(meta_dat)]		
	names(meta_dat)[!names(meta_dat) %in% names(meta_dat9)]
	meta<-plyr::rbind.fill(meta_dat9,meta_dat)
	meta$b<-as.numeric(meta$b)
	meta$se<-as.numeric(meta$se)
	meta$shape<-15
	meta$shape[meta$study == "Overall effect"]<-23		
	meta$weight<-1/meta$se/20
	meta$outcome_name<-paste0(meta$study,"\nN. cases=",meta$ncase.outcome)
	meta$outcome_name<-gsub("Esophageal squamous cell carcinoma","Esophageal\nsquamous\ncell carcinoma",meta$outcome_name)
	meta<-meta[order(as.numeric(meta$ncase.outcome),decreasing=TRUE),]	
	return(meta)
}

# load("~/fatty-acids/mr/results/results_rs3734398_meta_analysis_europeans.Rdata")
format_meta_analysis_rs3734398_wald_ratio<-function(){
	meta_dat<-meta_dat[meta_dat$SNP=="rs3734398",] #seems like rs3734398 not present in FinnGen
	meta_dat<-meta_dat[meta_dat$population=="European",]

	# meta_dat2[meta_dat2$outcome2 == "Esophageal squamous cell carcinoma",]
	meta_dat1<-meta_dat[meta_dat$study=="UKB" & meta_dat$outcome2 == "Overall cancer" ,]
	# meta_dat3<-meta_dat2[meta_dat2$study!="UKB",]
	meta_dat<-meta_dat[meta_dat$outcome2 != "Basal cell carcinoma" ,] #weak evidence that ara/D5D causes bcc 	
	meta_dat<-meta_dat[meta_dat$study != "UKB",]
	meta_dat<-rbind(meta_dat1,meta_dat)
	# meta_dat<-ara_increasing_effect_allele(Plot_dat=meta_dat)		
	meta_dat<-elovl2_function(Plot_dat=meta_dat)
	meta_dat$se.elovl2
	meta_dat$study2<-meta_dat$study
	meta_dat$study<-meta_dat$outcome2
	meta_dat9<-meta_analysis_1snp_moutcomes(dat=meta_dat,wr=TRUE)
	names(meta_dat9)[names(meta_dat9) == "b"]<-"b.elovl2"
	names(meta_dat9)[names(meta_dat9) == "se"]<-"se.elovl2"
	meta_dat9$b.elovl2<-as.numeric(meta_dat9$b.elovl2)
	meta_dat9$se.elovl2<-as.numeric(meta_dat9$se.elovl2)
	# names(meta_dat9)[!names(meta_dat9) %in% names(meta_dat)]		
	names(meta_dat)[!names(meta_dat) %in% names(meta_dat9)]
	meta<-plyr::rbind.fill(meta_dat9,meta_dat)
	meta$shape<-15
	meta$shape[meta$study == "Overall effect"]<-23		
	meta$weight<-1/meta$se.elovl2/20
	meta$outcome_name<-paste0(meta$study,"\nN. cases=",meta$ncase.outcome)
	meta$outcome_name<-gsub("Esophageal squamous cell carcinoma","Esophageal\nsquamous\ncell carcinoma",meta$outcome_name)
	meta<-meta[order(as.numeric(meta$ncase.outcome),decreasing=TRUE),]	
	return(meta)
}


format_meta_analysis_rs3734398_wald_ratio2<-function(meta_dat=NULL){
	
	
	# meta_dat<-meta_dat[meta_dat$SNP=="rs3734398",] #seems like rs3734398 not present in FinnGen
	meta_dat<-meta_dat[meta_dat$population=="European",]

	# meta_dat2[meta_dat2$outcome2 == "Esophageal squamous cell carcinoma",]
	meta_dat1<-meta_dat[which(meta_dat$study=="UK Biobank" & meta_dat$outcome2 == "Overall cancer") ,]


	# meta_dat3<-meta_dat2[meta_dat2$study!="UKB",]
	meta_dat<-meta_dat[which(meta_dat$outcome2 != "Basal cell carcinoma") ,] #weak evidence that ara/D5D causes bcc 	
	meta_dat<-meta_dat[which(meta_dat$study.abbreviation!="FinnGen") ,]
	meta_dat$outcome2
	meta_dat<-meta_dat[which(meta_dat$study.abbreviation != "UKB"),]
	meta_dat<-rbind(meta_dat1,meta_dat)
	ID<-paste(meta_dat$outcome2,meta_dat$SNP)
	meta_dat<-meta_dat[!duplicated(ID),]

	meta_dat$study2<-meta_dat$study
	meta_dat$study<-meta_dat$outcome2
	
	snps<-unique(meta_dat$SNP)
	meta_list<-NULL
	for(i in 1:length(snps)){
		test_dat<-meta_dat[meta_dat$SNP==snps[i],]
		meta_list[[i]]<-meta_analysis_1wr_moutcomes(dat=test_dat)
	}
	meta_dat9<-do.call(rbind,meta_list)	
	# names(meta_dat9)[!names(meta_dat9) %in% names(meta_dat)]		
	names(meta_dat)[!names(meta_dat) %in% names(meta_dat9)]
	meta_dat9$outcome2<-"Overall effect"
	names(meta_dat9)[names(meta_dat9) == "pval"]<-"p"
	# names(meta_dat)[names(meta_dat) == "pval"]

	meta<-plyr::rbind.fill(meta_dat9,meta_dat)
	meta$b<-as.numeric(meta$b)
	meta$se<-as.numeric(meta$se)
	meta$p<-as.numeric(meta$p)
	meta$shape<-15
	meta$shape[meta$study == "Overall effect"]<-23		
	meta$weight<-1/meta$se/5
	meta$outcome_name<-paste0(meta$study,"\nN. cases=",meta$ncase.outcome)
	meta<-meta[order(as.numeric(meta$ncase.outcome),decreasing=TRUE),]	
	meta$or<-round(exp(meta$b),3)
	meta$uci<-round(exp(meta$b+1.96*meta$se),3)
	meta$lci<-round(exp(meta$b-1.96*meta$se),3)
	return(meta)
}





elovl2_function<-function(Plot_dat=NULL){
	Plot_dat$b<-as.numeric(Plot_dat$b)
	Plot_dat$se<-as.numeric(Plot_dat$se)
	Plot_dat$eaf<-as.numeric(Plot_dat$eaf)
	b.elovl2<-0.23
	se.elvol2<- 0.015
	ea.elvol2<-"T"
	oa.elvol2<-"C"
	eaf.elovl2<-0.57
	ea<-Plot_dat$effect_allele
	oa<-Plot_dat$other_allele
	Pos<-ea!=ea.elvol2	
	b<-Plot_dat$b
	se<-Plot_dat$se
	eaf<-Plot_dat$eaf
	Plot_dat$b[Pos]<-b[Pos]*-1
	Plot_dat$eaf[Pos]<-1-eaf[Pos]
	Plot_dat$effect_allele[Pos]<-oa[Pos]
	Plot_dat$other_allele[Pos]<-ea[Pos]
	Plot_dat$b.elovl2<-Plot_dat$b/b.elovl2
	Plot_dat$se.elovl2<-Plot_dat$se/b.elovl2
	return(Plot_dat)
}


ara_increasing_effect_allele<-function(Plot_dat=NULL){
	# Plot_dat[,c("SNP","b","effect_allele")]
	beta_direction<-data.frame(matrix(c("rs4985155", "0.09248261","a","g","rs3734398", 0.05844874,"c","t","rs174546", 0.88347054 ,"c","t"),nrow=3,ncol=4,byrow=TRUE))
	names(beta_direction)<-c("SNP","b","effect_allele","other_allele")
	beta_direction$effect_allele<-toupper(beta_direction$effect_allele)
	beta_direction$other_allele<-toupper(beta_direction$other_allele)
	names(beta_direction)<-c("SNP","b","effect_allele","other_allele")
	Plot_dat<-merge(Plot_dat,	beta_direction,by="SNP")
	Pos<-Plot_dat$effect_allele.x != Plot_dat$effect_allele.y
	Plot_dat$b.x[Pos]<-as.numeric(Plot_dat$b.x[Pos])*-1
	Plot_dat$eaf[Pos]<-1-as.numeric(Plot_dat$eaf[Pos])
	ea<-Plot_dat$effect_allele.x[Pos]
	oa<-Plot_dat$other_allele.x[Pos]
	Plot_dat$effect_allele.x[Pos]<-oa
	Plot_dat$other_allele.x[Pos]<-ea
	return(Plot_dat)
}



# assume overlap only amongst controls only (e.g. applicable to biobank studies)
correlation_between_cancers1<-function(test=NULL,ncontrols_ukb_gecco=21407,ukb_controls=372016,finngen_controls=96357,bj_controls=211217,interlymph_controls=8107,biobanks=c("UKB","FinnGen","BJ","InterLymph"),male_cancers=c("Prostate cancer","Male genital cancer"),female_cancers=c("Breast cancer","Cervical cancer","Ovarian cancer","Endometrial cancer", "Female genital cancer")){
	
	# consortia or biobanks that overlap between cancer sites 
	# uk biobank
	# finngen
	# biobank japan
	# InterLymph 

	# consortia with overlap amongst cancer subtypes within sites (cases and controls) 
	# OCAC 
	# InterLymph 

	# consortia with overlap between cancer subtypes and overall cancer (overlap in cases and controls) 
	# OCAC 
	# BCAC 
	# PRACTICAL
	# ILCCO 
	# GECCO
	# INHANCE


	# consortia that overlap between cancer subtypes within site (controls) 
	# BCAC
	# OCAC 
	# PRACTICAL
	# ILCCO 
	# GECCO
	# 23NMSC
	# INHANCE
	# HNMSC (nonmelanoma skin cancer subtypes in nhs and phs)


	# consortia or biobanks with overlap between cancer sites (overlap in controls not cases)
	# uk biobank
	# finngen
	# biobank japan
	# InterLymph 
	# GECCO (ukb)
	

	 #ncontrols in biobanks in order of biobanks argument

	
	# nk1 = number of cases in kth study 
	# nk0 = number of controls in kth study 
	# nk = total number of subjects in kth study
	# nl1 = number of cases in lth study 
	# nl0 = number of controls in lth study 
	# nl = total number of subjects in the lth study 
	# nkl0 = number of controls that overlap between kth and lth studies 
	# nkl1 = number of cases that overlap between the kth and lth studies.
	
	print(test$study.abbreviation[1])
	studies1<-trimws(unlist(strsplit(test$study.abbreviation[1],split=";")))
	studies2<-trimws(unlist(strsplit(test$study.abbreviation[2],split=";")))


	if(any(studies1 %in% "GECCO/CORECT/CCFR")) {
		studies1<-unique(c(studies1,"UKB"))
		ukb_controls<-ncontrols_ukb_gecco
	}
	
	overlapping_controls<-c(ukb_controls,finngen_controls,bj_controls,interlymph_controls)

	if(any(studies2 %in% "GECCO/CORECT/CCFR")) {
		studies2<-unique(c(studies2,"UKB"))
	}

	Pos1<-which(biobanks %in% studies1 )
	Pos2<-which(biobanks %in% studies2 )
	Pos<-c(Pos1,Pos2)
	Pos<-Pos[duplicated(Pos)]
	nkl0_biobanks<-sum(overlapping_controls[Pos])
	# meta.tab$controls[meta.tab$study.abbreviation=="InterLymph"]
	
	Cancer1<-test$cancer[1]
	Cancer2<-test$cancer[2]
	if(any(male_cancers %in% Cancer1 | male_cancers %in% Cancer2) &
	any(female_cancers %in% Cancer1 | female_cancers %in% Cancer2)) nkl0_biobanks<-0 
	#if both conditions true means that one cancer is male cancer and other cancer is female cancer. therefore there can't be any overlap 

	# if(any(male_cancers %in% Cancer1 | male_cancers %in% Cancer2) |
	# 	any(female_cancers %in% Cancer1 | female_cancers %in% Cancer2)) nkl0_biobanks<-nkl0_biobanks/2 #if one cancer has been estimated without stratification on sex and the second has been estimated in either males or females, then the number of overlapping controls will be overestimated. Assuming a balanced sex ratio, we can divide the overlap estimate by half. This will however likely overestimate overlap with female cancers and underestimate overlap with male cancers (because most cases tend to be male). 

	nk1<-as.numeric(test$cases[1])
	nk0<-as.numeric(test$controls[1])
	nk=nk1+nk0
	nl1<-as.numeric(test$cases[2])
	nl0<-as.numeric(test$controls[2])
	# test[2,]
	nl<-nl1+nl0
	# nkl0<-min(nk0,nl0)
	nkl0<-nkl0_biobanks
	nkl1<-0 
	
	corr<-(nkl0*sqrt((nk1*nl1)/(nk0*nl0))+nkl1*sqrt((nk0*nl0)/(nk1*nl1)))/sqrt(nk*nl)

	# jack method
	# k=min(test$controls[1],test$controls[2])
	# rho=1
	# # k = number of overlapping participants
	# # N1 size of study 1
	# # N2 size of study 2	
	# correlation <- (k^2/(N1*N2))*rho
	return(corr)
}


# assume overlap amongst cases and amongst controls
correlation_between_cancers2<-function(){
	# nk1 = number of cases in kth study 
	# nk0 = number of controls in kth study 
	# nk = total number of subjects in kth study
	# nl1 = number of cases in lth study 
	# nl0 = number of controls in lth study 
	# nl = total number of subjects in the lth study 
	# nkl0 = number of controls that overlap between kth and lth studies 
	# nkl1 = number of cases that overlap between the kth and lth studies.
	nk1<-as.numeric(test$cases[1])
	nk0<-as.numeric(test$controls[1])
	nk=nk1+nk0
	nl1<-as.numeric(test$cases[2])
	nl0<-as.numeric(test$controls[2])	
	nl<-nl1+nl0
	nkl0<-min(nk0,nl0)
	nkl1<-min(nk1,nl1) 	
	corr<-(nkl0*sqrt((nk1*nl1)/(nk0*nl0))+nkl1*sqrt((nk0*nl0)/(nk1*nl1)))/sqrt(nk*nl)

	# jack method
	# k=min(test$controls[1],test$controls[2])
	# rho=1
	# # k = number of overlapping participants
	# # N1 size of study 1
	# # N2 size of study 2	
	# correlation <- (k^2/(N1*N2))*rho
	return(corr)
}


decoupling <- function(s, C, C.inv=NULL) 
{
	i <- !is.na(s)
	if (is.null(C.inv)) {
	Omega.inv <- solve(diag(s[i]) %*% C[i,i] %*% diag(s[i]))
	} else {
	Omega.inv <- diag(1/s[i]) %*% C.inv[i,i] %*% diag(1/s[i])
	}
	s.new <- sqrt(1/rowSums(Omega.inv))
	s[i] <- s.new
	s
}

# do i need to test this on an example in the paper?

# decouple_se<-function(lci=NULL,uci=NULL){
# 	s<-unlist(lapply(1:length(lci),FUN=function(i) 
# 		get_se(lci=lci[i],uci=uci[i])
# 		))
# 	C<-make_cov_matric()

# 	s2<-decoupling(s=s,C=C)
# 	return(s2)
# }

# mr_res1$outcome[grep("breast",mr_res1$outcome)]

# 1. Is it a unique cancer site? does it overlap with other cancer sites? e.g. high risk breast cancer is same site as breast cancer. subtypes at same site
# 2. Are the cases independent of cases for other cancer types? (e.g. BCAC cases across overall breast cancer and breast cancer subtypes)
# restrict_to_european_studies=TRUE
format_metareg_v2<-function(restrict_to_european_studies=FALSE,vogelstein=FALSE,seer=FALSE){
	# load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata") #in version 3 we included one additinoal outcome study that was missing  "5; 49; 86" and also restrict the analyses to be exclusivel AA:DGLA in Eur and GLA:LA in EAS. In the meta-analyses of Eur+Eas, was using GLA:LA in Eur. but this leads to overestimates of the effect compared to Eur only analyses because effect on GLA:LA is much weaker than for AA:DGLA

	# this section is redundant because in mr_results_rep_v3 the MR analysis is restrict to AA:DGLA/Europeans and GLA:LA/East Asians. When combining MR results across Eur and Eas, we model effect in Eur using AA:DGLA and effect in Eas using GLA:LA
	# mr_res2<-mr_res1[mr_res1$population == "European",]
	# mr_res2<-mr_res2[mr_res2$exposure == "AA:DGLA",]
	# mr_res3<-mr_res1[mr_res1$population != "European",]
	# mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	# mr_res1<-rbind(mr_res2,mr_res3)
	load("~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
	# mr_res1$outcome[grep("49", mr_res1$outcome)]
	mr_res1<-fix_ids(Dat=mr_res1,id_col="id.outcome")
	# mr_res1$cancer
	mr_res1$Cancer.Group[mr_res1$Cancer.Group %in% 	c("Kidney cancer in females","Kidney cancer in males")]<-"Kidney cancer"

	disc.tab9<-fix_ids(Dat=disc.tab9,id_col="ID")
	meta1<-disc.tab9[,c("ID2","study.abbreviation")]
	meta2<-	meta.tab9[c("ID","study.abbreviation")]
	names(meta1)[names(meta1) == "ID2"]<-"ID"
	meta_dat<-unique(rbind(meta1,meta2))
		
	mr_res1<-merge(mr_res1,meta_dat,by.x="ID2",by.y="ID")
	# mr_res$ID2[is.na(mr_res1$study.abbreviation.y)]
	unique(mr_res1[,c("exposure","population")])

	mr_res1$cases<-as.numeric(mr_res1$cases)	
	mr_res1$cancer[is.na(mr_res1$cancer)]<-mr_res1$outcome[is.na(mr_res1$cancer)]
	mr_res1<-mr_res1[!mr_res1$id.outcome %in% c("993","994","995","996","997","998","999","1499"),]

	mr_res1<-mr_res1[!mr_res1$id.outcome %in% "5; 49; 86",] #decided to keep individual studies separate because were actually studies of different blood cancers. # "5; 49; 86" were combined for purpose of discovery analysis because individually these studies have fewer than 1000 cases
	#ID=86: InterLymph, Marginal zone lymphoma         
	# ID=5: BC-NHL      B cell non-hodgkin lymphoma         
	#ID=49: FinnGen Non-hodgkin lymphoma unspecified         
	# Marginal zone lymphoma is classified as inflammatory cancer 
	
	# Temp<-unique(mr_res1[,c("cancer","site","Cancer.Group","system","cases","study.abbreviation.y")])
	# Temp<-Temp[order(Temp$cases,decreasing=TRUE),]
	# Temp<-Temp[!duplicated(Temp$cancer),]
	# write.table(Temp[order(Temp$site,Temp$Cancer.Group),],"~/temp.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
	


	mr_res1<-mr_res1[!mr_res1$cancer %in% c(
		"Cancer (all cause)","Cancer (excluding non-melanoma skin cancer)","Cancer of digestive organs",#cannot be clearly assigned to smoking/non smoking or inflammatory/non-inflammatory groups; also, cases will overlap with other digestive system cancers
		"Respiratory and intrathoracic cancer",#cases  overlap with lung cancer 
		"Lymphoid leukaemia", #analysed in UKB and nested within leukaemia analysis in UKB
		"Brain cancer", #analysysed in UKB and nested within Central nervous system and eye cancer
		"Female genital cancer", #cases overlap with Endometrial cancer analysis because of FinnGen
		"Malignant non-melanoma skin cancer","Malignant skin cancer", #retain melanoma, squamous cell and basal cell instead (no overlapping cases) and exclude the other skin cancers.  
		"Urinary tract cancer" #overlap in cases with bladder cancer
		),]


	# comment on Lymphoid leukaemias:
	# it's unlikely there are overlapping cases between the lymphoid leukaemia studies, either because they are distinct subtypes or because the studies were independent
		# Acute lymphoblastic leukaemia (BC-ALL; C-ALL; SJ-COG)
		# Chronic lymphocytic leukaemia (InterLymph)
		# Lymphoid leukaemia (UKB)

	# comment on non-hodgkin lymphomas:
	# it's unlikely there are overlapping cases between the non-hodkin lymphoma studies, either because they are distinct subtypes or because the studies were independent
	# the following non-hodgkin lymphomas are not likely to have overlapping cases: Diffuse large b cell,Follicular lymphoma, zone lymphoma. 
	# B cell non-hodgkin lymphoma likely includes Diffuse large b cell,Follicular lymphoma, Marginal zone lymphoma but comes from an independent study.
	# Non-hodgkin lymphoma unspecified (FinnGen) does not likely include cases from Follicular lymphoma (FinnGen). I interpret unspecified to mean exclusion of Follicular lymphoma 

	# Diffuse large b cell lymphoma	Blood	(InterLymph)
	# Follicular lymphoma	(InterLymph; FinnGen)
	# Marginal zone lymphoma	InterLymph
	# B cell non-hodgkin lymphoma	BC-NHL
	# Non-hodgkin lymphoma unspecified	FinnGen

	#comment on Myeloid leukemias
	# Overlap in cases and controls unlikely because studies are independent
	# Myeloid leukaemia (UKB)
	# Chronic myeloid leukaemia KCML

	# Retain the single largest cancer analysis within each of the following cancer groups
	cancer.group<-c(
		"Breast cancer",
		"Colorectal cancer",
		"Gastric adenocarcinoma",
		"Kidney cancer",
		"Liver & biliary tract cancer",
		"Lung cancer",
		"Mouth & throat cancer",
		"Ovarian cancer",
		"Prostate cancer",
		"Liver & biliary tract cancer")
		
	# do not drop this cancer from the Mouth & throat cancer group
	cancer.group_excl<-c(
		"Nasopharyngeal carcinoma")#nested within Mouth and throat cancer group but is independent to INHANCE 
		
	smoking_cancers<-c(
		"Bladder cancer",
		"Cervical cancer",
		"Colorectal cancer", 
		"Esophageal adenocarcinoma","Esophageal squamous cell carcinoma", 
		"Kidney cancer", 
		"Myeloid leukaemia", #acute myeloid leaukemia is smoking related but FAMRC has myeloid leukemia. it is expected that ~2/3 of Myeloid leukemia is AML (see note below)
		"Liver cancer", 
		"Lung cancer", 
		"Oral cavity and pharyngeal cancer",
		"Nasopharyngeal carcinoma",
		"Pancreatic cancer",
		"Gastric adenocarcinoma"
		)
		

	# "Myeloid leukaemia"
	# according to the American cancer society 2021:
	# leuk<-61090
	# aml<-20240 #new cases of aml in 2021
	#   #new cases of leukemia 
	# cml<-9110
	# ml<-cml+aml
	# aml/ml*100 #69% of ml cases are aml
	# use ml as proxy for aml in smoking meta regression analysis; or exclude from analysis

	# https://www.cancer.org/cancer/acute-myeloid-leukemia/about/key-statistics.html
	# The American Cancer Society’s estimates for leukemia in the United States for 2021 are:

	# About 61,090 new cases of leukemia (all kinds) and 23,660 deaths from leukemia (all kinds)
	# About 20,240 new cases of acute myeloid leukemia (AML). Most will be in adults.
	# About 11,400 deaths from AML. Almost all will be in adults.

	# https://www.cancer.org/cancer/chronic-myeloid-leukemia/about/statistics.html
	# The American Cancer Society's estimates for chronic myeloid leukemia (CML) in the United States for 2021 are:

	# About 9,110 new cases will be diagnosed with CML (5,150 in men and 3,960 in women).
	# About 1,220 people will die of CML (680 men and 540 women).

	Infl_cancers<-c(
		"Pleural mesothelioma",
		"Lung cancer",
		"Bladder cancer",
		"Oral cavity and pharyngeal cancer",
		"Colorectal cancer",
		"Pancreatic cancer",
		"Esophageal adenocarcinoma","Esophageal squamous cell carcinoma",
		"Marginal zone lymphoma",
		"Melanoma") 

		#exclude "Uveal melanoma" 

	# Chronic inflammatory cancers in ref41	(Coussens doi:10.1038/nature01322) and the matched cancer in FAMRC. 
	Infl_cancers_coussens<-c("mesothelioma", #"Pleural mesothelioma"
		"lung cancer", #Lung cancer
		"bladder cancer", #Bladder cancer
		"oral squamous cell carcinoma", #"Oral cavity and pharyngeal cancer
		"colorectal cancer", #Colorectal cancer
		"vulvar squamous cell carcinoma", #NA
		"pancreatic cancer", #Pancreatic cancer
		"esophageal cancer", #"Esophageal adenocarcinoma" & "Esophageal squamous cell carcinoma"
		"salivary gland carcinoma", #NA
		"mucosa-associated lymphoid tissue (MALT) lymphoma", #Marginal zone lymphoma. Extranodal MZL or Mucosa-Associated Lymphoid Tissue (MALT) is the most common form of MZL, accounting for about two-thirds of all MZL cases per year.https://lymphoma.org/aboutlymphoma/nhl/mzl/; https://lymphoma-action.org.uk/types-lymphoma-non-hodgkin-lymphoma/malt-lymphoma-gastric-and-non-gastric
		"melanoma" #"Melanoma"
		) 

	#Marginal zone lymphoma = MALT lymphoma Marginal zone lymphoma (MZL) is a group of indolent (slow growing) NHL B-cell lymphomas, which account for approximately eight percent of all NHL cases. The average age at diagnosis is 60 years, and it is slightly more common in women than in men. 

	# with name of cancer in coussens et al after the hash
	agents<-c("Biliary tract cancer",# Cholangiosarcoma, 
			"Colon cancer",# colon carcnoma, 
			"Gastric adenocarcinoma",# Gastric adenocarcinoma, 
			 "Liver cancer",# Hepatocellular carcinoma
			 "B cell non-hodgkin lymphoma",   # B-cell non-Hodgkin’s lymphoma,     
			 "Non-hodgkin lymphoma unspecified",# Non-Hodgkin’s lymphoma, 
			 "Ovarian cancer",     # Ovarian carcinoma, 
			 "Cervical cancer",   # cervical/anal carcinoma   
			 "Bladder cancer",# Bladder, 
			  "Rectal cancer") # rectal carcinoma, 

	mr_res1$infl<-0
	mr_res1$infl[mr_res1$cancer %in% Infl_cancers]<-1
	mr_res1$infl_agent<-mr_res1$infl
	mr_res1$infl_agent[mr_res1$cancer %in% agents]<-1
	mr_res1$smoking1<-0 #all cancers with no evidence or unclear evidence
	mr_res1$smoking1[mr_res1$cancer %in% smoking_cancers]<-1
	
	mr_res1$system2<-mr_res1$system	
	mr_res1$system2[mr_res1$Cancer.Group =="Mouth & throat cancer"]<-"Digestive"
	mr_res1$system_num<-mr_res1$system
	mr_res1$system_num[mr_res1$system=="Digestive"]<-"1"
	mr_res1$system_num[mr_res1$system!="Digestive"]<-"0"
	mr_res1$system_num<-as.numeric(mr_res1$system_num)
	mr_res1$system2_num<-mr_res1$system2
	mr_res1$system2_num[mr_res1$system2=="Digestive"]<-"1"
	mr_res1$system2_num[mr_res1$system2!="Digestive"]<-"0"
	mr_res1$system2_num<-as.numeric(mr_res1$system2_num)

	mr_res1$external<-mr_res1$system
	mr_res1$external[mr_res1$system %in% c("Integumentary","Respiratory","Digestive")]<-"external"
	mr_res1$external[!mr_res1$system %in% c("Integumentary","Respiratory","Digestive")]<-"non-external"
	mr_res1$external_num<-mr_res1$external
	mr_res1$external_num[mr_res1$external=="external"]<-"1"
	mr_res1$external_num[mr_res1$external=="non-external"]<-"0"
	mr_res1$external_num<-as.numeric(mr_res1$external_num)
	# table(mr_res1$system,mr_res1$external)
	# mr_res1$smoking2<-0 #all cancers with unclear evidence
	# mr_res1$smoking2[mr_res1$Cancer.Group %in% cancer.group]<-1
	# mr_res1$smoking2[mr_res1$cancer == "Cervical cancer"]<-1
	# mr_res1$smoking2[which(mr_res1$site == "Brain")]<-2
	# mr_res1$smoking2[which(mr_res1$Cancer.Group %in% c("Breast cancer","Prostate cancer"))]<-2 #no evidence

	# mr_res1$smoking3<-0 #all cancers with unclear evidence
	# mr_res1$smoking3[mr_res1$Cancer.Group %in% cancer.group]<-1
	# mr_res1$smoking3[mr_res1$cancer == "Endometrial cancer"]<-1 #1 = cancers with causal association but endometrial cancer shows protective effect 
	# mr_res1$smoking3[which(mr_res1$site == "Brain")]<-2
	# mr_res1$smoking3[which(mr_res1$Cancer.Group %in% c("Breast cancer","Prostate cancer"))]<-2

	if(restrict_to_european_studies){	
		mr_res1<-mr_res1[mr_res1$population == "European",]
	}

	if(seer | vogelstein){

	
		# add in cancer level characteristics from SEERS and Tomasetti and Vogelstein
				# Add in cancer level characteristics from SEER and Tomasetti and Vogelstein	
		url<-"https://docs.google.com/spreadsheets/d/1sTliqObb08y3rdfB-vsYecjveud1pPk4FTWSYcpx3Xk/edit?usp=drive_open&ouid=103571138679663288585"
		meta.tab<-data.frame(gsheet2tbl(url))
		cols_keep_meta.tab1<-c("cancer","incidence","survival_time","name_in_seer","median.age.diagnosis","Cumulativ.e.number.of.divisions.of.all.stem.cells.per.lifetime..lscd.","name_in_Tomasetti_Vogelstein","ID")
		meta.tab<-meta.tab[,cols_keep_meta.tab1]
		# meta.tab<-meta.tab[meta.tab$population=="European",]
		meta.tab<-meta.tab[!is.na(meta.tab$ID),]
		meta.tab<-meta.tab[!meta.tab$ID %in% c("993","994","995","996","997","998","999","1499"),]
		meta.tab$merge_index<-tolower(meta.tab$cancer)
		mr_res1$merge_index<-tolower(mr_res1$cancer)
		meta.tab$merge_index[meta.tab$merge_index=="non-melanoma skin cancer"]<-"malignant non-melanoma skin cancer"
		meta.tab$merge_index[meta.tab$merge_index=="small cell lung cancer"  ]<-"small cell lung carcinoma" 
		meta.tab<-meta.tab[!duplicated(meta.tab$cancer),]

		
		mr_res1<-merge(mr_res1,meta.tab,by="merge_index",all.x=FALSE)
		names(mr_res1)[names(mr_res1)=="cancer.x"]<-"cancer"
		names(mr_res1)[names(mr_res1)=="Cumulativ.e.number.of.divisions.of.all.stem.cells.per.lifetime..lscd."]<-"stem_cell_divisions"
		
		if(vogelstein){
			mr_res1<-mr_res1[!is.na(mr_res1$stem_cell_divisions),]

		}

		if(seer){
			mr_res1<-mr_res1[mr_res1$population == "European",]
			mr_res1<-mr_res1[!is.na(mr_res1$incidence),] #incidence, survival time and median age diagnosis are avaialable for the same set of cancers so can just drop cancers missing incidence

			# unique(mr_res1[order(mr_res1$incidence),c("incidence","survival_time","median.age.diagnosis","merge_index")])
		}
	}

	# Prune cancers. 
	mr_res1<-mr_res1[order(mr_res1$cases,decreasing=T),]
	mr_res1<-mr_res1[!duplicated(mr_res1$cancer),]
	# Prune cancer groups
	mr_res2<-mr_res1[mr_res1$Cancer.Group %in% cancer.group,]
	mr_res3<-mr_res1[mr_res1$cancer == cancer.group_excl,]
	mr_res4<-mr_res1[!mr_res1$Cancer.Group %in% cancer.group,]	
	mr_res2<-mr_res2[!duplicated(mr_res2$Cancer.Group),]
	mr_res2<-rbind(mr_res2,mr_res3)
	mr_res<-rbind(mr_res2,mr_res4)
	
	mr_res$b<-as.numeric(mr_res$b)
	mr_res$or<-exp(mr_res$b)
	mr_res$se<-as.numeric(mr_res$se)
	names(mr_res)[names(mr_res)=="study.abbreviation.y"]<-"study.abbreviation"

	# formating for plots
	mr_res$Colour<-mr_res$smoking1
	mr_res$Colour[mr_res$cancer %in% c("Meningioma","Glioma","Central nervous system and eye cancer","Breast cancer","Prostate cancer")]<-2
	mr_res$Colour[mr_res$cancer == "Endometrial cancer"]<-3
	# mr_res$Colour[mr_res$Colour == 0]<-4
	mr_res$Colour[mr_res$Colour == 0]<-2 
	mr_res$weight<-1/mr_res$se/5
	mr_res$Colour1<-NA
	mr_res$Colour1[mr_res$Colour == 1]<-"Increases risk"
	mr_res$Colour1[mr_res$Colour ==2]<-"Uncertain"
	# mr_res$Colour1[mr_res$Colour == 2]<-"No relationship"
	mr_res$Colour1[mr_res$Colour == 3]<-"Decreases risk"
	
	# mr_res$Colour1[mr_res$Colour == 4]<-"Unclear relationship"

	mr_res<-mr_res[order(mr_res$cases,decreasing=T),]
	mr_res1<-mr_res[mr_res$Colour == 1,]
	mr_res2<-mr_res[mr_res$Colour != 1,]
	# mr_res2<-mr_res[mr_res$Colour == 2,]
	# mr_res3<-mr_res[mr_res$Colour == 3,]
	# mr_res4<-mr_res[mr_res$Colour == 4,]
	mr_res<-rbind(mr_res1,mr_res2)
	# mr_res<-do.call(rbind,list(mr_res1,mr_res2,mr_res3,mr_res4))
	# mr_res[,c("cases","Colour1")]
	mr_res$weight<-1/mr_res$se^2
	mr_res$plot_name<-paste0(mr_res$cancer,"\n",mr_res$cases)
	# mr_res$smoking <- mr_res$Colour1
	mr_res$Shape<-"square"
	return(mr_res)
}


create_correlation_matrix<-function(study=NULL){
	corr_results_list1<-NULL
	corr_results_list2<-NULL
	M<-rep(0,nrow(study)^2)
	M<-matrix(M,nrow=nrow(study),ncol=nrow(study))
	j<-0
	# study$study.abbreviation

	while(j<nrow(study)-1){
		j<-j+1
		print(paste("j=",j))	
		# for(i in 1:nrow(study)){
		j2<-j+1
		print(paste("j2=",j2))
		for(i in j2:nrow(study)){

			print(paste("i=",i))
			test<-study[c(j,i),]
			corr_results_list1[[i]]<-correlation_between_cancers1(test=test)
			M[j,i]<-correlation_between_cancers1(test=test)
			M[i,j]<-correlation_between_cancers1(test=test)
		}
		corr_results_list2[[j]]<-corr_results_list1
		# M[j,i]<-unlist((corr_results_list2[[j]][i]))
		# M[i,j]<-unlist((corr_results_list2[[j]][i]))
		corr_results_list1<-NULL
	}
	for(i in 1:nrow(M)){
		M[i,i]<-1
	}
	out<-list("matrix"=M,"corr_results_list2"=corr_results_list2)
	# out<-M
	return(out)
}

# assess results across each independent study
# FinnGen typically has much fewer sNPs available compared to other studies. 
individual_studies_sensitvity_analysis<-function(dat=NULL,exclude_finngen=FALSE,excl_fads=FALSE){

	# if(exclude_finngen){
	# 	dat3<-dat3[dat3$study.abbreviation!="FinnGen",]
	# }
	dat3<-dat
	dat.meta_list<-NULL
	exposudat<-unique(dat3$exposure)

	for(i in 1:length(exposudat)){
		print(exposudat[i])
		input_dat<-dat3[dat3$exposure == exposudat[i],]
		dups<-unique(input_dat$outcome2[duplicated(input_dat$outcome2)])	
		print(input_dat[order(input_dat$outcome2),c("outcome2","nsnp","study.abbreviation")])
		dat.meta_list[[i]]<-meta_analysis3(dat=input_dat,exclude_finngen=exclude_finngen)
	} 
	dat.meta<-do.call(rbind,dat.meta_list)	
	dat4<-format_results4(dat=dat3,dat.meta=dat.meta,exclude_finngen=exclude_finngen)
	dat5<-format_results5(dat=dat3,dat.meta=dat.meta,exclude_finngen=exclude_finngen)
	Temp<-unique(dat5[!is.na(dat5$nsnp),c("exposure","nsnp","study.abbreviation")])
	Temp[order(Temp$exposure),]
	dat4$study.abbreviation
	dat5$b<-as.numeric(dat5$b)
	dat5$se<-as.numeric(dat5$se)

	# outcomes<-unique(dat5$outcome2[grep("FinnGen",dat5$study.abbreviation)])
	# dat5<-dat5[dat5$outcome2 %in% outcomes, ]
	exposures<-unique(dat5$exposure)

	for(i in 1:length(exposures)){
		plot_dat<-dat5[dat5$exposure == exposures[i] ,]
		print(exposures[i])
		plot_dat$name<-paste0(plot_dat$study.abbreviation," nsnps=",plot_dat$nsnp)
		P1<-forestplot(df = plot_dat,logodds = TRUE,name=name,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",weight=NULL)+theme(plot.title = element_text(size = ""),text = element_text(size=10))+ggtitle(exposures[i])
		P1<-P1+ggforce::facet_col(
		    facets = ~outcome2,
		    scales = "free_y",
		    space = "free" )
		exposure<-exposures[i]
		if(exclude_finngen){
			exposure<-paste0("excl_finngen_",exposure)
		}
		if(excl_fads){
			exposure<-paste0("excl_fads_",exposure)
		}

		png(paste0("~/fatty-acids/mr/results/plots/ggforest_secondary_pufas_",exposure,"by_individual_study_sensitivity_analysis.png"), width = 600, height = 600)
			print(P1) 
		dev.off()
	}
}

meta_p<-function(dat=NULL)
{
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
	meta_dat<-meta.tab9[,c("ID","study.abbreviation","population")]
	dat<-merge(dat,meta_dat,by.x="id.outcome",by.y="ID")
	outcomes<-unlist(strsplit(dat$outcome,split=" \\| "))
	dat$outcome2<-outcomes[seq(1,length(outcomes),2)]
	dat<-dat[dat$study.abbreviation != "FinnGen",]
	dat<-dat[dat$population == "European",]
	dat$outcome2<-gsub("Basal cell skin carcinoma" ,"Basal cell carcinoma",dat$outcome2 )
	
	outcomes<-unique(dat$outcome2)
	exposures<-unique(dat$exposure)
	methods<-unique(dat$method)
	b<-NULL
	for(i in 1:length(outcomes))
	{
		# i<-4
		print(outcomes[i])
		dat1<-dat[dat$outcome2 == outcomes[i],]
		# j<-8
		for(j in 1:length(exposures))
		{
			print(exposures[j])
			dat2<-dat1[dat1$exposure == exposures[j],]
			
			for(k in 1:length(methods))
			{
				# k<-2
				dat3<-dat2[dat2$method == methods[k],]		
				if(nrow(dat3)>1)
				{
					meta_p<-metap::sumlog(dat3$Q_pval)
					N_studies<-nrow(dat3)
					Q_pval<-paste0(round(dat3$Q_pval,5), collapse="; ")
					Q_df<-paste0(dat3$Q_df, collapse="; ")
					Q<-paste0(round(dat3$Q,5), collapse="; ")
					studies<-paste0(dat3$study.abbreviation, collapse="; ")	
					chisq_fisher<-meta_p$chisq
					df_fisher<-meta_p$df

					a<-data.frame(matrix(c(meta_p$p,outcomes[i],exposures[j],methods[k],N_studies,Q_pval,Q_df,Q,studies,chisq_fisher,df_fisher),nrow=1,ncol=11))

					names(a)<-c("Q_pval_fisher","outcome","exposure","method","nstudies","Q_pval","Q_df","Q","studies","chisq_fisher","df_fisher")

					b[[paste0(i,j,k)]]<-a
					c<-do.call(rbind,b)
				}
			}
		}
	}
	return(c)
}

exclude_ukb135_ala<-function()
{#exclude because there is only one SNP ("rs10215255") for ALA outside the FADS region, which is mising in UKB but present in 23andme (r2 proxies did not help retrieve the SNP)
	dat1<-dat[dat$id.outcome == 135, ]
	dat<-dat[dat$id.outcome != 135, ]
	dat1<-dat1[dat1$exposure!="Alpha-linolenic acid (18:3n3)",]
	dat<-rbind(dat,dat1)
	return(dat)
}



format_outcomes_csi<-function(dat=NULL,all_cols.keep=FALSE){

	names(dat)[names(dat) =="BETA"]<-"beta.outcome"
	names(dat)[names(dat) =="SE"]<-"se.outcome"
	names(dat)[names(dat) =="P"]<-"pval.outcome"
	names(dat)[names(dat) =="EFFECT_ALLELE"]<-"effect_allele.outcome"
	names(dat)[names(dat) =="OTHER_ALLELE"]<-"other_allele.outcome"
	names(dat)[names(dat) =="EAF"]<-"eaf.outcome"
	names(dat)[names(dat) =="SNP"]<-"SNP"
	dat$id.outcome<-"csi"
	dat$population<-"European"
	dat$study <- "UK Biobank" 
	dat$outcome<- "lifetime smoking index" 
	
	dat$se.outcome<-as.numeric(dat$se.outcome)
	dat$beta.outcome<-as.numeric(dat$beta.outcome)
	dat$eaf.outcome<-as.numeric(dat$eaf.outcome)
	dat$pval.outcome<-as.numeric(dat$pval.outcome)
	# Fix TRUE alleles as T alleles
	
	Pos<-which(dat$other_allele) 
	dat$other_allele.outcome[Pos]<-"T"

	# names(dat)[names(dat) =="cancer"]<-"outcome"
	# Cols.keep<-c("SNP","outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","ncase.outcome","ncontrol.outcome","effect_allele.outcome","other_allele.outcome","id.outcome","population","proxy","study")

	# dat<-dat[, Cols.keep]
	# dat$outcome2<-dat$outcome
	# dat$outcome<-paste0(dat$outcome," | ",dat$id.outcome)	
	
	# dat$index<-paste0(dat$outcome,dat$SNP)
	# Dups<-unique(dat$id.outcome[which(duplicated(dat$index))])
	# if(any(duplicated(dat$index))) warning(paste0("duplicate SNPs present, IDs=: ",paste(Dups,collapse=" | ")))

	return(dat)
}


format_snp_csi<-function()
{
	dat<-read.table("~/fatty-acids/mr/data/rs174546_rs2524299_lookup_csi_ukb.txt",sep="\t",head=T,stringsAsFactors=F)
	names(dat)[names(dat) =="BETA"]<-"beta.csi"
	names(dat)[names(dat) =="SE"]<-"se.csi"
	names(dat)[names(dat) =="P"]<-"pval.csi"
	names(dat)[names(dat) =="EFFECT_ALLELE"]<-"effect_allele.csi"
	names(dat)[names(dat) =="OTHER_ALLELE"]<-"other_allele.csi"
	names(dat)[names(dat) =="EAF"]<-"eaf.csi"
	names(dat)[names(dat) =="INFO"]<-"info.csi"
	names(dat)[names(dat) =="SNP"]<-"rsid"
	dat$population.csi<-"European"
	dat$study.csi <- "UK Biobank" 
	dat$exposure<- "lifetime smoking index" 
	
	dat$se.csi<-as.numeric(dat$se.csi)
	dat$beta.csi<-as.numeric(dat$beta.csi)
	dat$eaf.csi<-as.numeric(dat$eaf.csi)
	dat$pval.csi<-as.numeric(dat$pval.csi)
	# Fix TRUE alleles as T alleles
	
	Pos<-which(dat$other_allele.csi) 
	dat$other_allele.csi[Pos]<-"T"

	# dat<-format_outcomes_csi(dat=dat)
	# flip effect to reflect allele associated with higher CSI
	Pos<-dat$beta.csi<0
	beta<-dat$beta.csi[Pos]*-1
	eaf<-1-dat$eaf.csi[Pos]
	ea<-dat$other_allele.csi[Pos]
	oa<-dat$effect_allele.csi[Pos]
	dat$beta.csi[Pos]<-beta
	dat$eaf.csi[Pos]<-eaf
	dat$effect_allele.csi[Pos]<-ea
	dat$other_allele.csi[Pos]<-oa
	dat$lci.csi<-round(dat$beta.csi-1.96*dat$se.csi,3)
	dat$uci.csi<-round(dat$beta.csi+1.96*dat$se.csi,3)
	dat<-dat[,!names(dat) %in% c("CHR","BP")]
	return(dat)
}

rs174546_cancer<-function()
{
	load("~/fatty-acids/mr/data/dat_outcomes_final.Rdata")
	dat_outcomes_final$lnor<-as.numeric(dat_outcomes_final$lnor)
	dat_outcomes_final$se<-as.numeric(dat_outcomes_final$se)
	# dat_outcomes_final$or_total<-exp(dat_outcomes_final$lnor)
	# dat_outcomes_final$lci_total<-exp(dat_outcomes_final$lnor-1.96*dat_outcomes_final$se)
	# dat_outcomes_final$uci_total<-exp(dat_outcomes_final$lnor+1.96*dat_outcomes_final$se)
	dat_outcomes_final<-dat_outcomes_final[dat_outcomes_final$cancer %in% c("Lung cancer in never smokers","Lung cancer in ever smokers","Lung cancer","Colorectal cancer","Basal cell carcinoma","Cancer (all cause)"),c("rsid","cancer","population","ID","Effect.Allele","Other.Allele","eaf","lnor","se", "cases"  )]
	dat_outcomes_final<-dat_outcomes_final[dat_outcomes_final$population == "European",]
	dat_outcomes_final$cases<-as.numeric(dat_outcomes_final$cases)
	dat_outcomes_final<-dat_outcomes_final[order(dat_outcomes_final$cases,decreasing=TRUE),]
	dat_outcomes_final<-dat_outcomes_final[!duplicated(dat_outcomes_final$cancer),]
	dat_outcomes_final$cancer[dat_outcomes_final$cancer=="Cancer (all cause)"]<-"Overall cancer"
	names(dat_outcomes_final)[names(dat_outcomes_final)=="cancer"]<-"outcome"
	return(dat_outcomes_final)
}



rs2524299_cancer<-function()
{	
	load("~/fatty-acids/outcome_data/data/harmonised_data_postqc.Rdata")

	snp_cancer<-Dat_fa1[which(Dat_fa1$rsid=="rs2524299"),]
	snp_cancer<-snp_cancer[snp_cancer$outcome %in% c("Lung cancer in ever smokers","Lung cancer in never smokers", "Lung cancer","Colorectal cancer","Basal cell carcinoma","Overall cancer"),]
	snp_cancer<-snp_cancer[,c("rsid","outcome","population","ID","Effect.Allele","Other.Allele","eaf","lnor","se","ncase")]
	Pos<-snp_cancer$eaf>0.5
	eaf<-1-snp_cancer$eaf[Pos]
	lnor<-snp_cancer$lnor[Pos]*-1
	effect_allele<-snp_cancer$Other.Allele[Pos]
	other_allele<-snp_cancer$Effect.Allele[Pos]
	snp_cancer$eaf[Pos]<-eaf
	snp_cancer$lnor[Pos]<-lnor
	snp_cancer$Effect.Allele[Pos]<-effect_allele
	snp_cancer$Other.Allele[Pos]<-other_allele
	names(snp_cancer)[names(snp_cancer)=="ncase"]<-"cases"
	# snp_cancer$or_total<-exp(snp_cancer$lnor)
	# snp_cancer$lci_total<-exp(snp_cancer$lnor-1.96*snp_cancer$se)
	# snp_cancer$uci_total<-exp(snp_cancer$lnor+1.96*snp_cancer$se)
	return(snp_cancer)
}



snp_cancer<-function(ea_is_major_allele=TRUE){
	snp1_cancer<-rs174546_cancer()
	snp2_cancer<-rs2524299_cancer()
	snps_cancer<-rbind(snp1_cancer,snp2_cancer)
	snps_cancer<-snps_cancer[snps_cancer$population=="European",]
	names(snps_cancer)[names(snps_cancer)=="eaf"]<-"eaf.cancer"	
	names(snps_cancer)[names(snps_cancer)=="population"]<-"population.cancer"	
	snps_cancer$eaf.cancer<-as.numeric(snps_cancer$eaf.cancer)
	# express SNP-cancer effect to reflect allele associated with higher cig smk per day. this is the major allele
	if(ea_is_major_allele){	
		Pos.change<-snps_cancer$eaf.cancer<0.5		
	}
	if(!ea_is_major_allele)
	{
		Pos.change<-snps_cancer$eaf.cancer>0.5		
	}

	eaf<-1-snps_cancer$eaf.cancer[Pos.change]
	lnor<-snps_cancer$lnor[Pos.change]*-1
	oa<-snps_cancer$Effect.Allele[Pos.change]
	ea<-snps_cancer$Other.Allele[Pos.change]
	snps_cancer$eaf.cancer[Pos.change]<-eaf
	snps_cancer$lnor[Pos.change]<-lnor
	snps_cancer$Effect.Allele[Pos.change]<-ea
	snps_cancer$Other.Allele[Pos.change]<-oa
	snps_cancer$or_total<-exp(snps_cancer$lnor)
	snps_cancer$lci_total<-exp(snps_cancer$lnor-1.96*snps_cancer$se)
	snps_cancer$uci_total<-exp(snps_cancer$lnor+1.96*snps_cancer$se)
	names(snps_cancer)[names(snps_cancer)=="lnor"]<-"lnor_total"
	names(snps_cancer)[names(snps_cancer)=="se"]<-"se_total"
	return(snps_cancer)
}

# expressed SNP effects to reflect allele associated with higher smoking
format_snps_gscan<-function()
{
	snps_gscan$exposure<-snps_gscan$trait
	snps_gscan<-snps_gscan[,!names(snps_gscan) %in%  c("or","lci","uci","beta2","lci2","uci2")]
	names(snps_gscan)[names(snps_gscan)=="eaf"]<-"eaf.gscan"	
	Pos<-snps_gscan$beta<0
	beta<-snps_gscan$beta[Pos]*-1
	eaf<-1-snps_gscan$eaf.gscan[Pos]
	effect_allele<-snps_gscan$nea[Pos]
	other_allele<-snps_gscan$ea[Pos]
	snps_gscan$beta[Pos]<-beta
	snps_gscan$eaf.gscan[Pos]<-eaf
	snps_gscan$ea[Pos]<-effect_allele
	snps_gscan$nea[Pos]<-other_allele
	names(snps_gscan)[names(snps_gscan)=="beta"]<-"beta.gscan"
	names(snps_gscan)[names(snps_gscan)=="se"]<-"se.gscan"
	return(snps_gscan)
}

z_test_diff<-function(dat=NULL)
{
	# gscan_cancer2$lnor_decon	
	# gscan_cancer2$se_decon
	# gscan_cancer2$lnor_total	
	# gscan_cancer2$se_total	
	diff<-dat$lnor_total-dat$lnor_decon
	diff_se = sqrt(dat$se_total^2 + dat$se_decon^2)
	Z <- diff/diff_se
	# pnorm(abs(1.96),lower.tail=FALSE)*2
	dat$diff_p<-pnorm(abs(Z),lower.tail=FALSE)*2
	return(dat)
}


z_test_diff2<-function(lnor_total=NULL,se_total=NULL,lnor_decon=NULL,se_decon=NULL)
{
	# gscan_cancer2$lnor_decon	
	# gscan_cancer2$se_decon
	# gscan_cancer2$lnor_total	
	# gscan_cancer2$se_total	
	diff<-lnor_total-lnor_decon
	diff_se = sqrt(se_total^2 + se_decon^2)
	Z <- diff/diff_se
	# pnorm(abs(1.96),lower.tail=FALSE)*2
	diff_p<-pnorm(abs(Z),lower.tail=FALSE)*2
	return(diff_p)
}

gscan_decon<-function(method=NULL)
{#smk_can * snp_smk
	gscan_cancer<-merge(gscan_cancer,ao[,c("id","trait")],by.x="id.exposure",by.y="id")
	names(gscan_cancer)[names(gscan_cancer) == "exposure"]<-"exposure_plus_id"
	names(gscan_cancer)[names(gscan_cancer) == "trait"]<-"exposure"
	gscan_cancer<-merge(gscan_cancer,ao[,c("id","trait","note")],by.x="id.outcome",by.y="id")
	gscan_cancer$outcome[!is.na(gscan_cancer$note)]<-gscan_cancer$note[!is.na(gscan_cancer$note)]
	gscan_cancer$outcome[is.na(gscan_cancer$note)]<-"Lung cancer"
	gscan_cancer2<-gscan_cancer[gscan_cancer$method %in% method,]
	gscan_cancer2$or_smk_can<-exp(gscan_cancer2$b)
	gscan_cancer2$lci_smk_can<-exp(gscan_cancer2$b-1.96*gscan_cancer2$se)
	gscan_cancer2$uci_smk_can<-exp(gscan_cancer2$b+1.96*gscan_cancer2$se)
	snps_gscan<-format_snps_gscan()	
	gscan_cancer2<-merge(gscan_cancer2,snps_gscan,by="exposure")

	# calculate indirect or deconvoluated effect (effect of snp on cancer mediated by smoking)
	gscan_cancer2$lnor_decon<-gscan_cancer2$b*gscan_cancer2$beta.gscan
	gscan_cancer2$se_decon<-gscan_cancer2$se*gscan_cancer2$beta.gscan
	gscan_cancer2$or_decon<-exp(gscan_cancer2$lnor_decon)
	gscan_cancer2$lci_decon<-exp(gscan_cancer2$lnor_decon-1.96*gscan_cancer2$se_decon)
	gscan_cancer2$uci_decon<-exp(gscan_cancer2$lnor_decon+1.96*gscan_cancer2$se_decon)
	snps_cancer<-snp_cancer()
	snps_cancer<-snps_cancer[snps_cancer$ID %in% c("75","76","77"),] #restrict to lung cancer, lung cancer in ever smokers and lung cancer in never smokers
	# snps_cancer[,c("rsid","outcome","or_total","lci_total","uci_total")]
	gscan_cancer2<-merge(gscan_cancer2,snps_cancer,by=c("rsid","outcome"))

	gscan_cancer2$lci.gscan<-gscan_cancer2$beta.gscan-1.96*gscan_cancer2$se.gscan
	gscan_cancer2$uci.gscan<-gscan_cancer2$beta.gscan+1.96*gscan_cancer2$se.gscan
	gscan_cancer2<-z_test_diff(dat=gscan_cancer2)
	
	gscan_cancer2<-gscan_cancer2[,c("outcome","exposure","beta.gscan","se.gscan","lci.gscan","uci.gscan","or_total","lci_total","uci_total","or_decon", "lci_decon", "uci_decon","diff_p","or_smk_can", "lci_smk_can","uci_smk_can","rsid", "cases","Effect.Allele","Other.Allele","eaf.gscan","eaf.cancer","population.cancer","ID")]
	gscan_cancer2<-gscan_cancer2[gscan_cancer2$exposure == "Cigarettes smoked per day",]
	gscan_cancer2<-gscan_cancer2[order(gscan_cancer2$outcome,gscan_cancer2$rsid),]
	gscan_cancer2$or_total<-round(gscan_cancer2$or_total,4)
	gscan_cancer2$lci_total<-round(gscan_cancer2$lci_total,4)
	gscan_cancer2$uci_total<-round(gscan_cancer2$uci_total,4)
	gscan_cancer2$or_decon<-round(gscan_cancer2$or_decon,4)
	gscan_cancer2$lci_decon<-round(gscan_cancer2$lci_decon,4)
	gscan_cancer2$uci_decon<-round(gscan_cancer2$uci_decon,4)
	gscan_cancer2$or_smk_can<-round(gscan_cancer2$or_smk_can,4)
	gscan_cancer2$lci_smk_can<-round(gscan_cancer2$lci_smk_can,4)
	gscan_cancer2$uci_smk_can<-round(gscan_cancer2$uci_smk_can,4)
	names(gscan_cancer2)[names(gscan_cancer2) == "population.cancer"]<-"population"
	gscan_cancer2<-gscan_cancer2[order(gscan_cancer2$rsid),]
	return(gscan_cancer2)
}



format_res_csi<-function()
{
	load("~/fatty-acids/mr/results/res_csi.Rdata")
	res_csi<-res_csi[res_csi$outcome2 %in% c("Lung cancer","Colorectal cancer","Basal cell carcinoma","Overall cancer"),]
	# res_csi[res_csi$outcome2 %in% c("Lung cancer","Colorectal cancer","Basal cell carcinoma","Overall cancer"), ]


	res_csi$b<-as.numeric(res_csi$b)
	res_csi$se<-as.numeric(res_csi$se)
	res_csi$or_smk_can<-round(exp(res_csi$b),3)
	res_csi$lci_smk_can<-round(exp(res_csi$b-1.96*res_csi$se),3)
	res_csi$uci_smk_can<-round(exp(res_csi$b+1.96*res_csi$se),3)
	res_csi$exposure<-"lifetime smoking index"
	res_csi$outcome<-res_csi$outcome2
	res_csi<-res_csi[,names(res_csi)!="outcome2"]
	return(res_csi)
}

csi_decon<-function()
{	
	res_csi<-merge(res_csi,snp_csi,by="exposure") 	
	res_csi$lnor_decon<-res_csi$b*res_csi$beta.csi
	res_csi$se_decon<-res_csi$se*res_csi$beta.csi
	res_csi$or_decon<-exp(res_csi$lnor_decon)
	res_csi$lci_decon<-exp(res_csi$lnor_decon-1.96*res_csi$se_decon)
	res_csi$uci_decon<-exp(res_csi$lnor_decon+1.96*res_csi$se_decon)
	
	snps_cancer<-snp_cancer(ea_is_major_allele=FALSE)#make sure effect allele is minor allele, which is the allele associated with higher CSI (albet not sig and very uncertain and close to zero

	res_csi<-merge(res_csi,snps_cancer,by=c("rsid","outcome"))
	res_csi$or_total<-round(res_csi$or_total,3)
	res_csi$lci_total<-round(res_csi$lci_total,3)
	res_csi$uci_total<-round(res_csi$uci_total,3)

	res_csi<-z_test_diff(dat=res_csi)
	res_csi<-res_csi[,c("outcome","exposure","rsid","beta.csi","se.csi","lci.csi","uci.csi","or_total","lci_total","uci_total","or_decon", "lci_decon", "uci_decon","diff_p","or_smk_can", "lci_smk_can","uci_smk_can", "cases","eaf.csi","eaf.cancer","Effect.Allele","Other.Allele","population","ID")]
	res_csi<-res_csi[order(res_csi$cases,decreasing=TRUE),]
	res_csi<-res_csi[!duplicated(paste(res_csi$outcome,res_csi$rsid)),]
	res_csi<-res_csi[order(res_csi$rsid,res_csi$outcome),]

	return(res_csi)
}


make_instab_csi<-function()
{
	Pos<-grep("exposure",names(dat))
	Pos1<-which(names(dat)=="SNP")
	Pos2<-which(names(dat)=="outcome2")
	Pos3<-which(names(dat)=="id.outcome")
	Pos4<-which(names(dat)=="ncase.outcome")
	Pos<-c(Pos1,Pos,Pos2,Pos3,Pos4)
	exposure_dataset<-dat[,Pos]
	exposure_dataset<-exposure_dataset[exposure_dataset$outcome2 %in% c("Basal cell carcinoma","Colorectal cancer","Lung cancer","Overall cancer"),]
	exposure_dataset<-exposure_dataset[order(exposure_dataset$ncase.outcome,decreasing=TRUE),]
	exposure_dataset<-exposure_dataset[!duplicated(paste(exposure_dataset$outcome2,exposure_dataset$SNP)),]
	snps<-unique(exposure_dataset$SNP)
	dat.list<-NULL
	for(i in 1:length(snps))
	{
		print(snps[i])
		temp<-exposure_dataset[exposure_dataset$SNP == snps[i],]
		temp$outcome<-paste(temp$outcome2,collapse="; ")
		temp$id.outcome<-paste(temp$id.outcome,collapse="; ")
		temp<-temp[!duplicated(temp$SNP),]
		dat.list[[i]]<-temp
	}
	dat.snp<-do.call(rbind,dat.list)
	dat.snp<-dat.snp[,!names(dat.snp) %in% c("ncase.outcome","outcome2")]

	write.table(dat.snp,"~/fatty-acids/mr/data/csi_instrument.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}


including_laryngeal_cancer<-function()
{
	laryngeal_cancer_lnor <-(log(0.71)*-1)/0.72 #lnor per 1 SD increase FADS1/2 activity 25194280, rs174549 r2=0.9786 with rs174546 in CHB+JPT (north east asians)
	 
	lnor<-log(0.71) #originally reported odds ratio
	lci<-log(0.63)
	uci<-log(0.80) 
	se<-(uci-lci)/(1.96*2)
	lnor_d6d<-lnor/0.72*-1 #originally reported odds ratio refers to allele that lowers d6d
	se_d6d<-se/0.72 #this refers to effect of rs174546 on DGLA:LA
	laryngeal_cancer_se <-se/0.72 #0.72 is estimated effect of allele on D6D in East Asians 
	p<-pnorm(laryngeal_cancer_lnor/laryngeal_cancer_se,lower.tail=F)*2
	dat<-data.frame(matrix(c(laryngeal_cancer_lnor,laryngeal_cancer_se,p),nrow=1,ncol=3))
	names(dat)<-c("lnor","se","p")
	return(dat)
}


extract_meta_reg_results<-function(dat=NULL)
{
	
	r2<-dat$R2
	k<-Model$k
	n1<-"moderator"
	b1<-dat$beta[2]
	se1<-dat$se[2]
	z1<-dat$zval[2]
	p1<-dat$pval[2]
	or1<-exp(b1)
	lci1<-exp(b1-1.96*se1)
	uci1<-exp(b1+1.96*se1)

	n0<-"intercept"
	b0<-dat$beta[1]
	se0<-dat$se[1]
	z0<-dat$zval[1]
	p0<-dat$pval[1]	
	or0<-exp(b0)
	lci0<-exp(b0-1.96*se0)
	uci0<-exp(b0+1.96*se0)

	a<-c(n1,b1,se1,z1,p1,or1,lci1,uci1,n0,b0,se0,z0,p0,or0,lci0,uci0)
	b<-data.frame(matrix(a,nrow=2,ncol=length(a)/2,byrow=TRUE))
	b$r2<-r2
	b$k<-k
	names(b)<-c("coef","lnor","se","z","p","or","lci","uci","r2","k")	
	return(b)
}

extract_meta_results<-function(dat=NULL)
{
	fe<-dat$TE.fixed
	se.fe<-dat$seTE.fixed
	or_fe<-exp(fe)
	or_fe_lci<-exp(dat$lower.fixed)
	or_fe_uci<-exp(dat$upper.fixed)
	re<-dat$TE.random
	se.re<-dat$seTE.random
	or_re<-exp(re)
	or_re_lci<-exp(dat$lower.random)
	or_re_uci<-exp(dat$upper.random)
	fe.pval<-dat$pval.fixed
	re.pval<-dat$pval.random
	Q<-dat$Q
	df.Q<-dat$df.Q
	pval.Q<-dat$pval.Q
	K<-dat$k
	I2<-dat$I2
	I2_lower<-dat$lower.I2
	I2_upper<-dat$upper.I2
	
	a<-c(fe,se.fe,or_fe,or_fe_lci,or_fe_uci,re,se.re,or_re,or_re_lci,or_re_uci,fe.pval,re.pval,Q,df.Q,pval.Q,K,I2,I2_lower,I2_upper)
	b<-data.frame(matrix(a,nrow=1,ncol=length(a)))
	names(b)<-c("fe","se.fe","or_fe","or_fe_lci","or_fe_uci","re","se.re","or_re","or_re_lci","or_re_uci","fe.pval","re.pval","Q","df.Q","pval.Q","K","I2","I2_lower","I2_upper")	
	return(b)
}