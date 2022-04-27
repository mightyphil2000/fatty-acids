

preformat_cer_id129_v2<-function(){
	File1<-"~/MR_FattyAcids/data/summary_data/cervical cancer/smoking_interaction/data_gwasLookupDataInclA2_complete_inclMAF_forSmoking.csv"
	File2<-"~/MR_FattyAcids/data/summary_data/cervical cancer/smoking_interaction/data_gwasLookupDataInclA2_complete_inclMAF_forSmoking_2.csv"
	cer1<-read.table(File1,sep=";",head=T,stringsAsFactors=F)
	cer2<-read.table(File2,sep=";",head=T,stringsAsFactors=F)
	cer<-rbind(cer1,cer2)
	cer<-cer[!duplicated(cer$SNP),]
	cer$OR<-as.numeric(gsub(",",".",cer$OR))
	cer$SE<-as.numeric(gsub(",",".",cer$SE))
	cer$MAF<-as.numeric(gsub(",",".",cer$MAF))
	return(cer)
}

preformat_bla_id105<-function(){
	AA<-bla$controls_AA 
	AB<-bla$controls_AB 
	BB<-bla$controls_BB
	Nchr<-(AA+AB+BB)*2
	N_b_alleles<-AB+BB*2
	bla$eaf<-N_b_alleles/Nchr
	return(bla)
}


preformat_kid_90<-function(dat=NULL){
	Pos<-dat$Effect_allele != dat$Allele_2 & dat$Effect_allele != dat$Allele_1 
	dat<-dat[!Pos,]
	Pos<-which(dat$Effect_allele != dat$Allele_2)
	ea<-dat$Allele_2
	oa<-dat$Allele_1
	oa[Pos]<-dat$Allele_2[Pos]
	ea[Pos]<-dat$Allele_1[Pos]
	dat$Effect_allele<-ea
	dat$Other_allele<-oa	
	return(dat)
}

preformat_ova_120<-function(dat=NULL){
	ref<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",sep=" ",head=F,stringsAsFactors=F)
	ref1<-format_ref(dat=ref)
	dat<-merge(dat,ref1,by.x=c("Chromosome","Position"),by.y=c("chr","V2"))
	return(dat)
}

format_ref<-function(dat=NULL){
	chr<-unlist(strsplit(dat$V1,split="chr")) 
	dat$chr<-as.numeric(chr[seq(2,length(chr),2)])	
	return(dat)
}


extract_snps_and_format_bcc_23andMe<-function(snplist=NULL){
	all<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/all_snp_info-4.1.txt",exact_match=TRUE,path_to_target_file_sep="\t")
	gt<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/gt_snp_stat-4.1.txt",exact_match=TRUE,path_to_target_file_sep="\t")
	im<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/im_snp_stat-4.1.txt",exact_match=TRUE,path_to_target_file_sep="\t")
	bcc<-extract_snps(snplist=all$all.data.id,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_bcc/Chahal_2016_basal_cell_carcinoma-4.1/basal_cell_carcinoma-4.1.dat",exact_match=TRUE,path_to_target_file_sep="\t")
	bcc2<-format_23andme(dat=bcc,ncase=12945, ncontrol=274252,outcome="Basal cell skin carcinoma",pmid=27539887,ID=1,efo="basal cell carcinoma",all=all,gt=gt,im=im)
	return(bcc2)
}


extract_snps_and_format_scc_23andMe<-function(snplist=NULL){
	all<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/all_snp_info-4.1.txt",exact_match=TRUE,path_to_target_file_sep="\t")
	gt<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/gt_snp_stat-4.1.txt",exact_match=TRUE,path_to_target_file_sep="\t")
	im<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/im_snp_stat-4.1.txt",exact_match=TRUE,path_to_target_file_sep="\t")
	scc<-extract_snps(snplist=all$all.data.id,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_scc/squamous_cell_carcinoma-4.1.dat",exact_match=TRUE,path_to_target_file_sep="\t")
	scc2<-format_23andme(dat=scc,ncase=6579, ncontrol=280558,outcome="Squamous cell skin carcinoma",pmid=27424798,ID=2,efo="cutaneous squamous cell carcinoma",all=all,gt=gt,im=im)
	return(scc2)
}



# scc2[scc2$rsid %in% c("rs2236212", "rs174576",  "rs603424",  "rs174546" ),]
format_23andme<-function(dat=NULL,outcome=NULL,ncase=NULL,ncontrol=NULL,pmid=NULL,ID=NULL,efo=NULL,all=NULL,gt=NULL,im=NULL){
    dat2<-merge(dat,all,by="all.data.id")
    datI<-dat2[dat2$src == "I",]
    datI2<-merge(datI,im,by="assay.name",)
    datG<-dat2[dat2$src == "G",]
    datG2<-merge(datG,gt,by="assay.name",)

    if(nrow(datG2)!=0){
        dat3<-plyr::rbind.fill(datG2,datI2)
    }else{
        dat3<-datI2
    }
    names(dat3)[names(dat3) == "file.x"]<-"file_results"
    names(dat3)[names(dat3) == "file.y"]<-"all_annotations"

   	dat3<-dat3[, names(dat3)[!names(dat3) %in% c("gt.data.id.x", "gt.data.id.y")]]
	dat4<-merge(dat3,gt,by="assay.name",all.x=T)

    # dat4[dat4$assay.name %in% c("rs2236212", "rs174576",  "rs603424",  "rs174546" ),]

    # min(dat3$avg.call)
    # dat1<-dat[dat$src == "I",]
    # dat2<-dat[dat$src == "G",]

    # dat3[dat3$assay.name %in% c("rs2236212", "rs174576",  "rs603424",  "rs174546" ),]
    # alleles The two SNP alleles, A/B, in alphabetical order
    Alleles<-unlist(strsplit(dat4$alleles,split="/"))
    dat4$other_allele<-Alleles[seq(1,length(Alleles),2)]
    dat4$effect_allele<-Alleles[seq(2,length(Alleles),2)]
    dat<-format_data(dat=dat4,outcome=outcome,population="European",pmid=pmid,study="23NMSC",ncase=ncase,ncontrol=ncontrol,UKbiobank=FALSE,rsid="assay.name",effect_allele="effect_allele",other_allele="other_allele",eaf="freq.b.x",lnor="effect",lnor_se="stderr",p="pvalue",effect_allele_confirmed=TRUE,ID=ID,all_summary_stats=TRUE,info1="avg.rsqr",HWEp="hw.p.value",efo=efo)

    # names(dat4)[names(dat4) == "avg.rsqr" ]<-"info"
    # names(dat4)[names(dat4) == "freq.b.x" ]<-"eaf"    
    # names(dat4)[names(dat4) == "effect" ]<-"lnor"
    # names(dat4)[names(dat4) == "stderr" ]<-"se"
    # names(dat4)[names(dat4) == "pvalue" ]<-"p"
    # names(dat4)[names(dat4) == "assay.name"]<-"rsid"
    # names(dat4)[names(dat4) == "hw.p.value"]<-"HWEp"
    # dat4$effect_allele_confirmed<-TRUE
    # dat4$UKbiobank<-FALSE
    # dat4$study<-"23andMe"
    # dat4$ncase<-ncase
    # dat4$ncontrol<-ncontrol
    # dat4$population<-"European"
    # dat4$outcome<-outcome
    # dat4$pmid<-pmid
	
    return(dat)
}



preformat_ewi_27<-function(snplist=NULL){
	Ewi<-read.csv("/projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Postel-Vinay22327514/data/Ewing_association_Apr2017.csv",stringsAsFactors=F,head=F,skip=1)

	Ewi<-Ewi[,names(Ewi) != "V1"]
	names(Ewi)<-c("initial_order","CHR","SNP","BP","A1..effect_allele.","TEST","NMISS","OR","ORs.allele.freqs.calculations","relative.error.OR","SE","L99","U99","STAT","P","A1.Ewing","A2.Ewing","MAF.Ewing","A1.controls","A2.controls","MAF.controls","maf.control.inferred.adj.with.ewing.allele","relative.error.maf.control","A1.HWE","O.HET.ALL","E.HET.ALL","HWE.Pval.ALL","O.HET.AFF","E.HET.AFF","HWE.Pval.AFF","O.HET.UNAFF","E.HET.UNAFF","HWE.Pval.UNAFF")
	Ewi$SNP<-trimws(Ewi$SNP)
	# write.table(Ewi,"/projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Postel-Vinay22327514/data/Ewing_association_Apr2017.tab",col.names=T,row.names=F,quote=F)
	# Dat<-lapply(names(Ewi), FUN=function(x) trimws(Ewi[,x]))
	# Dat2<-do.call(cbind,Dat)

	# snplist<-unique(readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"))
	Ewi2<-Ewi[Ewi$SNP %in% snplist,]
	Ewi2$A1.Ewing<-trimws(Ewi2$A1.Ewing)
	Ewi2$A1.controls<-trimws(Ewi2$A1.controls)
	Ewi2$A1..effect_allele.<-trimws(Ewi2$A1..effect_allele.)
	Ewi2<-Ewi2[Ewi2$A1..effect_allele. == Ewi2$A1.controls,]
	return(Ewi2)
}



preformat_snps_interlymph<-function(dat=NULL){
	# snplist<-unique(readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"))
	# Res<-extract_snps(snplist=snplist,File=File,exact_match=TRUE,file_sep="\t",fill=TRUE)
	
	dat<-dat[order(dat$NUM_CASE,decreasing=T),]
	snps<-unique(dat$SNP)
	# i<-which(snps=="rs175133")
	Dat_list<-NULL
	for(i in 1:length(snps)){
		# print(snps[i])
		Dat<-dat[dat$SNP == snps[i],]
		Dat1<-Dat[!is.na(Dat$INFO),]
		Dat2<-Dat[is.na(Dat$INFO),]

		INFO<-sum(Dat1$INFO*Dat1$NUM_CASE)/sum(Dat1$NUM_CASE)
		eaf <-Dat1$EFFECT_ALLELE_FREQ_CONTROL
		Test<-all(eaf < 0.5 | eaf > 0.5)
		eaf_test<-FALSE
		if(!Test) eaf_test<-TRUE #one example observed where eaf was 0.5 in controls and was a bit less than 0.5 in controls in other studies. 
		# stop("allele frequency not consistent across studies")
		
		EAF<-sum(Dat1$EFFECT_ALLELE_FREQ_CONTROL*Dat1$NUM_CONTROL)/sum(Dat1$NUM_CONTROL)
		oa<-unique(Dat1$REFERENCE_ALLELE)
		ea<-unique(Dat1$EFFECT_ALLELE)
		
		allele_test<-FALSE
		if(length(oa) != 1 | length(ea) != 1  ) allele_test<-TRUE

		# stop("different effect alleles across studies")
		
		lnor<-log(Dat2$OR)
		CI<-unlist(strsplit(Dat2$CI,split="-"))
		CI<-unlist(strsplit(CI,split="\\("))
		CI<-as.numeric(unlist(strsplit(CI,split=")")))
		lci<-CI[1]
		uci<-CI[2]
		se<-(log(uci) - log(lci))/(1.96*2)
		Dat1$Direction<-NA
		Dat1$Direction[Dat1$OR>1]<-"+"
		Dat1$Direction[Dat1$OR<1]<-"-"
		Direction<-paste(Dat1$Direction,collapse="")
		
		Dat3<-data.frame(matrix(c(INFO,EAF,oa,ea,lnor,lci,uci,se,Direction,eaf_test,allele_test),ncol=11,nrow=1),stringsAsFactors=F)

		names(Dat3)<-c("info","eaf","other_allele","effect_allele","lnor","lci","uci","se","Direction","eaf_test","allele_test")
		Dat4<-cbind(Dat2,Dat3)

		# names(Dat4)[names(Dat4)=="SNP"]<-"rsid"
		# names(Dat4)[names(Dat4)=="P"]<-"p"
		# names(Dat4)[names(Dat4)=="Phet"]<-"phet"
		# names(Dat4)[names(Dat4)=="NUM_CASE"]<-"ncase"
		# names(Dat4)[names(Dat4)=="NUM_CONTROL"]<-"ncontrol"
		# Dat4$outcome<-outcome
		# Dat4$study<-"InterLymph"
		# Dat4$pmid<-pmid
		# Dat4$population<-"European"
		# Dat4$UKbiobank<-FALSE
		# Dat4$effect_allele_confirmed<-TRUE #summary data contains a column called EFFECT_ALLELE
		Dat_list[[i]]<-Dat4
	}
	Dat5<-do.call(rbind,Dat_list)
	return(Dat5)
}

preformat_gli_67<-function(dat=NULL){
	N<-unlist(strsplit(dat$Subjects,split="\\|"))
	dat$controls<-as.numeric(N[seq(1,length(N),by=2)])
	dat$cases<-as.numeric(N[seq(2,length(N),by=2)])
	N<-unlist(strsplit(dat$MAF,split="\\|"))
	dat$eaf.controls<-as.numeric(N[seq(1,length(N),by=2)])
	return(dat)
}





transform_betas<-function(dat=NULL,effect="lnor",effect.se="lnor_se"){
	# formula: log OR = beta / (u(1-u)); where u=ncases/(ncases + ncontrol) REPEAT with SE 	
	beta<-dat[,effect]
	se<-dat[,effect.se]
	u<-dat$ncase/(dat$ncase+dat$ncontrol)
	dat[,effect] <- beta / (u * (1 - u))
	dat[,effect.se]<-se / (u * (1 - u)) 	
	return(dat)
}


