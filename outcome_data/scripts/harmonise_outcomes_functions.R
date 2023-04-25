library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

format_results<-function(Dat=NULL){
	Dat<-Dat[Dat$nstudies==max(Dat$nstudies),]
	Dat$study<-Dat$study[1]
	Dat2<-Dat[!Dat$alleles %in% c("GC","CG","AT","TA"),]
	Dat2$study<-paste0(Dat2$study," excluding palindromic SNPs")	
	Dat<-rbind(Dat,Dat2)
	Dat$fisherp<-as.numeric(Dat$fisherp)
	Dat$Z_fisher<-qnorm(Dat$fisherp/2,lower.tail=F)
	return(Dat)
}

clean_and_combine_data2<-function(Dat=NULL){
	Dat<-do.call(plyr::rbind.fill,Dat)

	if(any("info" == names(Dat))){		
		Dat<-Dat[Dat$info>=0.80 | is.na(Dat$info),]		
	}
	
	Dat<-Dat[Dat$eaf>0.01,]

	Dat$alleles<-paste(Dat$Effect.Allele,Dat$Other.Allele,sep="")
	Dat<-Dat[which(nchar(Dat$alleles) == 2),]

	# Pos<-Dat$alleles %in% c("GC","CG","TA","AT")
	# Dat1<-Dat[Pos,]
	# Dat1<-Dat1[which(Dat1$eaf <0.40 | Dat1$eaf > 0.60), ]
	# Dat2<-Dat[!Pos,]
	# Dat<-rbind(Dat1,Dat2)		
	return(Dat)
}


clean_and_combine_data<-function(Dat=NULL){
	Dat<-do.call(plyr::rbind.fill,Dat)

	if(any("info" == names(Dat))){		
		Dat<-Dat[Dat$info>=0.80 | is.na(Dat$info),]		
	}
	
	Dat<-Dat[Dat$eaf>0.01 | is.na(Dat$eaf),]

	Dat$alleles<-paste(Dat$Effect.Allele,Dat$Other.Allele,sep="")
	Dat<-Dat[which(nchar(Dat$alleles) == 2),]

	Pos<-Dat$alleles %in% c("GC","CG","TA","AT")
	Dat1<-Dat[Pos,]
	Dat1<-Dat1[which(Dat1$eaf <0.40 | Dat1$eaf > 0.60), ]
	Dat2<-Dat[!Pos,]
	Dat<-rbind(Dat1,Dat2)		
	return(Dat)
}


check_eaf<-function(eaf=NULL){
	if( all(eaf<0.5) | all(eaf>0.5)) message("eaf consistent across SNPs")
	if(!all(eaf<0.5)) {
		if(!all(eaf>0.5)){
			stop("eaf inconsistent across SNP")
		}
	}
}


Combinep_fisher<-function(Datp=NULL){
	Chi<--2*sum(log(Datp))
	df<-length(Datp)*2
	p_fisher<-pchisq(Chi, df,lower.tail=F)
	return(p_fisher)
}



meta_analysis2<-function(Dat=NULL){
	snps<-unique(Dat$rsid)
	Res_list<-NULL	
	for(i in 1:length(snps)){
		print(i)
		Dat1<-Dat[Dat$rsid ==snps[i],]
		Dat1<-harmonise_effect_allele(Dat=Dat1)
		EA<-unique(Dat1$Effect.Allele)
		OA<-unique(Dat1$Other.Allele)
		alleles<-paste0(EA,OA)
		if(any(!alleles %in% c("GC","CG","TA","AT"))){
			if(length(EA)!=1) stop("effect allele not harmonised")
			if(length(OA)!=1) stop("effect allele not harmonised")
		}
		EA<-EA[1]
		OA<-OA[1]
		alleles<-alleles[1]

		eaf<-Dat1$eaf		
		if(any(alleles %in% c("GC","CG","TA","AT"))){
			check_eaf(eaf=eaf)
		}

		b<-Dat1$lnor
		se<-Dat1$se
		w<-1/se^2
		EAF<-round(sum((Dat1$eaf*w))/sum(w),3)
		b.fixed<-sum(b*w)/(sum(w))
		se.fixed<-sqrt(sum(w)^-1)
		z<-abs(b.fixed/se.fixed)
		p.fixed<-pnorm(z,lower.tail=F)*2
		nstudies.fixed<-length(b)
		cancer<-unique(Dat1$outcome)		
		if(length(cancer)!=1) stop("length of cancer not 1")
		ids.fixed<-paste(Dat1$ID,collapse="; ")
		cases<-sum(Dat1$ncase)
		controls<-sum(Dat1$ncontrol)
		study<-paste(Dat1$study,collapse="/")
		Q<-sum((b.fixed-b)^2*w)
		df.Q<-length(b)-1		
		Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)

		P_fisher<-Combinep_fisher(Datp=Dat1$p)
		rsid<-snps[i]
		
		if(length(EA)>1) stop("effect allele not consistent across studies")	
		# Cols<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
		# for(j in 1:length(Cols)){
		# 	Cols[j]
		# 	Cols[j]<-paste(unique(Dat1[,Cols[j]]),collapse="; ")
		# }
		# Col.dat<-data.frame(matrix(Cols,ncol=length(Cols),nrow=1),stringsAsFactors=F)
		# names(Col.dat)<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
		
		dat.matrix<-c(rsid,alleles,cancer,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,cases,controls,study,Q.p,EA,OA,EAF,P_fisher)
		Res<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
		
		names(Res)<-c("rsid","alleles","outcome","lnor","se","p","nstudies","ID","ncase","ncontrol","study","Q.p","Effect.Allele","Other.Allele","eaf","fisherp")

		# Res.c<-cbind(Res,Col.dat)		
		Res_list[[i]]<-Res
	}
	
	Res1<-do.call(rbind,Res_list)
	Res1$p<-as.numeric(Res1$p)

	return(Res1)
}

extract_snps_and_format_bcc_23andMe<-function(snplist=NULL){
	all<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/all_snp_info-4.1.txt",exact_match=TRUE,file_sep="\t")
	gt<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/gt_snp_stat-4.1.txt",exact_match=TRUE,file_sep="\t")
	im<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/im_snp_stat-4.1.txt",exact_match=TRUE,file_sep="\t")
	bcc<-extract_snps(snplist=all$all.data.id,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_bcc/Chahal_2016_basal_cell_carcinoma-4.1/basal_cell_carcinoma-4.1.dat",exact_match=TRUE,file_sep="\t")
	bcc2<-format_23andme(dat=bcc,ncase=12945, ncontrol=274252,outcome="Basal cell skin carcinoma",pmid=27539887,ID=1,efo="basal cell carcinoma",all=all,gt=gt,im=im)
	return(bcc2)
}


extract_snps_and_format_scc_23andMe<-function(snplist=NULL){
	all<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/all_snp_info-4.1.txt",exact_match=TRUE,file_sep="\t")
	gt<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/gt_snp_stat-4.1.txt",exact_match=TRUE,file_sep="\t")
	im<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/im_snp_stat-4.1.txt",exact_match=TRUE,file_sep="\t")
	scc<-extract_snps(snplist=all$all.data.id,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_scc/squamous_cell_carcinoma-4.1.dat",exact_match=TRUE,file_sep="\t")
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
    dat4$Other.Allele<-Alleles[seq(1,length(Alleles),2)]
    dat4$Effect.Allele<-Alleles[seq(2,length(Alleles),2)]
    dat<-format_data(Dat=dat4,outcome=outcome,population="European",pmid=pmid,study="23NMSC",ncase=ncase,ncontrol=ncontrol,UKbiobank=FALSE,rsid="assay.name",Effect.Allele="Effect.Allele",Other.Allele="Other.Allele",eaf="freq.b.x",lnor="effect",se="stderr",p="pvalue",effect_allele_confirmed=TRUE,ID=ID,all_summary_stats=TRUE,info="avg.rsqr",HWEp="hw.p.value",efo=efo)

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



# format_23andme2<-function(dat=NULL,outcome=NULL,ncase=NULL,ncontrol=NULL,pmid=NULL,ID=NULL){
#     dat2<-merge(dat,all,by="all.data.id")
#     datI<-dat2[dat2$src == "I",]
#     datI2<-merge(datI,im,by="assay.name",)
#     datG<-dat2[dat2$src == "G",]
#     datG2<-merge(datG,gt,by="assay.name",)

#     if(nrow(datG2)!=0){
#         dat3<-plyr::rbind.fill(datG2,datI2)
#     }else{
#         dat3<-datI2
#     }
    
#    # dat3[dat3$assay.name %in% c("rs2236212", "rs174576",  "rs603424",  "rs174546" ),]
#     # alleles The two SNP alleles, A/B, in alphabetical order
#     Alleles<-unlist(strsplit(dat3$alleles,split="/"))
#     dat3$Other.Allele<-Alleles[seq(1,length(Alleles),2)]
#     dat3$Effect.Allele<-Alleles[seq(2,length(Alleles),2)]
#     Chr<-unlist(strsplit(dat3$scaffold,split="chr"))
#     Chr<-Chr[Chr!=""]
#     dat3$chr<-Chr
    

#     dat<-format_data(Dat=dat3,outcome=outcome,population="European",pmid=pmid,study="23NMSC",ncase=ncase,ncontrol=ncontrol,UKbiobank=FALSE,rsid="assay.name",Effect.Allele="Effect.Allele",Other.Allele="Other.Allele",eaf="freq.b",lnor="effect",se="stderr",p="pvalue",effect_allele_confirmed=TRUE,ID=ID,all_summary_stats=TRUE,info="avg.rsqr",HWEp="hw.p.value")

#     # names(dat3)[names(dat3) == "avg.rsqr" ]<-"info"
#     # if(all(!names(dat3) == "freq.b")) stop("freq.b not present")
#     # names(dat3)[names(dat3) == "freq.b" ]<-"eaf"
#     # names(dat3)[names(dat3) == "effect" ]<-"lnor"
#     # names(dat3)[names(dat3) == "stderr" ]<-"se"
#     # names(dat3)[names(dat3) == "pvalue" ]<-"p"
#     # names(dat3)[names(dat3) == "assay.name"]<-"rsid"
#     # names(dat3)[names(dat3) == "hw.p.value"]<-"HWEp"
#     # dat3$effect_allele_confirmed<-TRUE
#     # dat3$UKbiobank<-FALSE
#     # dat3$study<-"23andMe"
#     # dat3$ncase<-ncase
#     # dat3$ncontrol<-ncontrol
#     # dat3$population<-"European"
#     # dat3$outcome<-outcome
#     # dat3$pmid<-pmid
    
#     # if(!is.null(ID)){
#     # 	dat3$ID<-ID
#     # }
#     return(dat)
# }


find_rsids<-function(dat=NULL,ref_dat=FALSE){	
	if(ref_dat){
		ref<-read.table("~/fatty-acids/colocalisation/data/UKBB_10K_bed_hg19.txt",head=F,stringsAsFactors=F,sep=" ")
	}
	dat2<-merge(dat,ref,by.x=c("chr","bp"),by.y=c("V1","V2"))
	names(dat2)[names(dat2) == "V4"]<-"rsid"
	return(dat2)
}



# preformat_snps_mzl_interlymph<-function(dat=NULL){
# 	# snplist<-unique(readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"))
# 	# Res<-extract_snps(snplist=snplist,File=File,exact_match=TRUE,file_sep="\t",fill=TRUE)
# 	names(dat)[names(dat) == "Locus"]<-"rsid"
# 	dat$outcome<-outcome
# 	dat$population<-"European"
# 	dat$pmid<-pmid
# 	dat$study<-"InterLymph"
# 	names(dat)[names(dat)=="Num_Case"]<-"ncase"
# 	names(dat)[names(dat)=="Num_Control"]<-"ncontrol"
# 	dat$UKbiobank<-FALSE
# 	names(dat)[names(dat)=="Effect_Allele"]<-"Effect.Allele"
# 	names(dat)[names(dat)=="Reference_Allele"]<-"Other.Allele"
# 	names(dat)[names(dat)=="Beta"]<-"lnor"
# 	names(dat)[names(dat)=="Standard_error_of_beta"]<-"se"
# 	names(dat)[names(dat)=="Effect_Allele_Freq_Control"]<-"eaf"
# 	names(dat)[names(dat)=="P_value"]<-"p"
# 	names(dat)[names(dat)=="Info"]<-"info"
# 	dat$effect_allele_confirmed<-TRUE #from original column names
# 	return(dat)
# }


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

		names(Dat3)<-c("info","eaf","Other.Allele","Effect.Allele","lnor","lci","uci","se","Direction","eaf_test","allele_test")
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


format_data<-function(Dat=NULL,outcome=NA,population=NA,pmid=NA,study=NA,ncase=NA,ncontrol=NA,UKbiobank=NA,rsid=NA,Effect.Allele=NA,Other.Allele=NA,lnor=NA,se=NA,eaf=NA,p=NA,info=NA,info1=NA,info2=NA,info3=NA,HWEp=NA,phet=NA,I2=NA,Q=NA,Direction=NA,effect_allele_confirmed=FALSE,or=NA,lci=NA,uci=NA,ref=NULL,chr=NA,pos=NA,ID=NULL,test_statistic=NA,all_summary_stats=FALSE,summary_set="FAsnps",open_gwas=FALSE,efo=NA,efo_id=NA){

	# if(any(is.na(Dat[,rsid]) |  Dat[,rsid] ==".")) stop("rsid missing")
	if(rsid=="ID") { 
		names(Dat)[names(Dat) == rsid]<-"rsid"
		rsid<-"rsid"
	}
	if(!is.null(ref)){
		ref_dat<-read.table(ref,sep=" ",head=F,stringsAsFactors=F)		
		Chr<-unlist(strsplit(ref_dat$V1,split="chr"))
		Chr<-Chr[Chr!=""]
		ref_dat$chr<-Chr
		if(is.na(chr)){
			Chr<-unlist(strsplit(Dat[,rsid],split=":"))
			Dat$Chr<-as.numeric(Chr[seq(1,by=2,length(Chr))])
			Dat$Pos<-as.numeric(Chr[seq(2,by=2,length(Chr))])
			chr<-"Chr"
			pos<-"Pos"
		}

		Dat<-merge(Dat,ref_dat,by.x=c(chr,pos),by.y=c("chr","V2"),all.x=T)

		if(!is.na(rsid)){
			Dat$V4[is.na(Dat$V4)]<-Dat[is.na(Dat$V4),rsid] 
		}

		Dat<-Dat[!is.na(Dat$V4),]
		rsid<-"V4"		
		# Dat<-Dat[grep("rs",Dat$V4),]				
	}

	if(!is.null(ID)){
		Dat$ID <- ID
	}

	# sometimes odds ratio and confidence intervals are reported
	if(!is.na(or) & is.na(se) & !is.na(uci)){
		Dat$lnor<-log(Dat[,or])
		Dat$se<-(log(Dat[,uci])-log(Dat[,lci]))/(1.96*2)
		lnor<-"lnor"
	}

	# odds ratio and p value but no standard error or confidence intervals
	if(!is.na(or) & is.na(se) & is.na(uci)){
		Dat$lnor<-log(Dat[,or])
		Dat$z<-qnorm(Dat[,p]/2,lower.tail=F)
		Dat$se<-abs(Dat$lnor)/Dat$z
	}

	if(!is.na(lnor) & is.na(se) & is.na(uci)){
		Dat$z<-qnorm(Dat[,p]/2,lower.tail=F)
		Dat$se<-abs(Dat[,lnor])/Dat$z
	}

	# sometimes Odds ratio and standard error of log odds ratio are reported
	if(!is.na(or) & !is.na(se)){
		Dat$lnor<-log(as.numeric(Dat[,or]))
		lnor<-"lnor"
	}

	if(is.na(p)  & is.na(test_statistic)){
		Dat$test_statistic<-abs(as.numeric(Dat[,lnor])/as.numeric(Dat[,se]))
		Dat$p<-pnorm(Dat$test_statistic ,lower.tail=F)*2
	}	

	if(is.na(p) & !is.na(test_statistic)){
		Dat$p<-pnorm(Dat[,test_statistic] ,lower.tail=F)*2
    }



	Name_cols<-c("rsid","Effect.Allele","Other.Allele","lnor","se","eaf","p","info","info1","info2","info3","HWEp","phet","I2","Q","Direction","chr","pos","test_statistic")


	if(is.numeric(ncase)){
		Dat$ncase<-ncase
		Dat$ncontrol<-ncontrol
	}else{
		Name_cols<-c(Name_cols,"ncase","ncontrol")
	}
	Dat$pmid<-pmid
	Dat$outcome<-outcome
	Dat$population<-population
	Dat$study<-study
	Dat$UKbiobank<-UKbiobank
	Dat$effect_allele_confirmed<-effect_allele_confirmed

	for(i in 1:length(Name_cols)){
		# print(i)
		if(!is.na(eval(parse(text=Name_cols[i])))){
			names(Dat)[names(Dat) == eval(parse(text=Name_cols[i]))]<-Name_cols[i]
		}
	}
	

	Dat$p<-as.numeric(Dat$p)
	Dat$se<-as.numeric(Dat$se)
	Dat$lnor<-as.numeric(Dat$lnor)
	
	if(!is.na(eaf)){
		Dat$eaf<-as.numeric(Dat$eaf)
	}

	Dat$Effect.Allele<-toupper(Dat$Effect.Allele)
	if(!is.na(Other.Allele)){		
		Dat$Other.Allele<-toupper(Dat$Other.Allele)
	}
	# Dat<-Dat[!duplicated(Dat$rsid),]
	# drop duplicated SNPs

	# remove duplicate rows. 
	# There are two types of duplicates. Those where the rsid and the results for the rsid are duplicated and those where only the rsid is duplicated (i.e. results vary across duplicate rsids). We first deal with the duplicates where results are also duplicated, retaining one of the duplicate rsids. Then we deal with the duplicates where only the rsid is duplicated. for the latter we drop the rsid entirely, i.e. we don't retain one of the duplicates. 
	study_id_temp<-paste(Dat$rsid,Dat$Effect.Allele,Dat$Other.Allele,Dat$lnor,Dat$se)
	Dat<-Dat[!duplicated(study_id_temp),]

	# Dups<-unique(study_id_temp[duplicated(study_id_temp)])
	# Pos.dups<-which(study_id_temp %in% Dups)
	# Temp<-Dat[Pos.dups,]
	# head(Temp[order(Temp$rsid),])
	
	# Drop duplicates rsids. These seem to usually correspond to trialleic SNPS
	Dups<-unique(Dat$rsid[duplicated(Dat$rsid)])	
	Dat<-Dat[!Dat$rsid %in% Dups,]
	
	Dat<-Dat[!is.na(Dat$lnor) & !is.na(Dat$se),]
	if(any(is.na(Dat$se)) & is.na(p)){
		Dat<-Dat[!is.na(Dat$se), ]	
	}

	if(any(is.na(Dat$se)) & !is.na(p)){		
		Dat1<-Dat[is.na(Dat$se),]
		if(sum(Dat1$p)!=0) stop("infer missing SE from p value")
	}
	Dat<-Dat[!is.na(Dat$se),]
	Dat<-Dat[which(Dat$se != "Inf"),]
	Dat<-Dat[which(Dat$se != 0), ]

	Dat$summary_set<-summary_set
	Dat$all_summary_stats<-FALSE
	if(all_summary_stats){
		Dat$all_summary_stats<-TRUE
		Dat$summary_set <- "full summary stats"
	}

	Dat$open_gwas<-FALSE
	if(open_gwas){
		Dat$open_gwas<-TRUE
	}
	Dat$efo<-paste(efo,collapse="; ")
	Dat$efo_id<-paste(efo_id,collapse="; ")
	return(Dat)
}	



# remove duplicate rows (duplicated across all columns) and return vector of duplicate rsids (only rsid is duplicated) 
check_duplicates<-function(Dat=NULL){
	Dat<-Dat[order(Dat$rsid),]
	
	Dups2<-Dat$rsid[duplicated(Dat)]
	Dat<-Dat[!duplicated(Dat),] #remove duplicate rows
	Dups<-unique(Dat$rsid[duplicated(Dat$rsid)])
	return(list(Dups,Dat))
}


format_ref<-function(dat=NULL){
	chr<-unlist(strsplit(dat$V1,split="chr")) 
	dat$chr<-as.numeric(chr[seq(2,length(chr),2)])	
	return(dat)
}


load_plinkfrq<-function(File=NULL,population=NULL,study=NULL){
	ref<-readLines(File)
	Dat<-NULL
	for(i in 1:length(ref)){
		# print(i)
		A<-unlist(strsplit(ref[i],split=" "))
		A<-A[A!=""]
		Dat[[i]]<-A
	}
	Dat2<-data.frame(do.call(rbind,Dat),stringsAsFactors=F)
	names(Dat2)<-Dat2[1,]
	Dat2<-Dat2[2:nrow(Dat2),]
	Dat2$population<-population
	Dat2$study<-study
	Dat2$file<-File
	Dat2$MAF<-as.numeric(Dat2$MAF)
	names(Dat2)[names(Dat2)=="A1"]<-"minor_allele"
	names(Dat2)[names(Dat2)=="A2"]<-"major_allele"
	names(Dat2)[names(Dat2) == "MAF"]<-"maf"
	names(Dat2)[names(Dat2) == "SNP"]<-"rsid"
	
	Dat2<-make_refdat(Dat2)
	Dat2<-Dat2[,c("rsid","maf","minor_allele","major_allele","minor_allele2","major_allele2","population","study")]
	return(Dat2)
}


load_ref<-function(Dat=NULL,population=NULL,study=NULL){
	Dat$effect_allele_freq<-as.numeric(Dat$effect_allele_freq)
	Dat<-Dat[which(!is.na(Dat$effect_allele_freq)),]
	Pos<-which(Dat$effect_allele_freq>0.5)
	# Beta<-Dat$beta[Pos]*-1
	# Dat$beta[Pos]<-Beta
	Dat$maf <- Dat$effect_allele_freq
	Dat$maf[Pos]<-1-Dat$effect_allele_freq[Pos]
	Dat$minor_allele<-Dat$effect_allele
	Dat$major_allele<-Dat$other_allele
	maj<-Dat$minor_allele[Pos]
	min<-Dat$major_allele[Pos]
	Dat$major_allele[Pos]<-maj
	Dat$minor_allele[Pos]<-min	
	Dat$population <- population
	Dat$study<-study

	# names(Dat)[names(Dat) == "effect_allele"]<-"Effect.Allele"
	# names(Dat)[names(Dat) == "other_allele"]<-"Other.Allele"
	# names(Dat)[names(Dat) == "effect_allele_freq"]<-"eaf"
	names(Dat)[names(Dat) == "snp"]<-"rsid"

	Dat<-make_refdat(Dat)
	Dat<-Dat[,c("rsid","maf","minor_allele","major_allele","minor_allele2","major_allele2","population","study")]
	return(Dat)
}

# load_charge_ref<-function(File=NULL){
# 	load(File)	
# 	Charge$effect_allele_freq<-as.numeric(Charge$effect_allele_freq)
# 	Charge<-Charge[which(!is.na(Charge$effect_allele_freq)),]
# 	Pos<-which(Charge$effect_allele_freq>0.5)
# 	# Beta<-Charge$beta[Pos]*-1
# 	# Charge$beta[Pos]<-Beta
# 	Charge$maf <- Charge$effect_allele_freq
# 	Charge$maf[Pos]<-1-Charge$effect_allele_freq[Pos]
# 	Charge$minor_allele<-Charge$effect_allele
# 	Charge$major_allele<-Charge$other_allele
# 	maj<-Charge$minor_allele[Pos]
# 	min<-Charge$major_allele[Pos]
# 	Charge$major_allele[Pos]<-maj
# 	Charge$minor_allele[Pos]<-min	

# 	# names(Charge)[names(Charge) == "effect_allele"]<-"Effect.Allele"
# 	# names(Charge)[names(Charge) == "other_allele"]<-"Other.Allele"
# 	# names(Charge)[names(Charge) == "effect_allele_freq"]<-"eaf"
# 	names(Charge)[names(Charge) == "snp"]<-"rsid"

# 	Charge<-make_refdat(Charge)
# 	Charge<-Charge[,c("rsid","maf","minor_allele","major_allele","minor_allele2","major_allele2")]
# 	return(Charge)
# }


make_refdat<-function(ref.dat=NULL){	
	strand1<-c("C","G","A","T")
	strand2<-c("G","C","T","A")
	ref.dat$minor_allele2<-strand2[match(ref.dat$minor_allele,strand1)] 
	ref.dat$major_allele2<-strand2[match(ref.dat$major_allele,strand1)] 	
	return(ref.dat)
}

flip_strand2<-function(dat=NULL,allele1=NULL,allele2=NULL){	
	strand1<-c("C","G","A","T")
	strand2<-c("G","C","T","A")
	dat[,allele1]<-strand2[match(dat[,allele1],strand1)]
	dat[,allele2]<-strand2[match(dat[,allele2],strand1)]
	return(dat)
}



format_refdat_1000G_superpops<-function(refdat_1000G_superpops=NULL,study=NULL){
		names(refdat_1000G_superpops)[names(refdat_1000G_superpops) == "SNP"]<-"rsid"
		names(refdat_1000G_superpops)[names(refdat_1000G_superpops) == "MAF"]<-"maf"
		refdat_1000G_superpops$study<-study
		return(refdat_1000G_superpops)
}



make_snplist2<-function(trait=NULL,efo_id=NULL,efo=NULL,population=NULL,Dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/",bed_37=FALSE){

	if(!is.null(efo_id)){
		top_hits<-gwas_catalog_hits(efo_id=efo_id)	
	}

	if(!is.null(trait)){
		top_hits<-gwas_catalog_hits(trait=trait)	
	}
	if(!is.null(efo)){
		top_hits<-gwas_catalog_hits(efo=efo)
	}
	
	top_hits_rsids<-top_hits$rsid	
	
	# snplist<-c(snplist,top_hits_rsids)
	load(paste0(Dir,"refdat_1000G_superpops.Rdata"))
	snplist<-c(top_hits_rsids,unique(refdat_1000G_superpops$SNP))	
	snplist<-unique(snplist)		
	
	# if(bed_37){
	# 	Ref<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",file_sep=" ",exact_match=TRUE,Head=FALSE)
	# 	Chr<-Ref$V1
	# 	Chr<-gsub("chr","",Chr)
	# 	bp<-Ref$V2
	# 	snplist1<-Ref$V4

	# 	snpbed37<-paste(Chr,bp,sep=":")			
	# 	snplist_dat<-data.frame(matrix(c(snplist1,Chr,bp,snpbed37),ncol=4,nrow=length(snplist1)),stringsAsFactors=F)
	# 	names(snplist_dat)<-c("rsid","chr","pos37","bed")

	# 	# read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",sep=" ",head=F,stringsAsFactors=F)
	# 	return(snplist_dat)
	# }else{
	return(snplist)
	# }	
}

make_snplist<-function(trait=NULL,efo=NULL,population=NULL,Dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/",bed_37=FALSE){

	
	if(population == "European"){
		snplist<-readLines(paste0(Dir,"snplist_Europeans_rsidsonly2.txt"))
	}
	if(population == "East Asian"){
		snplist<-readLines(paste0(Dir,"snplist_East_Asians_rsidsonly2_nodups.txt"))
		snplist<-snplist[snplist!="x"]
	}
	
	
	# if(population == "East Asian" & bed_37){
	# 	snplist<-readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_East_Asians_rsidsonly2_nodups.txt")
	# 	snplist<-snplist[snplist!="x"]
	# }

	if(!is.null(trait)){
		top_hits<-gwas_catalog_hits(trait=trait)	
	}
	if(!is.null(efo)){
		top_hits<-gwas_catalog_hits(efo=efo)
	}
	
	top_hits_rsids<-top_hits$rsid	
	
	snplist<-c(snplist,top_hits_rsids)
	load(paste0(Dir,"refdat_1000G_superpops.Rdata"))
	snplist<-c(snplist,unique(refdat_1000G_superpops$SNP))	
	snplist<-unique(snplist)		
	
	if(bed_37){
		Ref<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",file_sep=" ",exact_match=TRUE,Head=FALSE)
		Chr<-Ref$V1
		Chr<-gsub("chr","",Chr)
		bp<-Ref$V2
		snplist1<-Ref$V4

		snpbed37<-paste(Chr,bp,sep=":")			
		snplist_dat<-data.frame(matrix(c(snplist1,Chr,bp,snpbed37),ncol=4,nrow=length(snplist1)),stringsAsFactors=F)
		names(snplist_dat)<-c("rsid","chr","pos37","bed")

		# read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",sep=" ",head=F,stringsAsFactors=F)
		return(snplist_dat)
	}else{
		return(snplist)
	}	
}

# gwas_catalog_hits2<-function(){
#     library(gwasrapidd)
#     library(purrr)
#     library(dplyr)

# 	gwas_studies<-gwasrapidd::get_studies(efo_trait = efo)
# 	if(nrow(gwas_studies@studies)==0){
# 		warning(paste("search for efo -",efo,"- returned 0 studies from the GWAS catalog"))
# 	}
# 	# studies <- get_studies(reported_trait = "colorectal cancer")
# 	if(nrow(gwas_studies@studies)!=0){
# 		# ancestry_tab<-make_ancestry_table(gwas_studies=gwas_studies)		
	
# 	study_ids<-gwas_studies@studies$study_id	
		


#     # study_ids <- studies@studies$study_id
#     names(study_ids) <- study_ids
#     associations <-
#       purrr::map(study_ids, ~ get_associations(study_id = .x))
#     study2association <-
#       purrr::imap_dfr(
#         associations,
#         ~ tibble::tibble(
#           study_id = .y,
#           association_id = .x@associations$association_id,

#         )
#       )

#     ancestries <-
#       dplyr::left_join(gwas_studies@ancestries,
#                        gwas_studies@ancestral_groups,
#                        by = c('study_id', 'ancestry_id')) %>%
#       dplyr::left_join(gwas_studies@countries_of_origin, by = c('study_id', 'ancestry_id')) %>%
#       dplyr::rename(
#         co_country_name = country_name,
#         co_major_area = major_area,
#         co_region = region
#       ) %>%
#       dplyr::left_join(gwas_studies@countries_of_recruitment,
#                        by = c('study_id', 'ancestry_id')) %>%
#       dplyr::rename(
#         cr_country_name = country_name,
#         cr_major_area = major_area,
#         cr_region = region
#       )

#     (study_assoc_ancestry <-
#         dplyr::left_join(study2association, ancestries, by = c('study_id')))

gwas_catalog_hits<-function(trait=NULL,efo=NULL,efo_id=NULL){
		
	if(!is.null(efo)){
		efo<-trimws(unlist(strsplit(efo,split=";")))	
		gwas_studies<-gwasrapidd::get_studies(efo_trait = efo)		
		# unique(gwas_studies@studies$reported_trait)
		if(nrow(gwas_studies@studies)==0){
			warning(paste("search for efo -",efo,"- returned 0 studies from the GWAS catalog"))
		}
	}

	if(!is.null(efo_id)){
		efo_id<-trimws(unlist(strsplit(efo_id,split=";")))	
		gwas_studies<-gwasrapidd::get_studies(efo_id = efo_id)		
		# unique(gwas_studies@studies$reported_trait)
		if(nrow(gwas_studies@studies)==0){
			warning(paste("search for efo -",efo_id,"- returned 0 studies from the GWAS catalog"))
		}
	}
	
	if(!is.null(trait)){
		gwas_studies<-gwasrapidd::get_studies(reported_trait = trait)
		if(nrow(gwas_studies@studies)==0){
			warning(paste("search for trait -",trait,"- returned 0 studies from the GWAS catalog"))
		}
	}
	
	if(nrow(gwas_studies@studies)!=0){
		ancestry_tab<-make_ancestry_table(gwas_studies=gwas_studies)		
		study_ids<-gwas_studies@studies$study_id	
		Dat<-NULL	
		for(i in 1:length(study_ids)){		
		# print(i)			
			gwas_associations<-gwasrapidd::get_associations(study_id = study_ids[i])
			if(nrow(gwas_associations@associations)!=0){
				associations<-data.frame(gwas_associations@associations,stringsAsFactors=F)
				risk_alleles<-data.frame(gwas_associations@risk_alleles)
				gwas_results<-merge(associations,risk_alleles,by="association_id")
				# p_values<-gwas_associations@associations$pvalue			
				# odds_ratios<-gwas_associations@associations$or_per_copy_number
				# risk_alleles<-gwas_associations@risk_alleles$risk_allele
				# rsids<-gwas_associations@risk_alleles$variant_id
				# eaf<-gwas_associations@risk_alleles$risk_frequency
			
				gwas_results$z_scores<-qnorm(gwas_results$pvalue/2,lower.tail=F)
				gwas_results$log_odds_ratios<-log(gwas_results$or_per_copy_number)
				all(is.na(gwas_results$log_odds_ratios))
				gwas_results$log_odds_ratios<-gwas_results$beta_number
				gwas_results$standard_errors<-gwas_results$log_odds_ratios/gwas_results$z_scores
				gwas_results$study_id<-study_ids[i]								
				Dat[[i]]<-gwas_results	
			}			
		}
		if(!is.null(trait)) trait_efo<-trait
		if(!is.null(efo)) trait_efo<-efo

		if(is.null(Dat)){
			warning(paste0("no results found in GWAS catalog for ",trait_efo))
		}
		if(!is.null(Dat)){
			Dat<-do.call(rbind,Dat)
			Dat<-Dat[,c("variant_id","risk_allele","or_per_copy_number", "log_odds_ratios","standard_errors","risk_frequency","pvalue","z_scores","study_id")]
			names(Dat)<-c("rsid","Effect.Allele","or","lnor","se","eaf","p","test_statistic","study_id")	
			Pos<-is.na(Dat$lnor)
			Dat$lnor[Pos]<-log(Dat$or[Pos])
			se<-Dat$lnor/Dat$test_statistic
			Pos<-is.na(Dat$se)
			Dat$se[Pos]<-se[Pos]
			Dat<-merge(Dat,ancestry_tab,by="study_id")				
			return(Dat)
		table(Dat$study_id,Dat$ancestral_group)
		}
	}			
}

make_ancestry_table<-function(gwas_studies=NULL){
	ancestry_tab<-data.frame(gwas_studies@ancestral_groups,stringsAsFactors=F)		
	ancestry_tab<-unique(ancestry_tab[,c("study_id","ancestral_group")])
	# ancestry_tab<-ancestry_tab[!is.na(ancestry_tab$ancestral_group),]
	study_ids<-unique(ancestry_tab$study_id)
	anc2<-NULL
	for(i in 1:length(study_ids)){
		anc1<-ancestry_tab[ancestry_tab$study_id == study_ids[i],]
		anc1$ancestral_group<-paste(anc1$ancestral_group,collapse="; ")		
		anc2[[i]]<-anc1[1,]
	}
	anc2<-do.call(rbind,anc2)
	return(anc2)
}

# library(gwasrapidd) 
	# gwasrapidd::get_variants())
	# gwasres<-get_studies(efo_trait = "glioma")
	# gwas_studies<-gwasrapidd::get_studies(efo_trait = efo)
	# gwas_studies@studies$reported_trait
	
# 	gwas_studies@ancestral_groups$ancestral_group
# 	gwas_studies@ancestral_groups$study_id
# 	gwas_studies@publications$pubmed_id
# 	gwas_associations@associations$study_id
# 	gwas_studies@studies$study_id
# # gwas_associations<-gwasrapidd::get_associations(efo_trait = efo)	
# gwas_associations@associations$association_id


# gwasres<-get_studies(efo_trait = "glioma")
# gwasres@studies$reported_trait
# gwasres@ancestries
# gwasres@ancestral_groups
# gwasres@countries_of_origin
# gwasres@countries_of_recruitment
# gwasres@publications


# gwasres<-get_variants(efo_trait = 'glioma')
# gwasres@variants
# gwasres@ensembl_ids
# gwasres@genomic_contexts
# gwasres@entrez_ids
# head(gwasres)





# find_allele_errors_eastasian<-function(dat=NULL,population=NULL){	
# 	ref.dat_eas<-load_plinkfrq(File<-"~/fatty-acids/outcome_data/data/fatty_acid_snps_eas.frq",population="EAS1")
# 	load("~/fatty-acids/outcome_data/data/ref_dat_schs.RData")
# 	ref.dat_schs<-load_ref(Dat=SCHS,population="EAS2")
# 	ref.dat<-plyr::rbind.fill(ref.dat_eas,ref.dat_schs)

# 	Dat_list<-find_allele_error(target_dat=dat,ref_dat=ref.dat,snp="rsid",eaf="eaf",ref_dat_maf="maf",target_dat_effect_allele = "Effect.Allele",target_dat_other_allele= "Other.Allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",target_study="study",ID="ID",ref_study="study")
# 	# Dat_list_schs<-find_allele_error(target_dat=dat,ref_dat=ref.dat_schs,snp="rsid",eaf="eaf",ref_dat_maf="maf",target_dat_effect_allele = "Effect.Allele",target_dat_other_allele= "Other.Allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",target_study="study",ID="ID",ref_study="study")

# 	png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/fasnps_",unlist(Dat_list[3]),".png"))
# 		print(Dat_list[1])
# 	dev.off()	

# 	# png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/schs_",unlist(Dat_list_schs[3]),".png"))
# 		# print(Dat_list_schs[1])
# 	# dev.off()

# 	allele_errors<-unlist(Dat_list[2])

# 	# write.table(allele_errors,paste0("~/fatty-acids/outcome_data/results/snps_allele_errors_",Dat_list_eas[3],".txt"),col.names=T,row.names=F,quote=F,sep="\t")

# 	return(list(Dat_list_eur[2],Dat_list_charge[2]))
# }

# find_allele_errors_european<-function(dat=NULL,population=NULL){	
# 	ref.dat_eur<-load_plinkfrq(File<-"~/fatty-acids/outcome_data/data/fatty_acid_snps_eur.frq",population="EUR1",study="1000G")	
# 	load("~/fatty-acids/outcome_data/data/ref_dat_charge_imputed.RData")
# 	ref.dat_charge<-load_ref(Dat=Charge,population="EUR2",study="CHARGE")
# 	ref.dat<-plyr::rbind.fill(ref.dat_eur,ref.dat_charge)
	
# 	Dat_list<-find_allele_error(target_dat=dat,ref_dat=ref.dat,snp="rsid",eaf="eaf",ref_dat_maf="maf",target_dat_effect_allele = "Effect.Allele",target_dat_other_allele= "Other.Allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",target_study="study",ID="ID",ref_study="study")
	
# 	out_file<-unlist(Dat_list[3])
# 	out_file<-gsub("/","_",	out_file)

# 	png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/fasnps_",out_file,".png"))
# 		print(Dat_list[1])
# 	dev.off()	

	
# 	allele_errors<-unlist(Dat_list[2])
		
# 	# write.table(allele_errors,paste0("~/fatty-acids/outcome_data/results/snps_allele_errors_",out_file,".txt"),col.names=T,row.names=F,quote=F,sep="\t")

# 	return(list(Dat_list_eur[2],Dat_list_charge[2]))
# }


# find_allele_errors_superpops<-function(dat=NULL){	
# 	load("~/fatty-acids/outcome_data/data/refdat_1000G_superpops.Rdata")
# 	ref.dat<-format_refdat_1000G_superpops(study="1000G")
# 	Dat_list<-find_allele_error(target_dat=dat,ref_dat=ref.dat,snp="rsid",eaf="eaf",ref_dat_maf="MAF",target_dat_effect_allele = "Effect.Allele",target_dat_other_allele= "Other.Allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",study="study",ID="ID",target_dat_population="population",ref_dat_population="population",target_study="study",ref_study="study")
		
# 	png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/superpops_",unlist(Dat_list[3]),".png"))
# 			print(Dat_list[1])
# 		dev.off()

# 	# png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/1000GEUR_",out_file,".png"))
# 	# 	print(Dat_list_eur[1])
# 	# dev.off()	

# 	# png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/charge_",out_file,".png"))
# 	# 	print(Dat_list_charge[1])
# 	# dev.off()
# 	allele_errors<-unlist(Dat_list[2])
		
# 	# write.table(allele_errors,paste0("~/fatty-acids/outcome_data/results/snps_allele_errors_",out_file,".txt"),col.names=T,row.names=F,quote=F,sep="\t")
# 	return(Dat_list[2])
# }



preformat_all_id4<-function(){ 
	SNPlist<-read.table("~/fatty-acids/outcome_data/data/snplist_Europeans.txt",sep="\t",head=T,stringsAsFactors=F)
	SNPlist1<-SNPlist[SNPlist$type=="lead variant",]
	SNPlist2<-SNPlist[SNPlist$type!="lead variant",]
	SNPlist<-rbind(SNPlist1,SNPlist2)
	SNPlist<-SNPlist[!duplicated(SNPlist$SNP),]

	Lead<-unique(SNPlist$SNP[SNPlist$type=="lead variant"])
	Proxies<-unique(SNPlist$SNP[SNPlist$type=="proxy"])
	Alias<-unique(SNPlist$SNP[SNPlist$type=="alias rsid"])
	SNPs_all<-unique(SNPlist$SNP)
	Bed_38<-paste(SNPlist$Chr,SNPlist$SNP.GRCh38.p12,sep="-")
	Bed_37<-paste(SNPlist$Chr,SNPlist$SNP.GRCh37,sep="-")
	# Bed_37<-Bed_37[grep("NA",Bed_37,invert=T)]
	SNPlist$Bed_37<-Bed_37
	SNPlist$Bed_37[grep("NA",SNPlist$Bed_37)]<-NA
	SNPlist<-SNPlist[!duplicated(SNPlist$Bed_37),]
	# dups<-SNPlist$Bed_37[which(duplicated(SNPlist$Bed_37))]
	# dups<-dups[!is.na(dups)]
	# Temp<-SNPlist[SNPlist$Bed_37 %in% dups,]
	# Temp[order(Temp$Bed_37),]
	# work out effect allele information
	Files<-c("UKGWASI.txt","UKGWASII.txt","GERGWAS.txt")
	L_dat<-NULL
	for(i in 1:length(Files)){
		print(i)
		print(Files[i])
		Eaf<-read.table(paste("bcell_ALL_29632299/",Files[i],sep=""),sep="\t",head=T,stringsAsFactors=F)
		Eaf<-Eaf[!duplicated(Eaf$rsid),]
		Eaf$bed<-paste(Eaf$chromosome,Eaf$pos,sep="-")
		Eaf<-Eaf[!duplicated(Eaf$bed),] #exclude duplicated chromosome positions. For example "10-102315709" corresponds to rsid rs71891071 and rsid 10-102315709 . The minor allele seems to have opposite directions of effect
		# Eaf[Eaf$rsid== "10-102315709",]
		# Eaf[Eaf$bed== "10-102315709",]
		# dups<-Eaf$bed[duplicated(Eaf$bed)]
		# Temp<-Eaf[Eaf$bed %in% dups,]
		# # Temp<-Temp[order(Temp$bed),c("rsid","bed","pheno_caco_frequentist_add_cov_01_cov_02_score_beta_1")]
		# Names<-names(table(Temp$bed))
		# N<-as.numeric(table(Temp$bed))
		# Dat<-data.frame(matrix(c(Names,N),nrow=length(Names),ncol=2,byrow=FALSE),stringsAsFactors=F)
		# Temp<-merge(Temp,Dat,by.x="bed",by.y="X1")
		# which(Eaf$bed == "10-65104500")
		# bcell_all[bcell_all$rsid == "rs7896518",]
		Eaf1<-merge(Eaf,unique(SNPlist[,c("SNP","Bed_37")]),by.x="bed",by.y="Bed_37",all.x=T)
		Eaf2<-Eaf1[!is.na(Eaf1$SNP),]
		Eaf3<-Eaf1[is.na(Eaf1$SNP),]
		Eaf4<-Eaf3[grep("rs",Eaf3$rsid),]
		Eaf5<-rbind(Eaf2,Eaf4)
		Eaf5$SNP[is.na(Eaf5$SNP)]<-Eaf5$rsid[is.na(Eaf5$SNP)]
		AA<-Eaf5$controls_AA
		AB<-Eaf5$controls_AB
		BB<-Eaf5$controls_BB
		Eaf5$freqB<-(AB+2*BB)/((AA+AB+BB)*2)
		Eaf5$Ncontrols<-AA+AB+BB
		Eaf5$cohort<-Files[i]
		head(Eaf5)
		# Eaf5[Eaf5$bed == ""10-65104500"",]
		L_dat[[i]]<-Eaf5
	}

	Eaf<-do.call(plyr::rbind.fill,L_dat)

	ukgwas1<-Eaf[Eaf$cohort=="UKGWASI.txt",c("SNP","freqB","allele_B","Ncontrols")]
	ukgwas2<-Eaf[Eaf$cohort=="UKGWASII.txt",c("SNP","freqB","allele_B","Ncontrols")]
	gergwas<-Eaf[Eaf$cohort=="GERGWAS.txt",c("SNP","freqB","allele_B","Ncontrols")]
	names(ukgwas1)<-c("SNP","freqB_ukgwas1","allele_B_ukgwas1","controls_ukgwas1")
	names(ukgwas2)<-c("SNP","freqB_ukgwas2","allele_B_ukgwas2","controls_ukgwas2")
	names(gergwas)<-c("SNP","freqB_gergwas","allele_B_gergwas","controls_gergwas")

	bcell_all<-read.table("bcell_ALL_29632299/bcell_all.txt",sep="\t",head=T,stringsAsFactors=F)
	head(bcell_all)
	# bcell_all[bcell_all$rsid == "rs10056811",]
	bcell_all<-bcell_all[!duplicated(bcell_all$rsid),]
	head(bcell_all)
	bcell_all$bed<-paste(bcell_all$chr,bcell_all$pos,sep="-")
	bcell_all<-bcell_all[!duplicated(bcell_all$bed),] 

	bcell_all1<-merge(bcell_all,unique(SNPlist[,c("SNP","Bed_37")]),by.x="bed",by.y="Bed_37",all.x=T)
	head(SNPlist)
	bcell_all2<-bcell_all1[!is.na(bcell_all1$SNP),]
	bcell_all3<-bcell_all1[is.na(bcell_all1$SNP),]
	bcell_all4<-bcell_all3[grep("rs",bcell_all3$rsid),]
	bcell_all5<-rbind(bcell_all2,bcell_all4)
	bcell_all5$SNP[is.na(bcell_all5$SNP)]<-bcell_all5$rsid[is.na(bcell_all5$SNP)]

	bcell_all5<-merge(bcell_all5,ukgwas1,by="SNP",all.x=T)
	bcell_all5<-merge(bcell_all5,ukgwas2,by="SNP",all.x=T)
	bcell_all5<-merge(bcell_all5,gergwas,by="SNP",all.x=T)
	# weighted average of eaf 
	bcell_all5$eafw_g<-bcell_all5$freqB_gergwas * bcell_all5$controls_gergwas
	bcell_all5$eafw_uk1<-bcell_all5$freqB_ukgwas1 * bcell_all5$controls_ukgwas1
	bcell_all5$eafw_uk2<-bcell_all5$freqB_ukgwas2 * bcell_all5$controls_ukgwas2
	bcell_all5$Total<-rowSums(bcell_all5[,c("controls_gergwas","controls_ukgwas1","controls_ukgwas2")],na.rm=T)
	bcell_all5$eafw_total<-rowSums(bcell_all5[,c("eafw_g","eafw_uk1","eafw_uk2")],na.rm=T)
	bcell_all5$eaf<-bcell_all5$eafw_total/bcell_all5$Total
	return(bcell_all5)
}


# Re. 1, 2: The effect refers to the minor allele, as does the frequency. The minor allele could be either A1 or A2. I have added a new column A1frq here, to help determine which of A1 and A2 is the minor. There can be uncertainties when the MAF is close to 0.5 so you might want to exclude those variants.
# Re. 3: (Table 1 in Swaminathan) Re. 4: Yes, that should be the case.
# Re. 5: No, I can only share the . Re 6: Yes, Swaminathan is the original ref.
# All the best,

preformat_mma_id96<-function(){
	mma<-read.table("~/MR_FattyAcids/data/summary_data/multiple_myeloma_26007630/PhilipHaycock_2019-06-03.txt",sep="\t",head=T,stringsAsFactors=F)
	Pos1<-mma$A1frq<0.5
	# mma[mma$rsid=="rs174546",]

	mma$Effect.Allele<-mma$A2
	mma$Effect.Allele[Pos1]<-mma$A1[Pos1]
	mma<-mma[mma$frq<0.46,]
	mma$Other.Allele<-mma$A1 
	mma$Other.Allele[mma$Other.Allele == mma$Effect.Allele]<-mma$A2[mma$Other.Allele == mma$Effect.Allele] 
	return(mma)
}

preformat_cer_id129<-function(){
	cer<-read.table("~/MR_FattyAcids/data/summary_data/cervical cancer/cervicalcancer.txt",sep="\t",head=T,stringsAsFactors=F)
	maf<-read.table("~/MR_FattyAcids/data/summary_data/cervical cancer/data_gwasLookupDataInclA2_complete_inclMAF.csv",head=T,stringsAsFactors=F,sep=";")
	cer.m<-merge(Cer,maf[,c("SNP","MAF")],by="SNP")
	cer.m$MAF<-as.numeric(gsub(",",".",cer.m$MAF))
	return(cer.m)
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


preformat_neu_id107<-function(){
	neu<-read.table("~/MR_FattyAcids/data/summary_data/neuroblastoma_28545128/NBL_fatty_acid_cancer_McDaniel_etal.txt",sep="\t",head=T,stringsAsFactors=F)
	AA<-neu$controls_AA 
	AB<-neu$controls_AB 
	BB<-neu$controls_BB
	# summary(AA+AB+BB)
	Nchr<-(AA+AB+BB)*2
	N_b_alleles<-AB+BB*2
	neu$eaf<-N_b_alleles/Nchr
	return(neu)
}

# format_data_ukb<-function(ukb=NULL){

# 	Disease<-c("Brain cancer",
# 		"Breast cancer",	
# 		"Colorectal cancer",
# 		"Blood cancer",
# 		"Oral cavity and pharyngeal cancer",
# 		"Leukaemia",
# 		"Liver & bile duct cancer",
# 		"Liver cancer",
# 		"Lung cancer",
# 		"Lung cancer unadjusted for chip",	
# 		"Lymphoid leukaemia",
# 		"Multiple myeloma",
# 		"Myeloid Leukaemia",
# 		"Non-melanoma skin cancer",
# 		"Esophageal cancer",
# 		"Ovarian cancer",
# 		"Cancer (excluding non-melanoma skin cancer)",
# 		"Cancer (all cause)",
# 		"Prostate cancer",	
# 		"Melanoma") #file is called overall skin cancer but Kim says this is melanoma only		


# 	IDS<-c(138,
# 		139,
# 		143,
# 		137,
# 		157,
# 		146,
# 		147,
# 		148,
# 		149,
# 		1499,
# 		150,
# 		154,
# 		155,
# 		156,
# 		144,
# 		158,
# 		141,
# 		140,
# 		159,
# 		153)

# 	EFOS<-list(
# 		"central nervous system cancer",
# 		"breast carcinoma",
# 		"colorectal cancer",
# 		c("lymphoma","multiple myeloma","lymphoid leukemia"),
# 		c("head and neck squamous cell carcinoma","oropharynx cancer","nasopharyngeal neoplasm","hypopharynx cancer","oral cavity cancer","mouth neoplasm","pharynx cancer"),
# 		"acute lymphoblastic leukemia",
# 		"hepatocellular carcinoma",
# 		"hepatocellular carcinoma",
# 		"lung carcinoma",
# 		"lung carcinoma",
# 		"acute lymphoblastic leukemia",
# 		"multiple myeloma",
# 		"acute myeloid leukemia",
# 		"non-melanoma skin carcinoma",
# 		"esophageal adenocarcinoma",
# 		"ovarian carcinoma",
# 		"cancer",
# 		"cancer",
# 		"prostate carcinoma",
# 		"melanoma")
# 	dat_list<-NULL	
# 	for(i in 1:length(ukb)){
# 		dat<-data.frame(ukb[i],stringsAsFactors=F)
# 		outcome=Disease[i]
# 		ID<-IDS[i]
# 		EFO<-unlist(EFOS[i])
# 		dat<-format_data(Dat=dat,outcome=outcome,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",ID=ID,all_summary_stats=TRUE,efo=EFO)
# 		dat_list[[i]]<-dat
# 	}
# 	dat<-do.call(rbind,dat_list)
# 	return(dat)
# }

# # zcat /projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_haem_cancer_imputed.txt.gz | head
# # cd /projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results

# extract_snps_ukb_cancer<-function(snplist=NULL){
# 	setwd("/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results")
# 	Files<-dir()
# 	Files<-Files[grep("maf_overall",Files)]
# 	ukb<-NULL
# 	for(i in 1:length(snplist)){
# 	# for(i in 4){
# 		print(i)
# 		snps<-unlist(snplist[i])
# 		ukb[[i]]<-extract_snps(snplist=snps,File=paste0("/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/",Files[i]),exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
# 	}
# 	return(ukb)
# }

# make_snplist_ukb_cancer<-function(){

# 	# Disease<-c(
# 	# 	"Brain cancer",
# 	# 	"Breast carcinoma",	
# 	# 	"Colorectal cancer",
# 	# 	"Blood cancer",
# 	# 	"Oral cavity and pharyngeal cancer",
# 	# 	"Leukaemia",
# 	# 	"Liver & bile duct cancer",
# 	# 	"Liver cancer",
# 	# 	"Lung cancer",
# 	# 	"Lung cancer unadjusted for chip",	
# 	# 	"Lymphoid leukaemia",
# 	# 	"Multiple myeloma",
# 	# 	"Myeloid Leukaemia",
# 	# 	"Non-melanoma skin cancer",
# 	# 	"Esophageal cancer",
# 	# 	"Ovarian cancer",
# 	# 	"Cancer (excluding non-melanoma skin cancer)",
# 	# 	"Cancer (all cause)",
# 	# 	"Prostate cancer",	
# 	# 	"Melanoma") #file is called overall skin cancer but Kim says this is melanoma only		

# 	snplist1<-make_snplist(efo= "central nervous system cancer",population="European")	
# 	snplist2<-make_snplist(efo="breast carcinoma",population="European")	
# 	snplist3<-make_snplist(efo="colorectal cancer",population="European")	
# 	snplist4<-make_snplist_blood(population="European",Dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/")
# 	snplist5<-make_snplist_hnc(population="European",Dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/")
# 	# snplist6<-make_snplist(trait=Disease[6],population="European")	
# 	snplist6<-make_snplist(efo="acute lymphoblastic leukemia",population="European")	
# 	# snplist7<-make_snplist(trait=Disease[7],population="European")
# 	snplist7<-make_snplist(efo="hepatocellular carcinoma",population="European")	
# 	# snplist8<-make_snplist(trait=Disease[8],population="European")	
# 	snplist8<-make_snplist(efo="hepatocellular carcinoma",population="European")	
# 	snplist9<-make_snplist(efo="lung carcinoma",population="European")	
# 	snplist10<-make_snplist(efo="Lung carcinoma",population="European")	
# 	snplist11<-make_snplist(efo="acute lymphoblastic leukemia",population="European")	
# 	snplist12<-make_snplist(efo="multiple myeloma",population="European")	
# 	snplist13<-make_snplist(efo="acute myeloid leukemia",population="European")	
# 	snplist14<-make_snplist(efo="non-melanoma skin carcinoma",population="European")	
# 	snplist15<-make_snplist(efo="esophageal adenocarcinoma",population="European")	
# 	snplist16<-make_snplist(efo="ovarian carcinoma",population="European")	
# 	snplist17<-make_snplist(efo="cancer",population="European")	
# 	snplist18<-make_snplist(efo="cancer",population="European")	
# 	snplist19<-make_snplist(efo="prostate carcinoma",population="European")	
# 	snplist20<-make_snplist(efo="melanoma",population="European")	
	
# 	snp_list<-ls()[grep("snplist[0-9]",ls())] 

# 	# Test<-unlist(lapply(1:length(snp_list),FUN=function(x) 
# 		# !is.null(nrow(eval(parse(text=snp_list[x]))))))
# 	# rm(list=snp_list[!Test]) #remove the objects with no data
# 	# ,envir=.GlobalEnv
	
# 	snplist<-lapply(1:length(snp_list),FUN=function(x) eval(parse(text=snp_list[x])))
# 	return(snplist)
# }


# format_data_ugi<-function(ugi=NULL){
# 	pmids<-c(26129866,25129146,26129866,26129866)
# 	Disease<-c("Gastric cardia adenocarcinoma","Esophageal squamous cell carcinoma","Gastric adenocarcinoma","Noncardia gastric adenocarcinoma")
# 	IDS<-c(102,99,101,103)
# 	cases<-c(1189,2013,2350,1027)
# 	controls<-c(2708,2701,2708,2708)
# 	EFOS<-c("Gastric adenocarcinoma","Esophageal squamous cell carcinoma","Gastric adenocarcinoma","Gastric adenocarcinoma")
# 	dat<-NULL
# 	for(i in 1:length(ugi)){
# 		print(i)
# 		dat[[i]]<-format_data(Dat=data.frame(ugi[i],stringsAsFactors=F),outcome=Disease[i],population="East Asian",pmid=pmids[i],ncase=cases[i],ncontrol=controls[i],study="N-UGC",UKbiobank=FALSE,rsid="rs",Effect.Allele="risk_allele",Other.Allele="reference_allele",lnor="beta",se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info="info",ID=IDS[i],all_summary_stats=TRUE,efo=EFOS[i])  
# 	}
# 	dat<-do.call(rbind,dat)
# 	return(dat)
# }


# extract_snps_ugi<-function(snplist=NULL){
# 	Files<-c("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_CC/summary_chr_all.txt","/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/summary_chr_all.txt","/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_gastric/summary_chr_all.txt","/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_NC/summary_chr_all.txt")

# 	ugi<-NULL
# 	for(i in 1:length(snplist)){
# 		print(i)
# 		snps<-unlist(snplist[i])
# 		ugi[[i]]<-extract_snps(snplist=snps,File=Files[i],exact_match=TRUE,file_sep="\t",Comment="")
# 	}
# 	return(ugi)
# }

# make_snplist_ugi<-function(){
# 	# Disease<-c("Gastric cardia adenocarcinoma","Esophageal squamous cell carcinoma","Gastric adenocarcinoma","Noncardia gastric adenocarcinoma")
# 	snplist1<-make_snplist(efo="Gastric adenocarcinoma",population="East Asian")	
# 	snplist2<-make_snplist(efo="Esophageal squamous cell carcinoma",population="East Asian")	
# 	snplist3<-make_snplist(efo="Gastric adenocarcinoma",population="East Asian")	
# 	snplist4<-make_snplist(efo="Gastric adenocarcinoma",population="East Asian")			
# 	snp_list<-ls()[grep("snplist[0-9]",ls())] 	
# 	snplist<-lapply(1:length(snp_list),FUN=function(x) eval(parse(text=snp_list[x])))
# 	return(snplist)
# }


make_snplist_eastasian_bedgrch37<-function(){
	# gunzip /projects/MRC-IEU/users/ph14916/fatty_acids_summary/data_maf0.01_rs.bim.gz
	ref<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",sep=" ",head=F,stringsAsFactors=F)

	ref<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/data_maf0.01_rs.bim",sep="\t",head=F,stringsAsFactors=F)
	ref<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_East_Asians_rsidsonly2.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/data_maf0.01_rs.bim",exact_match=TRUE,file_sep="\t",Head=FALSE)
	ref$bed<-paste(ref$V1,"_",ref$V4,sep="")
	# ref1<-convert_to_bed(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt")
	write.table(ref$bed,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_East_Asians_rsidsonly2_bed_grch37.txt",quote=F,row.names=F,col.names=F)	
}


make_snplist_bj<-function(){

	Disease<-c(
	"cervical carcinoma",
	"cholangiocarcinoma",#"Biliary tract cancer"
	"colorectal cancer",
	"endometrial carcinoma",
	"esophageal squamous cell carcinoma",#"Esophageal cancer"
	"gastric carcinoma",
	"hematological malignancy",#"Hematological malignancy"
	"hepatocellular carcinoma",
	"lung carcinoma",
	"ovarian carcinoma",
	"pancreatic carcinoma",
	"prostate carcinoma")

	snplist1<-make_snplist(efo=Disease[1],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist2<-make_snplist(efo=Disease[2],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist3<-make_snplist(efo=Disease[3],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist4<-make_snplist(efo=Disease[4],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist5<-make_snplist(efo=Disease[5],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist6<-make_snplist(efo=Disease[6],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist7_1<-make_snplist(efo="lymphoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")	
	snplist7_2<-make_snplist(efo="multiple myeloma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")	
	snplist7_3<-make_snplist(efo="lymphoid leukemia",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist7<-unique(c(snplist7_1,snplist7_2,snplist7_3))
	rm(list=c("snplist7_1","snplist7_2","snplist7_3"))
	snplist8<-make_snplist(efo=Disease[8],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist9<-make_snplist(efo=Disease[9],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist10<-make_snplist(efo=Disease[10],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist11<-make_snplist(efo=Disease[11],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
	snplist12<-make_snplist(efo=Disease[12],population="East Asian",Dir="~/fatty-acids/outcome_data/data/")

	snp_list<-ls()[grep("snplist[0-9]",ls())] 
	snplist<-lapply(1:length(snp_list),FUN=function(x) eval(parse(text=snp_list[x])))
}

# i<-4
# IDS[i]
# Disease[i]
# library(TwoSampleMR)
# cd /mnt/storage/private/mrcieu/research/scratch/IGD/data/public/bbj-a-113

extract_snps_bbj_cancer<-function(snplist=NULL){	
	dat1<-ieugwasr::associations(id= "bbj-a-98", variants=unlist(snplist[1])) 
	dat2<-ieugwasr::associations(id= "bbj-a-92", variants=unlist(snplist[2])) 
	dat3<-ieugwasr::associations(id= "bbj-a-107", variants=unlist(snplist[3])) 
	dat4<-ieugwasr::associations(id= "bbj-a-113", variants=unlist(snplist[4])) 
	dat5<-ieugwasr::associations(id= "bbj-a-117", variants=unlist(snplist[5])) 
	dat6<-ieugwasr::associations(id= "bbj-a-119", variants=unlist(snplist[6])) 
	dat7<-ieugwasr::associations(id= "bbj-a-125", variants=unlist(snplist[7])) 
	dat8<-ieugwasr::associations(id= "bbj-a-158", variants=unlist(snplist[8])) 
	dat9<-ieugwasr::associations(id= "bbj-a-133", variants=unlist(snplist[9])) 
	dat10<-ieugwasr::associations(id= "bbj-a-139", variants=unlist(snplist[10])) 
	dat11<-ieugwasr::associations(id= "bbj-a-140", variants=unlist(snplist[11])) 
	dat12<-ieugwasr::associations(id= "bbj-a-148", variants=unlist(snplist[12])) 
	dat_list<-ls()[grep("dat[0-9]",ls())] 
	dat<-lapply(1:length(dat_list),FUN=function(x) eval(parse(text=dat_list[x])))
	return(dat)
}


format_data_bbj<-function(bbj=NULL){
	IDS<-c(
	11,
	9,
	12,
	13,
	14,
	15,
	10,
	16,
	17,
	18,
	19,
	20)

	id_open_gwas<-unique(do.call(rbind,bbj)$id)			
	# ao<-TwoSampleMR::available_outcomes()
	# ao1<-ao[ao$id %in% id_open_gwas,c("trait","id","ncase","ncontrol")]
	# save(ao1,file="~/fatty-acids/outcome_data/data/ao1_bbj.Rdata")	
	load("~/fatty-acids/outcome_data/data/ao1_bbj.Rdata")
	dat_list<-NULL	
	for(i in 1:length(bbj)){
		dat<-data.frame(bbj[i],stringsAsFactors=F)
		dat<-merge(dat,ao1[,c("id","ncase","ncontrol")],by="id")
		ID<-IDS[i]
		outcome<-unique(dat$trait)
		outcome<-gsub("Esophageal cancer",	 "Esophageal squamous cell carcinoma",outcome)
		outcome<-gsub("Gastric cancer","Gastric adenocarcinoma",outcome)
		outcome<-gsub("Hematological malignancy","Blood cancer",outcome)
		outcome<-gsub("hepatocellular carcinoma","Liver cancer",outcome)
		dat<-format_data(Dat=dat,outcome=outcome,population="East Asian",pmid=32514122,study="BJ",ncase="ncase",ncontrol="ncontrol",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",ID=ID,all_summary_stats=TRUE)
		dat_list[[i]]<-dat
	}
	dat<-do.call(rbind,dat_list)	
	return(dat)
}


preformat_accc_3<-function(){
	Crc_accc1<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt",sep="\t",stringsAsFactors=F,head=T)	
	Crc_accc2<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_batch2_regions.txt",sep="\t",stringsAsFactors=F,head=T)	
	accc<-rbind(Crc_accc1,Crc_accc2)
	SNP<-gregexpr(":",accc$snp)
	Test<-unlist(lapply(1:length(SNP),FUN=function(x)
		length(unlist(SNP[x]))))
	Pos<-which(Test==3)
	accc2<-accc[Pos,]
	SNP<-unlist(strsplit(accc2$snp,split=":"))
	accc2$chr<-SNP[seq(1,length(SNP),by=4)]
	accc2$bp<-SNP[seq(2,length(SNP),by=4)]
	accc2$chr<-paste0("chr",accc2$chr)
	acc<-find_rsids(dat=accc2,ref_dat=TRUE)
}


preformat_ewi_27<-function(snplist=NULL){
	Ewi<-read.csv("/projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Postel-Vinay22327514/data/Ewing_association_Apr2017.csv",stringsAsFactors=F,head=F,skip=1)

	Ewi<-Ewi[,names(Ewi) != "V1"]
	names(Ewi)<-c("initial_order","CHR","SNP","BP","A1..Effect.allele.","TEST","NMISS","OR","ORs.allele.freqs.calculations","relative.error.OR","SE","L99","U99","STAT","P","A1.Ewing","A2.Ewing","MAF.Ewing","A1.controls","A2.controls","MAF.controls","maf.control.inferred.adj.with.ewing.allele","relative.error.maf.control","A1.HWE","O.HET.ALL","E.HET.ALL","HWE.Pval.ALL","O.HET.AFF","E.HET.AFF","HWE.Pval.AFF","O.HET.UNAFF","E.HET.UNAFF","HWE.Pval.UNAFF")
	Ewi$SNP<-trimws(Ewi$SNP)
	# write.table(Ewi,"/projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Postel-Vinay22327514/data/Ewing_association_Apr2017.tab",col.names=T,row.names=F,quote=F)
	# Dat<-lapply(names(Ewi), FUN=function(x) trimws(Ewi[,x]))
	# Dat2<-do.call(cbind,Dat)

	# snplist<-unique(readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"))
	Ewi2<-Ewi[Ewi$SNP %in% snplist,]
	Ewi2$A1.Ewing<-trimws(Ewi2$A1.Ewing)
	Ewi2$A1.controls<-trimws(Ewi2$A1.controls)
	Ewi2$A1..Effect.allele.<-trimws(Ewi2$A1..Effect.allele.)
	Ewi2<-Ewi2[Ewi2$A1..Effect.allele. == Ewi2$A1.controls,]
	return(Ewi2)
}

preformat_mes_98<-function(){
	Dat2<-read.table("~/MR_FattyAcids/data/summary_data/mesothelioma_23626673/Results_GWAS_Cugliari_SNP.csv",sep=";",head=T,stringsAsFactors=F)
	Dat2$beta<-as.numeric(gsub(",",".",Dat2$beta))
	Dat2$se<-as.numeric(gsub(",",".",Dat2$se))
	Dat2$p<-as.numeric(gsub(",",".",Dat2$p))
	Dat2$info<-as.numeric(gsub(",",".",Dat2$info))	
	Dat2<-Dat2[,!names(Dat2) %in% c("ProxySNP","trait","pval","pmid","type","alias_rs","author","id" )]
}

preformat_pan_121<-function(){
	pan<-read.table("~/MR_FattyAcids/data/summary_data/panc4/Summary Data - snplist.txt",sep="\t",head=T,stringsAsFactors=F)
	ref<-read.table("~/fatty-acids/colocalisation/data/UKBB_10K_bed_hg19.txt",sep=" ",head=F,stringsAsFactors=F)
	ref1<-format_ref(dat=ref)
	pan<-merge(pan,ref1,by.x=c("chromosome","position"),by.y=c("chr","V2"))
	return(pan)
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

preformat_gli_67<-function(dat=NULL){
	N<-unlist(strsplit(dat$Subjects,split="\\|"))
	dat$controls<-as.numeric(N[seq(1,length(N),by=2)])
	dat$cases<-as.numeric(N[seq(2,length(N),by=2)])
	N<-unlist(strsplit(dat$MAF,split="\\|"))
	dat$eaf.controls<-as.numeric(N[seq(1,length(N),by=2)])
	return(dat)
}


preformat_pro_128<-function(){
	prac<-read.table("/projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Practical_data_app_262/App262_MRbase/MetaOverall/meta_v3_onco_euro_overall_1_fixed.txt",sep="\t",head=T,stringsAsFactors=F)
	ref<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",sep=" ",head=F,stringsAsFactors=F)
	ref1<-format_ref(dat=ref)	
	dat<-merge(prac,ref1,by.x=c("Chr","bp"),by.y=c("chr","V2"))
	write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/practical/meta_v3_onco_euro_overall_1_fixed_plusrsids.txt",sep="\t",col.names=T,row.names=F,quote=F)	
}


find_cancerstudies_opengwas<-function(){	
	Cancers_selected_ukb<-c("Cancer code, self-reported: basal cell carcinoma","Diagnoses - main ICD10: C67 Malignant neoplasm of bladder","Malignant neoplasm of digestive organs","Type of cancer: ICD10: C64 Malignant neoplasm of kidney, except renal pelvis", "Lymphomas","Malignant neoplasm of skin","Malignant neoplasm of respiratory system and intrathoracic organs","Cancer code  self-reported: small intestine/small bowel cancer","Cancer code  self-reported: squamous cell carcinoma","Malignant neoplasm of urinary organs")
	 
	 Cancers_selected_fin<-c("Malignant neoplasm of bladder","Primary_lymphoid and hematopoietic malignant neoplasms","Malignant neoplasm of brain","Malignant neoplasm of breast","Any event in cancer register","Neoplasms","Malignant neoplasm of digestive organs","Malignant neoplasm of eye, brain and central nervous system","Colorectal cancer","Malignant neoplasm of endocrine gland","Malignant neoplasm of corpus uteri","malignant neoplasm of female genital organs","Follicular lymphoma","Malignant neoplasm of kidney, except renal pelvis","Lung cancer and mesothelioma","Malignant neoplasm of bronchus and lung","Lymphoid leukaemia","malignant neoplasm of male genital organs","Malignant neoplasm of skin","Multiple myeloma and malignant plasma cell neoplasms","Non-follicular lymphoma","Other and unspecified types of non-Hodgkin lymphoma","Malignant neoplasm of lip, oral cavity and pharynx","Malignant neoplasm of ovary","Malignant neoplasm of pancreas","Malignant neoplasm of prostate","Malignant neoplasm of respiratory system and intrathoracic organs","Malignant neoplasm of thyroid gland","Malignant neoplasm of urinary organs")
	
	Token<-ieugwasr::get_access_token()
	# Token<-ieugwasr::check_access_token()
	ao<-ieugwasr::gwasinfo(access_token = Token)	
	# Disease<-Disease
	# Disease<-c("ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)","ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)","Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)")
	Cancers<-c("Hematological malignancy", "neuroblastoma","glioma","cancer","carcinoma","adenocarcinoma")		
	Pos1<-unlist(lapply(Cancers,FUN=function(x) grep(x,ao$trait,ignore.case=T)))
	Pos2<-which(ao$trait %in% Cancers_selected_ukb)
	Pos3<-which(ao$trait %in% Cancers_selected_fin)
	Pos<-unique(c(Pos1,Pos2,Pos3))
	Dat<-ao[Pos,c("trait","consortium","id","ncase","ncontrol","pmid")]
	Dat<-Dat[order(Dat$ncase,decreasing=T),]
	cancerstudies_opengwas<-Dat
	save(cancerstudies_opengwas,file="~/fatty-acids/outcome_data/data/cancerstudies_opengwas.RData")
	return(Dat)
}


# https://risteys.finngen.fi/
make_snplist_malignant_skin_cancer<-function(){
	malignant_skin_cancers<-c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma")
	# "kaposi sarcoma"
	snplist<-NULL
	for(i in 1:length(malignant_skin_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_skin_cancers[i],population="European",Dir="~/fatty-acids/outcome_data/data/")
	}
	snplist<-unique(unlist(snplist))
	return(snplist)
}

# C3_NASAL_CAVITY_MIDDLE_EAR, C3_ACCESSORY_SINUS, C3_LARYNX, C3_TRACHEA, C3_BRONCHUS_LUNG, C3_THYMUS, C3_HEART_MEDIASTINUM_PLEURA, C3_RESPIRATORY_INTRATHORACIC3_NAS
make_snplist_respiratory<-function(){
	malignant_respiratory_cancers<-c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer")
	snplist<-NULL
	for(i in 1:length(malignant_respiratory_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_respiratory_cancers[i],population="European",Dir="~/fatty-acids/outcome_data/data/")
	}
	snplist<-unique(unlist(snplist))
	return(snplist)
}


# C3_KIDNEY_NOTRENALPELVIS, C3_RENAL_PELVIS, C3_URETER, C3_BLADDER, C3_URINARY_TRACT_NAS
make_snplist_urinary<-function(){
	malignant_urinary_cancers<-c("kidney cancer","nephroblastoma","renal cell carcinoma","bladder carcinoma")
	snplist<-NULL
	for(i in 1:length(malignant_urinary_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_urinary_cancers[i],population="European",Dir="~/fatty-acids/outcome_data/data/")
	}
	snplist<-unique(unlist(snplist))
	return(snplist)
}


make_snplist_blood<-function(population="European",Dir="~/fatty-acids/outcome_data/data/"){
	malignant_blood_cancers<-c("lymphoma","multiple myeloma","lymphoid leukemia")
	snplist<-NULL
	for(i in 1:length(malignant_blood_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_blood_cancers[i],population=population,Dir=Dir)
	}
	snplist<-unique(unlist(snplist))		
	return(snplist)
}


# endocrine neoplasm: A benign or malignant neoplasm arising from the epithelial cells of an endocrine organ. Representative examples include pituitary gland adenoma, pituitary gland carcinoma, thyroid gland carcinoma, carcinoid tumor, and neuroendocrine carcinoma.

# C3_THYROID_GLAND, C3_ADRENAL_GLAND, C3_ENDOCRINE_NAS

make_snplist_end<-function(){
	malignant_end_cancers<-c("pituitary gland adenoma","thyroid carcinoma","carcinoid tumor","neuroendocrine neoplasm")
	snplist<-NULL
	for(i in 1:length(malignant_end_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_end_cancers[i],population="European",Dir="~/fatty-acids/outcome_data/data/")
	}
	snplist<-unique(unlist(snplist))	
	return(snplist)
}



# malignant neoplasm of female genital organs
# Genital neoplasm, female: Tumor or cancer of the female reproductive tract (GENITALIA, FEMALE).
# C3_VULVA, C3_VAGINA, C3_CERVIX_UTERI, C3_CORPUS_UTERI, C3_UTERUS_NAS, C3_OVARY, C3_FEMALE_GENITAL_NAS, C3_PLACENTA
make_snplist_fgen<-function(){
	malignant_fgen_cancers<-c("cervical carcinoma","endometrial carcinoma","ovarian carcinoma")
	snplist<-NULL
	for(i in 1:length(malignant_fgen_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_fgen_cancers[i],population="European",Dir="~/fatty-acids/outcome_data/data/")
	}
	snplist<-unique(unlist(snplist))	
	return(snplist)
}



# malignant neoplasm of male genital organs
# male reproductive organ cancer: A reproductive organ cancer that is manifested in the male genital system. This includes organs such as the penis and scrotum.
# C3_PENIS, C3_PROSTATE, C3_TESTIS, C3_MALE_GENITAL_NAS
make_snplist_mgen<-function(){
	malignant_mgen_cancers<-c("prostate carcinoma","testicular carcinoma")
	snplist<-NULL
	for(i in 1:length(malignant_mgen_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_mgen_cancers[i],population="European",Dir="~/fatty-acids/outcome_data/data/")
	}
	snplist<-unique(unlist(snplist))	
	return(snplist)
}


preformat_hpc_23<-function(){
	hpc<-readLines("~/MR_FattyAcids/data/summary_data/Hepatocellular carcinoma_22807686/Reply/fatty_acid_instruments_SNPtable_20190920/fatty_acid_instruments_SNPtable_results.assoc.logistic")
	Line<-lapply(1:length(hpc), FUN=function(i) unlist(strsplit(hpc[i],split=" ")))
	Lines<-lapply(1:length(Line),FUN=function(i) unlist(Line[i])[unlist(Line[i])!=""])
	hpc<-data.frame(do.call(rbind,Lines),stringsAsFactors=F)
	names(hpc)<-hpc[1,]
	hpc<-hpc[2:nrow(hpc),]
	meta.dat<-readLines("~/MR_FattyAcids/data/summary_data/Hepatocellular carcinoma_22807686/Reply/fatty_acid_instruments_SNPtable_20190920/fatty_acid_instruments_SNPtable_freq-case-control.frq.cc")
	Line<-lapply(1:length(meta.dat), FUN=function(i) unlist(strsplit(meta.dat[i],split=" ")))
	Lines<-lapply(1:length(Line),FUN=function(i) unlist(Line[i])[unlist(Line[i])!=""])	
	meta.dat<-data.frame(do.call(rbind,Lines),stringsAsFactors=F)	
	names(meta.dat)<-meta.dat[1,]
	meta.dat<-meta.dat[2:nrow(meta.dat),]
	names(meta.dat)<-gsub("\t","",names(meta.dat))
	hpc2<-merge(hpc,meta.dat,by="SNP")
	hpc2$ncase<-as.numeric(hpc2$NCHROBS_A)/2
	hpc2$ncontrol<-as.numeric(hpc2$NCHROBS_U)/2
	Info<-read.table("~/MR_FattyAcids/data/summary_data/Hepatocellular carcinoma_22807686/Reply/175_snps_rsid_info.txt",sep="\t",stringsAsFactors=F,head=T)
	# hpc2$SNP[!hpc2$SNP %in% Info$rsid]
	
	hpc3<-merge(hpc2,Info,by.x="SNP",by.y="rsid",all.x=F)
	return(hpc3)
}

# C3_LIP_ORAL_PHARYNX
# Malignant neoplasm of lip, oral cavity and pharynx
# C3_LIP, C3_BASEOFTONGUE, C3_TONGUENAS, C3_GUM, C3_FLOOROFMOUTH, C3_PALATE, C3_MOUTHNAS, C3_PAROTIDGLAND, C3_SALIVARYNAS, C3_TONSIL, C3_OROPHARYNG, C3_NASOPHARYNX, C3_PIRIFORMSINUS, C3_HYPOPHARYNX, C3_LIP_ORAL_PHARYNX_NAS
make_snplist_hnc<-function(population="European",Dir="~/fatty-acids/outcome_data/data/"){
	malignant_hnc_cancers<-c("head and neck squamous cell carcinoma","oropharynx cancer","nasopharyngeal neoplasm","hypopharynx cancer","oral cavity cancer","mouth neoplasm","pharynx cancer")
	snplist<-NULL
	for(i in 1:length(malignant_hnc_cancers)){
		snplist[[i]]<-make_snplist(efo=malignant_hnc_cancers[i],population=population,Dir=Dir)
	}
	snplist<-unique(unlist(snplist))
	return(snplist)
}


transform_betas<-function(dat=NULL,effect="lnor",effect.se="se"){
	# formula: log OR = beta / (u(1-u)); where u=ncases/(ncases + ncontrol) REPEAT with SE 	
	beta<-dat[,effect]
	se<-dat[,effect.se]
	u<-dat$ncase/(dat$ncase+dat$ncontrol)
	dat[,effect] <- beta / (u * (1 - u))
	dat[,effect.se]<-se / (u * (1 - u)) 	
	return(dat)
}




# transform_betas<-function(dat=NULL){
# 	# formula: log OR = beta / (u(1-u)); where u=ncases/(ncases + ncontrol) REPEAT with SE 	
# 	beta<-dat$lnor
# 	se<-dat$se
# 	u<-dat$ncase/(dat$ncase+dat$ncontrol)
# 	dat$lnor <- beta / (u * (1 - u))
# 	dat$se<-se / (u * (1 - u)) 	
# 	return(dat)
# }

harmonise_effect_allele2<-function(Dat=NULL){
	alleles<-paste0(unique(c(Dat1$Effect.Allele,Dat1$Other.Allele)),collapse="")
	
	if(any(alleles %in% c("GC","CG","TA","AT"))){
		warning("palindromic SNP detected")
	}
	EA<-Dat$Effect.Allele[1]
	Pos<-Dat$Effect.Allele!=EA
	lnor<-Dat$lnor[Pos]*-1
	Dat$lnor[Pos]<-lnor
	oa<-Dat$Effect.Allele[Pos]
	ea<-Dat$Other.Allele[Pos]
	Dat$Effect.Allele[Pos]<-ea
	Dat$Other.Allele[Pos]<-oa
	eaf<-1-Dat$eaf[Pos]
	Dat$eaf[Pos]<-eaf	
	return(Dat)
}
