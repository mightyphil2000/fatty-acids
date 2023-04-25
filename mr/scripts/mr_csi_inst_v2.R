library(TwoSampleMR)
source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/outcome_data/scripts/functions_combine_and_format_outcomes.R")
# Exposure
Csi<-read.table("~/fatty-acids/mr/data/Csi_sig_clumped.txt",sep=" ",head=TRUE,stringsAsFactors=FALSE)
Csi<-Csi[Csi$clump_r2==0.001,]
dim(Csi)
exposure_dat<-format_exposure2(dat=Csi,standardise_beta=FALSE,beta="BETA",se="SE",pval="pval",effect_allele="EFFECT_ALLELE",other_allele="OTHER_ALLELE",eaf="EAF",rsid="SNP",ID="id",snps=NULL,exposure="csi")
exposure_dat$id.exposure<-exposure_dat$id
exposure_dat$exposure<-"csi"
exposure_dat$beta.exposure<-exposure_dat$beta.exposure/0.6940093
exposure_dat$se.exposure<-exposure_dat$se.exposure/0.6940093

# Outcome
# Dat2
load("~/fatty-acids/mr/data/cancer_csi_inst.Rdata")
Temp<-Dat2[Dat2$ID == 91,]
unique(Temp$file.outcome)
outcome_dat<-format_outcomes3(dat=Dat2)

# harmonise data
dat <- harmonise_data(exposure_dat = exposure_dat,outcome_dat =outcome_dat)
# dat$outcome<-"csi"
head(dat)
make_instab_csi()

res<-mr(dat,method_list="mr_ivw")
res2<-format_results(res=res)

load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
# Pos<-grep("East Asian",disc.tab9$population,invert=TRUE)
# ID<-disc.tab9$ID[Pos]
dat.meta<-meta_analysis2(dat=res2,IDS=disc.tab9$ID)

res_csi<-format_results2()
# res_csi$b<-as.numeric(res_csi$b)*0.6940093 #transform to SD scale
# res_csi$se<-as.numeric(res_csi$se)*0.6940093 #transform to SD scale

# exp(as.numeric(res_csi$b[res_csi$outcome2 == "Colorectal cancer"]))

save(res_csi,file="~/fatty-acids/mr/results/res_csi_v2.Rdata") #transformed to SD scale


# format_results2<-function(){
# 	ids<-res2$id.outcome	
# 	ID<-dat.meta$ID
# 	ID<-trimws(unlist(strsplit(IDS,split=";")))
# 	ids_keep<-ids[!ids %in% ID]
# 	res3<-res2[res2$id.outcome %in% ids_keep,]
# 	res3$nstudies <- 1
# 	res3$Q.p<-NA 
# 	# names(res3)[!names(res3) %in% names(dat.meta)]
# 	res4<-plyr::rbind.fill(res3,dat.meta)
# 	res4<-res4[order(as.numeric(res4$ncase.outcome),decreasing=TRUE),]
# 	res4<-res4[!duplicated(res4$outcome2),]
# 	return(res4)
# }


# format_results<-function()
# {
		
# 		out<-unique(outcome_dat[,c("id.outcome","population","ncase.outcome","ncontrol.outcome","outcome2")])		
# 		out$ncase.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
# 			median(out$ncase.outcome[out$id.outcome==x])))
# 		out$ncontrol.outcome<-unlist(lapply(out$id.outcome,FUN=function(x)
# 			median(out$ncontrol.outcome[out$id.outcome==x])))
# 		out<-unique(out)

# 		res.m<-merge(res,out,by="id.outcome")
# 		res.m<-merge(res.m,meta.tab9[,c("ID","Cancer.Group","system","study","study.abbreviation")],by.x="id.outcome",by.y="ID")

# 		# if(exclude_east_asians)
# 		# {
# 		# 	res.m<-res.m[res.m$population!="East Asian",]
# 		# }
# 		return(res.m)
# } 

# format_outcomes3<-function(dat=NULL,all_cols.keep=FALSE){

# 	names(dat)[names(dat) =="lnor"]<-"beta.outcome"
# 	names(dat)[names(dat) =="lnor_se"]<-"se.outcome"
# 	names(dat)[names(dat) =="p"]<-"pval.outcome"
# 	names(dat)[names(dat) %in% c("cases","ncase")]<-"ncase.outcome"
# 	names(dat)[names(dat) %in% c("controls","ncontrol")]<-"ncontrol.outcome"
# 	names(dat)[names(dat) =="effect_allele"]<-"effect_allele.outcome"
# 	names(dat)[names(dat) =="other_allele"]<-"other_allele.outcome"
# 	names(dat)[names(dat) =="eaf"]<-"eaf.outcome"
# 	names(dat)[names(dat) =="rsid"]<-"SNP"
# 	names(dat)[names(dat) =="ID"]<-"id.outcome"

# 	# dat[dat$population == "European; East Asian",c("id.outcome")]
# 	# exclude meta-analysed outcome data, to avoid MR of mixed European and East Asian populations. Need to MR East Asian and European studies separately prior to meta analysis
# 	dat<-dat[grep(";",dat$id.outcome,invert=T),]	
# 	dat$population[dat$id.outcome %in% 993:999]<-"European"
# 	# dat$study.abbreviation[is.na(dat$study.abbreviation)]<-"survival"
# 	dat$outcome2<-paste(dat$outcome,dat$id.outcome)
# 	if(any(names(dat) == "cancer")){
# 		dat$outcome2<-paste(dat$cancer,dat$id.outcome)
# 	}

# 	# Dups<-unique(dat$outcome2[duplicated(dat$outcome2)])
# 	# Pos<-dat$outcome2 %in% Dups
# 	# if(any(Pos)) stop("duplicate outcomes when none expected")
# 	# dat$outcome2[Pos]<-paste(dat$outcome2[Pos],dat$id.outcome[Pos])
# 	# outcome<-dat$outcome2
# 	# dat$outcome2<-dat$cancer
# 	# dat$outcome<-outcome

# 	# if(all(c("cancer","outcome") %in% names(dat))) warning("two column names for cancer outcome when only one expected")
# 	# names(dat)[names(dat) %in% c("cancer","outcome")]<-"outcome"
# 	dat$se.outcome<-as.numeric(dat$se.outcome)
# 	dat$beta.outcome<-as.numeric(dat$beta.outcome)
# 	dat$eaf.outcome<-as.numeric(dat$eaf.outcome)
# 	dat$pval.outcome<-as.numeric(dat$pval.outcome)
# 	# names(dat)[names(dat) =="cancer"]<-"outcome"
# 	Cols.keep<-c("SNP","outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","ncase.outcome","ncontrol.outcome","effect_allele.outcome","other_allele.outcome","id.outcome","population")

# 	dat<-dat[, Cols.keep]
# 	dat$outcome2<-dat$outcome
# 	dat$outcome<-paste0(dat$outcome," | ",dat$id.outcome)	
	
# 	dat$index<-paste0(dat$outcome,dat$SNP)
# 	Dups<-unique(dat$id.outcome[which(duplicated(dat$index))])
# 	if(any(duplicated(dat$index))) warning(paste0("duplicate SNPs present, IDs=: ",paste(Dups,collapse=" | ")))

# 	return(dat)
# }


# meta_analysis2<-function(dat=NULL,IDS=IDS){
# 	IDS<-IDS[grep(";",IDS)]
# 	Res_list<-NULL	
# 	for(i in 1:length(IDS)){
# 		print(IDS[i])
# 		ID<-trimws(unlist(strsplit(IDS[i],split=";")))
# 		if(length(ID)==1) stop("ID=1")
# 		dat1<-dat[dat$id.outcome %in% ID,]
# 		dat1<-dat1[dat1$population=="European",]
# 		if(nrow(dat1)>0)
# 		{

			
# 			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
# 			# if(sum(grep("FinnGen",dat1$study.abbreviation))!=0) stop("fingen")
# 			b<-dat1$b
# 			se<-dat1$se
# 			# p<-temp$p
# 			w<-1/se^2
# 			b.fixed<-sum(b*w)/(sum(w))
# 			se.fixed<-sqrt(sum(w)^-1)
# 			z<-abs(b.fixed/se.fixed)
# 			p.fixed<-pnorm(z,lower.tail=F)*2
# 			nstudies.fixed<-length(b)
# 			cancer<-unique(dat1$outcome2)
# 			if(length(cancer)!=1) cancer<-unique(dat1$Cancer.Group)
# 			if(length(cancer)!=1) 	cancer<-unique(paste(dat1$system,"system cancers"))
# 			if(length(cancer)!=1) stop("length of cancer not 1")
# 			if(length(unique(dat1$population))!=1) stop("population not same across studies")
# 			ids.fixed<-paste(dat1$id.outcome,collapse="; ")
# 			cases<-sum(dat1$ncase.outcome)
# 			controls<-sum(dat1$ncontrol.outcome)
# 			study<-"Overall fixed effect"
# 			Q<-sum((b.fixed-b)^2*w)
# 			df.Q<-length(b)-1		
# 			Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
# 			# names(Meta)
# 			# Meta$Q
# 			# Meta$pval.Q
# 			# Meta<-metagen(TE=b,seTE=se,comb.fixed=T,sm="MD")
# 			# # Meta$TE.fixed
# 			# Meta$seTE.fixed
# 			# Meta$pval.fixed
# 			# Q.p
# 			# dat1[,c("cancer","Effect.Allele","Other.Allele","eaf","lnor","se","cases","study.abbreviation")]
# 			# EA<-unique(dat1$Effect.Allele)
# 			# OA<-unique(dat1$Other.Allele)
# 			# EAF<-round(sum((dat1$eaf*w))/sum(w),3)
# 			# if(length(EA)>1) stop("effect allele not consistent across studies")

# 			# studies<-paste(dat1$study.abbreviation,collapse="; ")
# 			# Cols<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
# 			# for(j in 1:length(Cols)){
# 			# 	Cols[j]
# 			# 	Cols[j]<-paste(unique(dat1[,Cols[j]]),collapse="; ")
# 			# }
# 			dat.matrix<-c(cancer,b.fixed,se.fixed,p.fixed,nstudies.fixed,ids.fixed,cases,controls,study,Q.p)
# 			Res<-data.frame(matrix(dat.matrix,nrow=length(cancer),ncol=length(dat.matrix)),stringsAsFactors=F)
# 			names(Res)<-c("outcome2","b","se","pval","nstudies","id.outcome","ncase.outcome","ncontrol.outcome","study","Q.p")
# 			# Col.dat<-data.frame(matrix(Cols,ncol=length(Cols),nrow=1),stringsAsFactors=F)
# 			# names(Col.dat)<-c("study.abbreviation","pmid","outcome","UKbiobank","rsid","effect_allele_confirmed","site","include","system","cell","Cancer.Group","original_name","population","subpopulation","total","Note","overlap","MAC_rs174546", "MAC100rs174546","set","id.mrbase","proxy","proxy_snp","proxy_snp.outcome","ProxySNP","proxy_outcome_snp")
# 			# Res.c<-cbind(Res,Col.dat)		
# 			round(as.numeric(Res$lnor),3)
# 			round(as.numeric(Res$se),3)
# 			# round(as.numeric(Res$p),5)
# 			# options("scipen"=2, "digits"=3)
# 			# format(as.numeric(Res$p), digits=3)		
# 			Res_list[[i]]<-Res
# 		}
# 	}
# 	Res1<-do.call(rbind,Res_list)
# 	Res1$pval<-as.numeric(Res1$pval)
# 	Res1$population<-"European"
# 	return(Res1)
# 	# i4<-which(as.numeric(Res1$Q.p)<0.10 & as.numeric(Res1$Q.p)>0.05)
# }
