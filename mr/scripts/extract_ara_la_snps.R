library(CheckSumStats)

snplist <- readLines("~/fatty-acids/mr/data/ara_la_snplist_independent.txt")


out_dir<-"~/fatty-acids/outcome_data/ara_la/"

##########################
# colorectal cancer GECCO#
###########################
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions_v2.R")
# source("~/fatty-acids/colocalisation/scripts/extract_SNPs_functions.R")

out_file<-paste0(out_dir,"crc60","_ara_la.txt")	
Crc_all<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer/crc_all.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=Crc_all,outcome="Colorectal cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=58131,ncontrol=67347,UKbiobank=TRUE,rsid="rs",effect_allele="Allele1",other_allele="Allele2",lnor="Effect",lnor_se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=60,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",z_score="LRT.Z")

dat<-dat[dat$rsid %in% snplist,]

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############################
# colorectal cancer in ACCC###
#################################
out_file<-paste0(out_dir,"crc3","_ara_la.txt")	
crc<-preformat_accc_3()
dat<-format_data(dat=crc,outcome="Colorectal cancer",population="East Asian",pmid=31826910,study="ACCC",ncase=23572,ncontrol=48700,UKbiobank=FALSE,rsid="rsid",effect_allele="Allele1",other_allele="Allele2",eaf="Freq1",lnor="Effect",lnor_se="StdErr",p="P.value",effect_allele_confirmed=TRUE,ID=3,all_summary_stats=FALSE,Direction="Direction",phet="HetPVal",I2="HetISq",Q="HetChiSq")
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
# lung cancer TRICL#
#####################

# out_file<-paste0(out_dir,"luc73","_ara_la.txt")
# luc<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/TRICL_LungCancer/TRICL_Meta_032114_Overall.csv",exact_match=TRUE,path_to_target_file_sep=",")
# dat<-format_data(dat=luc,outcome="Lung cancer",population="European",pmid=27488534,ncase="N_Cases",ncontrol="N_Controls",study="ILCCO",UKbiobank=FALSE,rsid="rs_number",effect_allele="effect_allele",other_allele="reference_allele",or="OR_fixed",lnor_se="StdError_fixed",eaf="EAF",p="Pvalue_fixed",I2="I2",phet="phete.Q._Pvalue",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=73,open_gwas=FALSE,efo = "lung carcinoma")

#########################
# Lung cancer UK Biobank#
#########################
snplist <- readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/ara_la_snplist_independent.txt")
out_dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/ara_la/"
out_file<-paste0(out_dir,"luc149","_ara_la.txt")	
ukb<-CheckSumStats::extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lung_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-CheckSumStats::format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Lung cancer",efo="lung carcinoma",ID=149)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################
#Lung cancer TRICL#
###################
ieugwasr::get_access_token()
out_file<-paste0(out_dir,"luc75","_ara_la.txt")
luc <- ieugwasr::associations(id="ieu-a-987", variants=snplist)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=28604730,ncase=29266,ncontrol=56450,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=75,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

############################
# Lung cancer biobank japan#
############################

out_file<-paste0(out_dir,"luc17","_ara_la.txt")
bbj<-ieugwasr::associations(id= "bbj-a-133", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Lung cancer",population="East Asian",pmid=32514122,study="BJ",ncase=4050,ncontrol=208403,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=17,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#####################
#Lung cancer FinnGen#
#####################

out_file<-paste0(out_dir,"luc42","_ara_la.txt")
fin <- ieugwasr::associations(id="finn-a-LUNG_CANCER_MESOT", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid="finn-a-LUNG_CANCER_MESOT",ncase=673,ncontrol=95826,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=42,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##################################
#Malignant skin cancer UK Biobank#
##################################

# Representative examples of malignant skin neoplasms include basal cell carcinoma, squamous cell carcinoma, melanoma, and Kaposi sarcoma.

out_file<-paste0(out_dir,"msc152_corrected","_ara_la.txt")
ukb <- ieugwasr::associations(id="ukb-d-C3_SKIN", variants=snplist)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="ukb-d-C3_SKIN",ncase=16531,ncontrol=344663,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=152,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###############################
#Malignant skin cancer FinnGen#
###############################

out_file<-paste0(out_dir,"msc46","_ara_la.txt")
fin <- ieugwasr::associations(id="finn-a-C3_SKIN", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="finn-a-C3_SKIN",ncase=895,ncontrol=95604,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=46,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################################################
#Respiratory and intrathoracic cancer UK Biobank#
#################################################

# C3_NASAL_CAVITY_MIDDLE_EAR, C3_ACCESSORY_SINUS, C3_LARYNX, C3_TRACHEA, C3_BRONCHUS_LUNG, C3_THYMUS, C3_HEART_MEDIASTINUM_PLEURA, C3_RESPIRATORY_INTRATHORACIC3_NAS

out_file<-paste0(out_dir,"ric160_corrected","_ara_la.txt")
ukb <- ieugwasr::associations(id="ukb-d-C3_RESPIRATORY_INTRATHORACIC", variants=snplist)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="ukb-d-C3_RESPIRATORY_INTRATHORACIC",ncase=1944,ncontrol=359250,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=160,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)
##############################################
#Respiratory and intrathoracic cancer FinnGen#
#############################################

out_file<-paste0(out_dir,"ric54","_ara_la.txt")
fin <- ieugwasr::associations(id="finn-a-C3_RESPIRATORY_INTRATHORACIC", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="finn-a-C3_RESPIRATORY_INTRATHORACIC",ncase=615,ncontrol=95884,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=54,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########################
#Overall cancer FinnGen#
########################

out_file<-paste0(out_dir,"can32","_ara_la.txt")
fin <- ieugwasr::associations(id="finn-a-ANY_CANC", variants=snplist,proxies=0)  
dat<-CheckSumStats::format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-ANY_CANC",ncase=9792,ncontrol=86707,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=32,open_gwas=TRUE,efo= "cancer")

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########################
#Overall cancer FinnGen#
########################

out_file<-paste0(out_dir,"can33","_ara_la.txt")

fin <- ieugwasr::associations(id="finn-a-II_NEOPLASM", variants=snplist,proxies=0)  
dat<-CheckSumStats::format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-II_NEOPLASM",ncase=31217,ncontrol=65282,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=33,open_gwas=TRUE,efo= "cancer")

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###############################################################
# overall cancer UK Biobank excluding non-melanoma skin cancer#
###############################################################

out_file<-paste0(out_dir,"oac141","_ara_la.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Overall cancer (excluding non-melanoma skin cancer)",efo="cancer",ID=141)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############################
# overall cancer UK Biobank#
############################

out_file<-paste0(out_dir,"oac140","_ara_la.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_inclc44_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Overall cancer",efo="cancer",ID=140)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#######################
# esophageal scc N-UGC#
########################

out_file<-paste0(out_dir,"esc99","_ara_la.txt")	
ugi<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/summary_chr_all.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
dat<-format_data(dat=ugi,outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=25129146,ncase=2013,ncontrol=2701,study="N-UGC",UKbiobank=FALSE,rsid="rs",effect_allele="risk_allele",other_allele="reference_allele",lnor="beta",lnor_se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info1="info",ID=99,all_summary_stats=TRUE,efo="esophageal squamous cell carcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################################################
# Esophageal squamous cell carcinoma biobank japan#
###################################################

out_file<-paste0(out_dir,"esc14","_ara_la.txt")
bbj<-ieugwasr::associations(id= "bbj-a-117", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=1300,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=14,open_gwas=TRUE,efo="esophageal squamous cell carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


######################################
# non melanoma skin cancer UK biobank#
######################################

out_file<-paste0(out_dir,"nmc156","_ara_la.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_nm_skin_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Non-melanoma skin cancer",efo="non-melanoma skin carcinoma",ID=156)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###############################
# basal cell carcinoma 23andMe#
###############################

out_file<-paste0(out_dir,"bcc1","_ara_la.txt")	
dat<-extract_snps_and_format_bcc_23andMe(snplist=snplist)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################################
#Basal cell carcinoma UK biobank#
##################################

out_file<-paste0(out_dir,"bcc135","_ara_la.txt")
ukb <- ieugwasr::associations(id="ukb-b-8837", variants=snplist)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Basal cell carcinoma",population="European",pmid="ukb-b-8837",ncase=4290,ncontrol=458643,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=135,open_gwas=TRUE,efo = "basal cell carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


# preformat_accc_3<-function(){
# 	Crc_accc1<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt",sep="\t",stringsAsFactors=F,head=T)	
# 	Crc_accc2<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_batch2_regions.txt",sep="\t",stringsAsFactors=F,head=T)	
# 	accc<-rbind(Crc_accc1,Crc_accc2)
# 	SNP<-gregexpr(":",accc$snp)
# 	Test<-unlist(lapply(1:length(SNP),FUN=function(x)
# 		length(unlist(SNP[x]))))
# 	Pos<-which(Test==3)
# 	accc2<-accc[Pos,]
# 	SNP<-unlist(strsplit(accc2$snp,split=":"))
# 	accc2$chr<-SNP[seq(1,length(SNP),by=4)]
# 	accc2$bp<-SNP[seq(2,length(SNP),by=4)]
# 	accc2$chr<-paste0("chr",accc2$chr)
# 	acc<-find_rsids(dat=accc2,ref_dat=TRUE)
# }


# find_rsids<-function(dat=NULL,ref_dat=FALSE){	
# 	if(ref_dat){
# 		ref<-read.table("~/fatty-acids/colocalisation/data/UKBB_10K_bed_hg19.txt",head=F,stringsAsFactors=F,sep=" ")
# 	}
# 	dat2<-merge(dat,ref,by.x=c("chr","bp"),by.y=c("V1","V2"))
# 	names(dat2)[names(dat2) == "V4"]<-"rsid"
# 	return(dat2)
# }


# transform_betas<-function(dat=NULL,effect="lnor",effect.se="lnor_se"){
# 	# formula: log OR = beta / (u(1-u)); where u=ncases/(ncases + ncontrol) REPEAT with SE 	
# 	beta<-dat[,effect]
# 	se<-dat[,effect.se]
# 	u<-dat$ncase/(dat$ncase+dat$ncontrol)
# 	dat[,effect] <- beta / (u * (1 - u))
# 	dat[,effect.se]<-se / (u * (1 - u)) 	
# 	return(dat)
# }


