library(CheckSumStats)

# secondary_pufa_instruments/instruments.Rdata created using clump_secondary_pufas.R

# system("scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_instruments/instruments.Rdata ~/fatty-acids/mr/data") 

# epi franklin server
load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_instruments/instruments.Rdata")
out_dir_server<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_outcomes/"

# local machine
load("~/fatty-acids/mr/data/instruments.Rdata")
out_dir<-"~/fatty-acids/outcome_data/secondary_pufas/"


snplist_eur<-unique(eur1$snp) #eur1 is eur restricted to single largest study for each trait
snplist_eas<-unique(eas$snp) #all traits in east asians are from the SCHS
suffix<-"_secondary_pufas.txt"

##########################
# colorectal cancer GECCO#
##########################
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions_v2.R")
# source("~/fatty-acids/colocalisation/scripts/extract_SNPs_functions.R")
out_dir_server<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/secondary_pufas/"
out_file<-paste0(out_dir,"crc60",suffix)	
# Crc_all<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer/crc_all.txt",sep="\t",head=T,stringsAsFactors=F)

Crc_all<-read.csv("~/MR_FattyAcids/data/summary_data/colorectal_cancer/18112021/1255_MarginalResults_HRC125K_20211118.csv",head=TRUE,stringsAsFactors=FALSE)

dat<-format_data(dat=Crc_all,outcome="Colorectal cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=58131,ncontrol=67347,UKbiobank=TRUE,rsid="rs",effect_allele="Allele1",other_allele="Allele2",lnor="Effect",lnor_se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=60,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",z_score="LRT.Z")

dat<-dat[dat$rsid %in% snplist_eur,]
dim(dat)

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############################
# colorectal cancer in ACCC###
#################################
out_file<-paste0(out_dir,"crc3",suffix)	
crc<-preformat_accc_3()
dat<-format_data(dat=crc,outcome="Colorectal cancer",population="East Asian",pmid=31826910,study="ACCC",ncase=23572,ncontrol=48700,UKbiobank=FALSE,rsid="rsid",effect_allele="Allele1",other_allele="Allele2",eaf="Freq1",lnor="Effect",lnor_se="StdErr",p="P.value",effect_allele_confirmed=TRUE,ID=3,all_summary_stats=FALSE,Direction="Direction",phet="HetPVal",I2="HetISq",Q="HetChiSq")
dat1<-dat[dat$rsid %in% snplist_eur,]
dat2<-dat[dat$rsid %in% snplist_eas,]
dat1$note<-"snp-trait sig in Europeans"
dat2$note<-"snp-trait sig in East Asians"
dat<-rbind(dat1,dat2)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################
#Lung cancer TRICL#
###################
# ieugwasr::get_access_token()
# out_file<-paste0(out_dir,"luc75",suffix)
# luc <- ieugwasr::associations(id="ieu-a-987", variants=snplist_eur)  
# ieugwasr::associations(id="ieu-a-987", variants="rs1741")  
# luc[luc$rsid == "rs1741",]
# exp_dat[exp_dat$SNP == "rs1741",]
# head(exp_dat)
# exp_dat[exp_dat$exposure == "Arachidonic acid (20:4n6)",]
# dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=28604730,ncase=29266,ncontrol=56450,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=75,open_gwas=TRUE,efo = "lung carcinoma")
# write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

############################
# Lung cancer biobank japan#
############################

out_file<-paste0(out_dir,"luc17",suffix)
bbj1<-ieugwasr::associations(id= "bbj-a-133", variants=snplist_eur) 
bbj2<-ieugwasr::associations(id= "bbj-a-133", variants=snplist_eas) 
bbj1$note<-"snp-trait sig in Europeans"
bbj2$note<-"snp-trait sig in East Asians"
bbj<-rbind(bbj1,bbj2)
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Lung cancer",population="East Asian",pmid=32514122,study="BJ",ncase=4050,ncontrol=208403,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=17,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#####################
#Lung cancer FinnGen#
#####################

out_file<-paste0(out_dir,"luc42",suffix)
fin <- ieugwasr::associations(id="finn-a-LUNG_CANCER_MESOT", variants=snplist_eur,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid="finn-a-LUNG_CANCER_MESOT",ncase=673,ncontrol=95826,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=42,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##################################
#Malignant skin cancer UK Biobank#
##################################

# Representative examples of malignant skin neoplasms include basal cell carcinoma, squamous cell carcinoma, melanoma, and Kaposi sarcoma.

out_file<-paste0(out_dir,"msc152_corrected",suffix)
ukb <- ieugwasr::associations(id="ukb-d-C3_SKIN", variants=snplist_eur)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="ukb-d-C3_SKIN",ncase=16531,ncontrol=344663,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=152,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###############################
#Malignant skin cancer FinnGen#
###############################

out_file<-paste0(out_dir,"msc46",suffix)
fin <- ieugwasr::associations(id="finn-a-C3_SKIN", variants=snplist_eur,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="finn-a-C3_SKIN",ncase=895,ncontrol=95604,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=46,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################################################
#Respiratory and intrathoracic cancer UK Biobank#
#################################################

# C3_NASAL_CAVITY_MIDDLE_EAR, C3_ACCESSORY_SINUS, C3_LARYNX, C3_TRACHEA, C3_BRONCHUS_LUNG, C3_THYMUS, C3_HEART_MEDIASTINUM_PLEURA, C3_RESPIRATORY_INTRATHORACIC3_NAS

out_file<-paste0(out_dir,"ric160_corrected",suffix)
ukb <- ieugwasr::associations(id="ukb-d-C3_RESPIRATORY_INTRATHORACIC", variants=snplist_eur)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="ukb-d-C3_RESPIRATORY_INTRATHORACIC",ncase=1944,ncontrol=359250,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=160,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)
##############################################
#Respiratory and intrathoracic cancer FinnGen#
#############################################

out_file<-paste0(out_dir,"ric54",suffix)
fin <- ieugwasr::associations(id="finn-a-C3_RESPIRATORY_INTRATHORACIC", variants=snplist_eur,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="finn-a-C3_RESPIRATORY_INTRATHORACIC",ncase=615,ncontrol=95884,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=54,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########################
#Overall cancer FinnGen#
########################

out_file<-paste0(out_dir,"can32",suffix)
fin <- ieugwasr::associations(id="finn-a-ANY_CANC", variants=snplist_eur,proxies=0)  
dat<-CheckSumStats::format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-ANY_CANC",ncase=9792,ncontrol=86707,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=32,open_gwas=TRUE,efo= "cancer")

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########################
#Overall cancer FinnGen#
########################

out_file<-paste0(out_dir,"can33",suffix)
fin <- ieugwasr::associations(id="finn-a-II_NEOPLASM", variants=snplist_eur,proxies=0)  
dat<-CheckSumStats::format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-II_NEOPLASM",ncase=31217,ncontrol=65282,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=33,open_gwas=TRUE,efo= "cancer")

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################################################
# Esophageal squamous cell carcinoma biobank japan#
###################################################

out_file<-paste0(out_dir,"esc14",suffix)
bbj1<-ieugwasr::associations(id= "bbj-a-117", variants=snplist_eur) 
bbj2<-ieugwasr::associations(id= "bbj-a-117", variants=snplist_eas) 
bbj1$note<-"snp-trait sig in Europeans"
bbj2$note<-"snp-trait sig in East Asians"
bbj<-rbind(bbj1,bbj2)
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=1300,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=14,open_gwas=TRUE,efo="esophageal squamous cell carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################################################
# colorectal cancer biobank japan#
###################################################

out_file<-paste0(out_dir,"crc12",suffix)
bbj1<-ieugwasr::associations(id= "bbj-a-107", variants=snplist_eur)
bbj2<-ieugwasr::associations(id= "bbj-a-107", variants=snplist_eas)
bbj1$note<-"snp-trait sig in Europeans"
bbj2$note<-"snp-trait sig in East Asians"
bbj<-rbind(bbj1,bbj2)

dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Colorectal cancer",population="East Asian",pmid=32514122,study="BJ",ncase=7062,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=12,open_gwas=TRUE,efo="colorectal cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################################
#Basal cell carcinoma UK biobank#
##################################
# "rs10215255" missing SNP for ALA from UKB
# snplist_eur =="rs10215255"
out_file<-paste0(out_dir,"bcc135",suffix)
ukb <- ieugwasr::associations(id="ukb-b-8837", variants=snplist_eur)  
# which(ukb$rsid =="rs10215255")
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Basal cell carcinoma",population="European",pmid="ukb-b-8837",ncase=4290,ncontrol=458643,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=135,open_gwas=TRUE,efo = "basal cell carcinoma")
dat<-transform_betas(dat=dat)
dim(dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###############################
# basal cell carcinoma HNMSC#
###############################
# out_file<-paste0(out_dir,"bcc70",suffix)
# bcc<-read.table("~/MR_FattyAcids/data/summary_data/BCC/fattyacidSNP_BCC_HarvardGWAS_4242cases_Meta_Final_Results.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)
# dat<-format_data(dat=bcc,outcome="Basal cell carcinoma",population="European",pmid=23548203,study="HNMSC",ncase=4242,ncontrol=12802,UKbiobank=FALSE,rsid="SNP",effect_allele="Allele1",other_allele="Allele2",lnor="Effect",lnor_se="StdErr",eaf="Freq1",p="Pvalue",effect_allele_confirmed=TRUE,ID=70,phet="HetPVal",I2="HetISq",Direction="Direction")

# dat<-dat[dat$rsid %in% snplist_eur,]
# dim(dat)
# too many missing SNPs. wrote to Jiali for additional SNPs 28/1/2022
# write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#########################
# associations on server#
#########################


####################
# lung cancer TRICL#
#####################

out_file<-paste0(out_dir_server,"luc75",suffix)
luc<-extract_snps(snplist=snplist_eur,path_to_target_file="/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc/Onco_TRICL_032116_Overall.csv.tab",exact_match=TRUE,path_to_target_file_sep="\t")
dat<-format_data(dat=luc,outcome="Lung cancer",population="European",pmid=28604730,ncase=29266,ncontrol=56450,study="ILCCO",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",lnor="beta",lnor_se="se",eaf="effect_allele_freq",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=75,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###############################################################
# overall cancer UK Biobank excluding non-melanoma skin cancer#
###############################################################

out_file<-paste0(out_dir_server,"oac141",suffix)	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_cancer_imputed.txt.gz",snplist=snplist_eur,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Overall cancer (excluding non-melanoma skin cancer)",efo="cancer",ID=141)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############################
# overall cancer UK Biobank#
############################

out_file<-paste0(out_dir_server,"oac140",suffix)	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_inclc44_cancer_imputed.txt.gz",snplist=snplist_eur,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Overall cancer",efo="cancer",ID=140)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#########################
# Lung cancer UK Biobank#
#########################
out_file<-paste0(out_dir_server,"luc149",suffix)	
ukb<-CheckSumStats::extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lung_cancer_imputed.txt.gz",snplist=snplist_eur,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-CheckSumStats::format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Lung cancer",efo="lung carcinoma",ID=149)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###############################
# Colorectal cancer UK Biobank#
###############################

out_file<-paste0(out_dir_server,"crc143",suffix)	
ukb<-CheckSumStats::extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_colorectal_cancer_imputed.txt.gz",snplist=snplist_eur,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)

dat<-CheckSumStats::format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Colorectal cancer",efo="colorectal cancer",ID=143)

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)



#######################
# esophageal scc N-UGC#
########################

out_file<-paste0(out_dir_server,"esc99",suffix)	
ugi1<-extract_snps(snplist=snplist_eur,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/summary_chr_all.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
ugi2<-extract_snps(snplist=snplist_eas,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/summary_chr_all.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
ugi1$note<-"snp-trait sig in Europeans"
ugi2$note<-"snp-trait sig in East Asians"
ugi<-rbind(ugi1,ugi2)
dat<-format_data(dat=ugi,outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=25129146,ncase=2013,ncontrol=2701,study="N-UGC",UKbiobank=FALSE,rsid="rs",effect_allele="risk_allele",other_allele="reference_allele",lnor="beta",lnor_se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info1="info",ID=99,all_summary_stats=TRUE,efo="esophageal squamous cell carcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


######################################
# non melanoma skin cancer UK biobank#
######################################

out_file<-paste0(out_dir_server,"nmc156",suffix)	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_nm_skin_cancer_imputed.txt.gz",snplist=snplist_eur,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Non-melanoma skin cancer",efo="non-melanoma skin carcinoma",ID=156)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###############################
# basal cell carcinoma 23andMe#
###############################

out_file<-paste0(out_dir_server,"bcc1",suffix)	
dat<-extract_snps_and_format_bcc_23andMe(snplist=snplist_eur)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

system("scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_outcomes/*.txt ~/fatty-acids/outcome_data/secondary_pufas/") 

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


find_rsids<-function(dat=NULL,ref_dat=FALSE){	
	if(ref_dat){
		ref<-read.table("~/fatty-acids/colocalisation/data/UKBB_10K_bed_hg19.txt",head=F,stringsAsFactors=F,sep=" ")
	}
	dat2<-merge(dat,ref,by.x=c("chr","bp"),by.y=c("V1","V2"))
	names(dat2)[names(dat2) == "V4"]<-"rsid"
	return(dat2)
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


