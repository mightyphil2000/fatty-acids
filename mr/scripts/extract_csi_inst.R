devtools::install_github("MRCIEU/CheckSumStats")

# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/

# head /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_unique_rsidonly_clumped.txt
setwd("~/mrQC")
library(CheckSumStats)
library(devtools)
load_all()
document()
# check()




out_dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/csi_inst/"
snplist<-readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_unique_rsidonly_clumped.txt")

# ####################
# Glioma 28346443#####
######################
out_file<-paste0(out_dir,"gli66","_csi_inst.txt")
Gli<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/research/projects/icep1/bioinformatics/013/working/data/new_data/meta_gicc_mda_allglioma.TBL",exact_match=TRUE)
dat<-format_data(dat=Gli,outcome="Glioma",population="European",pmid=28346443,study="GICC/MDA",ncase=5747,ncontrol=5522,UKbiobank=FALSE,rsid="MarkerName",effect_allele="Allele1",other_allele="Allele2",lnor="Effect",lnor_se="StdErr",eaf="Freq1",p="P.value",info=NA,info1=NA,info2=NA,info3=NA,HWEp=NA,phet=NA,I2=NA,Q=NA,Direction=NA,effect_allele_confirmed=TRUE,or=NA,or_lci=NA,or_uci=NA,chr="CHR",pos="POS",z_score=NA,ID=66,all_summary_stats=TRUE,efo="glioma")

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


########################
# Endometrial cancer####
########################
out_file<-paste0(out_dir,"enc25","_csi_inst.txt")
End<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/endometrial_cancer/ecac_gwas_data.txt",exact_match=TRUE,path_to_target_file_sep=" ")
dat<-format_data(dat=End,outcome="Endometrial cancer",population="European",pmid=30093612,study="ECAC",ncase=12906,ncontrol=108979,UKbiobank=TRUE,rsid="SNPID",effect_allele="EA",other_allele="OA",lnor="ALL_BETA",lnor_se="ALL_SE",eaf="MEAN_EAF",p="ALL_PVALUE",info1="OA_INFO",effect_allele_confirmed=TRUE,ID=25,all_summary_stats=TRUE,efo="endometrial carcinoma")
any(duplicated(dat$rsid))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


# ############################################
# acute lymphoblastic leukemia 22076464/####
############################################
out_file<-paste0(out_dir,"all21","_csi_inst.txt")
All<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/acute_lymphoblastic_leukemia_22076464/Assoc_ALL_Affy5.assoc",exact_match=TRUE)
dat<-format_data(dat=All,outcome="Acute lymphoblastic leukaemia",population="European",pmid=22076464,study="C-ALL",ncase=419,ncontrol=474,UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele="A2",or="OR",lnor_se="SE",or_lci="L95",or_uci="U95",eaf="F_U",p="P",effect_allele_confirmed=TRUE,ID=21,all_summary_stats=TRUE,efo="acute lymphoblastic leukemia")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


################################
# esophageal adenocarcinoma####
################################
# no eaf or maf reported
out_file<-paste0(out_dir,"esa24","_csi_inst.txt")
Eso<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/esophageal_adenocarcinoma/Oesophageal_adenocarcinoma_Lancet_Oncol2016_Gharahkhani_et_al.txt",exact_match=TRUE,path_to_target_file_sep=" ")

dat<-format_data(dat=Eso,outcome="Esophageal adenocarcinoma",population="European",pmid=27527254,study="EAS",ncase="N_cases",ncontrol="N_controls",UKbiobank=FALSE,rsid="MarkerName",effect_allele="effect_allele",other_allele="non.effect_allele",lnor="Effect",lnor_se="StdErr",p="P.value",effect_allele_confirmed=TRUE,ID=24,all_summary_stats=TRUE,I2="HetISq",Q="HetChiSq",phet="HetPVal",Direction="Direction",efo="esophageal adenocarcinoma")

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


##################
# Bladder cancer###
#####################
# chrALL* created using bash script in pre_harmonise_outcomes.sh
out_file<-paste0(out_dir,"blc105","_csi_inst.txt")
bla<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/bladder_cancer/Haycock/chrALL_GWAS_NBCS_risico_jan2017_methodscore_4mds_info.out",exact_match=TRUE,path_to_target_file_sep=" ")
bla<-preformat_bla_id105()
dat<-format_data(dat=bla,outcome="Bladder cancer",population="European",pmid=21750109,study="NBCS",ncase=1799,ncontrol=4745,UKbiobank=FALSE,rsid="rsid",effect_allele="allele_B",other_allele="allele_A",eaf="eaf",lnor="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_beta_1",lnor_se="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_se_1",p="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_pvalue",effect_allele_confirmed=TRUE,info1="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_info",ID=105,all_summary_stats=TRUE,efo="bladder carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


##################
# uveal melanoma###
#####################
out_file<-paste0(out_dir,"uvm165","_csi_inst.txt")
uve<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/uveal_melanoma/GWAS_results_17122018.txt",exact_match=TRUE,path_to_target_file_sep="\t")
dat<-format_data(dat=uve,outcome="Uveal melanoma",population="European",pmid=28781888,study="UMS",ncase=259,ncontrol=401,UKbiobank=FALSE,rsid="SNP",effect_allele="MINOR",eaf="MAF",or="OR",lnor_se="SE",or_lci="L95",or_uci="U95",p="P",effect_allele_confirmed=TRUE,HWEp="P_HWE",z_score="STAT",ID=165,all_summary_stats=TRUE,efo="uveal melanoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

################
# Melanoma  ####
################
# Melanoma_metaanalysis_chrALL* creared using bash script in pre_harmonise_outcomes.sh
out_file<-paste0(out_dir,"mel95","_csi_inst.txt")
mel<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/melanoma/Melanoma_meta_single_files/Melanoma_metaanalysis_chrALL_15052019_EAF_RSQ.txt",exact_match=TRUE,path_to_target_file_sep=" ")
dat<-format_data(dat=mel,outcome="Melanoma",population="European",pmid=26237428,study="MMAC",ncase="CASE",ncontrol="CONTROL",UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele="A2",or="OR",lnor_se="SE_fixed_qnorm",p="P",effect_allele_confirmed=TRUE,Q="Q",z_score="Z_fixed_sign",eaf="eaf",I2="I",info1="RSQ_median",ID=95,all_summary_stats=TRUE,efo="melanoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###############
# UK Biobank###
###############

#####
#CNS#
#####
out_file<-paste0(out_dir,"brc138","_csi_inst.txt")

ukb<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_brain_cancer_imputed.txt.gz",exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,outcome="Brain cancer",population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",ID=138,all_summary_stats=TRUE,efo="central nervous system cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################
# breast cancer#
################
out_file<-paste0(out_dir,"brc139","_csi_inst.txt")
ukb<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_breast_cancer_imputed.txt.gz",exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,outcome="Breast cancer",population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",ID=139,all_summary_stats=TRUE,efo="breast carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#####################
# colorectal cancer#
####################
out_file<-paste0(out_dir,"crc143","_csi_inst.txt")
ukb<-extract_snps(snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE,path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_colorectal_cancer_imputed.txt.gz")
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Colorectal cancer",efo="colorectal cancer",ID=143)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###############
# blood cancer#
###############
out_file<-paste0(out_dir,"blc137","_csi_inst.txt")
ukb<-extract_snps(snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE,path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_haem_cancer_imputed.txt.gz")
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Blood cancer",efo=c("lymphoma","multiple myeloma","lymphoid leukemia"),ID=137)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

############
#h&n cancer#
############
out_file<-paste0(out_dir,"opc157","_csi_inst.txt")
ukb<-extract_snps(snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE,path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_headneck_cancer_imputed.txt.gz")
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Oral cavity and pharyngeal cancer",efo=c("head and neck squamous cell carcinoma","oropharynx cancer","nasopharyngeal neoplasm","hypopharynx cancer","oral cavity cancer","mouth neoplasm","pharynx cancer"),ID=157)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######
# ALL#
#######
out_file<-paste0(out_dir,"leu146","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_leuk_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Leukaemia",efo="acute lymphoblastic leukemia",ID=146)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##########################
#Liver & bile duct cancer#
##########################
out_file<-paste0(out_dir,"lbc147","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_liver_bile_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Liver & bile duct cancer",efo="hepatocellular carcinoma",ID=147)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###############
# Liver cancer#
################
out_file<-paste0(out_dir,"lic148","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_liver_cell_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Liver cancer",efo="hepatocellular carcinoma",ID=148)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############
# Lung cancer#
##############
out_file<-paste0(out_dir,"luc149","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lung_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Lung cancer",efo="lung carcinoma",ID=149)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################
# lung carcinoma 1499#
######################

out_file<-paste0(out_dir,"luc1499","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lung_cancer_unadj_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Lung cancer unadjusted for chip",efo="lung carcinoma",ID=1499
	)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#####
# ALL#
####

out_file<-paste0(out_dir,"lle150","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lymph_leuk_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Lymphoid leukaemia",efo="acute lymphoblastic leukemia",ID=150)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################
# multiple myeloma#
###################
out_file<-paste0(out_dir,"mum154","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_mult_myel_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Multiple myeloma",efo="multiple myeloma",ID=154)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################
# Myeloid leukaemia#
###############
out_file<-paste0(out_dir,"myl155","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_myel_leuk_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Myeloid Leukaemia",efo="acute myeloid leukemia",ID=155)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################
# non melanoma skin cancer
#########
out_file<-paste0(out_dir,"nmc156","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_nm_skin_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Non-melanoma skin cancer",efo="non-melanoma skin carcinoma",ID=156)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#####################
# esophageal cancer#
####################
out_file<-paste0(out_dir,"esa144","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_oesoph_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Esophageal adenocarcinoma",efo="esophageal adenocarcinoma",ID=144)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

################
# ovarian cancer#
#################

out_file<-paste0(out_dir,"ovc158","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_ovarian_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Ovarian cancer",efo="ovarian carcinoma",ID=158)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##################
# overall cancer
################
out_file<-paste0(out_dir,"oac141","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,
	outcome="Overall cancer (excluding non-melanoma skin cancer)",efo="cancer",ID=141)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

################
# overall cancer2
#####

out_file<-paste0(out_dir,"oac140","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_inclc44_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Overall cancer",efo="cancer",ID=140)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##################
# prostate cancer#
##################

out_file<-paste0(out_dir,"pro159","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_prostate_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,outcome="Prostate cancer",efo="prostate carcinoma",ID=159)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###########
# melanoma
###########

out_file<-paste0(out_dir,"mel153","_csi_inst.txt")	
ukb<-extract_snps(path_to_target_file="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_skin_cancer_imputed.txt.gz"
	,snplist=snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=TRUE)
dat<-format_data(dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",effect_allele="ALLELE1",other_allele="ALLELE0",lnor="beta",lnor_se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info1="INFO",all_summary_stats=TRUE,	outcome="Melanoma",	efo="melanoma",ID=153)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)




###################################
# Upper gastrointestinal cancers#####
########################################
# created files using pre_harmonise_outcomes.sh
out_file<-paste0(out_dir,"gca102","_csi_inst.txt")	
ugi<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_CC/summary_chr_all.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
dat<-format_data(dat=ugi,outcome="Gastric cardia adenocarcinoma",population="East Asian",pmid=26129866,ncase=1189,ncontrol=2708,study="N-UGC",UKbiobank=FALSE,rsid="rs",effect_allele="risk_allele",other_allele="reference_allele",lnor="beta",lnor_se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info1="info",ID=102,all_summary_stats=TRUE,efo="gastric adenocarcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########
# esophageal scc
####

out_file<-paste0(out_dir,"esc99","_csi_inst.txt")	
ugi<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/summary_chr_all.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
dat<-format_data(dat=ugi,outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=25129146,ncase=2013,ncontrol=2701,study="N-UGC",UKbiobank=FALSE,rsid="rs",effect_allele="risk_allele",other_allele="reference_allele",lnor="beta",lnor_se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info1="info",ID=99,all_summary_stats=TRUE,efo="esophageal squamous cell carcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############
# gastric adenocarcinoma
##############
out_file<-paste0(out_dir,"gac101","_csi_inst.txt")	
ugi<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_gastric/summary_chr_all.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
dat<-format_data(dat=ugi,outcome="Gastric adenocarcinoma",population="East Asian",pmid=26129866,ncase=2350,ncontrol=2708,study="N-UGC",UKbiobank=FALSE,rsid="rs",effect_allele="risk_allele",other_allele="reference_allele",lnor="beta",lnor_se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info1="info",ID=101,all_summary_stats=TRUE,efo="gastric adenocarcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##########
# gastric adenocarcinoma
#####
out_file<-paste0(out_dir,"nga103","_csi_inst.txt")	
ugi<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_NC/summary_chr_all.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
dat<-format_data(dat=ugi,outcome="Noncardia gastric adenocarcinoma",population="East Asian",pmid=26129866,ncase=1027,ncontrol=2708,study="N-UGC",UKbiobank=FALSE,rsid="rs",effect_allele="risk_allele",other_allele="reference_allele",lnor="beta",lnor_se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info1="info",ID=103,all_summary_stats=TRUE,efo="gastric adenocarcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###################################################
# Additional datasets identified in GWAS catalog###
######################################################

####################
# cervical cancer####
####################
# eaf/ maf not reported
out_file<-paste0(out_dir,"cec92","_csi_inst.txt")	
cec<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Cervical_cancer_28806749/28806749-GCST004833-EFO_0001061.h.tsv",exact_match=TRUE,path_to_target_file_sep="\t")
dat<-format_data(dat=cec,outcome="Cervical cancer",population="European",pmid=28806749,study="MCCS",ncase=2866,ncontrol=6481,UKbiobank=FALSE,rsid="variant_id",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_frequency",or="odds_ratio",lnor_se="standard_error",p="p_value",effect_allele_confirmed=TRUE,z_score="stat",ID=92,all_summary_stats=TRUE,efo="cervical carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


##########################################
# BRCA negative high risk breast cancer 30323354#
##########################################
# eaf/ maf not reported
out_file<-paste0(out_dir,"hrb88","_csi_inst.txt")	
# snplist<-make_snplist(efo="breast carcinoma")
brc<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/BRCA1_2_negative_high_risk_breast_cancer_30323354/30323354-GCST006719-EFO_0009443-build37.f.tsv",exact_match=TRUE,path_to_target_file_sep="\t")
dat<-format_data(dat=brc,outcome="BRCA negative high risk breast cancer",population="East Asian",pmid=30323354,study="KHBC",ncase=1469,ncontrol=5979,UKbiobank=FALSE,rsid="variant_id",effect_allele="effect_allele",other_allele=NA,or="odds_ratio",eaf=NA,p="gc",effect_allele_confirmed=TRUE,z_score="stat",ID=88,all_summary_stats=TRUE,efo="breast carcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Dat$Z<-qnorm(Dat$gc/2,lower.tail=F)#gc corrected p value. I looked up the top SNP from the paper and the gc value corresponds exactly to the reported P value for this SNP


##########################
# thyroid cancer #
##########################

out_file<-paste0(out_dir,"thc163","_csi_inst.txt")	
thy<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Thyroid_cancer_30104761/PheCode_193_SAIGE_MACge20.txt.vcf",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
dat<-format_data(dat=thy,outcome="Thyroid cancer",population="European",pmid=30104761,study="UKB",ncase="num_cases",ncontrol="num_controls",UKbiobank=TRUE,rsid="ID",effect_allele="ALT",other_allele="REF",lnor="beta",lnor_se="sebeta",eaf="af",p="pval",effect_allele_confirmed=FALSE,all_summary_stats=TRUE,ID=163,z_score="Tstat",efo="thyroid carcinoma")  # GAME-ON had same names for allele columns.  
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


####################################################
#Kidney cancer in females
####################################################
out_file<-paste0(out_dir,"kif90","_csi_inst.txt")	
kid<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Kidney_cancer_31231134/Laskar_31231134_Females",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
kid<-preformat_kid_90(dat=kid)
dat<-format_data(dat=kid,outcome="Kidney cancer in females",population="European",pmid=31231134,study="KidRISK",ncase="Number_cases",ncontrol="Number_controls",UKbiobank=FALSE,rsid="Variant_ID",effect_allele="Effect_allele",other_allele="Other_allele",or="Odds_ratio",lnor_se="Standard_error",eaf="effect_allele_frequency_controls",p="P_Value",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=90,efo="kidney cancer")  
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


####################################################
#Kidney cancer in males
####################################################
out_file<-paste0(out_dir,"kim91","_csi_inst.txt")	
kid<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Kidney_cancer_31231134/Laskar_31231134_Males",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
kid<-preformat_kid_90(dat=kid)
dat<-format_data(dat=kid,outcome="Kidney cancer in males",population="European",pmid=31231134,study="KidRISK",ncase="Number_cases",ncontrol="Number_controls",UKbiobank=FALSE,rsid="Variant_ID",effect_allele="Effect_allele",other_allele="Other_allele",or="Odds_ratio",lnor_se="Standard_error",eaf="effect_allele_frequency_controls",p="P_Value",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=91,efo="kidney cancer") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################################################### 
# Ovarian_cancer_EastAsians_30898391 OCAC
####################################################
# input file created using pre_harmonise_outcomes.sh
out_file<-paste0(out_dir,"ovc120","_csi_inst.txt")	
ova<-extract_snps(snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Ovarian_cancer_EastAsians_30898391/SummaryResults_Asian_chr_all.txt",exact_match=TRUE,path_to_target_file_sep=",",Comment="")
ova2<-preformat_ova_120(dat=ova)

# summary(ova2$R2_oncoarray[ova2$R2_oncoarray!=-99])
dat<-format_data(dat=ova2,outcome="Ovarian cancer",population="East Asian",pmid=30898391,study="OCAC (EAS)",ncase=3238,	ncontrol=4083,UKbiobank=FALSE,rsid="V4",effect_allele="Effect",other_allele="Baseline",lnor="overall_OR",lnor_se="overall_SE",eaf="EAF",p="overall_pvalue",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=120,efo="ovarian carcinoma")  
# dat[dat$lnor < -2,c("lnor","eaf")]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###########################################
# Bcell_nonHodgkinlymphoma_Chinese_23749188#
###########################################

# ambiguous effect allele column. drop
# out_file<-paste0(out_dir,"bnh5","_csi_inst.txt")	
# nhl<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Bcell_nonHodgkinlymphoma_Chinese_23749188/NHL_253cases_1438controls_536555snps_summary_stats.txt",exact_match=TRUE,path_to_target_file_sep="\t")
# dat<-format_data(dat=nhl,outcome="B cell non-Hodgkin lymphoma",population="East Asian",pmid=23749188,study="BC-NHL",ncase=253,ncontrol=1438,UKbiobank=FALSE,rsid="SNP",effect_allele="TestAllele",other_allele="MajorAlelle",eaf="TAF_Controls",or="OR",lnor_se="SE",p="P",effect_allele_confirmed=TRUE,ID=5,all_summary_stats=TRUE,efo="neoplasm of mature b-cells")
# write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###########################
# Non melanoma skin cancer#
###########################

# source("~/fatty-acids-mr/instruments/Extract_SNPs_function.R")
# bcc<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt",path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_bcc/Chahal_2016_basal_cell_carcinoma-4.1/basal_cell_carcinoma-4.1.dat",exact_match=TRUE,path_to_target_file_sep="\t")

out_file<-paste0(out_dir,"bcc1","_csi_inst.txt")	
dat<-extract_snps_and_format_bcc_23andMe(snplist=snplist)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

out_file<-paste0(out_dir,"scc2","_csi_inst.txt")	
dat<-extract_snps_and_format_scc_23andMe(snplist=snplist)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


# ############
#Ewing sarcoma#
##############
# ewi<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/ewing_sarcoma/Postel-Vinay22327514/data/Ewing_association_Apr2017.csv",exact_match=TRUE,path_to_target_file_sep=",")

# MAF reported but unknown whether corresponds to EAF 
out_file<-paste0(out_dir,"ews27","_csi_inst.txt")	
ewi<-preformat_ewi_27(snplist=snplist)
dat<-format_data(dat=ewi,outcome="Ewing's sarcoma",population="European",pmid=22327514,study="ESS",ncase=427,ncontrol=684,UKbiobank=FALSE,rsid="SNP",effect_allele="A1.controls",other_allele="A2.controls",or="OR",lnor_se="SE",p="P",effect_allele_confirmed=TRUE,ID=27,all_summary_stats=TRUE,z_score="STAT",HWEp="HWE.Pval.UNAFF",efo="Ewing sarcoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


############################################
# Chronic lymphocytic leukemia Interlymph###
############################################
# input file prepared using pre_harmonise_outcomes.sh

out_file<-paste0(out_dir,"cll83","_csi_inst.txt")	
cll<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/CLL/cll.out",exact_match=TRUE,path_to_target_file_sep="\t",fill=TRUE)
dat<-preformat_snps_interlymph(dat=cll)
dat<-format_data(dat=dat,outcome="Chronic lymphocytic leukaemia",population="European",pmid=26956414,study="InterLymph",ncase="NUM_CASE",ncontrol="NUM_CONTROL",eaf="eaf",UKbiobank=FALSE,rsid="SNP",effect_allele="effect_allele",other_allele="other_allele",lnor="lnor",lnor_se="se",p="P",effect_allele_confirmed=TRUE,ID=83,all_summary_stats=TRUE,info1="info",Direction="Direction",phet="Phet",I2="I2",efo="chronic lymphocytic leukemia")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)



#####################################
# Follicular lymphoma Interlymph ####
#####################################
# input file prepared using pre_harmonise_outcomes.sh
out_file<-paste0(out_dir,"fll85","_csi_inst.txt")	
fl<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/FL/fl.out",exact_match=TRUE,path_to_target_file_sep="\t",fill=TRUE)
dat<-preformat_snps_interlymph(dat=fl)
dat<-format_data(dat=dat,outcome="Follicular lymphoma",population="European",pmid=25279986,study="InterLymph",ncase="NUM_CASE",ncontrol="NUM_CONTROL",eaf="eaf",UKbiobank=FALSE,rsid="SNP",effect_allele="effect_allele",other_allele="other_allele",lnor="lnor",lnor_se="se",p="P",effect_allele_confirmed=TRUE,ID=85,all_summary_stats=TRUE,info1="info",Direction="Direction",phet="Phet",I2="I2",efo="follicular lymphoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###########################################
# Diffuse large B cell lymphoma interlymph####
###########################################
# input file prepared using pre_harmonise_outcomes.sh
out_file<-paste0(out_dir,"dlb84","_csi_inst.txt")	
dlbcl<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/DLBCL/dlbcl.out",exact_match=TRUE,path_to_target_file_sep="\t",fill=TRUE)
dat<-preformat_snps_interlymph(dat=dlbcl)
dat<-format_data(dat=dat,outcome="Diffuse large B cell lymphoma",population="European",pmid=25261932,study="InterLymph",ncase="NUM_CASE",ncontrol="NUM_CONTROL",eaf="eaf",UKbiobank=FALSE,rsid="SNP",effect_allele="effect_allele",other_allele="other_allele",lnor="lnor",lnor_se="se",p="P",effect_allele_confirmed=TRUE,ID=84,all_summary_stats=TRUE,info1="info",Direction="Direction",phet="Phet",I2="I2",efo="diffuse large b-cell lymphoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#####################################
# Marginal zone lymphoma        ####
#####################################
out_file<-paste0(out_dir,"mzl86","_csi_inst.txt")	
mzl<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/MZL/mzl.out",exact_match=TRUE,path_to_target_file_sep="\t",fill=TRUE)
dat<-format_data(dat=mzl,outcome="Marginal zone lymphoma",population="European",pmid=25569183,study="InterLymph",ncase="Num_Case",ncontrol="Num_Control",eaf="Effect_Allele_Freq_Control",UKbiobank=FALSE,rsid="Locus",effect_allele="Effect_Allele",other_allele="Reference_Allele",lnor="Beta",lnor_se="Standard_error_of_beta",p="P_value",effect_allele_confirmed=TRUE,ID=86,all_summary_stats=TRUE,info1="Info",efo="marginal zone B-cell lymphoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# cp /projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Rajaraman22886559/for_Philip_data_delivery.txt .
##########################
# GliomaScan #
##########################
out_file<-paste0(out_dir,"gli67","_csi_inst.txt")
gli<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gliomascan/for_Philip_data_delivery.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
gli<-preformat_gli_67(dat=gli)
dat<-format_data(dat=gli,outcome="Glioma",population="European",pmid=22886559,study="GliomaScan",ncase="cases",ncontrol="controls",UKbiobank=FALSE,rsid="Locus",effect_allele="Allele2",other_allele="Allele1",or="OR",or_lci="OR_95._CI_l",or_uci="OR_95._CI_u",eaf=NA,p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=67,efo="glioma")   #confirmed by correspondence 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################################
# Prostate cancer PRACTICAL Build 37#
######################################
out_file<-paste0(out_dir,"pro128","_csi_inst.txt")
#input file created using preformat_pro_128()
pro<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/practical/meta_v3_onco_euro_overall_1_fixed_plusrsids.txt",exact_match=TRUE,path_to_target_file_sep="\t",Comment="")
dat<-format_data(dat=pro,outcome="Prostate cancer",population="European",pmid=29892016,ncase=79148,ncontrol=61106,study="PRACTICAL",UKbiobank=FALSE,rsid="V4",effect_allele="Allele1",other_allele="Allele2",lnor="Effect",lnor_se="StdErr",eaf="Freq1",p="P.value",Direction = "Direction",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=128,efo="prostate carcinoma") #effect allele = allele 1 as confirmed in documentation ~/OpenGWAS/PRACTICAL/documentations
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


####################################
# hepatocellular_carcinoma_22174901##
######################################

out_file<-paste0(out_dir,"hpc68","_csi_inst.txt")
hpc<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/livercancer/HKU_data_summary.txt",exact_match=TRUE)
dat<-format_data(dat=hpc,outcome="Liver cancer",population="East Asian",pmid=22174901,study="HKHC",ncase=95,ncontrol=97,UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele="A2",lnor="ln.OR.",lnor_se="SE",eaf=NA,effect_allele_confirmed=TRUE,ID=68,all_summary_stats=TRUE,efo="hepatocellular carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


############
# GAMEONE####
############
# no MAF/EAF reported
out_file<-paste0(out_dir,"can57","_csi_inst.txt")
gam<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/GAMEON/CrossCancer_MetaAnalysis_Results_2016_Release_20190222.txt",exact_match=TRUE)
dat<-format_data(dat=gam,outcome="Cancer (5 sites)",population="European",pmid=27197191,study="GAME-ON",ncase=61851,ncontrol=61820,UKbiobank=FALSE,rsid="ID",effect_allele="ALT",other_allele="REF",or="meta.OR",p="meta.Pvalue",effect_allele_confirmed=TRUE,ID=57,all_summary_stats=TRUE,efo=c("breast carcinoma","lung carcinoma","colorectal cancer","ovarian carcinoma","prostate carcinoma"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############
# Lung cancer#
##############
out_file<-paste0(out_dir,"luc73","_csi_inst.txt")
luc<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/TRICL_LungCancer/TRICL_Meta_032114_Overall.csv",exact_match=TRUE,path_to_target_file_sep=",")
dat<-format_data(dat=luc,outcome="Lung cancer",population="European",pmid=27488534,ncase="N_Cases",ncontrol="N_Controls",study="ILCCO",UKbiobank=FALSE,rsid="rs_number",effect_allele="effect_allele",other_allele="reference_allele",or="OR_fixed",lnor_se="StdError_fixed",eaf="EAF",p="Pvalue_fixed",I2="I2",phet="phete.Q._Pvalue",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=73,open_gwas=FALSE,efo = "lung carcinoma") 
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#############
# OpenGWAS###
############



################
# Breast cancer#
################
out_file<-paste0(out_dir,"brc6","_csi_inst.txt")
brc <- ieugwasr::associations(id="ieu-a-1126", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(brc,stringsAsFactors=F),outcome="Breast cancer",population="European",pmid=29059683,ncase=122977,ncontrol=105974,study="BCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=6,open_gwas=TRUE,efo="breast carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################
#ER- Breast cancer#
###################
out_file<-paste0(out_dir,"erneg7","_csi_inst.txt")
brc <- ieugwasr::associations(id="ieu-a-1128", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(brc,stringsAsFactors=F),outcome="ER- breast cancer",population="European",pmid=29059683,ncase=21468,ncontrol=105974,study="BCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=7,open_gwas=TRUE,efo="estrogen-receptor negative breast cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###################
#ER+ Breast cancer#
###################
out_file<-paste0(out_dir,"erpos8","_csi_inst.txt")
brc$trait <- ieugwasr::associations(id="ieu-a-1127", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(brc,stringsAsFactors=F),outcome="ER+ breast cancer",population="European",pmid=29059683,ncase=69501,ncontrol=105974,study="BCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=8,open_gwas=TRUE,efo="estrogen-receptor positive breast cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################
#Ovarian cancer#
###################

out_file<-paste0(out_dir,"ovc118","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1120", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Ovarian cancer",population="European",pmid=28346442,ncase=25509,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=118,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###################
#Serous ovarian cancer#
###################
out_file<-paste0(out_dir,"soc119","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1228", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Serous ovarian cancer",population="European",pmid=28346442,ncase=14049,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=119,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##################################
#High grade serous ovarian cancer#
##################################

out_file<-paste0(out_dir,"hso110","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1121", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="High grade serous ovarian cancer",population="European",pmid=28346442,ncase=13037,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=110,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##################################
#Low malignant potential ovarian cancer#
##################################

out_file<-paste0(out_dir,"lmp115","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1233", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Low malignant potential ovarian cancer",population="European",pmid=28346442,ncase=3103,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=115,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################################################################
#Low grade & low malignant potential serous ovarian cancer#
####################################################################
out_file<-paste0(out_dir,"llo112","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1229", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Low grade & low malignant potential serous ovarian cancer",population="European",pmid=28346442,ncase=2966,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=112,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############################
# Endometrioid ovarian cancer#
##############################

out_file<-paste0(out_dir,"eoc109","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1125", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Endometrioid ovarian cancer",population="European",pmid=28346442,ncase=2810,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=109,open_gwas=TRUE,efo="ovarian endometrioid carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#########################
#Mucinous ovarian cancer#
#########################
out_file<-paste0(out_dir,"moc117","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1231", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Mucinous ovarian cancer",population="European",pmid=28346442,ncase=2566,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=117,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###############################################
#Low malignant potential serous ovarian cancer#
###############################################

out_file<-paste0(out_dir,"lso116","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1230", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Low malignant potential serous ovarian cancer",population="European",pmid=28346442,ncase=1954,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=116,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##################################
#Invasive mucinous ovarian cancer#
##################################

out_file<-paste0(out_dir,"imo111","_csi_inst.txt")

ovc <- ieugwasr::associations(id="ieu-a-1123", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Invasive mucinous ovarian cancer",population="European",pmid=28346442,ncase=1417,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=111,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


###########################
#Clear cell ovarian cancer#
###########################

out_file<-paste0(out_dir,"cco108","_csi_inst.txt")

ovc <- ieugwasr::associations(id="ieu-a-1124", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Clear cell ovarian cancer",population="European",pmid=28346442,ncase=1366,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=108,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################################################
#Low malignant potential mucinous ovarian cancer#
################################################

out_file<-paste0(out_dir,"lmm114","_csi_inst.txt")

ovc <- ieugwasr::associations(id="ieu-a-1232", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Low malignant potential mucinous ovarian cancer",population="European",pmid=28346442,ncase=1149,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=114,open_gwas=TRUE,efo="ovarian carcinoma")
# rs4927355
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################################
#Low grade serous ovarian cancer#
#################################

out_file<-paste0(out_dir,"lso113","_csi_inst.txt")
ovc <- ieugwasr::associations(id="ieu-a-1122", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ovc,stringsAsFactors=F),outcome="Low grade serous ovarian cancer",population="European",pmid=28346442,ncase=1012,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=113,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###########################
#Advanced prostate cancer #
###########################

out_file<-paste0(out_dir,"apc126","_csi_inst.txt")
pro <- ieugwasr::associations(id="ieu-a-1238", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(pro,stringsAsFactors=F),outcome="Advanced prostate cancer",population="European",pmid=29892016,ncase=15167,ncontrol=58308,study="PRACTICAL",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=126,open_gwas=TRUE,efo = "prostate carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#############################
#Early-onset prostate cancer#
#############################

out_file<-paste0(out_dir,"eop127","_csi_inst.txt")

pro <- ieugwasr::associations(id="ieu-a-1240", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(pro,stringsAsFactors=F),outcome="Early-onset prostate cancer",population="European",pmid=29892016,ncase=6988,ncontrol=44256,study="PRACTICAL",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=127,open_gwas=TRUE,efo = "prostate carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)



#############
#Lung cancer#
#############

out_file<-paste0(out_dir,"luc74","_csi_inst.txt")

luc <- ieugwasr::associations(id="ieu-a-966", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=24880342,ncase=11348,ncontrol=15861,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=74,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer#
#############

out_file<-paste0(out_dir,"luc75","_csi_inst.txt")
luc <- ieugwasr::associations(id="ieu-a-987", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=28604730,ncase=29266,ncontrol=56450,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=75,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer in ever smokers#
#############
out_file<-paste0(out_dir,"lce76","_csi_inst.txt")

luc <- ieugwasr::associations(id="ieu-a-985", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer in ever smokers",population="European",pmid=28604730,ncase=23848,ncontrol=16605,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=76,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)
 
#############
#Lung cancer in never smokers#
#############
out_file<-paste0(out_dir,"lcn77","_csi_inst.txt")

luc <- ieugwasr::associations(id="ieu-a-986", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer in never smokers",population="European",pmid=28604730,ncase=2303,ncontrol=6995,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=77,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#######################
#Lung adenocarcinoma#
#######################

out_file<-paste0(out_dir,"lad72","_csi_inst.txt")

luc <- ieugwasr::associations(id="ieu-a-984", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Lung adenocarcinoma",population="European",pmid=28604730,ncase=11245,ncontrol=54619,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=72,open_gwas=TRUE,efo = "lung adenocarcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#######################
#Squamous cell lung cancer #
#######################

out_file<-paste0(out_dir,"scl79","_csi_inst.txt")

luc <- ieugwasr::associations(id="ieu-a-989", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Squamous cell lung cancer",population="European",pmid=28604730,ncase=7704,ncontrol=54763,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=79,open_gwas=TRUE,efo = "squamous cell lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########################
#Small cell lung cancer#
########################

out_file<-paste0(out_dir,"scl78","_csi_inst.txt")

luc <- ieugwasr::associations(id="ieu-a-988", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(luc,stringsAsFactors=F),outcome="Small cell lung cancer",population="European",pmid=28604730,ncase=2791,ncontrol=20580,study="ILCCO",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=78,open_gwas=TRUE,efo = "small cell lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################
# Japanese Biobank####
########################
# Cervical cancer
out_file<-paste0(out_dir,"cec11","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-98", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Cervical cancer",population="East Asian",pmid=32514122,study="BJ",ncase=605,ncontrol=89731,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=11,open_gwas=TRUE,efo="cervical carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Biliary tract cancer
out_file<-paste0(out_dir,"btc9","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-92", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Biliary tract cancer",population="East Asian",pmid=32514122,study="BJ",ncase=339,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=9,open_gwas=TRUE,efo="cholangiocarcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Colorectal cancer
out_file<-paste0(out_dir,"crc12","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-107", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Colorectal cancer",population="East Asian",pmid=32514122,study="BJ",ncase=7062,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=12,open_gwas=TRUE,efo="colorectal cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Endometrial cancer
out_file<-paste0(out_dir,"enc13","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-113", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Endometrial cancer",population="East Asian",pmid=32514122,study="BJ",ncase=999,ncontrol=89731,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=13,open_gwas=TRUE,efo="endometrial carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Esophageal squamous cell carcinoma
out_file<-paste0(out_dir,"esc14","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-117", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=1300,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=14,open_gwas=TRUE,efo="esophageal squamous cell carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Gastric adenocarcinoma
out_file<-paste0(out_dir,"gac15","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-119", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Gastric adenocarcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=6563,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=15,open_gwas=TRUE,efo="gastric carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Blood cancer
out_file<-paste0(out_dir,"blc10","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-125", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Blood cancer",population="East Asian",pmid=32514122,study="BJ",ncase=1236,ncontrol=211217,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=10,open_gwas=TRUE,efo=c("lymphoma","multiple myeloma","lymphoid leukemia"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Liver cancer
out_file<-paste0(out_dir,"lic16","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-158", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Liver cancer",population="East Asian",pmid=32514122,study="BJ",ncase=1866,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=16,open_gwas=TRUE,efo="hepatocellular carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Lung cancer
out_file<-paste0(out_dir,"luc17","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-133", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Lung cancer",population="East Asian",pmid=32514122,study="BJ",ncase=4050,ncontrol=208403,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=17,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Ovarian cancer
out_file<-paste0(out_dir,"ovc18","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-139", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Ovarian cancer",population="East Asian",pmid=32514122,study="BJ",ncase=720,ncontrol=89731,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=18,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Pancreatic cancer
out_file<-paste0(out_dir,"pac19","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-140", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Pancreatic cancer",population="East Asian",pmid=32514122,study="BJ",ncase=442,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=19,open_gwas=TRUE,efo="pancreatic carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# Prostate cancer
out_file<-paste0(out_dir,"pro20","_csi_inst.txt")

bbj<-ieugwasr::associations(id= "bbj-a-148", variants=snplist) 
dat<-format_data(dat=data.frame(bbj,stringsAsFactors=F),outcome="Prostate cancer",population="East Asian",pmid=32514122,study="BJ",ncase=5408,ncontrol=103939,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=20,open_gwas=TRUE,efo="prostate carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#######################
Gallbladder cancer#
#######################
# <50 cases
# out_file<-paste0(out_dir,"gbc58","_csi_inst.txt")

# gal <- ieugwasr::associations(id="ieu-a-1057", variants=snplist,proxies=0)  
# dat<-format_data(dat=data.frame(gal,stringsAsFactors=F),outcome="Gallbladder cancer",population="East Asian",pmid=22318345,ncase=41,ncontrol=866,study="GCS",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=58,open_gwas=TRUE,efo="gallbladder neoplasm")
# write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########################
#Prostate cancer CGEMS#
########################
# not in open gwas and not independent of PRACTICAL
# no effect allele, other allele or eaf
# snplist<-make_snplist(trait = "Prostate cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
# pro <- ieugwasr::associations(id="ieu-a-823", variants=snplist,proxies=0)  


########
#Glioma#
########

out_file<-paste0(out_dir,"gli967_corrected","_csi_inst.txt")
gli <- ieugwasr::associations(id="ieu-a-1013", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(gli,stringsAsFactors=F),outcome="Glioma",population="European",pmid=22886559,ncase=1856	,ncontrol=4955,study="GliomaScan",UKbiobank=FALSE,rsid="rsid",effect_allele="nea",other_allele="ea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=967,open_gwas=TRUE,efo = "glioma")
# effect_allele="ea",other_allele="nea" QC plots clearly indicate that these are wrong way around

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

################
#"Neuroblastoma#
################
# maf/eaf not reported
out_file<-paste0(out_dir,"neu106","_csi_inst.txt")

neu <- ieugwasr::associations(id="ieu-a-816", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(neu,stringsAsFactors=F),outcome="Neuroblastoma",population="European",pmid=23222812,ncase=1627,ncontrol=3254,study="NBS",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=106,open_gwas=TRUE,efo = "neuroblastoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


################
#Thyroid cancer#
################

out_file<-paste0(out_dir,"thy131","_csi_inst.txt")

thy <- ieugwasr::associations(id="ieu-a-1082", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(thy,stringsAsFactors=F),outcome="Thyroid cancer",population="European",pmid=23894154,ncase=690,ncontrol=497,study="TCS",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=131,open_gwas=TRUE,efo = "thyroid carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#################
#"Upper gastrointestinal cancers || id:825"   #
################
# not in Open GWAS and not independent of N-UGC
# snplist<-make_snplist(trait = "Gastric cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
# gac <- ieugwasr::associations(id="ieu-a-825", variants=snplist,proxies=0)  


#######################################################
# Cancers from the ukb and finn-a batches in Open GWAS#
#######################################################
# Explore FinnGen data at the phenotype level
# https://risteys.finngen.fi/

######################
#Basal cell carcinoma#
######################

out_file<-paste0(out_dir,"bcc135","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-b-8837", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Basal cell carcinoma",population="European",pmid="ukb-b-8837",ncase=4290,ncontrol=458643,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=135,open_gwas=TRUE,efo = "basal cell carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################
#Bladder cancer#
######################
out_file<-paste0(out_dir,"bla136_corrected","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-d-C67", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Bladder cancer",population="European",pmid="ukb-d-C67",ncase=1554,ncontrol=359640,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=136,open_gwas=TRUE,efo = "bladder carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################
#Cancer of digestive organs#
######################

out_file<-paste0(out_dir,"cdo142_corrected","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-d-C3_DIGESTIVE_ORGANS", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Cancer of digestive organs",population="European",pmid="ukb-d-C3_DIGESTIVE_ORGANS",ncase=5690,ncontrol=355504,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=142,open_gwas=TRUE,efo = "digestive system carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################
#Kidney cancer#
######################

out_file<-paste0(out_dir,"kic145_corrected","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-b-1316", variants=snplist,proxies=0)  

dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Kidney cancer",population="European",pmid="ukb-b-1316",ncase=1114,ncontrol=461896,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=145,open_gwas=TRUE,efo = "kidney cancer")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################
#Lymphoma
######################

out_file<-paste0(out_dir,"lym151_corrected","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-d-C_LYMPHOMA", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Lymphoma",population="European",pmid="ukb-d-C_LYMPHOMA",ncase=1752,ncontrol=359442,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=151,open_gwas=TRUE,efo = "lymphoma")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################
#Malignant skin cancer
######################
# Representative examples of malignant skin neoplasms include basal cell carcinoma, squamous cell carcinoma, melanoma, and Kaposi sarcoma.

out_file<-paste0(out_dir,"msc152_corrected","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-d-C3_SKIN", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="ukb-d-C3_SKIN",ncase=16531,ncontrol=344663,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=152,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

######################################
#Respiratory and intrathoracic cancer#
######################################

# C3_NASAL_CAVITY_MIDDLE_EAR, C3_ACCESSORY_SINUS, C3_LARYNX, C3_TRACHEA, C3_BRONCHUS_LUNG, C3_THYMUS, C3_HEART_MEDIASTINUM_PLEURA, C3_RESPIRATORY_INTRATHORACIC3_NAS

out_file<-paste0(out_dir,"ric160_corrected","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-d-C3_RESPIRATORY_INTRATHORACIC", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="ukb-d-C3_RESPIRATORY_INTRATHORACIC",ncase=1944,ncontrol=359250,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=160,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
#Small bowel cancer#
####################

out_file<-paste0(out_dir,"sbc161_corrected.txt","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-a-56", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Small bowel cancer",population="European",pmid="ukb-a-56",ncase=156,ncontrol=337003,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=161,open_gwas=TRUE,efo= "small intestine neuroendocrine tumor")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
#Squamous cell carcinoma#
####################

out_file<-paste0(out_dir,"scc162_corrected","_csi_inst.txt")

ukb <- ieugwasr::associations(id="ukb-a-60", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Squamous cell carcinoma",population="European",pmid="ukb-a-60",ncase=404,ncontrol=336755,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=162,open_gwas=TRUE,efo= "squamous cell carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
#Urinary tract cancer#
####################
out_file<-paste0(out_dir,"utc164_corrected.txt","_csi_inst.txt")
ukb <- ieugwasr::associations(id="ukb-d-C3_URINARY_TRACT", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(ukb,stringsAsFactors=F),outcome="Urinary tract cancer",population="European",pmid="ukb-d-C3_URINARY_TRACT",ncase=1841,ncontrol=359353,study="UKB",UKbiobank=TRUE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=164,open_gwas=TRUE,efo=c("kidney cancer","nephroblastoma","renal cell carcinoma","bladder carcinoma"))
dat<-transform_betas(dat=dat)
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
#Bladder cancer#
####################
out_file<-paste0(out_dir,"bla28","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_BLADDER", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Bladder cancer",population="European",pmid="finn-a-C3_BLADDER",ncase=366,ncontrol=96133,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=28,open_gwas=TRUE,efo= "bladder carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############
#Blood cancer#
##############

out_file<-paste0(out_dir,"blc29","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-CD2_PRIMARY_LYMPHOID_HEMATOPOIETIC", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Blood cancer",population="European",pmid="finn-a-CD2_PRIMARY_LYMPHOID_HEMATOPOIETIC",ncase=1001,ncontrol=95498,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=29,open_gwas=TRUE,efo=c("lymphoma","multiple myeloma","lymphoid leukemia"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############
#Brain cancer#
##############

out_file<-paste0(out_dir,"brc30","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_BRAIN", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Brain cancer",population="European",pmid="finn-a-C3_BRAIN",ncase=142,ncontrol=96357,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=30,open_gwas=TRUE,efo= "central nervous system cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

##############
#Breast cancer
##############

out_file<-paste0(out_dir,"brc31","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_BREAST_3", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Breast cancer",population="European",pmid="finn-a-C3_BREAST_3",ncase=2589,ncontrol=93910,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=31,open_gwas=TRUE,efo= "breast carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

################
#Overall cancer#
################

out_file<-paste0(out_dir,"can32","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-ANY_CANC", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-ANY_CANC",ncase=9792,ncontrol=86707,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=32,open_gwas=TRUE,efo= "cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

################
#Overall cancer#
################

out_file<-paste0(out_dir,"can33","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-II_NEOPLASM", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-II_NEOPLASM",ncase=31217,ncontrol=65282,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=33,open_gwas=TRUE,efo= "cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

################
#Cancer of digestive organs#
################

out_file<-paste0(out_dir,"cdo34","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_DIGESTIVE_ORGANS", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Cancer of digestive organs",population="European",pmid="finn-a-C3_DIGESTIVE_ORGANS",ncase=1582,ncontrol=94917,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=34,open_gwas=TRUE,efo= "digestive system carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#######################################
#Central nervous system and eye cancer#
#######################################

# Malignant neoplasm of eye, brain and central nervous system
# Include	C3_EYE_ADNEXA, C3_MENINGES, C3_BRAIN, C3_SPINAL_CORD_CRANIAL_AND_OTHER_CNS
out_file<-paste0(out_dir,"cns35","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-C3_EYE_BRAIN_NEURO", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Central nervous system and eye cancer",population="European",pmid="finn-a-C3_EYE_BRAIN_NEURO",ncase=207,ncontrol=96292,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=35,open_gwas=TRUE,efo= "central nervous system cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
# Colorectal cancer#
####################

out_file<-paste0(out_dir,"crc36","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-CUSTOM_COLORECTAL_CANCER_EXALLC", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Colorectal cancer",population="European",pmid="finn-a-CUSTOM_COLORECTAL_CANCER_EXALLC",ncase=843,ncontrol=95656,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=36,open_gwas=TRUE,efo= "colorectal cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
# Endocrine gland cancer#
####################

out_file<-paste0(out_dir,"egc37","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-C3_ENDOCRINE", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Endocrine gland cancer",population="European",pmid="finn-a-C3_ENDOCRINE",ncase=328,ncontrol=96171,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=37,open_gwas=TRUE,efo=c("pituitary gland adenoma","thyroid carcinoma","carcinoid tumor","neuroendocrine neoplasm"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
# Endometrial cancer#
####################

out_file<-paste0(out_dir,"enc38","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_CORPUS_UTERI", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Endometrial cancer",population="European",pmid="finn-a-C3_CORPUS_UTERI",ncase=366,ncontrol=53896,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=38,open_gwas=TRUE,efo= "endometrial carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
#Female genital cancer#
####################

out_file<-paste0(out_dir,"fgc39","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_FEMALE_GENITAL", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Female genital cancer",population="European",pmid="finn-a-C3_FEMALE_GENITAL",ncase=672,ncontrol=53590,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=39,open_gwas=TRUE,efo=c("cervical carcinoma","endometrial carcinoma","ovarian carcinoma"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
#Follicular lymphoma#
####################

out_file<-paste0(out_dir,"fol40","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-CD2_FOLLICULAR_LYMPHOMA", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Follicular lymphoma",population="European",pmid="finn-a-CD2_FOLLICULAR_LYMPHOMA",ncase=158,ncontrol=96341,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=40,open_gwas=TRUE,efo="follicular lymphoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

####################
#Kidney cancer#
####################

out_file<-paste0(out_dir,"kid41","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_KIDNEY_NOTRENALPELVIS", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Kidney cancer",population="European",pmid="finn-a-C3_KIDNEY_NOTRENALPELVIS",ncase=301,ncontrol=96198,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=41,open_gwas=TRUE,efo="kidney cancer")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer#
#############

out_file<-paste0(out_dir,"luc42","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-LUNG_CANCER_MESOT", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid="finn-a-LUNG_CANCER_MESOT",ncase=673,ncontrol=95826,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=42,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer#
#############

out_file<-paste0(out_dir,"luc43","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-C3_BRONCHUS_LUNG", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid="finn-a-C3_BRONCHUS_LUNG",ncase=516,ncontrol=95983,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=43,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lymphoid leukaemia#
#############

out_file<-paste0(out_dir,"lyl44","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-CD2_LYMPHOID_LEUKAEMIA", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Lymphoid leukaemia",population="European",pmid="finn-a-CD2_LYMPHOID_LEUKAEMIA",ncase=198,ncontrol=96301,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=44,open_gwas=TRUE,efo="lymphoid leukemia")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
# Male genital cancer
#############

out_file<-paste0(out_dir,"mgc45","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-C3_MALE_GENITAL", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Male genital cancer",population="European",pmid="finn-a-C3_MALE_GENITAL",ncase=1887,ncontrol=94612,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=45,open_gwas=TRUE,efo=c("prostate carcinoma","testicular carcinoma"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
# Malignant skin cancer
#############

out_file<-paste0(out_dir,"msc46","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_SKIN", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="finn-a-C3_SKIN",ncase=895,ncontrol=95604,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=46,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
# Multiple myeloma
#############
out_file<-paste0(out_dir,"mmm47","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-CD2_MULTIPLE_MYELOMA_PLASMA_CELL", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Multiple myeloma",population="European",pmid="finn-a-CD2_MULTIPLE_MYELOMA_PLASMA_CELL",ncase=180,ncontrol=96319,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=47,open_gwas=TRUE,efo="multiple myeloma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
# Non-follicular lymphoma
#############
out_file<-paste0(out_dir,"nfl48","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-CD2_NONFOLLICULAR_LYMPHOMA", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Non-follicular lymphoma",population="European",pmid="finn-a-CD2_NONFOLLICULAR_LYMPHOMA",ncase=344,ncontrol=96155,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=48,open_gwas=TRUE,efo="lymphoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#non-Hodgkin lymphoma unspecified
#############
out_file<-paste0(out_dir,"nhl49","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-CD2_NONHODGKIN_NAS", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="non-Hodgkin lymphoma unspecified",population="European",pmid="finn-a-CD2_NONHODGKIN_NAS",ncase=155,ncontrol=96344,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=49,open_gwas=TRUE,efo="non-Hodgkins lymphoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Oral cavity and pharyngeal cancer
#############


out_file<-paste0(out_dir,"opc50","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-C3_LIP_ORAL_PHARYNX", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Oral cavity and pharyngeal cancer",population="European",pmid="finn-a-C3_LIP_ORAL_PHARYNX",ncase=234,ncontrol=96265,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=50,open_gwas=TRUE,efo=c("head and neck squamous cell carcinoma","oropharynx cancer","nasopharyngeal neoplasm","hypopharynx cancer","oral cavity cancer","mouth neoplasm","pharynx cancer"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Ovarian cancer
#############
out_file<-paste0(out_dir,"ovc51","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_OVARY", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Ovarian cancer",population="European",pmid="finn-a-C3_OVARY",ncase=184,ncontrol=54078,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=51,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Pancreatic cancer
#############
out_file<-paste0(out_dir,"pac52","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_PANCREAS", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Pancreatic cancer",population="European",pmid="finn-a-C3_PANCREAS",ncase=229,ncontrol=96270,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=52,open_gwas=TRUE,efo="pancreatic carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Prostate cancer
#############
out_file<-paste0(out_dir,"pro53","_csi_inst.txt")

fin <- ieugwasr::associations(id="finn-a-C3_PROSTATE", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Prostate cancer",population="European",pmid="finn-a-C3_PROSTATE",ncase=1824,ncontrol=40413,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=53,open_gwas=TRUE,efo="prostate carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Respiratory and intrathoracic cancer
#############

out_file<-paste0(out_dir,"ric54","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_RESPIRATORY_INTRATHORACIC", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="finn-a-C3_RESPIRATORY_INTRATHORACIC",ncase=615,ncontrol=95884,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=54,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Thyroid cancer
#############

out_file<-paste0(out_dir,"thc55","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_THYROID_GLAND", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Thyroid cancer",population="European",pmid="finn-a-C3_THYROID_GLAND",ncase=321,ncontrol=96178,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=55,open_gwas=TRUE,efo="thyroid carcinoma")
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

#############
#Urinary tract cancer
#############

out_file<-paste0(out_dir,"utc56","_csi_inst.txt")
fin <- ieugwasr::associations(id="finn-a-C3_URINARY_TRACT", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(fin,stringsAsFactors=F),outcome="Urinary tract cancer",population="European",pmid="finn-a-C3_URINARY_TRACT",ncase=690,ncontrol=95809,study="FinnGen",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=56,open_gwas=TRUE,efo=c("kidney cancer","nephroblastoma","renal cell carcinoma","bladder carcinoma"))
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

# pancreatic cancer
# eaf/maf not reported
out_file<-paste0(out_dir,"pan122","_csi_inst.txt")
pan <- ieugwasr::associations(id="ieu-a-822", variants=snplist,proxies=0)  
dat<-format_data(dat=data.frame(pan,stringsAsFactors=F),outcome="Pancreatic cancer",population="European",pmid=19648918,ncase=1896,ncontrol=1939,study="PanScan I",UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",lnor="beta",lnor_se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=122,open_gwas=TRUE,efo="pancreatic carcinoma")
write.table(dat,out_file,sep="\t",col.name=T,row.name=F,quote=F)

######################################
# datasets obtained by correspondence#
######################################
library(CheckSumStats)
out_dir<-"~/fatty-acids/mr/data/csi_inst/"
snplist<-readLines("~/fatty-acids/mr/data/Csi_sig_unique_rsidonly_clumped.txt")


####################
# Colorectal cancer 
####################
out_file<-paste0(out_dir,"crc_60","_csi_inst.txt")
Crc<-read.csv("/Users/ph14916/MR_FattyAcids/data/summary_data/colorectal_cancer/smoking/1255_MarginalResults_HRC125K_20210728.csv",head=TRUE,stringsAsFactors=FALSE)
dat<-format_data(dat=Crc,outcome="Colorectal cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=58131,ncontrol=67347,UKbiobank=TRUE,rsid="rs",effect_allele="Allele1",other_allele="Allele2",lnor="Effect",lnor_se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=60,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",z_score="LRT.Z")
dat<-dat[dat$rsid %in% snplist,]
dat$ID<-dat$id
# head(dat$id)
write.table(dat,out_file,sep="\t",col.name=T,row.name=F,quote=F)

#########################
# cervical cancer ID=129#
##########################
out_file<-paste0(out_dir,"cervical_129","_csi_inst.txt")
cer<-preformat_cer_id129_v2()
dat<-format_data(dat=cer,outcome="Cervical cancer",population="European",pmid=23482656,study="SCCS",ncase=1034,ncontrol=3948,UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele="A2",or="OR",lnor_se="SE",eaf="MAF",effect_allele_confirmed=TRUE,ID=129)
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.name=T,row.name=F,quote=F)

######################
# glioma SFAGs ID=133#
######################

out_file<-paste0(out_dir,"glioma_133","_csi_inst.txt")
File<-"/Users/ph14916/MR_FattyAcids/data/summary_data/glioma_SFAGS_28346443/smoking_interaction/sfags_illumina_haycock__smokingsnps_results_v2.txt"
gli<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=gli,outcome="Glioma",population="European",pmid=28346443,study="UCSF_AGS + SFAGS",ncase=677,ncontrol=3940,UKbiobank=FALSE,rsid="rsid",effect_allele="allele_B",other_allele="allele_A",lnor="beta",lnor_se="se",p="P_value",info1="info",effect_allele_confirmed=TRUE,ID=133,all_summary_stats=FALSE,HWEp="controls_hwe",eaf=NA)
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.name=T,row.name=F,quote=F)

########################
#Hodgkin lymphoma ID=69#
########################

out_file<-paste0(out_dir,"hodgkin_69","_csi_inst.txt")
File<-"~/MR_FattyAcids/data/summary_data/Hodgkin_lymphoma_29196614/smoking_interactions/snp_data_haycock_HL.txt"
Hdl<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=Hdl,outcome="Hodgkin lymphoma",population="European",pmid=29196614,study="HLS",ncase=3077,ncontrol=13680,UKbiobank=FALSE,rsid="rsid",effect_allele="allele_B",other_allele="allele_A",lnor="beta",lnor_se="se",eaf="coded_af",p="P_value",info1="cohort_1_info",info2= "cohort_2_info",info3="cohort_3_info",effect_allele_confirmed=TRUE,ID=69,phet="P_heterogeneity",I2="I2",Q="Q")
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

########################
# Kidney cancer ID = 89#
########################
out_file<-paste0(out_dir,"kidney_89","_csi_inst.txt")
File<-"~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/smoking_interaction/RCC/Smoking_0321_Kidney_IARC_NCIs_UK_MDA_2021-03-15.txt"
kid<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=kid,outcome="Kidney cancer" ,population="European",pmid=28598434,study="KidRISK",UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele= "A2",lnor="Effect",lnor_se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info1="Quality",ncase=10784 ,ncontrol=20406,ID=89) 
dat<-dat[dat$rsid %in% snplist,]

write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

###########################
# multiple myeloma ID=96.2#
###########################
out_file<-paste0(out_dir,"multiplemyeloma_96_2.txt","_csi_inst.txt")
File<-"~/MR_FattyAcids/data/summary_data/multiple_myeloma_Duran_Lozano/smoking_results.txt"
mma<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=mma,outcome="Multiple myeloma",population="European",pmid="unpublished",study="MMS",ncase=2338,ncontrol=11971 ,UKbiobank=FALSE,rsid="rsName",effect_allele="EA",other_allele="OA",or="OR",lnor_se="SE",eaf="EAfrq",p="P",effect_allele_confirmed=TRUE,ID=96.2,info1="Info")
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)

############################## 
#head and neck cancer Oral ID=80###
##############################
out_file<-paste0(out_dir,"oral80","_csi_inst.txt")
File<-"~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/smoking_interaction/HNC/Smoking_0321_HeadNeck_OC_2021-03-15.txt"
hnc<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=hnc,outcome="Oral cancer",population="European",pmid=27749845,study="INHANCE",UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele= "A2",lnor="Effect",lnor_se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info1="Quality",ncase=2990 ,ncontrol=6585,ID=80) 
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


############################## 
#head and neck cancer overall ID=81#
##############################
out_file<-paste0(out_dir,"headneck_81","_csi_inst.txt")
File<-"~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/smoking_interaction/HNC/Smoking_0321_HeadNeck_HNC_2021-03-15.txt"
hnc<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=hnc,outcome="Oral cavity and pharyngeal cancer",population="European",pmid=27749845,study="INHANCE",UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele= "A2",lnor="Effect",lnor_se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info1="Quality",ncase=6034 ,ncontrol=6585,ID=81) 
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


############################## 
#head and neck cancer oropharyngeal ###
##############################
out_file<-paste0(out_dir,"oropharyngeal_82","_csi_inst.txt")
File<-"~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/smoking_interaction/HNC/Smoking_0321_HeadNeck_OPC_2021-03-15.txt"
hnc<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(dat=hnc,outcome="Oropharyngeal cancer",population="European",pmid=27749845,study="INHANCE",UKbiobank=FALSE,rsid="SNP",effect_allele="A1",other_allele= "A2",lnor="Effect",lnor_se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info1="Quality",ncase=2641 ,ncontrol=6585,ID=82) 
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)


#######################
# Thyroid cancer####
####################
out_file<-paste0(out_dir,"thyroid_26","_csi_inst.txt")
File<-"~/MR_FattyAcids/data/summary_data/thyroid_cancer/smoking_interaction/cpd_csi_snplist_v2_extraction.txt"

thy<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
Names<-names(thy)[names(thy)!="snp.id"]
names(thy)<-Names
thy$snp.id<-row.names(thy)

dat<-format_data(dat=thy,outcome="Thyroid cancer",population="European",pmid="unpublished",study="EPITHYR",ncase=1554,ncontrol=1973,UKbiobank=FALSE,rsid="snp.id",effect_allele="A2",other_allele="A1",lnor="beta.all",lnor_se="se.all",eaf="freq2",p="p.value.all",effect_allele_confirmed=TRUE,ID=26,info1="info")
dat<-dat[dat$rsid %in% snplist,]
write.table(dat,out_file,sep="\t",col.names=T,row.names=F,quote=F)
