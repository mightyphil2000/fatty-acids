
###################################################################
# studies with full GWAS summary statistics available not in OpenGWAS##############
###################################################################

# ####################
# Glioma 28346443#####
######################

snplist<-make_snplist(efo="glioma",population="European")
Gli<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/research/projects/icep1/bioinformatics/013/working/data/new_data/meta_gicc_mda_allglioma.TBL",exact_match=TRUE)
dat<-format_data(Dat=Gli,outcome="Glioma",population="European",pmid=28346443,study="GICC/MDA",ncase=5747,ncontrol=5522,UKbiobank=FALSE,rsid="MarkerName",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="P.value",info=NA,info1=NA,info2=NA,info3=NA,HWEp=NA,phet=NA,I2=NA,Q=NA,Direction=NA,effect_allele_confirmed=TRUE,or=NA,lci=NA,uci=NA,ref=NULL,chr="CHR",pos="POS",test_statistic=NA,ID=66,all_summary_stats=TRUE,efo="glioma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/glioma/gli66.txt",sep="\t",col.names=T,row.names=F,quote=F)


# cd ~/fatty-acids/outcome_data/data/harmonised
# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/glioma/glioma66.txt .   

# write.table(Gbm,"Final/glioblastoma.txt",sep="\t",col.names=T,row.names=F,quote=F)
# write.table(nongbm,"Final/nonglioblastoma.txt",sep="\t",col.names=T,row.names=F,quote=F)

########################
# Endometrial cancer####
########################

snplist<-make_snplist(efo="endometrial carcinoma",population="European")
End<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/endometrial_cancer/ecac_gwas_data.txt",exact_match=TRUE,file_sep=" ")
dat<-format_data(Dat=End,outcome="Endometrial cancer",population="European",pmid=30093612,study="ECAC",ncase=12906,ncontrol=108979,UKbiobank=TRUE,rsid="SNPID",Effect.Allele="EA",Other.Allele="OA",lnor="ALL_BETA",se="ALL_SE",eaf="MEAN_EAF",p="ALL_PVALUE",info="OA_INFO",effect_allele_confirmed=TRUE,ID=25,all_summary_stats=TRUE,efo="endometrial carcinoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/endometrial_cancer/enc25.txt",sep="\t",col.names=T,row.names=F,quote=F)


# cd ~/MR_FattyAcids/data/summary_data/Final
# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/endometrial_cancer/endometrial.txt . 

# ############################################
# acute lymphoblastic leukemia 22076464/####
############################################

snplist<-make_snplist(efo="acute lymphoblastic leukemia",population="European")
All<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/acute_lymphoblastic_leukemia_22076464/Assoc_ALL_Affy5.assoc",exact_match=TRUE)
dat<-format_data(Dat=All,outcome="Acute lymphoblastic leukaemia",population="European",pmid=22076464,study="C-ALL",ncase=419,ncontrol=474,UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele="A2",or="OR",se="SE",lci="L95",uci="U95",eaf="F_U",p="P",effect_allele_confirmed=TRUE,ID=21,all_summary_stats=TRUE,efo="acute lymphoblastic leukemia")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/acute_lymphoblastic_leukemia_22076464/All21.txt",sep="\t",col.names=T,row.names=F,quote=F)

################################
# esophageal adenocarcinoma####
################################
# no eaf or maf reported
snplist<-make_snplist(efo="esophageal adenocarcinoma",population="European")
Eso<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/esophageal_adenocarcinoma/Oesophageal_adenocarcinoma_Lancet_Oncol2016_Gharahkhani_et_al.txt",exact_match=TRUE,file_sep=" ")

dat<-format_data(Dat=Eso,outcome="Esophageal adenocarcinoma",population="European",pmid=27527254,study="EAS",ncase="N_cases",ncontrol="N_controls",UKbiobank=FALSE,rsid="MarkerName",Effect.Allele="effect_allele",Other.Allele="non.effect_allele",lnor="Effect",se="StdErr",p="P.value",effect_allele_confirmed=TRUE,ID=24,all_summary_stats=TRUE,I2="HetISq",Q="HetChiSq",phet="HetPVal",Direction="Direction",efo="esophageal adenocarcinoma")

write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/esophageal_adenocarcinoma/esa24.txt",sep="\t",col.names=T,row.names=F,quote=F)


##################
# Bladder cancer###
#####################
# chrALL* creared using bash script in pre_harmonise_outcomes.sh
snplist<-make_snplist(efo="bladder carcinoma",population="European")
bla<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/bladder_cancer/Haycock/chrALL_GWAS_NBCS_risico_jan2017_methodscore_4mds_info.out",exact_match=TRUE,file_sep=" ")
bla<-preformat_bla_id105()
dat<-format_data(Dat=bla,outcome="Bladder cancer",population="European",pmid=21750109,study="NBCS",ncase=1799,ncontrol=4745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="allele_B",Other.Allele="allele_A",eaf="eaf",lnor="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_beta_1",se="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_se_1",p="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_pvalue",effect_allele_confirmed=TRUE,info="UBC_yn_frequentist_add_mds1_mds2_mds3_mds4_score_info",ID=105,all_summary_stats=TRUE,efo="bladder carcinoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/bladder_cancer/Haycock/blc105.txt",sep="\t",col.names=T,row.names=F,quote=F)

##################
# uveal melanoma###
#####################

snplist<-make_snplist(efo="uveal melanoma",population="European")
uve<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/uveal_melanoma/GWAS_results_17122018.txt",exact_match=TRUE,file_sep="\t")
dat<-format_data(Dat=uve,outcome="Uveal melanoma",population="European",pmid=28781888,study="UMS",ncase=259,ncontrol=401,UKbiobank=FALSE,rsid="SNP",Effect.Allele="MINOR",eaf="MAF",or="OR",se="SE",lci="L95",uci="U95",p="P",effect_allele_confirmed=TRUE,HWEp="P_HWE",test_statistic="STAT",ID=165,all_summary_stats=TRUE,efo="uveal melanoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/uveal_melanoma/uvm165.txt",sep="\t",col.names=T,row.names=F,quote=F)


################
# Melanoma  ####
################

# A<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/melanoma/Melanoma_meta_single_files/Melanoma_metaanalysis_chrALL_15052019_EAF_RSQ.txt",sep="\t",head=T,stringsAsFactors=FALSE)

snplist<-make_snplist(efo="melanoma",population="European")
# Melanoma_metaanalysis_chrALL* creared using bash script in pre_harmonise_outcomes.sh
mel<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/melanoma/Melanoma_meta_single_files/Melanoma_metaanalysis_chrALL_15052019_EAF_RSQ.txt",exact_match=TRUE,file_sep=" ")
dat<-format_data(Dat=mel,outcome="Melanoma",population="European",pmid=26237428,study="MMAC",ncase="CASE",ncontrol="CONTROL",UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele="A2",or="OR",se="SE_fixed_qnorm",p="P",effect_allele_confirmed=TRUE,Q="Q",test_statistic="Z_fixed_sign",eaf="eaf",I2="I",info="RSQ_median",ID=95,all_summary_stats=TRUE,efo="melanoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/melanoma/mel95.txt",sep="\t",col.names=T,row.names=F,quote=F)


###############
# UK Biobank###
###############
snplist<-make_snplist(efo= "central nervous system cancer",population="European")	
ukb<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_brain_cancer_imputed.txt.gz",exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,outcome="Brain cancer",population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",ID=138,all_summary_stats=TRUE,efo="central nervous system cancer")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/brc138.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo= "breast carcinoma",population="European")	
ukb<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_breast_cancer_imputed.txt.gz",exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,outcome="Breast cancer",population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",ID=139,all_summary_stats=TRUE,efo="breast carcinoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/brc139.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo= "colorectal cancer",population="European")	
ukb<-extract_snps(snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE,
	File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_colorectal_cancer_imputed.txt.gz")
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Colorectal cancer",efo="colorectal cancer",ID=143
	)

write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/crc143.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist_blood(population="European",Dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/")
ukb<-extract_snps(snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE,
	File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_haem_cancer_imputed.txt.gz")
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Blood cancer",	efo=c("lymphoma","multiple myeloma","lymphoid leukemia"),ID=137)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/blc137.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist_hnc(population="European",Dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/")
ukb<-extract_snps(snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE,
	File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_headneck_cancer_imputed.txt.gz")
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Oral cavity and pharyngeal cancer",efo=c("head and neck squamous cell carcinoma","oropharynx cancer","nasopharyngeal neoplasm","hypopharynx cancer","oral cavity cancer","mouth neoplasm","pharynx cancer"),ID=157)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/opc157.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "acute lymphoblastic leukemia",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_leuk_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Leukaemia",efo="acute lymphoblastic leukemia",ID=146)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/leu146.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo= "hepatocellular carcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_liver_bile_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Liver & bile duct cancer",efo="hepatocellular carcinoma",ID=147)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/lbc147.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo= "hepatocellular carcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_liver_cell_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Liver cancer",efo="hepatocellular carcinoma",ID=148)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/lic148.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo= "lung carcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lung_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Lung cancer",efo="lung carcinoma",ID=149)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/luc149.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "lung carcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lung_cancer_unadj_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Lung cancer unadjusted for chip",efo="lung carcinoma",ID=1499
	)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/luc1499.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "acute lymphoblastic leukemia",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lymph_leuk_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Lymphoid leukaemia",efo="acute lymphoblastic leukemia",ID=150)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/lle150.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo= "multiple myeloma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_mult_myel_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Multiple myeloma",efo="multiple myeloma",ID=154)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/mum154.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "acute myeloid leukemia",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_myel_leuk_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Myeloid Leukaemia",efo="acute myeloid leukemia",ID=155)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/myl155.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "non-melanoma skin carcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_nm_skin_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Non-melanoma skin cancer",efo="non-melanoma skin carcinoma",ID=156)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/nmc156.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "esophageal adenocarcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_oesoph_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Esophageal adenocarcinoma",efo="esophageal adenocarcinoma",ID=144)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/esa144.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "ovarian carcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_ovarian_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Ovarian cancer",efo="ovarian carcinoma",ID=158)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/ovc158.txt",sep="\t",col.names=T,row.names=F,quote=F)



snplist<-make_snplist(efo= "cancer",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Overall cancer (excluding non-melanoma skin cancer)",efo="cancer",ID=141)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/oac141.txt",sep="\t",col.names=T,row.names=F,quote=F)



snplist<-make_snplist(efo= "cancer",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_pan_inclc44_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Overall cancer",efo="cancer",ID=140
	)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/oac140.txt"
	,sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo= "prostate carcinoma",population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_prostate_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Prostate cancer",efo="prostate carcinoma",ID=159
	)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/pro159.txt"
	,sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo= "melanoma"	,population="European")	
ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_skin_cancer_imputed.txt.gz"
	,snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Melanoma",	efo="melanoma",ID=153)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/mel153.txt",sep="\t",col.names=T,row.names=F,quote=F)

###################################
# Upper gastrointestinal cancers#####
########################################
# created files using pre_harmonise_outcomes.sh
snplist<-make_snplist(efo="gastric adenocarcinoma",population="East Asian")	
ugi<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_CC/summary_chr_all.txt",exact_match=TRUE,file_sep="\t",Comment="")
dat<-format_data(Dat=ugi,outcome="Gastric cardia adenocarcinoma",population="East Asian",pmid=26129866,ncase=1189,ncontrol=2708,study="N-UGC",UKbiobank=FALSE,rsid="rs",Effect.Allele="risk_allele",Other.Allele="reference_allele",lnor="beta",se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info="info",ID=102,all_summary_stats=TRUE,efo="gastric adenocarcinoma") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/gca102.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo="esophageal squamous cell carcinoma",population="East Asian")
ugi<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/summary_chr_all.txt",exact_match=TRUE,file_sep="\t",Comment="")
dat<-format_data(Dat=ugi,outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=25129146,ncase=2013,ncontrol=2701,study="N-UGC",UKbiobank=FALSE,rsid="rs",Effect.Allele="risk_allele",Other.Allele="reference_allele",lnor="beta",se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info="info",ID=99,all_summary_stats=TRUE,efo="esophageal squamous cell carcinoma") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/esc99.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo="Gastric adenocarcinoma",population="East Asian")	
ugi<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_gastric/summary_chr_all.txt",exact_match=TRUE,file_sep="\t",Comment="")
dat<-format_data(Dat=ugi,outcome="Gastric adenocarcinoma",population="East Asian",pmid=26129866,ncase=2350,ncontrol=2708,study="N-UGC",UKbiobank=FALSE,rsid="rs",Effect.Allele="risk_allele",Other.Allele="reference_allele",lnor="beta",se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info="info",ID=101,all_summary_stats=TRUE,efo="gastric adenocarcinoma") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/gac101.txt",sep="\t",col.names=T,row.names=F,quote=F)


snplist<-make_snplist(efo="Gastric adenocarcinoma",population="East Asian")		
ugi<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_NC/summary_chr_all.txt",exact_match=TRUE,file_sep="\t",Comment="")
dat<-format_data(Dat=ugi,outcome="Noncardia gastric adenocarcinoma",population="East Asian",pmid=26129866,ncase=1027,ncontrol=2708,study="N-UGC",UKbiobank=FALSE,rsid="rs",Effect.Allele="risk_allele",Other.Allele="reference_allele",lnor="beta",se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info="info",ID=103,all_summary_stats=TRUE,efo="gastric adenocarcinoma") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/nga103.txt",sep="\t",col.names=T,row.names=F,quote=F)


###################################################
# Additional datasets identified in GWAS catalog###
######################################################

####################
# cervical cancer####
####################
# eaf/ maf not reported
snplist<-make_snplist(efo="cervical carcinoma",population="European")
cec<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Cervical_cancer_28806749/28806749-GCST004833-EFO_0001061.h.tsv",exact_match=TRUE,file_sep="\t")
head(cec)
dat<-format_data(Dat=cec,outcome="Cervical cancer",population="European",pmid=28806749,study="MCCS",ncase=2866,ncontrol=6481,UKbiobank=FALSE,rsid="variant_id",Effect.Allele="effect_allele",Other.Allele="other_allele",eaf="effect_allele_frequency",or="odds_ratio",se="standard_error",p="p_value",effect_allele_confirmed=TRUE,test_statistic="stat",ID=92,all_summary_stats=TRUE,efo="cervical carcinoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Cervical_cancer_28806749/cec92.txt",sep="\t",col.names=T,row.names=F,quote=F)


##########################################
# BRCA negative high risk breast cancer 30323354#
##########################################
# eaf/ maf not reported

snplist<-make_snplist(efo="breast carcinoma",population="East Asian")
brc<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/BRCA1_2_negative_high_risk_breast_cancer_30323354/30323354-GCST006719-EFO_0009443-build37.f.tsv",exact_match=TRUE,file_sep="\t")

dat<-format_data(Dat=brc,outcome="BRCA negative high risk breast cancer",population="East Asian",pmid=30323354,study="KHBC",ncase=1469,ncontrol=5979,UKbiobank=FALSE,rsid="variant_id",Effect.Allele="effect_allele",Other.Allele=NA,or="odds_ratio",eaf=NA,p="gc",effect_allele_confirmed=TRUE,test_statistic="stat",ID=88,all_summary_stats=TRUE,efo="breast carcinoma") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/BRCA1_2_negative_high_risk_breast_cancer_30323354/hrb88.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Dat$Z<-qnorm(Dat$gc/2,lower.tail=F)#gc corrected p value. I looked up the top SNP from the paper and the gc value corresponds exactly to the reported P value for this SNP


##########################
# thyroid cancer #
##########################

snplist<-make_snplist(efo="thyroid carcinoma",population="European")
thy<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Thyroid_cancer_30104761/PheCode_193_SAIGE_MACge20.txt.vcf",exact_match=TRUE,file_sep="\t",Comment="")
dat<-format_data(Dat=thy,outcome="Thyroid cancer",population="European",pmid=30104761,study="UKB",ncase="num_cases",ncontrol="num_controls",UKbiobank=TRUE,rsid="ID",Effect.Allele="ALT",Other.Allele="REF",lnor="beta",se="sebeta",eaf="af",p="pval",effect_allele_confirmed=FALSE,all_summary_stats=TRUE,ID=163,test_statistic="Tstat",efo="thyroid carcinoma")  # GAME-ON had same names for allele columns.  
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Thyroid_cancer_30104761/thc163.txt",sep="\t",col.names=T,row.names=F,quote=F)



####################################################
#Kidney cancer in females
####################################################
snplist<-make_snplist(efo="kidney cancer",population="European")
kid<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Kidney_cancer_31231134/Laskar_31231134_Females",exact_match=TRUE,file_sep="\t",Comment="")
kid<-preformat_kid_90(dat=kid)
dat<-format_data(Dat=kid,outcome="Kidney cancer in females",population="European",pmid=31231134,study="KidRISK",ncase="Number_cases",ncontrol="Number_controls",UKbiobank=FALSE,rsid="Variant_ID",Effect.Allele="Effect_allele",Other.Allele="Other_allele",or="Odds_ratio",se="Standard_error",eaf="effect_allele_frequency_controls",p="P_Value",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=90,efo="kidney cancer")  
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Kidney_cancer_31231134/kif90",sep="\t",col.names=T,row.names=F,quote=F)


####################################################
#Kidney cancer in males
####################################################
snplist<-make_snplist(efo="kidney cancer",population="European")
kid<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Kidney_cancer_31231134/Laskar_31231134_Males",exact_match=TRUE,file_sep="\t",Comment="")
kid<-preformat_kid_90(dat=kid)
dat<-format_data(Dat=kid,outcome="Kidney cancer in males",population="European",pmid=31231134,study="KidRISK",ncase="Number_cases",ncontrol="Number_controls",UKbiobank=FALSE,rsid="Variant_ID",Effect.Allele="Effect_allele",Other.Allele="Other_allele",or="Odds_ratio",se="Standard_error",eaf="effect_allele_frequency_controls",p="P_Value",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=91,efo="kidney cancer") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Kidney_cancer_31231134/kim91",sep="\t",col.names=T,row.names=F,quote=F)

#################################################### 
# Ovarian_cancer_EastAsians_30898391 OCAC
####################################################
# input file created using pre_harmonise_outcomes.sh
snplist<-make_snplist(efo="ovarian carcinoma",population="East Asian")
ova<-extract_snps(snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Ovarian_cancer_EastAsians_30898391/SummaryResults_Asian_chr_all.txt",exact_match=TRUE,file_sep=",",Comment="")
ova<-preformat_ova_120(dat=ova)
summary(ova$R2_oncoarray[ova$R2_oncoarray!=-99])
dat<-format_data(Dat=ova,outcome="Ovarian cancer",population="East Asian",pmid=30898391,study="OCAC (EAS)",ncase=3238,	ncontrol=4083,UKbiobank=FALSE,rsid="V4",Effect.Allele="Effect",Other.Allele="Baseline",lnor="overall_OR",se="overall_SE",eaf="EAF",p="overall_pvalue",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=120,efo="ovarian carcinoma") 
# dat[dat$lnor < -2,c("lnor","eaf")]
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Ovarian_cancer_EastAsians_30898391/ovc120.txt",sep="\t",col.names=T,row.names=F,quote=F)


###########################################
# Bcell_nonHodgkinlymphoma_Chinese_23749188#
###########################################

snplist<-make_snplist(efo="neoplasm of mature b-cells",population="East Asian")
nhl<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Bcell_nonHodgkinlymphoma_Chinese_23749188/NHL_253cases_1438controls_536555snps_summary_stats.txt",exact_match=TRUE,file_sep="\t")
head(nhl)
dat<-format_data(Dat=nhl,outcome="B cell non-Hodgkin lymphoma",population="East Asian",pmid=23749188,study="BC-NHL",ncase=253,ncontrol=1438,UKbiobank=FALSE,rsid="SNP",Effect.Allele="TestAllele",Other.Allele="MajorAlelle",eaf="TAF_Controls",or="OR",se="SE",p="P",effect_allele_confirmed=TRUE,ID=5,all_summary_stats=TRUE,efo="neoplasm of mature b-cells")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Bcell_nonHodgkinlymphoma_Chinese_23749188/bnh5.txt",sep="\t",col.names=T,row.names=F,quote=F)

###########################
# Non melanoma skin cancer#
###########################

# source("~/fatty-acids-mr/instruments/Extract_SNPs_function.R")
# bcc<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_bcc/Chahal_2016_basal_cell_carcinoma-4.1/basal_cell_carcinoma-4.1.dat",exact_match=TRUE,file_sep="\t")

snplist<-make_snplist(efo="basal cell carcinoma",population="European")
dat<-extract_snps_and_format_bcc_23andMe(snplist=snplist)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/bcc1.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-make_snplist(efo="cutaneous squamous cell carcinoma",population="European")
dat<-extract_snps_and_format_scc_23andMe(snplist=snplist)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/scc2.txt",sep="\t",col.names=T,row.names=F,quote=F)


# ############
#Ewing sarcoma#
##############
# ewi<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/ewing_sarcoma/Postel-Vinay22327514/data/Ewing_association_Apr2017.csv",exact_match=TRUE,file_sep=",")

# MAF reported but unknown whether corresponds to EAF 
snplist<-make_snplist(efo="Ewing sarcoma",population="European")
ewi<-preformat_ewi_27(snplist=snplist)

dat<-format_data(Dat=ewi,outcome="Ewing's sarcoma",population="European",pmid=22327514,study="ESS",ncase=427,ncontrol=684,UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1.controls",Other.Allele="A2.controls",or="OR",se="SE",p="P",effect_allele_confirmed=TRUE,ID=27,all_summary_stats=TRUE,test_statistic="STAT",HWEp="HWE.Pval.UNAFF",efo="Ewing sarcoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/ewings_sarcoma/ews27.txt",sep="\t",col.names=T,row.names=F,quote=F)


############################################
# Chronic lymphocytic leukemia Interlymph###
############################################
# input file prepared using pre_harmonise_outcomes.sh

snplist<-make_snplist(efo="chronic lymphocytic leukemia",population="European")
cll<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/CLL/cll.out",exact_match=TRUE,file_sep="\t",fill=TRUE)
dat<-preformat_snps_interlymph(dat=cll)
dat<-format_data(Dat=dat,outcome="Chronic lymphocytic leukaemia",population="European",pmid=26956414,study="InterLymph",ncase="NUM_CASE",ncontrol="NUM_CONTROL",eaf="eaf",UKbiobank=FALSE,rsid="SNP",Effect.Allele="Effect.Allele",Other.Allele="Other.Allele",lnor="lnor",se="se",p="P",effect_allele_confirmed=TRUE,ID=83,all_summary_stats=TRUE,info="info",Direction="Direction",phet="Phet",I2="I2",efo="chronic lymphocytic leukemia")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/cll83.txt",sep="\t",col.names=T,row.names=F,quote=F)
head(dat)


#####################################
# Follicular lymphoma Interlymph ####
#####################################
# input file prepared using pre_harmonise_outcomes.sh
snplist<-make_snplist(efo="follicular lymphoma",population="European")
fl<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/FL/fl.out",exact_match=TRUE,file_sep="\t",fill=TRUE)
dat<-preformat_snps_interlymph(dat=fl)
dat<-format_data(Dat=dat,outcome="Follicular lymphoma",population="European",pmid=25279986,study="InterLymph",ncase="NUM_CASE",ncontrol="NUM_CONTROL",eaf="eaf",UKbiobank=FALSE,rsid="SNP",Effect.Allele="Effect.Allele",Other.Allele="Other.Allele",lnor="lnor",se="se",p="P",effect_allele_confirmed=TRUE,ID=85,all_summary_stats=TRUE,info="info",Direction="Direction",phet="Phet",I2="I2",efo="follicular lymphoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/fll85.txt",sep="\t",col.names=T,row.names=F,quote=F)


###########################################
# Diffuse large B cell lymphoma interlymph####
###########################################
# input file prepared using pre_harmonise_outcomes.sh
snplist<-make_snplist(efo="diffuse large b-cell lymphoma",population="European")
dlbcl<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/DLBCL/dlbcl.out",exact_match=TRUE,file_sep="\t",fill=TRUE)
dat<-preformat_snps_interlymph(dat=dlbcl)
dat<-format_data(Dat=dat,outcome="Diffuse large B cell lymphoma",population="European",pmid=25261932,study="InterLymph",ncase="NUM_CASE",ncontrol="NUM_CONTROL",eaf="eaf",UKbiobank=FALSE,rsid="SNP",Effect.Allele="Effect.Allele",Other.Allele="Other.Allele",lnor="lnor",se="se",p="P",effect_allele_confirmed=TRUE,ID=84,all_summary_stats=TRUE,info="info",Direction="Direction",phet="Phet",I2="I2",efo="diffuse large b-cell lymphoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/dlb84.txt",sep="\t",col.names=T,row.names=F,quote=F)
head(dat)


#####################################
# Marginal zone lymphoma        ####
#####################################
snplist<-make_snplist(efo="marginal zone B-cell lymphoma",population="European")
mzl<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/MZL/mzl.out",exact_match=TRUE,file_sep="\t",fill=TRUE)
dat<-format_data(Dat=mzl,outcome="Marginal zone lymphoma",population="European",pmid=25569183,study="InterLymph",ncase="Num_Case",ncontrol="Num_Control",eaf="Effect_Allele_Freq_Control",UKbiobank=FALSE,rsid="Locus",Effect.Allele="Effect_Allele",Other.Allele="Reference_Allele",lnor="Beta",se="Standard_error_of_beta",p="P_value",effect_allele_confirmed=TRUE,ID=86,all_summary_stats=TRUE,info="Info",efo="marginal zone B-cell lymphoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/mzl86.txt",sep="\t",col.names=T,row.names=F,quote=F)

# cp /projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Rajaraman22886559/for_Philip_data_delivery.txt .
##########################
# GliomaScan #
##########################
snplist<-make_snplist(efo="glioma",population="European")
gli<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gliomascan/for_Philip_data_delivery.txt",exact_match=TRUE,file_sep="\t",Comment="")
gli<-preformat_gli_67(dat=gli)
dat<-format_data(Dat=gli,outcome="Glioma",population="European",pmid=22886559,study="GliomaScan",ncase="cases",ncontrol="controls",UKbiobank=FALSE,rsid="Locus",Effect.Allele="Allele2",Other.Allele="Allele1",or="OR",lci="OR_95._CI_l",uci="OR_95._CI_u",eaf="eaf.controls",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=67,efo="glioma")   #confirmed by correspondence 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gliomascan/gli67.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################################
# Prostate cancer PRACTICAL Build 37#
######################################
snplist<-make_snplist(efo="prostate carcinoma",population="European",bed_37=FALSE)
#input file created using preformat_pro_128()
pro<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/practical/meta_v3_onco_euro_overall_1_fixed_plusrsids.txt",exact_match=TRUE,file_sep="\t",Comment="")
dat<-format_data(Dat=pro,outcome="Prostate cancer",population="European",pmid=29892016,ncase=79148,ncontrol=61106,study="PRACTICAL",UKbiobank=FALSE,rsid="V4",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="P.value",Direction = "Direction",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=128,efo="prostate carcinoma") #effect allele = allele 1 as confirmed in documentation ~/OpenGWAS/PRACTICAL/documentations
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/practical/pro128.txt",sep="\t",col.names=T,row.names=F,quote=F)


####################################
# hepatocellular_carcinoma_22174901##
######################################

snplist<-make_snplist(efo="hepatocellular carcinoma",population="European")
hpc<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/livercancer/HKU_data_summary.txt",exact_match=TRUE)
dat<-format_data(Dat=hpc,outcome="Liver cancer",population="East Asian",pmid=22174901,study="HKHC",ncase=95,ncontrol=97,UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele="A2",lnor="ln.OR.",se="SE",,eaf="Freq_A1",effect_allele_confirmed=TRUE,ID=68,all_summary_stats=TRUE,efo="hepatocellular carcinoma")
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/livercancer/hpc68.txt",sep="\t",col.names=T,row.names=F,quote=F)


############
# GAMEONE####
############
# no MAF/EAF reported
snplist<-make_snplist(efo= c("breast carcinoma","lung carcinoma","colorectal cancer","ovarian carcinoma","prostate carcinoma"),population="European")
gam<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/GAMEON/CrossCancer_MetaAnalysis_Results_2016_Release_20190222.txt",exact_match=TRUE)
dat<-format_data(Dat=gam,outcome="Cancer (5 sites)",population="European",pmid=27197191,study="GAME-ON",ncase=61851,ncontrol=61820,UKbiobank=FALSE,rsid="ID",Effect.Allele="ALT",Other.Allele="REF",or="meta.OR",p="meta.Pvalue",effect_allele_confirmed=TRUE,ID=57,all_summary_stats=TRUE,efo=c("breast carcinoma","lung carcinoma","colorectal cancer","ovarian carcinoma","prostate carcinoma"))
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/GAMEON/can57.txt",sep="\t",col.names=T,row.names=F,quote=F)

##############
# Lung cancer#
##############
snplist<-make_snplist(efo = "lung carcinoma",population="European")
luc<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/TRICL_LungCancer/TRICL_Meta_032114_Overall.csv",exact_match=TRUE,file_sep=",")
dat<-format_data(Dat=luc,outcome="Lung cancer",population="European",pmid=27488534,ncase="N_Cases",ncontrol="N_Controls",study="ILCCO",UKbiobank=FALSE,rsid="rs_number",Effect.Allele="effect_allele",Other.Allele="reference_allele",or="OR_fixed",se="StdError_fixed",eaf="EAF",p="Pvalue_fixed",I2="I2",phet="phete.Q._Pvalue",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=73,open_gwas=FALSE,efo = "lung carcinoma") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/TRICL_LungCancer/luc73.txt",sep="\t",col.names=T,row.names=F,quote=F)

