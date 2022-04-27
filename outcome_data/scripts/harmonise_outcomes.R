rm(list=ls())
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/colocalisation/scripts/extract_SNPs_functions.R")
# library(ggplot2)
# library(cowplot)
# library(plyr)
# library(devtools)
# remotes::install_github("ramiromagno/gwasrapidd")
# library(gwasrapidd)
# install_github("MRCIEU/TwoSampleMR")
# library(TwoSampleMR)

# devtools::install_github('MarkEdmondson1234/googleAuthR@v0.8.1')
# setwd("~/MR_FattyAcids/data/summary_data/")

#####################################
#Squamous cell carcinoma skin cancer#
#####################################

#The L95 and U95 columns don't seem to match the provided StdErr and P value columns. 
#Use the Effect and StdErr columns 
scc<-read.table("~/MR_FattyAcids/data/summary_data/SCC/Fatty_candidate_SNPs_Harvard_SCC_data.txt",sep="\t",head=T,stringsAsFactors=F )
dat<-format_data(Dat=scc,outcome="Squamous cell skin cancer",population="European",pmid=23548203,study="HNMSC",UKbiobank=FALSE,rsid="SNP",Effect.Allele="Effect.Allele",Other.Allele= "Non.Effect.Allele",lnor="Effect",se="StdErr",eaf="Effect.Allele.Freq",p="P.value",effect_allele_confirmed=TRUE,phet="HetPVal",I2="HetISq",ncase=825,ncontrol=11518,ID=71,Direction="Direction") 
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/skinscc71.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################################
##################################
#Basal cell carcinoma skin cancer#
##################################

bcc<-read.table("~/MR_FattyAcids/data/summary_data/BCC/fattyacidSNP_BCC_HarvardGWAS_4242cases_Meta_Final_Results.txt",
	sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=bcc,outcome="Basal cell skin cancer",population="European",pmid=23548203,study="HNMSC",UKbiobank=FALSE,rsid="SNP",Effect.Allele="Allele1",Other.Allele= "Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="Pvalue",effect_allele_confirmed=TRUE,phet="HetPVal",I2="HetISq",ncase=4242 ,ncontrol=12802,ID=70,Direction="Direction") 
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/skinbcc70.txt",sep="\t",col.names=T,row.names=F,quote=F)

############################## 
#head and neck cancer overall#
##############################

hnc<-read.table("~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/HeadNeck_14042020/FA_HeadNeck_HNC_2020-04-14.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=hnc,outcome="Oral cavity and pharyngeal cancer",population="European",pmid=27749845,study="INHANCE",UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele= "A2",lnor="Effect",se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info="Quality",ncase=6034 ,ncontrol=6585,ID=81) 
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/headneck81.txt",sep="\t",col.names=T,row.names=F,quote=F)

############################## 
#head and neck cancer Oral ###
##############################

hnc<-read.table("~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/HeadNeck_14042020/FA_HeadNeck_OC_2020-04-14.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=hnc,outcome="Oral cancer",population="European",pmid=27749845,study="INHANCE",UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele= "A2",lnor="Effect",se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info="Quality",ncase=2990 ,ncontrol=6585,ID=80) 
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/oral80.txt",sep="\t",col.names=T,row.names=F,quote=F)

############################## 
#head and neck cancer oropharyngeal ###
##############################

hnc<-read.table("~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/HeadNeck_14042020/FA_HeadNeck_OPC_2020-04-14.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=hnc,outcome="Oropharyngeal cancer",population="European",pmid=27749845,study="INHANCE",UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele= "A2",lnor="Effect",se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info="Quality",ncase=2641 ,ncontrol=6585,ID=82) 
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/oropharyngeal82.txt",sep="\t",col.names=T,row.names=F,quote=F)

# ###############
#Kidney cancer#
###############
kid<-read.table("~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/Kidney_14042020/FA_Kidney_IARC_NCIs_UK_MDA_2020-04-14.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=kid,outcome="Kidney cancer" ,population="European",pmid=28598434,study="KidRISK",UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele= "A2",lnor="Effect",se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info="Quality",ncase=10784 ,ncontrol=20406,ID=89) 
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/kidney89.txt",sep="\t",col.names=T,row.names=F,quote=F)

# #################
#Pancreatic cancer#
###################

pan<-read.table("~/MR_FattyAcids/data/summary_data/IARC/Summary_stats/Pancreas_09042020/FA_Pancreas_2020-04-09.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=pan,outcome="Pancreatic cancer" ,population="European",pmid=26098869,study="PanScan I+II+PanC4",UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele= "A2",lnor="Effect",se="StdErr",eaf="EAF",p="P.value",effect_allele_confirmed=TRUE,info="Quality",ncase=7110 ,ncontrol=7264,ID=124) 

write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pan124.txt",sep="\t",col.names=T,row.names=F,quote=F)



################################
#acute lymphoblastic leukemia #
################################

All<-read.table("~/MR_FattyAcids/data/summary_data/ALL_19684603/allrisk.txt",sep="\t",head=T,stringsAsFactors=F)
All.freq<-read.table("~/MR_FattyAcids/data/summary_data/ALL_19684603/fatty_acid_all_freq.txt",sep="\t",head=T,stringsAsFactors=F)
All<-merge(All,All.freq[,c("rs","case.b.freq","control.b.freq")],by.x="rs",by.y="rs")
dat<-format_data(Dat=All,outcome="Acute lymphoblastic leukemia",population="European",pmid=19684603,study="SJ-COG",ncase=317,ncontrol=17958,UKbiobank=FALSE,rsid="rs",Effect.Allele="B",Other.Allele="A",or="susceptibility_OR..B.",eaf="control.b.freq",p="p.value",effect_allele_confirmed=TRUE,ID=130)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/all130.txt",sep="\t",col.names=T,row.names=F,quote=F)

# ########################
# Meningioma_29762745######
##############################
Men<-read.table("~/MR_FattyAcids/data/summary_data/Meningioma_29762745/meningioma_fa_snps_data.tsv",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Men,outcome="Meningioma",population="European",pmid=29762745,study="MENC",ncase=2138,ncontrol=12081,UKbiobank=FALSE,rsid="snp",Effect.Allele="effect_allele",Other.Allele="other_allele",lnor="beta",se="se",eaf="effect_allele_frequency",p="pvalue",info=NA,info1="german_cohort_info",info2="usa_cohort_info",info3=NA,HWEp=NA,phet="p_heterogeneity",I2="I2",Q=NA,Direction=NA,effect_allele_confirmed=TRUE,or=NA,lci=NA,uci=NA,ID=94)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/men94.txt",sep="\t",col.names=T,row.names=F,quote=F)

#########################
# East Asian studies#####
#########################

########################################
# nasopharyngeal_carcinoma_19664746####
########################################
# 43 SNPs duplicated with incompatible effect information between the duplicates. Do not consider this study reliable. However the corresponding author wrote that something went wrong with the imputed data and that the GWAS data was correct. So presumably the genotyped SNPs are correct. 
 
Ncc.g<-read.table("~/MR_FattyAcids/data/summary_data/nasopharyngeal_carcinoma_19664746/ncc_gwas.txt",sep="\t",head=T,stringsAsFactors=F)
Ncc.g$File<-"gwas"
Ncc.i<-read.table("~/MR_FattyAcids/data/summary_data/nasopharyngeal_carcinoma_19664746/ncc_impute.txt",sep="\t",head=T,stringsAsFactors=F)
Ncc.i<-Ncc.i[,1:8]
Ncc.i$File <- "impute"
Ncc<-plyr::rbind.fill(Ncc.g,Ncc.i)
# Ncc<-Ncc[Ncc$File != "impute",]

dat<-format_data(Dat=Ncc,outcome="Nasopharyngeal carcinoma",population="East Asian",pmid=19664746,study="TNC",ncase=277,ncontrol=285,UKbiobank=FALSE,rsid="SNP",Effect.Allele="Effect.allele",Other.Allele="Non.effect.allele",lnor="Natural.log.odds.ratio",se="Standard.error.of.the.log.odds.ratio",info="info",HWEp="P.values.for.Hardy.Weinberg.equilibrium",effect_allele_confirmed=TRUE,ID=132,eaf="Effect.allele.frequency")

# eaf="Effect.allele.frequency". failed QC pipeline. set eaf to NA
# but need minor allele frequency to caclulate minor allele count. So keep and drop eaf later. 

write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/ncc132.txt",sep="\t",col.names=T,row.names=F,quote=F)

##########################################
# Hepatocellular carcinoma_22807686######
##########################################

# excluded because required information not provided
# hpc<-read.table("Hepatocellular carcinoma_22807686/fatty_acid_instruments_SNPtable_20190920/fatty_acid_instruments_SNPtable_results.assoc.logistic.txt",sep="\t",head=T,stringsAsFactors=F,fill=T)

hpc<-preformat_hpc_23()
head(hpc)
dat<-format_data(Dat=hpc,outcome="Hepatocellular carcinoma in chronic hepatitis B virus carriers",population="East Asian",pmid=22807686,study="CHC",ncase="ncase",ncontrol="ncontrol",UKbiobank=FALSE,rsid="SNP",Effect.Allele="minor_Allele.x",Other.Allele="major_Allele",se="SE",eaf="minor_Allele_frequency(controls)",p="P",info="info",effect_allele_confirmed=TRUE,or="OR",lci="L95",uci="U95",ID=23,test_statistic="STAT")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/hpc23.txt",sep="\t",col.names=T,row.names=F,quote=F)

# readLines("Hepatocellular carcinoma_22807686/Reply/LD_SNPtable_20190920/LD_SNPtable_results.assoc.logistic.logistic")


##########################################
# chronic_myeloid_leukemia_21540461######
##########################################

Cml<-read.table("~/MR_FattyAcids/data/summary_data/chronic_myeloid_leukemia_21540461/snplist_Kim_21540461.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Cml,outcome="Chronic myeloid leukemia",population="East Asian",pmid=21540461,study="KCML",ncase=201,ncontrol=497,UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele="A2",lnor=NA,se="SE",eaf="F_U",p="P",info=NA,info1=NA,info2=NA,info3=NA,HWEp=NA,phet=NA,I2=NA,Q=NA,Direction=NA,effect_allele_confirmed=FALSE,or="OR",lci= "L95",uci="U95",ID=87,test_statistic=NA)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/cmyeloidleukemia87.txt",sep="\t",col.names=T,row.names=F,quote=F)


################################
#glioma_UCSF_Mayo_28346443####
##############################
# genome build is in hg19 confirmed by correspondence
setwd("~/MR_FattyAcids/data/summary_data/glioma_UCSF_Mayo_28346443")
Files<-dir()
Files<-Files[grep(".pdf",Files,invert=T)]
gli2<-lapply(1:length(Files), FUN=function(i)
    read.table(Files[i],sep=" ",head=T,stringsAsFactors=F))
gli2<-do.call(rbind,gli2)

dat<-format_data(Dat=gli2,outcome="Glioma",population="European",pmid=28346443,study="UCSF_MAYO",ncase=1519,ncontrol=804,UKbiobank=FALSE,rsid="SNP",Effect.Allele="A2",Other.Allele="A1",se="SE",eaf="FRQ",p="P",info="INFO",effect_allele_confirmed=TRUE,or="OR",ref="~/fatty-acids/colocalisation/data/UKBB_10K_bed_hg19.txt",ID=134,all_summary_stats=FALSE,summary_set="FAregions")

write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/glioma134.txt",sep="\t",col.names=T,row.names=F,quote=F)

################################
#glioma_SFAGS_28346443####
##############################

gli<-read.table("~/MR_FattyAcids/data/summary_data/glioma_SFAGS_28346443/sfags_illumina_haycock_results v2.txt",sep="\t",head=T,stringsAsFactors=F)
# gli<-read.table("~/MR_FattyAcids/data/summary_data/glioma_SFAGS_28346443/sfags_illumina_haycock_results.txt",sep="\t",head=T,stringsAsFactors=F)

# gli<-gli[gli$rsid %in% c("rs114935910","rs7068072","rs16993510", "rs3848771"),]
# gli$z<-gli$beta/gli$se
# gli$z.p<-qnorm(round(gli$P_value,4)/2,lower.tail=F)
# gli[,c("rsid","P_value","beta","se","z","z.p")]

dat<-format_data(Dat=gli,outcome="Glioma",population="European",pmid=28346443,study="UCSF_AGS + SFAGS",ncase=677,ncontrol=3940,UKbiobank=FALSE,rsid="rsid",Effect.Allele="allele_B",Other.Allele="allele_A",lnor="beta",se="se",p="P_value",info="info",effect_allele_confirmed=TRUE,ID=133,all_summary_stats=FALSE,summary_set="FAregions",HWEp="controls_hwe",eaf="controls_maf") #eaf was coded as "controls_maf" by mistake in first harmonisation 

# write.table(dat,"~/fatty-acids/outcome_data/data/harmonised_preqc/glioma133.txt",sep="\t",col.names=T,row.names=F,quote=F)

# drop "rs114935910","rs7068072","rs16993510","rs3848771", which had very unsually large betas and Z scores given their p values
# dat<-dat[!dat$rsid %in% c("rs114935910","rs7068072","rs16993510","rs3848771"),]

# dat[abs(dat$lnor) > 2,]
# dat<-dat[dat$info>=0.8,]
# dat$z<-abs(dat$lnor/dat$se)
# dat$zp<-qnorm(dat$p/2,lower.tail=F)
# dat<-dat[dat$z<5,]
# plot(dat$zp,dat$z)
# dat[dat$z>5,c("rsid","p","lnor","se","")]

# we make the datasets in both harmonised and harmonise_preqc the same because we set eaf NA after first collating the summary statistics. The harmonised directory is basically same as harmonised_preqc, except that for a small number of studies with the effect allele column switched to the non-effect allele e.g. ID 967
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/glioma133.txt",sep="\t",col.names=T,row.names=F,quote=F)

write.table(dat,"~/fatty-acids/outcome_data/data/harmonised_preqc/glioma133.txt",sep="\t",col.names=T,row.names=F,quote=F)


#############################
# Hodgkin lymphoma 29196614####
##############################

Hdl<-read.table("~/MR_FattyAcids/data/summary_data/Hodgkin_lymphoma_29196614/HL_UKGERNSHLGMETA_FA_SNPS_INFO.txt",sep="\t",head=T,stringsAsFactors=F)
Hdl1<-Hdl[grep(";",Hdl$rsid),]
Hdl2<-Hdl[grep(";",Hdl$rsid,invert=T),]
for(i in 1:nrow(Hdl1)){
	Hdl1$rsid[i]<-unlist(strsplit(Hdl1$rsid[i],split=";"))[1]
}
Hdl<-rbind(Hdl1,Hdl2)

dat<-format_data(Dat=Hdl,outcome="Hodgkin lymphoma",population="European",pmid=29196614,study="HLS",ncase=3077,ncontrol=13680,UKbiobank=FALSE,rsid="rsid",Effect.Allele="allele_B",Other.Allele="allele_A",lnor="beta",se="se",eaf="coded_af",p="P_value",info1="cohort_1_info",info2= "cohort_2_info",info3="cohort_3_info",effect_allele_confirmed=TRUE,ID=69,phet="P_heterogeneity",I2="I2",Q="Q")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/hodgkin69.txt",sep="\t",col.names=T,row.names=F,quote=F)

# ############
# Colorectal####
################

Crc_all<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer/crc_all.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Crc_all,outcome="Colorectal cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=58131,ncontrol=67347,UKbiobank=TRUE,rsid="rs",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=60,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",test_statistic="LRT.Z")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/colorectal60.txt",sep="\t",col.names=T,row.names=F,quote=F)

# CRC males

Crc_males<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer/crc_males.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Crc_males,outcome="Colorectal cancer in males",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=31288,ncontrol=34527,UKbiobank=TRUE,rsid="rs",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=62,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",test_statistic="LRT.Z")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/colorectal_males62.txt",sep="\t",col.names=T,row.names=F,quote=F)

# CRC females

Crc_females<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer/crc_females.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Crc_females,outcome="Colorectal cancer in females",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=26843,ncontrol=32820,UKbiobank=TRUE,rsid="rs",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=61,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",test_statistic="LRT.Z")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/colorectal_females61.txt",sep="\t",col.names=T,row.names=F,quote=F)

# CRC rectal

Crc_rectal<-read.table("colorectal_cancer/crc_rectal.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Crc_rectal,outcome="Rectal cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=15775,ncontrol=67694,UKbiobank=TRUE,rsid="rs",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=65,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",test_statistic="LRT.Z")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/rectal65.txt",sep="\t",col.names=T,row.names=F,quote=F)

# CRC colon

Crc_colon<-read.table("colorectal_cancer/crc_colon.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Crc_colon,outcome="Colon cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=31083,ncontrol=67694,UKbiobank=TRUE,rsid="rs",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=59,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",test_statistic="LRT.Z")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/colon59.txt",sep="\t",col.names=T,row.names=F,quote=F)


# CRC proximal

Crc_proximal<-read.table("colorectal_cancer/crc_proximal.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Crc_proximal,outcome="Proximal colorectal cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=13857,ncontrol=67694,UKbiobank=TRUE,rsid="rs",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=64,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",test_statistic="LRT.Z")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/colorectal_proximal64.txt",sep="\t",col.names=T,row.names=F,quote=F)

# CRC distal

Crc_distal<-read.table("colorectal_cancer/crc_distal.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=Crc_distal,outcome="Distal colorectal cancer",population="European",pmid=30510241,study="GECCO/CORECT/CCFR",ncase=15306,ncontrol=67694,UKbiobank=TRUE,rsid="rs",Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="LRT.pval",effect_allele_confirmed=TRUE,ID=63,phet="HetPVal",I2="HetISq",Q="HetChiSq",Direction="Direction",test_statistic="LRT.Z")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/colorectal_distal63.txt",sep="\t",col.names=T,row.names=F,quote=F)

# ####################
# bcell_ALL_29632299# 
# ####################

bcell_all<-preformat_all_id4()
dat<-format_data(Dat=bcell_all,outcome="Acute lymphoblastic leukaemia",population="European",pmid=29632299,study="BC-ALL",ncase=2442,ncontrol=14609,UKbiobank=FALSE,rsid="SNP",Effect.Allele="Risk_allele",Other.Allele="Non_risk_allele",lnor="beta",se="se",eaf="eaf",p="P_value",effect_allele_confirmed=TRUE,ID=4,phet="P_heterogeneity",I2="I2",Q="Q")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/bcell_all4.txt",sep="\t",col.names=T,row.names=F,quote=F)

# ####################
#multiple_myeloma_26007630
# ####################
mma<-preformat_mma_id96()
dat<-format_data(Dat=mma,outcome="Multiple myeloma",population="European",pmid=26007630,study="MMS",ncase=1714,ncontrol=10391,UKbiobank=FALSE,rsid="SNP",Effect.Allele="Effect.Allele",Other.Allele="Other.Allele",lnor="beta",eaf="frq",p="p",effect_allele_confirmed=TRUE,ID=96,info="info")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/multiplemyeloma96.txt",sep="\t",col.names=T,row.names=F,quote=F)


#######################
# Thyroid cancer####
####################
# 8. age, sex, study and PC
# 9. inflation factor = 1.02
# 10. We excluded SNPs with a call rate <95% by study, duplicate SNPs, SNPs for which the cluster plots were judged to be not ideal,  SNPs not in HWE in each ethnic group (p<10-7 in controls, p<10-12 in cases), SNP with MAF <1%, imputed SNPs with R2<0.3

# thy<-readLines("thyroid_cancer/Fatty_acid_SNPlist_results.txt")
thy<-read.table("~/MR_FattyAcids/data/summary_data/thyroid_cancer/thy.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=thy,outcome="Thyroid cancer",population="European",pmid="unpublished",study="EPITHYR",ncase=1554,ncontrol=1973,UKbiobank=FALSE,rsid="snp.id",Effect.Allele="A2",Other.Allele="A1",lnor="beta.all",se="se.all",eaf="freq2",p="p.value.all",effect_allele_confirmed=TRUE,ID=26,info="info")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/thyroid26.txt",sep="\t",col.names=T,row.names=F,quote=F)


###################
# Cervical cancer#
###################

# Received update on 25 September 2019 with effect allele frequency information that need to incorporate
cer<-preformat_cer_id129()
dat<-format_data(Dat=cer,outcome="Cervical cancer",population="European",pmid=23482656,study="SCCS",ncase=1034,ncontrol=3948,UKbiobank=FALSE,rsid="SNP",Effect.Allele="A1",Other.Allele="A2",or="OR",se="SE",eaf="MAF",effect_allele_confirmed=TRUE,ID=129)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/cervical129.txt",sep="\t",col.names=T,row.names=F,quote=F)


##################
# Neuroblastoma######
##################

neu<-preformat_neu_id107()
dat<-format_data(Dat=neu,outcome="Neuroblastoma",population="European",pmid=28545128,study="NBS",ncase=2101,ncontrol=4202,UKbiobank=FALSE,rsid="rsid",Effect.Allele="allele_B",Other.Allele="allele_A",lnor="Phenotype_frequentist_add_C1_C2_C3_C4_C5_C6_C7_C8_C9_C10_C11_C12_C13_C14_C15_C16_C17_C18_C19_C20_score_beta_1",se="Phenotype_frequentist_add_C1_C2_C3_C4_C5_C6_C7_C8_C9_C10_C11_C12_C13_C14_C15_C16_C17_C18_C19_C20_score_se_1",eaf="eaf",effect_allele_confirmed=TRUE,ID=107,info="Phenotype_frequentist_add_C1_C2_C3_C4_C5_C6_C7_C8_C9_C10_C11_C12_C13_C14_C15_C16_C17_C18_C19_C20_score_info")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/neuroblastoma107.txt",sep="\t",col.names=T,row.names=F,quote=F)



##############################
# colorectal cancer in ACCC###
#################################


snplist<-make_snplist(trait="colorectal cancer",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
crc<-preformat_accc_3()
dat<-format_data(Dat=crc,outcome="Colorectal cancer",population="East Asian",pmid=31826910,study="ACCC",ncase=23572,ncontrol=48700,UKbiobank=FALSE,rsid="rsid",Effect.Allele="Allele1",Other.Allele="Allele2",eaf="Freq1",lnor="Effect",se="StdErr",p="P.value",effect_allele_confirmed=TRUE,ID=3,all_summary_stats=FALSE,summary_set="FAregions",Direction="Direction",phet="HetPVal",I2="HetISq",Q="HetChiSq")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/crc3.txt",sep="\t",col.names=T,row.names=F,quote=F)



########################################
# nasopharyngeal carcinoma	27236004 #####
#########################################
# previously pmid 19478819
# snplist<-make_snplist(trait="nasopharyngeal carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
npc<-read.table("~/MR_FattyAcids/data/summary_data/Nasopharyngealcancer_19478819_27236004/mendelianrandomizationproject/NPC Malaysia GWAS fatty acids SNP  shortlist.txt",sep="\t",head=T,stringsAsFactors=F)
head(npc)
dat<-format_data(Dat=npc,outcome="Nasopharyngeal carcinoma",population="East Asian",pmid=27236004,study="MNC",ncase=271,ncontrol=456,UKbiobank=FALSE,rsid="rsid",Effect.Allele="Effect.allele.alleleB",Other.Allele="Non.effect.allele.alleleA",lnor="frequentist_add_beta_1.genotype",se="frequentist_add_se_1.genotype",eaf="control_effect_allele_freq",p="frequentist_add_pvalue",effect_allele_confirmed=TRUE,ID=97,info="all_impute_info",HWEp="controls_hwe",Q="Cochrane.Q.test")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/npc97.txt",sep="\t",col.name=T,row.name=F,quote=F)



#########################################
#Malignant pleural mesothelioma 23626673#
#########################################
# maf/eaf not reported

# Dat1<-read.table("~/MR_FattyAcids/data/summary_data/mesothelioma_23626673/Results_GWAS_Cugliari_SNP_52_proxy.csv",sep=";",head=T,stringsAsFactors=F)

mes<-preformat_mes_98()
dat<-format_data(Dat=mes,ncase=407,ncontrol=389,pmid=23626673,outcome="Pleural mesothelioma",population="European",study="MPM",UKbiobank=FALSE,rsid="SNP",Effect.Allele="allele1",Other.Allele="allele2",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,ID=98,info="info")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/mes98.txt",sep="\t",col.name=T,row.name=F,quote=F)



##############################################
# Noncardia gastric adenocarcinoma	26701879##
##############################################
nga<-read.table("~/MR_FattyAcids/data/summary_data/gastric_cancer_26701879/gc_fatty_acid_MR_gc_new_1kgv3_snptest-20181129.txt",sep="\t",head=T,stringsAsFactors=F)

dat<-format_data(Dat=nga,outcome="Noncardia gastric adenocarcinoma",population="East Asian",pmid=26701879,study="NB-UGC",ncase=1006,ncontrol=2273,UKbiobank=FALSE,rsid="rsid",Effect.Allele="eff",Other.Allele="ref",lnor="beta" ,se="se",eaf="eff.allele.freq.in.controls",info="info",HWEp="HWE_Controls",effect_allele_confirmed=TRUE,ID=104) #effect allele confirmed from column names

write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/gastric104.txt",sep="\t",col.name=T,row.name=F,quote=F)



##########################################
# Pancreatic cancer	PanScan I+II 29422604#
##########################################

pan<-read.table("~/MR_FattyAcids/data/summary_data/panscan/PanScanI+II.txt",sep="\t",head=T,stringsAsFactors=F)
dat<-format_data(Dat=pan,outcome="Pancreatic cancer",population="European",pmid=29422604,study="PanScan I+II",ncase="Num_Case",ncontrol="Num_Control",UKbiobank=FALSE,rsid="Locus",Effect.Allele="Effect_Allele",Other.Allele="Reference_Allele",lnor="Beta",se="Standard_error_of_beta",eaf="Effect_Allele_Freq_Control",p="P_value",info="Info",effect_allele_confirmed=TRUE,ID=123) #effect allele confirmed from column names
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pan123.txt",sep="\t",col.name=T,row.name=F,quote=F)

##########################################
# Pancreatic cancer	PanScan III 29422604#
##########################################

pan<-read.table("~/MR_FattyAcids/data/summary_data/panscan/PanScanIII.txt",sep="\t",head=T,stringsAsFactors=F)

dat<-format_data(Dat=pan,outcome="Pancreatic cancer",population="European",pmid=29422604,study="PanScan III",ncase="Num_Case",ncontrol="Num_Control",UKbiobank=FALSE,rsid="Locus",Effect.Allele="Effect_Allele",Other.Allele="Reference_Allele",lnor="Beta",se="Standard_error_of_beta",eaf="Effect_Allele_Freq_Control",p="P_value",info="Info",effect_allele_confirmed=TRUE,ID=125) #effect allele confirmed from column names

write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pan125.txt",sep="\t",col.name=T,row.name=F,quote=F)


###########################
# Pancreatic cancer panc4###
##############################

pan<-preformat_pan_121()
dat<-format_data(Dat=pan,outcome="Pancreatic cancer",population="European",pmid=29422604,study="PanC4",ncase="cases_total",ncontrol="controls_total", UKbiobank=FALSE,rsid="V4",Effect.Allele="effect_al",Other.Allele="ref_al",lnor="frequentist_add_beta_1",se="frequentist_add_se_1",eaf="controls_EAF",p="frequentist_add_pvalue",effect_allele_confirmed=TRUE,info="info",ID=121) 
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pan121.txt",sep="\t",col.names=T,row.names=F,quote=F)


############################################################
# Full summary association statistics available in OpenGWAS#
############################################################

Dat<-find_cancerstudies_opengwas()

################
# Breast cancer#
################

snplist<-make_snplist(efo="breast carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
brc <- ieugwasr::associations(id="ieu-a-1126", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(brc,stringsAsFactors=F),outcome="Breast cancer",population="European",pmid=29059683,ncase=122977,ncontrol=105974,study="BCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=6,open_gwas=TRUE,efo="breast carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/brc6.txt",sep="\t",col.names=T,row.names=F,quote=F)

###################
#ER- Breast cancer#
###################

snplist<-make_snplist(efo="estrogen-receptor negative breast cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
brc <- ieugwasr::associations(id="ieu-a-1128", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(brc,stringsAsFactors=F),outcome="ER- breast cancer",population="European",pmid=29059683,ncase=21468,ncontrol=105974,study="BCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=7,open_gwas=TRUE,efo="estrogen-receptor negative breast cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/erneg7.txt",sep="\t",col.names=T,row.names=F,quote=F)


###################
#ER+ Breast cancer#
###################

snplist<-make_snplist(efo="estrogen-receptor positive breast cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
brc$trait <- ieugwasr::associations(id="ieu-a-1127", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(brc,stringsAsFactors=F),outcome="ER+ breast cancer",population="European",pmid=29059683,ncase=69501,ncontrol=105974,study="BCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=8,open_gwas=TRUE,efo="estrogen-receptor positive breast cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/erpos8.txt",sep="\t",col.names=T,row.names=F,quote=F)

###################
#Ovarian cancer#
###################

snplist<-make_snplist(efo="ovarian carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1120", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Ovarian cancer",population="European",pmid=28346442,ncase=25509,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=118,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/ovc118.txt",sep="\t",col.names=T,row.names=F,quote=F)

###################
#Serous ovarian cancer#
###################

snplist<-make_snplist(efo="ovarian serous carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1228", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Serous ovarian cancer",population="European",pmid=28346442,ncase=14049,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=119,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/soc119.txt",sep="\t",col.names=T,row.names=F,quote=F)

##################################
#High grade serous ovarian cancer#
##################################

snplist<-make_snplist(efo="ovarian serous carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1121", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="High grade serous ovarian cancer",population="European",pmid=28346442,ncase=13037,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=110,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/hso110.txt",sep="\t",col.names=T,row.names=F,quote=F)

##################################
#Low malignant potential ovarian cancer#
##################################

snplist<-make_snplist(efo="ovarian carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1233", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Low malignant potential ovarian cancer",population="European",pmid=28346442,ncase=3103,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=115,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lmp115.txt",sep="\t",col.names=T,row.names=F,quote=F)

####################################################################
#Low grade & low malignant potential serous ovarian cancer#
####################################################################

snplist<-make_snplist(efo="ovarian serous carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1229", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Low grade & low malignant potential serous ovarian cancer",population="European",pmid=28346442,ncase=2966,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=112,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/llo112.txt",sep="\t",col.names=T,row.names=F,quote=F)

##############################
# Endometrioid ovarian cancer#
##############################

snplist<-make_snplist(efo="ovarian endometrioid carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1125", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Endometrioid ovarian cancer",population="European",pmid=28346442,ncase=2810,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=109,open_gwas=TRUE,efo="ovarian endometrioid carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/eoc109.txt",sep="\t",col.names=T,row.names=F,quote=F)


#########################
#Mucinous ovarian cancer#
#########################
snplist<-make_snplist(efo="ovarian carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1231", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Mucinous ovarian cancer",population="European",pmid=28346442,ncase=2566,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=117,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/moc117.txt",sep="\t",col.names=T,row.names=F,quote=F)

###############################################
#Low malignant potential serous ovarian cancer#
###############################################

snplist<-make_snplist(efo="ovarian serous carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1230", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Low malignant potential serous ovarian cancer",population="European",pmid=28346442,ncase=1954,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=116,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lso116.txt",sep="\t",col.names=T,row.names=F,quote=F)

##################################
#Invasive mucinous ovarian cancer#
##################################

snplist<-make_snplist(efo="ovarian carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1123", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Invasive mucinous ovarian cancer",population="European",pmid=28346442,ncase=1417,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=111,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/imo111.txt",sep="\t",col.names=T,row.names=F,quote=F)


###########################
#Clear cell ovarian cancer#
###########################

snplist<-make_snplist(efo="ovarian carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1124", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Clear cell ovarian cancer",population="European",pmid=28346442,ncase=1366,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=108,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/cco108.txt",sep="\t",col.names=T,row.names=F,quote=F)

#################################################
#Low malignant potential mucinous ovarian cancer#
################################################

snplist<-make_snplist(efo="ovarian carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1232", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Low malignant potential mucinous ovarian cancer",population="European",pmid=28346442,ncase=1149,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=114,open_gwas=TRUE,efo="ovarian carcinoma")
# rs4927355
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lmm114.txt",sep="\t",col.names=T,row.names=F,quote=F)

#################################
#Low grade serous ovarian cancer#
#################################

snplist<-make_snplist(efo="ovarian serous carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ovc <- ieugwasr::associations(id="ieu-a-1122", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ovc,stringsAsFactors=F),outcome="Low grade serous ovarian cancer",population="European",pmid=28346442,ncase=1012,ncontrol=40941,study="OCAC",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=113,open_gwas=TRUE,efo="ovarian serous carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lso113.txt",sep="\t",col.names=T,row.names=F,quote=F)

###########################
#Advanced prostate cancer #
###########################

snplist<-make_snplist(efo = "prostate carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
pro <- ieugwasr::associations(id="ieu-a-1238", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(pro,stringsAsFactors=F),outcome="Advanced prostate cancer",population="European",pmid=29892016,ncase=15167,ncontrol=58308,study="PRACTICAL",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=126,open_gwas=TRUE,efo = "prostate carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/apc126.txt",sep="\t",col.names=T,row.names=F,quote=F)


#############################
#Early-onset prostate cancer#
#############################

snplist<-make_snplist(efo = "prostate carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
pro <- ieugwasr::associations(id="ieu-a-1240", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(pro,stringsAsFactors=F),outcome="Early-onset prostate cancer",population="European",pmid=29892016,ncase=6988,ncontrol=44256,study="PRACTICAL",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=127,open_gwas=TRUE,efo = "prostate carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/eop127.txt",sep="\t",col.names=T,row.names=F,quote=F)



#############
#Lung cancer#
#############

snplist<-make_snplist(efo = "lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
luc <- ieugwasr::associations(id="ieu-a-966", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=24880342,ncase=11348,ncontrol=15861,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=74,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/luc74.txt",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer#
#############

snplist<-make_snplist(efo = "lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
luc <- ieugwasr::associations(id="ieu-a-987", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=28604730,ncase=29266,ncontrol=56450,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=75,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/luc75.txt",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer in ever smokers#
#############
snplist<-make_snplist(efo = "lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
luc <- ieugwasr::associations(id="ieu-a-985", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer in ever smokers",population="European",pmid=28604730,ncase=23848,ncontrol=16605,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=76,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lce76.txt",sep="\t",col.names=T,row.names=F,quote=F)
 
#############
#Lung cancer in never smokers#
#############
snplist<-make_snplist(efo = "lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
luc <- ieugwasr::associations(id="ieu-a-986", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(luc,stringsAsFactors=F),outcome="Lung cancer in never smokers",population="European",pmid=28604730,ncase=2303,ncontrol=6995,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=77,open_gwas=TRUE,efo = "lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lcn77.txt",sep="\t",col.names=T,row.names=F,quote=F)


#######################
#Lung adenocarcinoma#
#######################

snplist<-make_snplist(efo = "lung adenocarcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
luc <- ieugwasr::associations(id="ieu-a-984", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(luc,stringsAsFactors=F),outcome="Lung adenocarcinoma",population="European",pmid=28604730,ncase=11245,ncontrol=54619,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=72,open_gwas=TRUE,efo = "lung adenocarcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lad72.txt",sep="\t",col.names=T,row.names=F,quote=F)

#######################
#Squamous cell lung cancer #
#######################

snplist<-make_snplist(efo = "squamous cell lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
luc <- ieugwasr::associations(id="ieu-a-989", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(luc,stringsAsFactors=F),outcome="Squamous cell lung cancer",population="European",pmid=28604730,ncase=7704,ncontrol=54763,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=79,open_gwas=TRUE,efo = "squamous cell lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/scl79.txt",sep="\t",col.names=T,row.names=F,quote=F)

########################
#Small cell lung cancer#
########################

snplist<-make_snplist(efo = "small cell lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
luc <- ieugwasr::associations(id="ieu-a-988", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(luc,stringsAsFactors=F),outcome="Small cell lung cancer",population="European",pmid=28604730,ncase=2791,ncontrol=20580,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=78,open_gwas=TRUE,efo = "small cell lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/scl78.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################
# Japanese Biobank####
########################
# Cervical cancer
snplist<-make_snplist(efo="cervical carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-98", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Cervical cancer",population="East Asian",pmid=32514122,study="BJ",ncase=605,ncontrol=89731,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=11,open_gwas=TRUE,efo="cervical carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/cec11.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Biliary tract cancer
snplist<-make_snplist(efo="cholangiocarcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-92", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Biliary tract cancer",population="East Asian",pmid=32514122,study="BJ",ncase=339,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=9,open_gwas=TRUE,efo="cholangiocarcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/btc9.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Colorectal cancer
snplist<-make_snplist(efo="colorectal cancer",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-107", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Colorectal cancer",population="East Asian",pmid=32514122,study="BJ",ncase=7062,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=12,open_gwas=TRUE,efo="colorectal cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/crc12.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Endometrial cancer
snplist<-make_snplist(efo="endometrial carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-113", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Endometrial cancer",population="East Asian",pmid=32514122,study="BJ",ncase=999,ncontrol=89731,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=13,open_gwas=TRUE,efo="endometrial carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/enc13.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Esophageal squamous cell carcinoma
snplist<-make_snplist(efo="esophageal squamous cell carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-117", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=1300,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=14,open_gwas=TRUE,efo="esophageal squamous cell carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/esc14.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Gastric adenocarcinoma
snplist<-make_snplist(efo="gastric carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-119", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Gastric adenocarcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=6563,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=15,open_gwas=TRUE,efo="gastric carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/gac15.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Blood cancer
snplist<-make_snplist_blood(population="East Asian")
bbj<-ieugwasr::associations(id= "bbj-a-125", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Blood cancer",population="East Asian",pmid=32514122,study="BJ",ncase=1236,ncontrol=211217,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=10,open_gwas=TRUE,efo=c("lymphoma","multiple myeloma","lymphoid leukemia"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/blc10.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Liver cancer
snplist<-make_snplist(efo="hepatocellular carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-158", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Liver cancer",population="East Asian",pmid=32514122,study="BJ",ncase=1866,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=16,open_gwas=TRUE,efo="hepatocellular carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lic16.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Lung cancer
snplist<-make_snplist(efo="lung carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-133", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Lung cancer",population="East Asian",pmid=32514122,study="BJ",ncase=4050,ncontrol=208403,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=17,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/luc17.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Ovarian cancer
snplist<-make_snplist(efo="ovarian carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-139", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Ovarian cancer",population="East Asian",pmid=32514122,study="BJ",ncase=720,ncontrol=89731,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=18,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/ovc18.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Pancreatic cancer
snplist<-make_snplist(efo="pancreatic carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-140", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Pancreatic cancer",population="East Asian",pmid=32514122,study="BJ",ncase=442,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=19,open_gwas=TRUE,efo="pancreatic carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pac19.txt",sep="\t",col.names=T,row.names=F,quote=F)

# Prostate cancer
snplist<-make_snplist(efo="prostate carcinoma",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
bbj<-ieugwasr::associations(id= "bbj-a-148", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Prostate cancer",population="East Asian",pmid=32514122,study="BJ",ncase=5408,ncontrol=103939,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=20,open_gwas=TRUE,efo="prostate carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pro20.txt",sep="\t",col.names=T,row.names=F,quote=F)


#######################
Gallbladder cancer#
#######################
# <100 cases
snplist<-make_snplist(efo = "gallbladder neoplasm",population="East Asian",Dir="~/fatty-acids/outcome_data/data/")
gal <- ieugwasr::associations(id="ieu-a-1057", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(gal,stringsAsFactors=F),outcome="Gallbladder cancer",population="East Asian",pmid=22318345,ncase=41,ncontrol=866,study="GCS",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=58,open_gwas=TRUE,efo="gallbladder neoplasm")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/gbc58.txt",sep="\t",col.names=T,row.names=F,quote=F)

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

snplist<-make_snplist(efo = "glioma",population="European",Dir="~/fatty-acids/outcome_data/data/")
gli <- ieugwasr::associations(id="ieu-a-1013", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(gli,stringsAsFactors=F),outcome="Glioma",population="European",pmid=22886559,ncase=1856	,ncontrol=4955,study="GliomaScan",UKbiobank=FALSE,rsid="rsid",Effect.Allele="nea",Other.Allele="ea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=967,open_gwas=TRUE,efo = "glioma")
# Effect.Allele="ea",Other.Allele="nea" QC plots clearly indicate that these are wrong way around

write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/gli967_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

################
#"Neuroblastoma#
################
# maf/eaf not reported
snplist<-make_snplist(efo = "neuroblastoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
neu <- ieugwasr::associations(id="ieu-a-816", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(neu,stringsAsFactors=F),outcome="Neuroblastoma",population="European",pmid=23222812,ncase=1627,ncontrol=3254,study="NBS",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=106,open_gwas=TRUE,efo = "neuroblastoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/neu106.txt",sep="\t",col.names=T,row.names=F,quote=F)


################
#Thyroid cancer#
################

snplist<-make_snplist(efo = "thyroid carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
thy <- ieugwasr::associations(id="ieu-a-1082", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(thy,stringsAsFactors=F),outcome="Thyroid cancer",population="European",pmid=23894154,ncase=690,ncontrol=497,study="TCS",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=131,open_gwas=TRUE,efo = "thyroid carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/thy131.txt",sep="\t",col.names=T,row.names=F,quote=F)

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

snplist<-make_snplist(efo = "basal cell carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ukb <- ieugwasr::associations(id="ukb-b-8837", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Basal cell carcinoma",population="European",pmid="ukb-b-8837",ncase=4290,ncontrol=458643,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=135,open_gwas=TRUE,efo = "basal cell carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/bcc135_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################
#Bladder cancer#
######################

snplist<-make_snplist(efo = "bladder carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ukb <- ieugwasr::associations(id="ukb-d-C67", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Bladder cancer",population="European",pmid="ukb-d-C67",ncase=1554,ncontrol=359640,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=136,open_gwas=TRUE,efo = "bladder carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/bla136_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################
#Cancer of digestive organs#
######################

snplist<-make_snplist(efo = "digestive system carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ukb <- ieugwasr::associations(id="ukb-d-C3_DIGESTIVE_ORGANS", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Cancer of digestive organs",population="European",pmid="ukb-d-C3_DIGESTIVE_ORGANS",ncase=5690,ncontrol=355504,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=142,open_gwas=TRUE,efo = "digestive system carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised_preqc/cdo142.txt",sep="\t",col.names=T,row.names=F,quote=F)
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/cdo142_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################
#Kidney cancer#
######################

snplist<-make_snplist(efo = "kidney cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
ukb <- ieugwasr::associations(id="ukb-b-1316", variants=snplist,proxies=0)  

dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Kidney cancer",population="European",pmid="ukb-b-1316",ncase=1114,ncontrol=461896,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=145,open_gwas=TRUE,efo = "kidney cancer")
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/kic145_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################
#Lymphoma
######################

snplist<-make_snplist(efo = "lymphoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ukb <- ieugwasr::associations(id="ukb-d-C_LYMPHOMA", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Lymphoma",population="European",pmid="ukb-d-C_LYMPHOMA",ncase=1752,ncontrol=359442,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=151,open_gwas=TRUE,efo = "lymphoma")
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lym151_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################
#Malignant skin cancer
######################
# Representative examples of malignant skin neoplasms include basal cell carcinoma, squamous cell carcinoma, melanoma, and Kaposi sarcoma.

snplist<-make_snplist_malignant_skin_cancer()
ukb <- ieugwasr::associations(id="ukb-d-C3_SKIN", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="ukb-d-C3_SKIN",ncase=16531,ncontrol=344663,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=152,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/msc152_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

######################################
#Respiratory and intrathoracic cancer#
######################################

# C3_NASAL_CAVITY_MIDDLE_EAR, C3_ACCESSORY_SINUS, C3_LARYNX, C3_TRACHEA, C3_BRONCHUS_LUNG, C3_THYMUS, C3_HEART_MEDIASTINUM_PLEURA, C3_RESPIRATORY_INTRATHORACIC3_NAS

snplist<-make_snplist_respiratory()
ukb <- ieugwasr::associations(id="ukb-d-C3_RESPIRATORY_INTRATHORACIC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="ukb-d-C3_RESPIRATORY_INTRATHORACIC",ncase=1944,ncontrol=359250,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=160,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/ric160_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

####################
#Small bowel cancer#
####################

snplist<-make_snplist(efo= "small intestine neuroendocrine tumor",population="European",Dir="~/fatty-acids/outcome_data/data/") #the closest match in the GWAS catalog I could find for small bowel cancer was small intestine neuroendocrine tumor
ukb <- ieugwasr::associations(id="ukb-a-56", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Small bowel cancer",population="European",pmid="ukb-a-56",ncase=156,ncontrol=337003,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=161,open_gwas=TRUE,efo= "small intestine neuroendocrine tumor")
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/sbc161_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

####################
#Squamous cell carcinoma#
####################

snplist<-make_snplist(efo= "squamous cell carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
ukb <- ieugwasr::associations(id="ukb-a-60", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Squamous cell carcinoma",population="European",pmid="ukb-a-60",ncase=404,ncontrol=336755,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=162,open_gwas=TRUE,efo= "squamous cell carcinoma")
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/scc162_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

####################
#Urinary tract cancer#
####################
snplist<-make_snplist_urinary()
ukb <- ieugwasr::associations(id="ukb-d-C3_URINARY_TRACT", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Urinary tract cancer",population="European",pmid="ukb-d-C3_URINARY_TRACT",ncase=1841,ncontrol=359353,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=164,open_gwas=TRUE,efo=c("kidney cancer","nephroblastoma","renal cell carcinoma","bladder carcinoma"))
dat<-transform_betas(dat=dat)
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/utc164_corrected.txt",sep="\t",col.names=T,row.names=F,quote=F)

####################
#Bladder cancer#
####################
snplist<-make_snplist(efo= "bladder carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_BLADDER", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Bladder cancer",population="European",pmid="finn-a-C3_BLADDER",ncase=366,ncontrol=96133,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=28,open_gwas=TRUE,efo= "bladder carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/bla28.txt",sep="\t",col.names=T,row.names=F,quote=F)

##############
#Blood cancer#
##############

snplist<-make_snplist_blood()
fin <- ieugwasr::associations(id="finn-a-CD2_PRIMARY_LYMPHOID_HEMATOPOIETIC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Blood cancer",population="European",pmid="finn-a-CD2_PRIMARY_LYMPHOID_HEMATOPOIETIC",ncase=1001,ncontrol=95498,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=29,open_gwas=TRUE,efo=c("lymphoma","multiple myeloma","lymphoid leukemia"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/blc29.txt",sep="\t",col.names=T,row.names=F,quote=F)

##############
#Brain cancer#
##############

snplist<-make_snplist(efo= "central nervous system cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_BRAIN", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Brain cancer",population="European",pmid="finn-a-C3_BRAIN",ncase=142,ncontrol=96357,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=30,open_gwas=TRUE,efo= "central nervous system cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/brc30",sep="\t",col.names=T,row.names=F,quote=F)

##############
#Breast cancer
##############

snplist<-make_snplist(efo= "breast carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_BREAST_3", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Breast cancer",population="European",pmid="finn-a-C3_BREAST_3",ncase=2589,ncontrol=93910,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=31,open_gwas=TRUE,efo= "breast carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/brc31",sep="\t",col.names=T,row.names=F,quote=F)

################
#Overall cancer#
################

snplist<-make_snplist(efo= "cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-ANY_CANC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-ANY_CANC",ncase=9792,ncontrol=86707,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=32,open_gwas=TRUE,efo= "cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/can32",sep="\t",col.names=T,row.names=F,quote=F)

################
#Overall cancer#
################

snplist<-make_snplist(efo= "cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-II_NEOPLASM", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Overall cancer",population="European",pmid="finn-a-II_NEOPLASM",ncase=31217,ncontrol=65282,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=33,open_gwas=TRUE,efo= "cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/can33",sep="\t",col.names=T,row.names=F,quote=F)

################
#Cancer of digestive organs#
################

snplist<-make_snplist(efo= "digestive system carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_DIGESTIVE_ORGANS", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Cancer of digestive organs",population="European",pmid="finn-a-C3_DIGESTIVE_ORGANS",ncase=1582,ncontrol=94917,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=34,open_gwas=TRUE,efo= "digestive system carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/cdo34",sep="\t",col.names=T,row.names=F,quote=F)

#######################################
#Central nervous system and eye cancer#
#######################################

# Malignant neoplasm of eye, brain and central nervous system
# Include	C3_EYE_ADNEXA, C3_MENINGES, C3_BRAIN, C3_SPINAL_CORD_CRANIAL_AND_OTHER_CNS
snplist<-make_snplist(efo= "central nervous system cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_EYE_BRAIN_NEURO", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Central nervous system and eye cancer",population="European",pmid="finn-a-C3_EYE_BRAIN_NEURO",ncase=207,ncontrol=96292,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=35,open_gwas=TRUE,efo= "central nervous system cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/cns35",sep="\t",col.names=T,row.names=F,quote=F)

####################
# Colorectal cancer#
####################

snplist<-make_snplist(efo= "colorectal cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-CUSTOM_COLORECTAL_CANCER_EXALLC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Colorectal cancer",population="European",pmid="finn-a-CUSTOM_COLORECTAL_CANCER_EXALLC",ncase=843,ncontrol=95656,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=36,open_gwas=TRUE,efo= "colorectal cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/crc36",sep="\t",col.names=T,row.names=F,quote=F)

####################
# Endocrine gland cancer#
####################

snplist<-make_snplist_end()
fin <- ieugwasr::associations(id="finn-a-C3_ENDOCRINE", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Endocrine gland cancer",population="European",pmid="finn-a-C3_ENDOCRINE",ncase=328,ncontrol=96171,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=37,open_gwas=TRUE,efo=c("pituitary gland adenoma","thyroid carcinoma","carcinoid tumor","neuroendocrine neoplasm"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/egc37",sep="\t",col.names=T,row.names=F,quote=F)

####################
# Endometrial cancer#
####################

snplist<-make_snplist(efo= "endometrial carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_CORPUS_UTERI", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Endometrial cancer",population="European",pmid="finn-a-C3_CORPUS_UTERI",ncase=366,ncontrol=53896,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=38,open_gwas=TRUE,efo= "endometrial carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/enc38",sep="\t",col.names=T,row.names=F,quote=F)

####################
#Female genital cancer#
####################

snplist<-make_snplist_fgen()
fin <- ieugwasr::associations(id="finn-a-C3_FEMALE_GENITAL", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Female genital cancer",population="European",pmid="finn-a-C3_FEMALE_GENITAL",ncase=672,ncontrol=53590,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=39,open_gwas=TRUE,efo=c("cervical carcinoma","endometrial carcinoma","ovarian carcinoma"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/fgc39",sep="\t",col.names=T,row.names=F,quote=F)

####################
#Follicular lymphoma#
####################

snplist<-make_snplist(efo="follicular lymphoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-CD2_FOLLICULAR_LYMPHOMA", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Follicular lymphoma",population="European",pmid="finn-a-CD2_FOLLICULAR_LYMPHOMA",ncase=158,ncontrol=96341,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=40,open_gwas=TRUE,efo="follicular lymphoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/fol40",sep="\t",col.names=T,row.names=F,quote=F)

####################
#Kidney cancer#
####################

snplist<-make_snplist(efo="kidney cancer",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_KIDNEY_NOTRENALPELVIS", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Kidney cancer",population="European",pmid="finn-a-C3_KIDNEY_NOTRENALPELVIS",ncase=301,ncontrol=96198,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=41,open_gwas=TRUE,efo="kidney cancer")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/kid41",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer#
#############

snplist<-make_snplist(efo="lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-LUNG_CANCER_MESOT", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid="finn-a-LUNG_CANCER_MESOT",ncase=673,ncontrol=95826,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=42,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/luc42",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lung cancer#
#############

snplist<-make_snplist(efo="lung carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_BRONCHUS_LUNG", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid="finn-a-C3_BRONCHUS_LUNG",ncase=516,ncontrol=95983,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=43,open_gwas=TRUE,efo="lung carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/luc43",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Lymphoid leukaemia#
#############

snplist<-make_snplist(efo="lymphoid leukemia",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-CD2_LYMPHOID_LEUKAEMIA", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Lymphoid leukaemia",population="European",pmid="finn-a-CD2_LYMPHOID_LEUKAEMIA",ncase=198,ncontrol=96301,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=44,open_gwas=TRUE,efo="lymphoid leukemia")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/lyl44",sep="\t",col.names=T,row.names=F,quote=F)

#############
# Male genital cancer
#############

snplist<-make_snplist_mgen()
fin <- ieugwasr::associations(id="finn-a-C3_MALE_GENITAL", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Male genital cancer",population="European",pmid="finn-a-C3_MALE_GENITAL",ncase=1887,ncontrol=94612,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=45,open_gwas=TRUE,efo=c("prostate carcinoma","testicular carcinoma"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/mgc45",sep="\t",col.names=T,row.names=F,quote=F)

#############
# Malignant skin cancer
#############

snplist<-make_snplist_malignant_skin_cancer()
fin <- ieugwasr::associations(id="finn-a-C3_SKIN", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="finn-a-C3_SKIN",ncase=895,ncontrol=95604,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=46,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/msc46",sep="\t",col.names=T,row.names=F,quote=F)

#############
# Multiple myeloma
#############
snplist<-make_snplist(efo="multiple myeloma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-CD2_MULTIPLE_MYELOMA_PLASMA_CELL", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Multiple myeloma",population="European",pmid="finn-a-CD2_MULTIPLE_MYELOMA_PLASMA_CELL",ncase=180,ncontrol=96319,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=47,open_gwas=TRUE,efo="multiple myeloma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/mmm47",sep="\t",col.names=T,row.names=F,quote=F)

#############
# Non-follicular lymphoma
#############
snplist<-make_snplist(efo="lymphoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-CD2_NONFOLLICULAR_LYMPHOMA", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Non-follicular lymphoma",population="European",pmid="finn-a-CD2_NONFOLLICULAR_LYMPHOMA",ncase=344,ncontrol=96155,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=48,open_gwas=TRUE,efo="lymphoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/nfl48",sep="\t",col.names=T,row.names=F,quote=F)

#############
#non-Hodgkin lymphoma unspecified
#############
snplist<-make_snplist(efo="non-Hodgkins lymphoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-CD2_NONHODGKIN_NAS", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="non-Hodgkin lymphoma unspecified",population="European",pmid="finn-a-CD2_NONHODGKIN_NAS",ncase=155,ncontrol=96344,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=49,open_gwas=TRUE,efo="non-Hodgkins lymphoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/nhl49",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Oral cavity and pharyngeal cancer
#############

snplist<-make_snplist_hnc()
fin <- ieugwasr::associations(id="finn-a-C3_LIP_ORAL_PHARYNX", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Oral cavity and pharyngeal cancer",population="European",pmid="finn-a-C3_LIP_ORAL_PHARYNX",ncase=234,ncontrol=96265,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=50,open_gwas=TRUE,efo=c("head and neck squamous cell carcinoma","oropharynx cancer","nasopharyngeal neoplasm","hypopharynx cancer","oral cavity cancer","mouth neoplasm","pharynx cancer"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/opc50",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Ovarian cancer
#############
snplist<-make_snplist(efo="ovarian carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_OVARY", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Ovarian cancer",population="European",pmid="finn-a-C3_OVARY",ncase=184,ncontrol=54078,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=51,open_gwas=TRUE,efo="ovarian carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/ovc51",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Pancreatic cancer
#############
snplist<-make_snplist(efo="pancreatic carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_PANCREAS", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Pancreatic cancer",population="European",pmid="finn-a-C3_PANCREAS",ncase=229,ncontrol=96270,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=52,open_gwas=TRUE,efo="pancreatic carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pac52",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Prostate cancer
#############
snplist<-make_snplist(efo="prostate carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_PROSTATE", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Prostate cancer",population="European",pmid="finn-a-C3_PROSTATE",ncase=1824,ncontrol=40413,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=53,open_gwas=TRUE,efo="prostate carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pro53",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Respiratory and intrathoracic cancer
#############

snplist<-make_snplist_respiratory()
fin <- ieugwasr::associations(id="finn-a-C3_RESPIRATORY_INTRATHORACIC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="finn-a-C3_RESPIRATORY_INTRATHORACIC",ncase=615,ncontrol=95884,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=54,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/ric54",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Thyroid cancer
#############

snplist<-make_snplist(efo="thyroid carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
fin <- ieugwasr::associations(id="finn-a-C3_THYROID_GLAND", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Thyroid cancer",population="European",pmid="finn-a-C3_THYROID_GLAND",ncase=321,ncontrol=96178,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=55,open_gwas=TRUE,efo="thyroid carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/thc55",sep="\t",col.names=T,row.names=F,quote=F)

#############
#Urinary tract cancer
#############

snplist<-make_snplist_urinary()
fin <- ieugwasr::associations(id="finn-a-C3_URINARY_TRACT", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Urinary tract cancer",population="European",pmid="finn-a-C3_URINARY_TRACT",ncase=690,ncontrol=95809,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=56,open_gwas=TRUE,efo=c("kidney cancer","nephroblastoma","renal cell carcinoma","bladder carcinoma"))
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/utc56",sep="\t",col.names=T,row.names=F,quote=F)

# pancreatic cancer
# eaf/maf not reported
snplist<-make_snplist(efo="pancreatic carcinoma",population="European",Dir="~/fatty-acids/outcome_data/data/")
pan <- ieugwasr::associations(id="ieu-a-822", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(pan,stringsAsFactors=F),outcome="Pancreatic cancer",population="European",pmid=19648918,ncase=1896,ncontrol=1939,study="PanScan I",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=122,open_gwas=TRUE,efo="pancreatic carcinoma")
write.table(dat,"~/fatty-acids/outcome_data/data/harmonised/pan122.txt",sep="\t",col.name=T,row.name=F,quote=F)
