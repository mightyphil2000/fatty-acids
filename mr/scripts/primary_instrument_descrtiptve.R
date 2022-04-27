# rm(list=ls())
source("~/fatty-acids/mr/scripts/mr_functions.R")
# source("~/fatty-acids/outcome_data/scripts/combine_and_format_outcomes.R")
library(TwoSampleMR)
library(tidyr)
library(data.table)
library(plyr)

exp<-read.table("~/fatty-acids/mr/data/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",head=T,stringsAsFactors=F)
exp[exp$exposure == "AA:DGLA" & exp$population=="East Asian",]

exposure_dat<-format_exposure2(dat=exp,standardise_beta=TRUE)
exposure_dat[exposure_dat$population=="East Asian",]
a<-exposure_dat[exposure_dat$SNP == "rs174546" & exposure_dat$exposure %in% c("AA:DGLA","GLA:LA"),]
z<-a$beta.exposure/a$se.exposure

a$z<-z
a$pval.exposure2<-pnorm(abs(a$z),lower.tail=FALSE)*2
beta<-a$beta.exposure*-1
a$beta.exposure<-beta
eaf<-1-a$eaf.exposure
a$eaf.exposure<-eaf
names(a)[names(a) == "effect_allele.exposure"]<-"other_allele"
names(a)[names(a) == "other_allele.exposure"]<-"effect_allele"
names(a)<-gsub(".exposure","",names(a))
a$pval2[a$pval2 == 0]<- pnorm(37.5,lower.tail=FALSE)*2
write.table(a,"~/fatty-acids/mr/results/primary_instrument_descriptive.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# head -1 score_lnD5D_control_allchr_qc.txt
#  grep -w rs174546 score_lnD5D_control_allchr_qc.txt
#  grep -w rs174546 score_lnD5D_case_allchr_qc.txt