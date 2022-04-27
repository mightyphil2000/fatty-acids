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

# load("~/fatty-acids/mr/data/clumped_snps_gwis_charge.Rdata")
# exposure_dat<-format_exposure(dat=FA2,exposure="AA_to_DGLA",snp="rs174546")

load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
load("~/fatty-acids/mr/data/dat_outcomes_final.Rdata")

# dat_outcomes_final[which(dat_outcomes_final$cancer %in% c("Lung cancer","Colorectal cancer")),c("lnor","study","outcome")]

# dat_outcomes_final$outcome[grep("breast",dat_outcomes_final$outcome,ignore.case=TRUE)]

# dat_list<-format_outcomes(dat=dat_outcomes_final,IDS=disc.tab9$ID)
# # dat_list<-format_outcomes(dat=dat_outcomes_final,IDS=meta.tab9$ID)
# outcome_disc<-data.frame(dat_list[1])
# # which(outcome_rep$id.outcome == "7")
# outcome_rep<-data.frame(dat_list[2])

outcome_dat<-format_outcomes2(dat=dat_outcomes_final)
# outcome_dat<-format_outcomes2(dat=outcome_rep)
# which(outcome_dat$id.outcome == "7")
# outcome_dat[grep("Esophageal",outcome_dat$outcome),]
# outcome_disc[outcome_disc$Cancer.Group == "Esophageal cancer","outcome"]

# table(outcome_disc$system)
# outcome_rep$id.outcome
# table(outcome_rep$system)
# sum(as.numeric(outcome_disc$ncase))
# sum(as.numeric(outcome_disc$ncontrol))

# names(outcome_dat)
# outcome_dat <- read_outcome_data(filename="~/fatty-acids/mr/data/outcome_discovery.txt", sep = "\t",snp_col="SNP")
# outcome_disc$id.outcome
# mr_res_disc<-tsmr(exposure_dat = exposure_dat,outcome_dat=outcome_disc,meta.dat=disc.tab9)
# mr_res_disc$id.outcome
# ids<-mr_res_disc$id.outcome[grep(";",mr_res_disc$id.outcome)]
# unique(ids[!ids %in% mr_res1$id.outcome ])
# unique(ids[!ids %in% mr_res1$id.outcome ])
# mr_res1$id.outcome[grep("99",mr_res1$id.outcome)]
# mr_res$id.outcome[grep("10§1",mr_res$id.outcome)]
# save( mr_res_disc,file="~/fatty-acids/mr/results/mr_results_discovery_v2.Rdata")
# mr_res$id.outcome

mr_res<-tsmr(exposure_dat = exposure_dat,outcome_dat=outcome_dat,meta.dat=meta.tab9)
# mr_res[grep("Lung cancer",mr_res$outcome),]
# mr_res[which(mr_res$id.outcome == "75"),c("AA:")]
 
mr_res1<-mr_res[mr_res$population == "European" & mr_res$exposure=="AA:DGLA",]
mr_res2<-mr_res[mr_res$population == "East Asian" & mr_res$exposure=="GLA:LA",]
mr_res1<-rbind(mr_res1,mr_res2)

# mr_res_meta$id.outcome
mr_res_meta<-mr_meta_analysis(Dat=mr_res1)
mr_res1<-rbind.fill(mr_res_meta,mr_res1)
# dim(mr_res1)
# mr_res1$id.outcome[mr_res1$id.outcome == "86"]
# mr_res_meta$id.outcome[grep("86",mr_res_meta$id.outcome )]
# mr_res1$id.outcome[grep("86",mr_res1$id.outcome )]
# "86; 5; 49"
save(mr_res1,file="~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")

# mr_res1<-mr_res1[!is.na(mr_res1$cases),]
# mr_res1$cases<-as.numeric(mr_res1$cases)
# mr_res3<-mr_res1[mr_res1$cases>=1000,]
# mr_res3$cancer[is.na(mr_res3$cancer)]<-mr_res3$outcome[is.na(mr_res3$cancer)]


# mr_res3[mr_res3$cancer == "Breast cancer",]
# unique(mr_res3$cancer)

# mr_res3<-mr_res3[order(mr_res3$system),]
# unique(mr_res3$cancer[mr_res3$system == "Reproductive"])
# table(mr_res3$cancer,mr_res3$system)
# combine MR results by meta analysis. Better to do this for outcomes from Europeans and East Asians


# load("~/fatty-acids/mr/results/mr_results_discovery.Rdata")

# table(mr_res2$cancer[mr_res2$pval<0.05],mr_res2$system[mr_res2$pval<0.05])

# mr_res3<-mr_res2[order(mr_res2$pval),c("outcome","pval","Cancer.Group") ]
# mr_res3<-mr_res3[mr_res3$pval<0.05,]
# length(table(mr_res3$Cancer.Group))
# dim(mr_res3)
# length(table(mr_res2$site))

# Temp<-mr_res[mr_res$p<0.05/68,]
# Temp[order(Temp$pval),]
# unique(mr_res2[mr_res2$pval<0.05,c("cancer","id.outcome")])


# source("~/TwoSampleMR/R/harmonise.R")
