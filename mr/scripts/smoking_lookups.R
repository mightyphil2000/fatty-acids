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
exposure_dat<-exposure_dat[exposure_dat$exposure=="AA:DGLA" & exposure_dat$population=="European" & exposure_dat$SNP == "rs174546",]

# load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")
# exposure_dat[which(exposure_dat$SNP=="rs174546"),]

ao<-ieugwasr::gwasinfo()

pos<-grep("smok",ao$trait,ignore.case=TRUE)
ao1<-ao[pos,]
ao1<-ao1[order(ao1$sample_size,decreasing=TRUE),]
# data.frame(ao1[1:20,c("trait","sample_size","consortium","id")])

data.frame(ao[which(ao$consortium == "GSCAN"),])
IDS<-ao$id[which(ao$consortium == "GSCAN")]

# ao1<-ao[grep("smok",ao$trait,ignore.case=TRUE),]
# ao1<-ao1[order(ao1$sample_size,decreasing=TRUE),]
# ao1<-ao1[!duplicated(ao1$trait),]
# ao1$trait
# ao1<-ao1[ao1$trait %in%  c("Ever smoked","Cigarettes smoked per day"),]
# ids1<-ao1$id
# data.frame(ao1[,c("consortium","trait","sample_size")])
# ids<-ao$id[grep("TAG",ao$consortium,ignore.case=TRUE)]
# a<-ieugwasr::associations(id=ids, variants="rs174546",proxies=0) 
b<-ieugwasr::associations(id=IDS, variants="rs174546",proxies=0) 
# rs2524299 
data.frame(b[,c("ea","nea","trait","p","beta")])
data.frame(b[b$p<0.05,c("trait","p","beta")])

out_dat <- TwoSampleMR::extract_outcome_data(
    snps = "rs174546",
    outcomes = ids1
)

# out_dat <- TwoSampleMR::extract_outcome_data(
#     snps = "rs174546",
#     outcomes = ids
# )

unlist(b_sd(dat=out_dat)[1])

out_dat$beta.outcome[out_dat$outcome ==  "Cigarettes smoked per day || id:ieu-b-142"] <- unlist(b_sd()[1])
out_dat$se.outcome[out_dat$outcome ==  "Cigarettes smoked per day || id:ieu-a-961"] <- unlist(b_sd()[2])
exposure_dat
out_dat[out_dat$outcome ==  "Ever vs never smoked || id:ieu-a-962",]
lnor<-out_dat$beta.outcome[out_dat$outcome ==  "Ever vs never smoked || id:ieu-a-962"]
lnor<-lnor*-1
se<-out_dat$se.outcome[out_dat$outcome ==  "Ever vs never smoked || id:ieu-a-962"]
or<-round(exp(lnor),2)
lci<-round(exp(lnor-1.96*se),2)
uci<-round(exp(lnor+1.96*se),2)
c(or,lci,uci)
 
dat <- harmonise_data(exposure_dat = exposure_dat,outcome_dat =out_dat)

res<-mr(dat,method_list="mr_wald_ratio")
lnor<-res$b[res$outcome ==  "Ever vs never smoked || id:ieu-a-962"]
se<-res$se[res$outcome ==  "Ever vs never smoked || id:ieu-a-962"]

or<-round(exp(lnor),2)
lci<-round(exp(lnor-1.96*se),2)
uci<-round(exp(lnor+1.96*se),2)
c(or,lci,uci)



b_sd<-function(dat=NULL){
	dat<-dat[dat$outcome == "Cigarettes smoked per day || id:ieu-b-142",]
	N<-dat$samplesize	
	z<-dat$beta.outcome/dat$se.outcome
	p<-exposure_dat$eaf
	n<-N
	b_sd<-z/sqrt(2*p*(1-p)*(n+z^2)) 
	se_sd<-b_sd/z
	return(list(b_sd,se_sd))
}

