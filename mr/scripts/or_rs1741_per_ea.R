# BiocManager::install("multtest")
library(ggforestplot)
library(ggplot2)
library(TwoSampleMR)
library(devtools)
setwd("~/TwoSampleMR")
load_all()
document()

source("~/fatty-acids/mr/scripts/mr_functions.R")
# exposure
source("~/fatty-acids/mr/scripts/mr_functions.R")
load("~/fatty-acids/mr/data/instruments.Rdata")

exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")

# exp_dat<-exp_dat[!exp_dat$FADS,]

#Confirmed one SNP per FADS region per exposure
# Temp<-exp_dat[exp_dat$FADS,c("chr","position","id.exposure","SNP","exposure")]
# Temp[order(Temp$id.exposure),]

# Outcome
# Dat2
load("~/fatty-acids/mr/data/outcome_dat_secondary_pufas.Rdata")
# ,"rs881803"


unique(exp_dat$exposure)
exp_dat[exp_dat$SNP=="rs3798713",]
exp_dat<-exp_dat[which(exp_dat$exposure=="Arachidonic acid (20:4n6)"),]

# harmonise data
dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat =outcome_dat,tolerance=0.05)
# "rs174546","rs1741","rs881803"
dat$beta.exposure
a<-dat[dat$SNP == "rs1741" & dat$outcome2 %in% c("Colorectal cancer", "Esophageal squamous cell carcinoma","Lung cancer","Basal cell skin carcinoma","Basal cell carcinoma","Overall cancer"),c("effect_allele.outcome","outcome2","beta.outcome","se.outcome","beta.exposure","ncase.outcome","study","id.outcome")]
a<-a[order(a$ncase.outcome,decreasing=TRUE),]
a$outcome2[a$outcome2=="Basal cell carcinoma" ]<-"Basal cell skin carcinoma" 
a<-a[!duplicated(a$outcome2),]
a$beta.exposure
a$beta.outcome<-a$beta.outcome*-1
a$OR<-exp(a$beta.outcome)
a$lci<-exp(a$beta.outcome-1.96*a$se.outcome)
a$uci<-exp(a$beta.outcome+1.96*a$se.outcome)
a[,c("outcome2","id.outcome","OR","lci","uci")]
