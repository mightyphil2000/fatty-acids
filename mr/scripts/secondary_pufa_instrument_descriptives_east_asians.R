load("~/fatty-acids/mr/data/instruments.Rdata")

source("~/fatty-acids/mr/scripts/mr_functions.R")

exp_dat<-format_exposure3(dat=eas,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")

exposures<-unique(exp_dat$exposure)

length(unique(exp_dat$SNP))
table(exp_dat$exposure[exp_dat$FADS])
table(exp_dat$exposure[!exp_dat$FADS])
table(exp_dat$exposure,exp_dat$FADS)
length(unique(exp_dat$SNP[!exp_dat$exposure %in% excl & !exp_dat$FADS]))


exp_dat[exp_dat$exposure==exposures[5],c("SNP","exposure","FADS","chr","position","r2")]
temp<-exp_dat[exp_dat$exposure==exposures[2],c("chr","r2")]
sum(temp$r2[!duplicated(temp$chr)])

var<-1
exp_dat$r2<-2*exp_dat$beta.exposure^2*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure)/var


