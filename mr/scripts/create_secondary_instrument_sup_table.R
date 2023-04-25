source("~/fatty-acids/mr/scripts/mr_functions.R")
# exposure
source("~/fatty-acids/mr/scripts/mr_functions.R")
load("~/fatty-acids/mr/data/instruments.Rdata")
exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")


exp_dat[exp_dat$exposure=="Eicosapentaenoic acid (20:5n3)",c("id","exposure","SNP","id.exposure","chr","eaf.exposure")]

head(exp_dat)
sup_tab<-exp_dat[,c("exposure","SNP","beta.exposure","se.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","consortium","FADS","chr","position","samplesize.exposure","pval.exposure","pmid","id")]