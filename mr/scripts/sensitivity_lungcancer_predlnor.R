#########################################
# Predicted log odds ratios versus reported effect sizes #
########################################
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
load("~/fatty-acids/outcome_data/data/harmonised_data.Rdata")
Datqc<-Dat
load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")
Dat<-rbind(dat1,dat2)
Dat<-estimate_bias(dat=Dat)
Dat<-fix_info(dat=Dat) 


illco<-Dat[which(Dat$ID==75),]
head(illco)


illco<-illco[which(illco$rsid %in% snps),]
illco$p_sh<-pnorm(abs(illco$lnor_sh/illco$se_sh),lower.tail=F)*2
Datqc<-Datqc[Datqc$ID %in% c(75,17,42,149),]
Datqc
