source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

load("~/fatty-acids/outcome_data/data/harmonised_data_preqc.Rdata")

Dat1<-Dat[Dat$ID == 68,]
Pred_lnor<-pred_lnor_sh(dat2=Dat1)