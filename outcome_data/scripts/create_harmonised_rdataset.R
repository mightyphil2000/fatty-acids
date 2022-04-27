source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

Dat<-collate_dat()
Dat<-basic_qc(dat=Dat)
save(Dat,file="~/fatty-acids/outcome_data/data/harmonised_data.Rdata")