source("~/fatty-acids/mr/scripts/mr_functions.R")
out_dat<-read.table("~/fatty-acids/mr/data/rs174546_rs2524299_lookup_csi_ukb.txt",sep="\t",head=T,stringsAsFactors=F)
out_dat<-format_outcomes_csi(dat=out_dat)


out_dat$beta<-round(out_dat$beta.outcome,3)
out_dat$lci<-round(out_dat$beta-1.96*out_dat$se,3)
out_dat$uci<-round(out_dat$beta+1.96*out_dat$se,3)
