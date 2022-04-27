# install.packages("data.table")
library(data.table)
setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/CSI_results")
a<-fread("2019.10.02 Lifetime Smoking GWAS Data Sheet 1.txt")
b<-a[a$SNP %in% c("rs174546","rs2524299") ,]
write.table(b,"rs174546_rs2524299_lookup_csi_ukb.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# cd ~/fatty-acids/mr/data
# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/CSI_results/rs174546_rs2524299_lookup_csi_ukb.txt .
source("~/fatty-acids/mr/scripts/mr_functions.R")
out_dat<-read.table("~/fatty-acids/mr/data/rs174546_rs2524299_lookup_csi_ukb.txt",sep="\t",head=T,stringsAsFactors=F)
out_dat<-format_outcomes_csi(dat=out_dat)



exp<-read.table("~/fatty-acids/mr/data/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",head=T,stringsAsFactors=F)
exp[exp$exposure == "AA:DGLA" & exp$population=="East Asian",]
exposure_dat<-format_exposure2(dat=exp,standardise_beta=TRUE)
exposure_dat<-exposure_dat[exposure_dat$exposure=="AA:DGLA" & exposure_dat$population=="European" & exposure_dat$SNP == "rs174546",]

library(TwoSampleMR)
dat <- harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = out_dat
)


res<-mr(dat,method_list = "mr_wald_ratio" )

write.table(res,"~/fatty-acids/mr/results/mr_results_aadgla_csi.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

out_dat$beta<-round(out_dat$beta.outcome,3)
out_dat$lci<-round(out_dat$beta-1.96*out_dat$se,3)
out_dat$uci<-round(out_dat$beta+1.96*out_dat$se,3)


# This study was a genome-wide association study (GWAS) of a lifetime smoking index (which combined smoking initiation, duration, heaviness and cessation), conducted in a sample of 462,690 current, former and never smokers in the UK Biobank. Full details of the methods are given in Wootton et al (2019). The GWAS was conducted via the University of Bristol MRC IEU GWAS pipeline (Elsworth et al, 2017) using BOLT LMM (Loh et al., 2015), which accounts for population stratification and relatedness using linear mixed modelling. 