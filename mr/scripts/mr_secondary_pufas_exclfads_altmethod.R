# Rerun analysis of lung cancer using alternative method suggested by review in second round of reviews - i.e. meta analyse the summary data across SNPs in outcome_dat before harmonising with the exposure data and running the MR analysis. Lung cancer is the only outcome analysed where multiple SNPs were used and where MR results were combined across separate studies. The analysis method suggested by the review will give same result as current analysis when only a single SNP is used - hence why we don't need to repeat the analysis for the primary analyses. 

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

exp_dat<-exp_dat[!exp_dat$FADS,]

#Confirmed one SNP per FADS region per exposure
# Temp<-exp_dat[exp_dat$FADS,c("chr","position","id.exposure","SNP","exposure")]
# Temp[order(Temp$id.exposure),]

# Outcome
# Dat2
load("~/fatty-acids/mr/data/outcome_dat_secondary_pufas.Rdata")
out<-outcome_dat[which(outcome_dat$outcome2 == "Lung cancer"),]
out<-meta_analysis_v3(dat2=out[out$study %in% c("UKB","ILCCO"),],beta.col="beta.outcome",se.col="se.outcome",outcome="outcome2",ncase="ncase.outcome",ncontrol="ncontrol.outcome")
out<-out[which(out$nstudies==2),]

# harmonise data
dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat =out,tolerance = 0.05) #relax tolerance from 0.08 to 0.05
res<-mr(dat,method_list=c("mr_ivw","mr_wald_ratio","mr_ivw_fe"))
res$id.analysis<-paste0(res$id.outcome,res$id.exposure)
# TwoSampleMR::mr_method_list()
pos1<-which(res$nsnp==1)
pos2<-which(res$nsnp==2)
pos3<-which(res$nsnp>2)
res1_1<-res[pos2,] #where nsnps == 2, use ivw fe
res1_1<-res1_1[res1_1$method == "Inverse variance weighted (fixed effects)",]
res1_2<-res[which(res$nsnp!=2),]

res1_2<-subset_on_method(mr_res=res1_2,multi_snp_method="Inverse variance weighted")
res1<-rbind(res1_1,res1_2)
dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat =out,tolerance=0.05)
het<-mr_heterogeneity(dat)
# note that fisherp method is used to combine Q P values across studies. In this altmethod,it is not needed because we combine summary data across studies at the SNP level before MR. however we are ignoring the between study heterogeneity between SNPs . 
#het_fisherp<-meta_p(dat=het)

res_single <- mr_singlesnp(dat)

save(res1,file="~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_altmethod.Rdata")
save(res_single,file="~/fatty-acids/mr/results/res_single_secondary_pufas_exclfads_altmethod.Rdata")
save(het,file="~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity_altmethod.Rdata")
#save(het_fisherp,file="~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity_fisher_p_altmethod.Rdata")
