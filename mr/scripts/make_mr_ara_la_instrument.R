library(TwoSampleMR)
library(ieugwasr)
library(CheckSumStats)

source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/outcome_data/scripts/functions_combine_and_format_outcomes.R")
# Exposure
snplist <- readLines("~/fatty-acids/mr/data/ara_la_snplist_independent.txt")
ara_snplist<-c("rs4985155","rs174546","rs3734398")
la_nmr<-ieugwasr::tophits(id="met-d-LA",clump=0)
la_nmr1<-la_nmr[la_nmr$rsid %in% snplist, ]
head(la_nmr1)
la_nmr1$n<-114999	
la_nmr1<-predict_beta_sd(dat=data.frame(la_nmr1),beta="beta",se="se",eaf="eaf",sample_size="n",pval="p") #looks like is already standardised
# head(la_nmr1)
# summary(la_nmr1$beta_sd)

ara_dat<-extract_snps(path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.fixed.tab",snplist=ara_snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=FALSE)

ara_dat_schs<-extract_snps(path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.fixed.tab",snplist=ara_snplist,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=FALSE)

ara_dat<-predict_beta_sd(dat=ara_dat,beta="beta",se="se",eaf="effect_allele_freq",sample_size="n",pval="p")
ara_dat$study<-"CHARGE"

write.table(ara_dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/ara_instrument.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

ara_dat<-read.table("~/fatty-acids/mr/data/ara_instrument.txt",sep="\t",head=TRUE,quote="")

head(ara_dat)
ara_dat$id<-"ara"
ara_dat$exposure<-"Arachidonic acid"
exposure_dat1<-format_exposure2(dat=ara_dat,standardise_beta=TRUE,beta="beta",se="se",pval="p",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="id",exposure="exposure",samplesize="n")
exposure_dat1$id.exposure<-"exp1"

la_nmr1$id<-"la"
exposure_dat2<-format_exposure2(dat=la_nmr1,standardise_beta=FALSE,beta="beta",se="se",pval="p",effect_allele="ea",other_allele="nea",eaf="eaf",rsid="rsid",ID="id",exposure="trait",samplesize="n")
exposure_dat2$id.exposure<-"exp2"
exposure_dat<-plyr::rbind.fill(exposure_dat2,exposure_dat1)
save(exposure_dat,file="~/fatty-acids/mr/data/exposure_dat_ara_la.RData")

