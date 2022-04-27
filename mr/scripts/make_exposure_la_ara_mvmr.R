# source("~/fatty-acids/mr/scripts/mr_functions.R")
# devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)
library(CheckSumStats)
library(TwoSampleMR)
# remotes::install_github("MRCIEU/TwoSampleMR")
# install.packages("devtools")

# devtools::install_github("MRCIEU/CheckSumStats")

# load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")
# load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/exposure_dat_ara_la.Rdata")
# snps.ara<-exposure_dat$SNP[exposure_dat$exposure=="Arachidonic acid"]
# snps.ara[snps.ara %in% exposure_dat$SNP[exposure_dat$exposure=="Linoleic acid"]]

# exposure_dat$SNP[exposure_dat$exposure=="Linoleic acid"]

snplist_ara_la<-readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/ara_la_snplist_independent.txt")

la_nmr_mv<-associations(id="met-d-LA",variants=unique(snplist_ara_la))
la_nmr_mv$id <-"LA"
la_nmr_mv$n<-114999
# la_nmr_mv$p #no 0 p values

la_clump<-ieugwasr::ld_clump( clump_r2 = 0.001,
    dplyr::tibble(rsid=la_nmr_mv$rsid, pval=la_nmr_mv$p, id=la_nmr_mv$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

la_nmr_mv1<-la_nmr_mv[la_nmr_mv$rsid %in% la_clump$rsid,]

# ara_dat<-extract_snps(path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.fixed.tab",snplist=exposure_dat$SNP,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=FALSE)
# ara_dat<-predict_beta_sd(dat=ara_dat,beta="beta",se="se",eaf="effect_allele_freq",sample_size="n",pval="p")
# ara_dat$study<-"CHARGE"
# dim(ara_dat)

ara_dat_imp<-extract_snps(path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/imputed/N6meta2041.tbl.fixed.tab",snplist=la_clump$rsid,exact_match=TRUE,path_to_target_file_sep="\t",Test.gz=FALSE)
ara_dat_imp<-predict_beta_sd(dat=ara_dat_imp,beta="beta",se="se",eaf="effect_allele_freq",sample_size="n",pval="p")
ara_dat_imp$study<-"CHARGE"
dim(ara_dat_imp)
ara_dat_imp$id<-"ara"
ara_dat_imp$exposure<-"Arachidonic acid"

exposure_dat1<-format_exposure2(dat=ara_dat_imp,standardise_beta=TRUE,beta="beta",se="se",pval="p",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="id",exposure="exposure",samplesize="n")
exposure_dat1$id.exposure<-"ara"

# la_nmr_mv$id<-"la"
la_nmr_mv1$exposure<-"Linoleic acid"
la_nmr_mv1$outcome<-"Linoleic acid"
exposure_dat2<-format_exposure2(dat=data.frame(la_nmr_mv1),standardise_beta=FALSE,beta="beta",se="se",pval="p",effect_allele="ea",other_allele="nea",eaf="eaf",rsid="rsid",ID="id",exposure="trait",samplesize="n")
exposure_dat2$id.exposure<-"la"
exposure_dat2$id.outcome<-"la"
Dat<-list(exposure_dat1,exposure_dat2)
save(Dat,file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/exposure_dat2_ara_la_mvmr_preharmonise.RData")
load("~/fatty-acids/mr/data/exposure_dat2_ara_la_mvmr_preharmonise.RData")  
exposure_dat1<-data.frame(Dat[1])
exposure_dat2<-data.frame(Dat[2])
exposure_dat2<-exposure_dat2[,!names(exposure_dat2) %in% c("exposure","id.exposure")]
# exposure_dat2$id.outcome
names(exposure_dat2)<-gsub("exposure","outcome",names(exposure_dat2))

exposure_dat<-harmonise_data(exposure_dat = exposure_dat1,outcome_dat =exposure_dat2,action=1)

names.outcome<-names(exposure_dat)[grep("outcome",names(exposure_dat))]
exposure_dat_la<-exposure_dat[,c("SNP",names.outcome)]
exposure_dat_la$exposure<-"Linoleic acid"
names(exposure_dat_la)<-gsub("outcome","exposure",names(exposure_dat_la))
names.exposure<-names(exposure_dat)[grep("exposure",names(exposure_dat))]
exposure_dat_ara<-exposure_dat[,c("SNP",names.exposure)]
exposure_dat<-rbind(exposure_dat_la,exposure_dat_ara)

save(exposure_dat,file="~/fatty-acids/mr/data/exposure_dat_ara_la_mvmr.RData")
# exposure_dat<-plyr::rbind.fill(exposure_dat2,exposure_dat1)

# all(exposure_dat1$SNP %in% exposure_dat2$SNP)

# # exposure_dat$id<-"id"
# la_clump<-ieugwasr::ld_clump( clump_r2 = 0.001,
#     dplyr::tibble(rsid=exposure_dat2$SNP, pval=exposure_dat2$pval, id=exposure_dat2$id),
#     plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
#     bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
# )

# LD_matrix<-ld_matrix(exposure_dat2$SNP)
# snp<-exposure_dat2$SNP[!exposure_dat2$SNP %in% la_clump$rsid]
# Pos<-grep(snp,row.names(LD_matrix))
# Pos2<-which(row.names(LD_matrix)=="rs174564_G_A")
# LD_matrix[Pos,Pos2]
# la_clump$rsid



# save(exposure_dat,file="~/fatty-acids/mr/data/exposure_dat_ara_la.RData")



# la_nmr_clump<-la_nmr[la_nmr$rsid %in% la_clump$rsid,]
# la_snplist<-la_nmr_clump$rsid


#  # rs780093 
# ara_snplist<-c("rs4985155","rs174546","rs3734398")
# snplist<-unique(c(ara_snplist,la_snplist))
# class(snplist)
# # snplist<-ara_snplist

# ara_snplist

# write.table(snplist,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/ara_la_snplist_independent.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)

# transfer data
# cp /newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/N6meta2041.tbl.fixed.tab  /projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/imputed/


