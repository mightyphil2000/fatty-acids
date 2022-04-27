exposure_dat<-read_exposure_data(
       "~/fatty-acids/mr/data/exposure_dat.txt",
       clump = FALSE,
       sep = "\t",
       phenotype_col = "exposure",
       snp_col = "SNP",
       beta_col = "beta",
       se_col = "se",
       eaf_col = "eaf",
       effect_allele_col = "effect_allele",
       other_allele_col = "other_allele",
       pval_col = "p",
       samplesize_col = "n")

outcome_dat <- read_outcome_data(filename="~/fatty-acids/mr/data/outcome_dat.txt", sep = "\t",
       phenotype_col = "outcome", snp_col = "rsid", beta_col = "lnor",
       se_col = "se", eaf_col = "eaf", effect_allele_col = "Effect.Allele",
       other_allele_col = "Other.Allele", pval_col = "p",
       ncase_col = "ncase", ncontrol_col = "ncontrol")

dat <- harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = outcome_dat,action=2)

mr_res<-mr(dat,method_list=c("mr_wald_ratio"))
dat$wr<-dat$beta.outcome/dat$beta.exposure
dat$sewr<-dat$se.outcome/abs(dat$beta.exposure)

dat<-dat[dat$SNP=="rs4985155",]
dat$z<-dat$wr/dat$sewr
dat[,c("SNP","wr","sewr","z","outcome","exposure")]
