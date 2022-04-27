setwd("C:/Users/stephen burgess/Dropbox (Cambridge University)/People/Fatima Batool")
table(Dat$trait)

load("fattyacids_mvmr/data_for_cismvmr_esophageal_scc.RData")
out = Dat[Dat$trait=="Esophageal squamous cell carcinoma",]
write.csv(out, "fattyacids_mvmr/esophageal_scc.csv")

load("fattyacids_mvmr/data_for_cismvmr_lung_cancer_bj.RData")
out = Dat[Dat$trait=="Lung cancer",]
write.csv(out, "fattyacids_mvmr/lung_cancer_bj.csv")

load("fattyacids_mvmr/data_for_cismvmr_colorectal_accc.RData")
out = Dat[Dat$trait=="Colorectal cancer",]
write.csv(out, "fattyacids_mvmr/colorectal_accc.csv")

load("fattyacids_mvmr/data_for_cismvmr_malignant_skin_cancer_ukb.RData")
out = Dat[Dat$trait=="Malignant skin cancer",]
write.csv(out, "fattyacids_mvmr/malignant_skin_cancer_ukb.csv")

load("fattyacids_mvmr/data_for_cismvmr_lung_cancer_ilcco_ukb.RData")
out = Dat[Dat$trait=="Lung cancer",]
write.csv(out, "fattyacids_mvmr/lung_cancer_ilcco_ukb.csv")

load("fattyacids_mvmr/data_for_cismvmr_lung_cancer_tricl.RData")
out = Dat[Dat$trait=="Lung cancer",]
write.csv(out, "fattyacids_mvmr/lung_cancer_tricl.csv")


