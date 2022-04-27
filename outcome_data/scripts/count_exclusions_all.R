source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

Dat<-collate_dat(postqc=FALSE) 
Dat<-basic_qc(dat=Dat)

nrow(Dat)
length(unique(Dat$study))
length(unique(Dat$ID))

snp1<-readLines("~/fatty-acids/outcome_data/data/snplist_Europeans_rsidsonly2.txt")
snp2<-readLines("~/fatty-acids/outcome_data/data/snplist_East_Asians_rsidsonly2_nodups.txt")

snps_fa<-unique(c(snp1,snp2))

Dat1<-Dat[Dat$summary_set == "FAregions",]
Dat1<-Dat1[Dat1$rsid %in% snps_fa,]

Dat2<-Dat[Dat$summary_set != "FAregions",]
# exclude FA region SNPs that anre't in the fatty acid snp set
Dat_xfar<-rbind(Dat1,Dat2)
dim(Dat_xfar)
dim(Dat)

# N SNPs that correspond to the fatty acid SNP set
dim(Dat[Dat$rsid %in% snps_fa,])
dim(Dat)

# datasets with incorrect effect allele 
N_incorrect_effect_allele<-nrow(Dat[Dat$ID == 967,])

# drop datasets Ambiguous effect allele column
Dat2<-Dat[!Dat$ID %in% c(106,5,74),]
nrow(Dat2)
length(unique(Dat2$study))
length(unique(Dat2$ID))

# drop datasets with info<0.8
Dat3<-prune_info80(Dat=Dat)
N<-nrow(Dat)-nrow(Dat3)
N/nrow(Dat)*100 #expected nrow(Dat)=401026
length(unique(Dat3$study))

# drop log odds ratios > 1 or <-1 (includes asociations with unusually large Zb scores (>100))
# glioma133
Dat4<-Dat[which(abs(Dat$lnor)<1),]
N<-nrow(Dat)-nrow(Dat4)
N/nrow(Dat)*100
length(unique(Dat4$study))

# # set eaf to missing
# Dat4$eaf[Dat4$ID==132]<-NA
# Dat4$eaf[Dat4$ID==133]<-NA
# Dat4$eaf[Dat4$ID==67]<-NA
# Dat4$eaf[Dat4$ID==68]<-NA

length(which(is.na(Dat$eaf[Dat$ID %in% c(132,133,67,68)])))
N_eaf_to_na<-length(which(!is.na(Dat$eaf[Dat$ID %in% c(132,133,67,68)])))

##################################################################
# drop SNPs that were contributed < median number of studies
##################################################################
Dat$Nstudies<-NA
Dat$Nstudies_median<-NA
# Dat3$Nstudies<-NA
# Dat3$Nstudies_median<-NA
Dat4<-Dat[!is.na(Dat$Direction),]
Dat5<-Dat[is.na(Dat$Direction),]

Dat4<-compute_nstudies(dat=Dat4)	

length(which(Dat4$Nstudies<Dat4$Nstudies_median))
length(which(Dat4$Nstudies>=Dat4$Nstudies_median))
length(unique(Dat4$ID))

Dat4<-Dat4[which(Dat4$Nstudies>=Dat4$Nstudies_median),]

Dat4<-rbind(Dat5,Dat4)
nrow(Dat3)-nrow(Dat4) 
length(unique(Dat4$study))

