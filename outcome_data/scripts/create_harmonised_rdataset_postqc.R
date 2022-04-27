source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

Dat<-collate_dat()
Dat<-basic_qc(dat=Dat)

nrow(Dat)
length(unique(Dat$study))
length(unique(Dat$ID))

# N<-length(which(Dat$MAC_case[Dat$ID == 133]<50 | Dat$MAC_control[Dat$ID == 133]<50))
# nrow(Dat)+N

dim(Dat[Dat$ID %in% c(106,5),])
# drop datasets Ambiguous effect allele column
Dat1<-Dat[!Dat$ID %in% c(106,5,74),]
nrow(Dat1)
length(unique(Dat1$study))
length(unique(Dat1$ID))

# drop datasets with info<0.8
Dat2<-prune_info80(Dat=Dat1)
nrow(Dat1)-nrow(Dat2)
length(unique(Dat2$study))

# drop log odds ratios > 1 or <-1 (includes asociations with unusually large Zb scores (>100))
# glioma133
Dat3<-Dat2[which(abs(Dat2$lnor)<1),]
nrow(Dat2)-nrow(Dat3)
length(unique(Dat3$study))

##################################################################
# drop SNPs that were contributed < median number of studies
##################################################################
Dat3$Nstudies<-NA
Dat3$Nstudies_median<-NA
Dat4<-Dat3[!is.na(Dat3$Direction),]
Dat5<-Dat3[is.na(Dat3$Direction),]

Dat4<-compute_nstudies(dat=Dat4)	

length(which(Dat4$Nstudies<Dat4$Nstudies_median))
length(which(Dat4$Nstudies>=Dat4$Nstudies_median))
length(unique(Dat4$ID))

Dat4<-Dat4[which(Dat4$Nstudies>=Dat4$Nstudies_median),]

Dat4<-rbind(Dat5,Dat4)
nrow(Dat3)-nrow(Dat4) 
length(unique(Dat4$study))

# set eaf to missing
N_eaf_to_na<-length(which(!is.na(Dat4$eaf[Dat4$ID %in% c(132,133,67,68)])))
Dat5<-Dat4
Dat5$eaf[Dat5$ID==132]<-NA
Dat5$eaf[Dat5$ID==133]<-NA
Dat5$eaf[Dat5$ID==67]<-NA
Dat5$eaf[Dat5$ID==68]<-NA

# datasets with effect allele column definitely wrong
N_967<-nrow(Dat4[Dat4$ID==967,])
##################################################
# total exclusions after applying the QC pipeline#
##################################################

Nassoc1<-nrow(Dat)
Nassoc2<-nrow(Dat4)
studies1<-length(unique(Dat$study))
studies2<-length(unique(Dat4$study))
Ndat1<-length(unique(Dat$ID)) 
Ndat2<-length(unique(Dat4$ID)) 

# total number of genetic associations with probable effect allele meta-data errors, low metrics of imputation quality (defined as info or imputation r2 score <0.8), high sample size conflict, unusually large log odds ratios or effect sizes that didn’t closely correspond to their reported P values. 

N<-Nassoc1-Nassoc2+N_eaf_to_na+N_967
N/nrow(Dat)*100

# total number of genetic associations with low metrics of imputation quality (defined as info or imputation r2 score <0.8), high sample size conflict, unusually large log odds ratios or effect sizes that didn’t closely correspond to their reported P values. 

N<-Nassoc1-Nassoc2

studies1-studies2
Ndat1-Ndat2

###################################################
# number that correspond to the fatty acid SNP set#
###################################################

snp1<-readLines("~/fatty-acids/outcome_data/data/snplist_Europeans_rsidsonly2.txt")
snp2<-readLines("~/fatty-acids/outcome_data/data/snplist_East_Asians_rsidsonly2_nodups.txt")

snps_fa<-unique(c(snp1,snp2))

Dat_fa1<-Dat4[which(Dat4$rsid %in% snps_fa),]

length(unique(Dat_fa1$ID))
N_fa1<-nrow(Dat_fa1)
length(unique(Dat_fa1$ID))
length(unique(Dat_fa1$study))

Dat_fa2<-Dat[which(Dat$rsid %in% snps_fa),]
length(unique(Dat_fa2$ID))
N_fa2<-nrow(Dat_fa2)

#Number of genetic associations with eaf set to NA that also corresponded to fatty acid SNPs, after the other exclusions.
N_fa_snps_eaf_to_na<-nrow(Dat4[!is.na(Dat4$eaf) & Dat4$rsid %in% snps_fa & Dat4$ID %in% c(132,133,67,68),])
# N_eaf_to_na<-length(which(!is.na(Dat3$eaf[Dat$ID %in% c(132,133,67,68)])))
N_fa_snps_967<-nrow(Dat4[Dat4$rsid %in% snps_fa & Dat4$ID == 967,])

# total number of genetic associations with errors or issues that corresponded to the fatty acid SNP set
N<-N_fa2-N_fa1+N_fa_snps_eaf_to_na+N_fa_snps_967
N/nrow(Dat_fa2)*100 #expected nrow(Dat_fa2) =223,970

# N associations in fatty acid SNP set after various exclusions
nrow(Dat_fa1)

###################################################
# number that correspond to the 1kg ref set#
###################################################

load("~/fatty-acids/outcome_data/data/refdat_1000G_superpops.Rdata")
snp1k<-unique(refdat_1000G_superpops$SNP)
Dat_1kg<-Dat4[which(Dat4$rsid %in% snp1k),]
length(unique(Dat_1kg$ID))
N_1k<-nrow(Dat_1kg)

###################################################
# number that correspond to the hits GWAS catalog
###################################################

Dat_gc<-Dat4[!Dat4$rsid %in% unique(c(snp1k,snps_fa)),]
Dat_gc<-Dat_gc[which(Dat_gc$all_summary_stats & !Dat_gc$rsid %in% unique(c(snp1k,snps_fa))),]
length(unique(Dat_gc$ID))
N_gc<-nrow(Dat_gc)

###################################################
# number that correspond to the FA regions
###################################################

Temp<-Dat[which(Dat$summary_set == "FAregions"),]
dim(Temp)
Temp1<-Temp[Temp$rsid %in% snps_fa,]
Temp2<-Temp[!Temp$rsid %in% snps_fa,]
unique(Temp$ID)
# number of SNPs in FAregions corresponding to the fatty acid SNP set (include these in the counts)
nrow(Temp1)
# number of SNPs in FAregions not corresponding to the fatty acid SNP set (exclude these from the counts)
nrow(Temp2)

Dat_fareg<-Dat4[Dat4$summary_set=="FAregions" & !Dat4$rsid %in% unique(c(snp1k,snps_fa)),]
length(unique(Dat_fareg$ID))
N_faregions<-nrow(Dat_fareg)


# studies with Nstudies less than median nstudies: c(59, 63,64,65,66,70,71)

# z<-Dat$lnor/Dat$se
# Dat[which(abs(z)>100),]

# Dat[which(Dat$rsid =="rs114935910"),]
# Dat[Dat$rsid %in% c("rs114935910","rs7068072","rs16993510","rs3848771") & Dat$ID ==133,]

save(Dat,file="~/fatty-acids/outcome_data/data/harmonised_data.Rdata")
save(Dat_fa1,file="~/fatty-acids/outcome_data/data/harmonised_data_postqc.Rdata")

