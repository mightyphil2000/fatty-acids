source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

Dat<-collate_dat(postqc=FALSE) 
Dat<-basic_qc(dat=Dat)

snp1<-readLines("~/fatty-acids/outcome_data/data/snplist_Europeans_rsidsonly2.txt")
snp2<-readLines("~/fatty-acids/outcome_data/data/snplist_East_Asians_rsidsonly2_nodups.txt")

snps_fa<-unique(c(snp1,snp2))
Dat_fa<-Dat[which(Dat$rsid %in% snps_fa),]
length(unique(Dat_fa$ID))
length(unique(Dat_fa$study))
N_fa<-nrow(Dat_fa)

Dat<-Dat_fa

nrow(Dat)
length(unique(Dat$study))

##################################################
# drop datasets Ambiguous effect allele column#
##################################################
Dat1<-Dat[!Dat$ID %in% c(106,5,74),]
nrow(Dat1)
length(unique(Dat1$study))
length(unique(Dat1$ID))

##############################
# drop datasets with info<0.8#
##############################

Dat2<-prune_info80(Dat=Dat)
N<-nrow(Dat)-nrow(Dat2)
N/nrow(Dat)*100
length(unique(Dat2$study))

# drop log odds ratios > 1 or <-1 (includes asociations with unusually large Zb scores (>100))
# glioma133
Dat3<-Dat[which(abs(Dat$lnor)<1),]
N<-nrow(Dat)-nrow(Dat3)
N/nrow(Dat)*100
length(unique(Dat3$study))

which(abs(Dat3$z)>98)

# set eaf to missing
Dat3$eaf[Dat3$ID==132]<-NA
Dat3$eaf[Dat3$ID==133]<-NA
Dat3$eaf[Dat3$ID==67]<-NA
Dat3$eaf[Dat3$ID==68]<-NA

#################################################################
# drop SNPs that were contributed < median number of studies
#################################################################

Dat$Nstudies<-NA
Dat$Nstudies_median<-NA
# Dat3$Nstudies<-NA
# Dat3$Nstudies_median<-NA
Dat4<-Dat[!is.na(Dat$Direction),]
Dat5<-Dat[is.na(Dat$Direction),]

Dat4<-compute_nstudies(dat=Dat4)	
length(unique(Dat4$ID[which(Dat4$Nstudies<Dat4$Nstudies_median)]))
length(which(Dat4$Nstudies<Dat4$Nstudies_median))
length(which(Dat4$Nstudies>=Dat4$Nstudies_median))
length(unique(Dat4$ID))

Dat4<-Dat4[which(Dat4$Nstudies>=Dat4$Nstudies_median),]
Dat4<-rbind(Dat5,Dat4)


nrow(Dat4)
length(unique(Dat4$study))
length(unique(Dat4$ID))

##################################################
# total exclusions after applying the QC pipeline#
##################################################

Nassoc1<-nrow(Dat)
Nassoc2<-nrow(Dat4)
studies1<-length(unique(Dat$study))
studies2<-length(unique(Dat4$study))
Ndat1<-length(unique(Dat$ID)) 
Ndat2<-length(unique(Dat4$ID)) 

Nassoc1-Nassoc2
studies1-studies2
Ndat1-Ndat2
