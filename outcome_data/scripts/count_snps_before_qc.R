source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

Dat<-collate_dat(postqc=FALSE)

# min(as.numeric(Dat$info),na.rm=TRUE)
# min(as.numeric(Dat$info),na.rm=TRUE)
# min(Dat$info,na.rm=TRUE)
length(unique(Dat$study))
length(unique(Dat$ID[Dat$all_summary_stats]))
length(unique(Dat$ID[!Dat$all_summary_stats]))

#N studies in Open GWAS
Dat1<-Dat[which(Dat$open_gwas),]
length(unique(Dat1$ID[Dat1$study == "FinnGen"]))+
length(unique(Dat1$ID[Dat1$study == "BJ"]))+
length(unique(Dat1$ID[Dat1$study == "UKB" ]))+
length(unique(Dat1$ID[Dat1$consortium]))

# N datasets with eaf, info
length(unique(Dat$ID[!is.na(Dat$eaf)]))
length(unique(Dat$ID[!is.na(Dat$info) | !is.na(Dat$info1) | !is.na(Dat$info2) | !is.na(Dat$info3)]))
length(unique(Dat$ID[!is.na(Dat$HWEp)]))
length(unique(Dat$ID[!is.na(Dat$phet)]))
length(unique(Dat$ID[!is.na(Dat$Direction)]))
length(unique(Dat$ID))
nrow(Dat[Dat$ID == 106,])+nrow(Dat[Dat$ID == 5,])

names(Dat)
# number of GWAS consortia in Open GWAS

length(unique(Dat$ID[which(Dat$correspondence)]))
length(unique(Dat$ID[which(Dat$gwas_catalog)]))

unique(Dat$ID[which(Dat$correspondence & Dat$open_gwas)])

unique(Dat$correspondence[which(Dat$open_gwas)])

length(unique(Dat$study[which(!Dat$open_gwas)]))

Dat<-basic_qc(dat=Dat)
nrow(Dat)
length(unique(Dat$study))
length(unique(Dat$ID))
# count number of fatty acid SNPs & excluding faregion SNPs 
# number that correspond to the fatty acid SNP set
snp1<-readLines("~/fatty-acids/outcome_data/data/snplist_Europeans_rsidsonly2.txt")
snp2<-readLines("~/fatty-acids/outcome_data/data/snplist_East_Asians_rsidsonly2_nodups.txt")

snps_fa<-unique(c(snp1,snp2))
Dat2<-Dat[which(Dat$rsid %in% snps_fa),]
nrow(Dat2)
Dat2_1<-Dat2[Dat2$summary_set !="FAregions",]
Dat2_2<-Dat2[Dat2$summary_set =="FAregions",]
Dat2_2<-Dat2_2[Dat2_2$rsid %in% snps_fa,]
Dat2<-rbind(Dat2_1,Dat2_2)
dim(Dat2)
dim(Dat2_1)
dim(Dat2_2)
length(unique(Dat2$study))
length(unique(Dat2$ID))

# N SNPs that correspnd to the 1000 genomes reference set
load("~/fatty-acids/outcome_data/data/refdat_1000G_superpops.Rdata")
snp1k<-unique(refdat_1000G_superpops$SNP)
Dat3<-Dat[which(Dat$rsid %in% snp1k),]
# Temp<-Dat[which(!Dat$rsid %in% snp1k),]
length(unique(Dat3$ID))
nrow(Dat3)
length(unique(Dat3$study))
length(unique(Dat3$ID))

# N SNPs that correspond to known associations in the GWAS catalog
# N_gc<-N_tot-N_fa-N_1k
# Dat4<-Dat[!Dat$rsid %in% unique(c(snp1k,snps_fa)),]
Dat4<-Dat[which(Dat$all_summary_stats & !Dat$rsid %in% unique(c(snp1k,snps_fa))),]
length(unique(Dat4$ID))
N_gc<-nrow(Dat4)
N_tot
N_fa
N_1k
N_gc

# N SNPs corresponding to the FAregions that do or don't overlap with fatty acid SNP set
Dat5<-Dat[Dat$summary_set == "FAregions",]
Dat5_1<-Dat5[!Dat5$rsid %in% snps_fa,]
Dat5_2<-Dat5[Dat5$rsid %in% snps_fa,]
nrow(Dat5_1)
nrow(Dat5_2)
unique(Dat5_1$ID)
unique(Dat5_2$ID)


# number of meta analyses
Dat$meta_analysis<-NA
Dat$meta_analysis[!is.na(Dat$Direction)]<-TRUE
Dat$meta_analysis[Dat$study %in% c("23NMSC","BJ","FinnGen","UKB","ESS","BC-NHL","C-ALL","HKHC","KCML","KHBC","MCCS","MMS","MNC","MPM","NBCS","NBS","SCCS","SJ-COG","TCS","TNC","UCSF_MAYO","UMS")]<-FALSE
Dat$meta_analysis[Dat$study %in% c("BCAC","ECAC","EPITHYR","GAME-ON","GliomaScan","ILCCO","INHANCE","InterLymph","KidRISK","OCAC","OCAC (EAS)","PanC4","PanScan I","PanScan I+II","PanScan I+II+PanC4","PanScan III","PRACTICAL","UCSF_AGS + SFAGS","BC-ALL","CHC","HLS","MENC","MMAC","NB-UGC","N-UGC")]<-TRUE
length(unique(Dat$study[Dat$meta_analysis]))
length(unique(Dat$ID[Dat$meta_analysis]))
length(unique(Dat$study[!Dat$meta_analysis]))
length(unique(Dat$ID[!Dat$meta_analysis]))



# number with full gwas summary statistics

nrow(Dat[Dat$all_summary_stats,])
length(unique(Dat$ID[Dat$all_summary_stats]))
length(unique(Dat$study[Dat$all_summary_stats]))

nrow(Dat[!Dat$all_summary_stats,])
length(unique(Dat$ID[!Dat$all_summary_stats]))
length(unique(Dat$study[!Dat$all_summary_stats]))
length(unique(Dat$ID))
N_tot<-nrow(Dat)

