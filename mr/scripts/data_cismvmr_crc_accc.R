source("~/fatty-acids/colocalisation/scripts/regional_association_plots_functions.R")
source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/outcome_data/scripts/functions_combine_and_format_outcomes.R")

library(plyr)
library(gassocplot)
library(biomaRt)

data_list<-load_data(
	cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata",
	cancer3="~/fatty-acids/colocalisation/data/cancer_data_colocalisation_v2.RData",
	bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	# gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",
	# eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_imputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata",
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
	ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",#ref data for fads region. much faster to work with
	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  
	ref_turn_off=FALSE) #don't load ref_dat to save time, e.g. because already loaded

gene_tab1<-data.frame(data_list[1],stringsAsFactors=F)
ref1<-data.frame(data_list[2],stringsAsFactors=F)
fa.tab1<-data.frame(data_list[3],stringsAsFactors=F)
gtex_data1<-data.frame(data_list[4],stringsAsFactors=F)
eqtlgen_data1<-data.frame(data_list[5],stringsAsFactors=F)
bbj_eqtl_data1<-data.frame(data_list[6],stringsAsFactors=F)
gwis_file1<-unlist(data_list[7])
lun_data1<-data.frame(data_list[8],stringsAsFactors=F)
crc_data1<-data.frame(data_list[9],stringsAsFactors=F)
cancer_data1<-data.frame(data_list[10],stringsAsFactors=F)

cancer_data2<-cancer_data1[cancer_data1$trait=="Colorectal cancer" & cancer_data1$study=="ACCC",]

Trait<-c("GLA:LA lnD6Dpooled","AA:DGLA lnD5Dpooled")

# tissues<-unique(bbj_eqtl_data1$tissue)
# tissues<-tissues[!tissues %in%c("artery coronary","colon sigmoid","colon transverse")]

fads1<-format_data3(lun_data=NULL,crc_data=NULL,cancer_data=cancer_data2,gtex_data=NULL,fa.tab=fa.tab1[fa.tab1$trait %in% Trait,],ref=ref1,eqtlgen_data=NULL,gene_tab=gene_tab1,ld_eas_pop="eas",gene="FADS1",region=500000,gtex_tissues=NULL,gwis_file=gwis_file1,studies="Dorajoo/SCHS",fix_log10=FALSE,bbj_eqtl_data=bbj_eqtl_data1,pre_formatted_data=TRUE) #fix_log10 rescales the Z scores when -log10Pvalue is greater than 1000. Should only be needed for the regional association plots 
fads2<-format_data3(lun_data=NULL,crc_data=NULL,cancer_data=cancer_data2,gtex_data=NULL,fa.tab=fa.tab1[fa.tab1$trait %in% Trait,],ref=ref1,eqtlgen_data=NULL,gene_tab=gene_tab1,ld_eas_pop="eas",gene="FADS2",region=500000,gtex_tissues=NULL,gwis_file=gwis_file1,studies="Dorajoo/SCHS",fix_log10=FALSE,bbj_eqtl_data=bbj_eqtl_data1,pre_formatted_data=TRUE) 

ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt"
ref <- read.table(ref_dat,stringsAsFactors=F,head=F)


data_list<-prep_data_eas2(Dat=fads1) #prep_data_eas2 uses the cancer as the reference. SNPs that are not present in the cancer dataset get dropped. prep_data_eas uses AA:DGLA as the reference and therefore SNPs not present in that dataset are dropped. Since cancer is always part of the MR analysis makes sense to use as reference. On other hand, doesn't seem to make if any difference to final number of SNPs which one is used
Dat<-data.frame(data_list[1],stringsAsFactors=F)
ref_dat<-data.frame(data_list[2],stringsAsFactors=F)
unique(Dat$trait2)
unique(ref_dat$trait2)
# Dat2<-harmonise_dat(Dat=Dat,ref_dat=ref_dat,marker="marker",assume_same_strand_palindromic_and_non_palindromic_SNPs=TRUE,drop_palindromic_SNPs=FALSE)
Dat2<-harmonise_dat(Dat=Dat,ref_dat=ref_dat,marker="marker",assume_same_strand_palindromic_and_non_palindromic_SNPs=TRUE)
fads1_data<-format_results(gtex=FALSE,drop_palindromic_SNPs=TRUE) #we're mainly interested in BJ eQTL effect on cancer. since BJ eQTL data don't have palindromic SNPS, makes sense to just exclude them from all studies. 

# fads1_data2<-format_results(gtex=FALSE,drop_palindromic_SNPs=FALSE)
# alleles<-paste0(Dat2$effect_allele,Dat2$other_allele)
# study<-Dat2$study == "BJ"
# dim(Dat2[which(alleles %in% c("TA","AT","GC","CG") & study),])
# unique(Dat2$trait2[which(alleles %in% c("TA","AT","GC","CG") & study)])

data_list<-prep_data_eas2(Dat=fads2)
Dat<-data.frame(data_list[1],stringsAsFactors=F)
ref_dat<-data.frame(data_list[2],stringsAsFactors=F)
Dat2<-harmonise_dat(Dat=Dat,ref_dat=ref_dat,marker="marker",assume_same_strand_palindromic_and_non_palindromic_SNPs=TRUE)
fads2_data<-format_results(gtex=FALSE,drop_palindromic_SNPs=TRUE)

fads2_data<-fads2_data[grep("FADS2",fads2_data$trait),]
Dat<-rbind(fads1_data,fads2_data)
Dat[Dat$marker=="rs174546",]
table(Dat$trait)
unique(Dat$trait)
dim(Dat)
head(Dat)
Dat<-Dat[order(Dat$trait2),]


save(Dat,file="~/fatty-acids/mr/data/data_for_cismvmr_colorectal_accc.RData")
# load("~/fatty-acids/mr/data/data_for_cismvmr_lung_cancer.RData")



min(Dat$eaf[which(!is.na(Dat$eaf))])

# fads1_dat<-format_data(Dat=fads1)
# fads2_dat<-format_data(Dat=fads2)

# # fads1_dat[fads1_dat$SNP=="rs174546",]
# fads2_dat<-fads2_dat<-fads2_dat[fads2_dat$trait == "FADS2 expression",]
# fads_dat<-rbind(fads1_dat,fads2_dat)

# save(list=c("fads_dat","ld.matrix"),file="~/fatty-acids/mr/data/data_for_cismvmr.RData")
# load(file="~/fatty-acids/mr/data/data_for_cismvmr.RData")
