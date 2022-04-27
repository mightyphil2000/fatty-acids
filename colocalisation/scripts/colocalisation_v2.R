source("~/fatty-acids/colocalisation/scripts/regional_association_plots_functions.R")
# install.packages("devtools")
# library(devtools)
# install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
# library(devtools)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(plyr)
library(gassocplot)
library(biomaRt)

# BiocManager::install("snpStats")


# Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
# Attr<-listAttributes(Mart)

data_list<-load_data(
	cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata",
	cancer3="~/fatty-acids/colocalisation/data/cancer_data_colocalisation.RData",
	# bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	# gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",
	gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_FADS_ukbsnps_v2.Rdata",
	eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_imputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
	ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",#ref data for fads region. much faster to work with
	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  
	ref_turn_off=FALSE) #don't load ref_dat to save time, e.g. because already loaded


#######################
# East Asian datasets#
#######################

data_list3<-load_data(
	# cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt",
	cancer3="~/fatty-acids/colocalisation/data/escc.RData",
	bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	# gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",
	# eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata",
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
	ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",#ref data for fads region. much faster to work with
	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  #this only used for mapping chromosome positions to rsids
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

fa.tab2<-data.frame(data_list2[3],stringsAsFactors=F)
gwis_file2<-unlist(data_list2[7])

fa.tab1$trait<-paste(fa.tab1$trait,"notimputed")
fa.tab2$trait<-paste(fa.tab2$trait,"imputed")
fa.tab3<-rbind(fa.tab1,fa.tab2)

gene_tab1<-data.frame(data_list3[1],stringsAsFactors=F)
ref1<-data.frame(data_list3[2],stringsAsFactors=F)
fa.tab1<-data.frame(data_list3[3],stringsAsFactors=F)
gtex_data1<-data.frame(data_list3[4],stringsAsFactors=F)
eqtlgen_data1<-data.frame(data_list3[5],stringsAsFactors=F)
bbj_eqtl_data1<-data.frame(data_list3[6],stringsAsFactors=F)
gwis_file1<-unlist(data_list3[7])
lun_data1<-data.frame(data_list3[8],stringsAsFactors=F)
crc_data1<-data.frame(data_list3[9],stringsAsFactors=F)
cancer_data1<-data.frame(data_list3[10],stringsAsFactors=F)

gene<-"FADS1"
region<-500000
# tissues<-unique(bbj_eqtl_data1$tissue)
tissues<-unique(gtex_data1$tissue)
tissues<-tissues[!tissues %in%c("artery coronary","colon sigmoid","colon transverse")]
 
tissues.crc<-tissues[tissues %in% c("adipose subcutaneous","adipose visceral omentum","colon sigmoid","colon transverse","liver", "whole blood"   )]
tissues.skin<-tissues[tissues %in% c("adipose subcutaneous","adipose visceral omentum","skin not sun exposed suprapubic" ,"skin sun exposed lower leg",  "transverse","liver", "whole blood"   )]

Studies<-c("Dorajoo/SCHS")
# Studies<-c("CHARGE","Framingham","Shin")
# Trait<-"AA:DGLA"
# Trait<-"GLA:LA"
Trait<-"AA:DGLA lnD5Dpooled"
# Trait<-"GLA:LA lnD6Dpooled" 
# Traits<-c("GLA:LA","AA:DGLA")

# escc_data2<-escc_data1[escc_data1$trait2 ==  "Esophageal squamous cell carcinoma BJ/N-UGC" ,] 

cancer_data2<-cancer_data1[cancer_data1$study == "ILCCO/UKB",]

cancer_data2<-cancer_data1[cancer_data1$trait !="Esophageal squamous cell carcinoma",]

fa.tab3<-fa.tab3[fa.tab3$trait %in% c( "AA:DGLA notimputed","AA:DGLA imputed"),]
cancer_data2<-cancer_data1[cancer_data1$outcome=="Squamous cell lung cancer" , ]
unique(escc_data2$outcome)
fa.tab[fa.tab$SNP == "rs174546",c("trait","pval")]
bbj_eqtl_data2<-bbj_eqtl_data1[bbj_eqtl_data1$gene %in% c("ENSG00000197977","ENSG00000134824"),]
bbj_eqtl_data2[bbj_eqtl_data2$SNP == "rs174546",c("gene","p.value")]
escc_data2[escc_data2$marker == "rs174546",c("outcome","pval")]

cancer_data2<-cancer_data1[cancer_data1$trait =="Malignant skin cancer" ,]
cancer_data2<-cancer_data1[cancer_data1$trait =="Malignant non-melanoma skin cancer"  ,]
cancer_data2<-cancer_data1[cancer_data1$trait =="Basal cell carcinoma"  ,]
cancer_data2<-cancer_data1[cancer_data1$study == "BJ",]
# cancer_data2<-cancer_data1[cancer_data1$trait =="Respiratory and intrathoracic cancer" ,]
Dat<-format_data3(lun_data=NULL,crc_data=NULL,cancer_data=cancer_data2,gtex_data=NULL,fa.tab=fa.tab1[fa.tab1$trait %in% Trait,],ref=ref1,eqtlgen_data=NULL,gene_tab=gene_tab1,ld_eas_pop="eas",gene="FADS1",region=500000,gtex_tissues=NULL,gwis_file=gwis_file1,studies=Studies[1],fix_log10=FALSE,bbj_eqtl_data=bbj_eqtl_data1) #fix_log10 rescales the Z scores when -log10Pvalue is greater than 1000. Should only be needed for the regional association plots 

# ld.matrix<-data.frame(Dat[1])
# snps<-c("rs174546","rs2727271","rs2524299","rs174570")
# Pos<-which(rownames(ld.matrix) %in% snps)
# ld.matrix[Pos,Pos]



Binary.outcomes<-c(1,rep(0,length(Dat$Traits)-1))
# Binary.outcomes<-c(rep(1,8),0,0)
# Binary.outcomes<-rep(0,length(Dat$Traits))
# length(Binary.outcomes)
hyprcoloc_res<-hyprcoloc_function(Dat=Dat,gene=gene,gwis_file=gwis_file1,binary.outcomes=Binary.outcomes)
File_name<-make_file_name(Dir="~/fatty-acids/colocalisation/results/",Method="hyprcoloc",Info=Dat$Traits)

write.table(hyprcoloc_res,File_name,sep="\t",col.names=T,row.name=F,quote=F)

# coloc_abf_res<-coloc_abf_function(Traits=Dat$Traits,type_1="cc",type_2="quant",s=0.28)

Dat$Traits
Dat$cases
coloc_abf_res<-call_coloc_abf_res(Dat=Dat,trait1_cancer=FALSE)
coloc_abf_res<-do.call(rbind,coloc_abf_res)
File_name<-make_file_name(Dir="~/fatty-acids/colocalisation/results/",Method="coloc",Info=Dat$Traits)
write.table(coloc_abf_res,File_name,sep="\t",col.names=T,row.names=F,quote=F)

# Names<-coloc_res$traits_nice_names
# coloc_res$dropped_trait_nice_names
# Names<-Names[!is.na(Names)]
# strsplit(Names,split=",")
# Res<-moloc_test(df.ls,prior_var = c(0.15, 0.2), 
# 	priors = c(1e-04, 1e-06),
# 	save.SNP.info = T)
# Moloc <- moloc_test(df.ls)

	
# H0: neither trait has a genetic association in the region
# H1: only trait 1 has a genetic association in the region
# H2: only trait 2 has a genetic association in the region
# H3: both traits are associated, but with different causal variants
# H4: both traits are associated and share a single causal variant


	



info_fail<-function(){
	# SEnsitivity of results to summary data imputation errors in AA or DGLA
	snp_info_fail_files<-c("N6meta2031.tbl.fixed.tab_snps_fail_info_score.txt","
		N6meta2041.tbl.fixed.tab_snps_fail_info_score.txt")
	snps_info_fail1<-readLines(paste0("~/fatty-acids/colocalisation/data/",snp_info_fail_files[1]))
	snps_info_fail2<-readLines(paste0("~/fatty-acids/colocalisation/data/",snp_info_fail_files[2]))

	snps_info_fail1<-snps_info_fail1[!duplicated(snps_info_fail1)]
	snps_info_fail2<-snps_info_fail2[!duplicated(snps_info_fail2)]
	snps_info_fail<-c(snps_info_fail1,snps_info_fail2)
	snps_info_fail<-snps_info_fail[!duplicated(snps_info_fail)]

	fa.tab2<-fa.tab1[!fa.tab1$SNP %in% snps_info_fail,]
	fa.tab3<-fa.tab2[which(fa.tab2$eaf>0.05),]
}