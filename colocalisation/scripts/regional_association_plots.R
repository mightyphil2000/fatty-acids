# rm(list=ls())
source("~/fatty-acids/colocalisation/scripts/regional_association_plots_functions.R")

# imputed fatty acid summary data
data_list<-load_data(
	cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata",
	# cancer3="~/fatty-acids/colocalisation/data/ilcco.RData",
	# cancer3="~/fatty-acids/colocalisation/data/bcc.RData",
	cancer3="~/fatty-acids/colocalisation/data/cancer_data_colocalisation.RData",
	# cancer3="~/fatty-acids/colocalisation/data/luca.RData",
	# bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_FADS_ukbsnps_v2.Rdata",
	eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_imputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  
	ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",
	ref_turn_off=FALSE) #don't load ref_dat to save time, e.g. because already loaded



data_list2<-load_data(
	cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata",
	cancer3="~/fatty-acids/colocalisation/data/luca.RData",
	# bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_FADS_ukbsnps_v2.Rdata",
	eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
	# ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",
	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  
	ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",
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


# # fa.tab2$samplesize[fa.tab2$file %in% c("POA_to_PA.tab","OA_to_SA.tab")]<-3521
# fa.tab2$z<-fa.tab2$beta/fa.tab2$se
# fa.tab2$b_sd<-b_sd(z=fa.tab2$z,maf = fa.tab2$eaf,n=fa.tab2$samplesize)
# fa.tab2$se_sd<-fa.tab2$b_sd/fa.tab2$z
# fa.tab2$r2<-2*fa.tab2$b_sd^2*fa.tab2$eaf*(1-fa.tab2$eaf)


# fa.tab2[fa.tab2$file %in% c("AA_to_DGLA.tab.gz","DHA_to_DPA_n3.tab.gz","POA_to_PA.tab.gz"),c("study","SNP","file","b_sd","se_sd","r2","eaf","pval","effect_allele","other_allele")]


gene_tab1<-data.frame(data_list2[1],stringsAsFactors=F)
ref1<-data.frame(data_list2[2],stringsAsFactors=F)
fa.tab1<-data.frame(data_list2[3],stringsAsFactors=F)
gtex_data1<-data.frame(data_list2[4],stringsAsFactors=F)
eqtlgen_data1<-data.frame(data_list2[5],stringsAsFactors=F)
bbj_eqtl_data1<-data.frame(data_list2[6],stringsAsFactors=F)
gwis_file1<-unlist(data_list2[7])
lun_data1<-data.frame(data_list2[8],stringsAsFactors=F)
crc_data1<-data.frame(data_list2[9],stringsAsFactors=F)
cancer_data1<-data.frame(data_list2[10],stringsAsFactors=F)


#######################
# East Asian datasets#
#######################

data_list3<-load_data(
	# cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt",
	cancer3="~/fatty-acids/colocalisation/data/cancer_data_colocalisation.RData",
	bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	# gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",
	# eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata",
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
	ref_dat="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",
	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  #this only used for mapping chromosome positions to rsids
	ref_turn_off=FALSE) #don't load ref_dat to save time, e.g. because already loaded

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

# fa.tab2<-fa.tab1[fa.tab1$SNP %in%  c("rs603424","rs174546","rs3734398" ),]
#
# fa.tab2$samplesize[fa.tab2$file %in% c("POA_to_PA.tab","OA_to_SA.tab")]<-3521
# fa.tab2$z<-fa.tab2$beta/fa.tab2$se
# fa.tab2$b_sd<-b_sd(z=fa.tab2$z,maf = fa.tab2$eaf,n=fa.tab2$samplesize)
# fa.tab2$se_sd<-fa.tab2$b_sd/fa.tab2$z
# fa.tab2$r2<-2*fa.tab2$b_sd^2*fa.tab2$eaf*(1-fa.tab2$eaf)

# fa.tab2[fa.tab2$file %in% c("GLA_to_LA_pooled.tab.gz","score_lnD6D_pooled_allchr_qc1.tab","POA_to_PA.tab"),c("SNP","file","b_sd","se_sd","r2","eaf","pval","effect_allele","other_allele")]


# rs603424
# fa.tab2<-fa.tab[fa.tab$file %in% c("AA_to_DGLA.tab", "GLA_to_LA.tab","AA_to_DGLA_adjSNP.tab", "GLA_to_LA_adjSNP.tab"),]

# length(which(bbj_eqtl_data1$gene == "ENSG00000075826"))

####################
# plots East Asians#
####################
Pos1<-grep("AA_to_DGLA",fa.tab1$file)
Pos2<-grep("score_lnD5D_pooled_allchr_qc1.tab",fa.tab1$file)
Pos<-c(Pos1,Pos2)
Pos<-unique(Pos)
Dat_fa.tab<-fa.tab1[Pos2,] 
unique(Dat_fa.tab$file)

Pos1<-grep("GLA_to_LA",fa.tab1$file)
Pos2<-grep("score_lnD6D_pooled_allchr_qc1.tab",fa.tab1$file)
Pos<-c(Pos1,Pos2)
Pos<-unique(Pos)
Dat_fa.tab<-fa.tab1[Pos2,]
# any(Dat_fa.tab$SNP %in% c("rs968567","rs61896141"))
Pos1<-grep("OA_to_SA",fa.tab1$file)
Dat_fa.tab<-fa.tab1[Pos1,]

Pos1<-grep("POA_to_PA",fa.tab1$file)
Dat_fa.tab<-fa.tab1[Pos1,]

Pos1<-grep("DGLA_to_GLA",fa.tab1$file)
Dat_fa.tab<-fa.tab1[Pos1,]

Pos1<-grep("DHA_to_EPA",fa.tab1$file)
Dat_fa.tab<-fa.tab1[Pos1,]

unique(fa.tab1[,c("file","study")])



i<-7
Dat_fa.tab[Dat_fa.tab$SNP=="rs603424",]
Studies<-c("Dorajoo/SCHS","ChargeEA")
for(i in 1:length(gene_tab1$gene)){
	print(i)
	regional_association_plot(gene=gene_tab1$gene[1],study=Studies[2],plot_strategy="alltraits_1study_alltissues",gene_tab=gene_tab1,ref=ref1,fa.tab=Dat_fa.tab,bbj_eqtl_data=bbj_eqtl_data1,gtex_data=NULL,eqtlgen_data=NULL,region=500000,gwis_file=gwis_file1,ld_eas_pop="eas",Top.marker="rs6584379")
}
# ,Plot_height=500
# ,Top.marker="rs603424"
# Top.marker="rs603424"
head(bbj_eqtl_data1)


####################
# plots East Asians#
####################
# highlight SNPs   #
####################
Studies<-c("Dorajoo/SCHS","ChargeEA")
Top.markers<-c("rs174584","rs174546","rs174572","rs174570","rs174575","rs603424","rs6584379","rs5792235")
gene<-"FADS1"
Top.markers<-"rs174546"

Trait<-"AA:DGLA lnD5Dpooled"
# Trait<-"GLA:LA lnD6Dpooled" 
# Traits<-"AA:DGLA"
# Traits<-c("GLA:LA","AA:DGLA")

# create duplicate blood panels
# bbj_eqtl_data2<-bbj_eqtl_data1[bbj_eqtl_data1$tissue == "Blood",]
# bbj_eqtl_data3<-bbj_eqtl_data2
# bbj_eqtl_data4<-bbj_eqtl_data2
# bbj_eqtl_data2$tissue<-paste0(bbj_eqtl_data2$tissue,"1")
# bbj_eqtl_data3$tissue<-paste0(bbj_eqtl_data3$tissue,"2")
# bbj_eqtl_data4$tissue<-paste0(bbj_eqtl_data4$tissue,"3")
# bbj_eqtl_data<-do.call(rbind,list(bbj_eqtl_data2,bbj_eqtl_data3,bbj_eqtl_data4))

cancer_data2<-cancer_data1[cancer_data1$study %in% c("BJ","BJ/N-UGC"),]

regional_association_plot(gene=gene,study=Studies[1],gene_tab=gene_tab1,ref=ref1,fa.tab=fa.tab1[fa.tab1$trait == Trait,],bbj_eqtl_data=bbj_eqtl_data1,gtex_data=NULL,eqtlgen_data=NULL,gwis_file=gwis_file1,region=500000,Top.marker=Top.markers,ld_eas_pop="eas",crc_data=NULL,cancer_data=NULL,gtex_tissues=NULL,fix_charge=FALSE,fix_log10=FALSE)

crc_data1[crc_data1$marker %in% c("rs174546","rs968567","rs61896141"),]


# rs6584379 top variant for POA:PA in East Asians
# rs5792235 causal variant for FADS1/FADS2 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124650/
# rs174575 causal variant for CRC/ FADS2 in https://cancerres.aacrjournals.org/content/early/2020/03/03/0008-5472.CAN-19-2389.short?casa_token=iiOiLi5UtFAAAAAA:yWqfRjn60Pq-cle2jviqq3m4fgCaMk8cd0TUrE7S7X5WrXjHc-kDETYmT72LzUnDKCLFR8ombfGPgsIE
# rs174584 this was top SNP in GLA:LA cases and FADS2 expression in CD8 T cells and has r2 ≥0.813 with all other top SNPs in regional association plots for FADS1, FADS2 and GLA:LA in Japanese, (r2 ≥ 0.958 in Han Chinese; r2 ≥0.907 in Southern Han Chinese; r2≥0.834 in Dai CHinese; r2 ≥0.924 in combined Japanese, Han Beinijing, Han Southern Chinese and Dai Chinese)
# rs174546 is selected FADS1 SNP in Europeans
# rs174572 is top SNP in AA:DGLA controls Dorajoo but SNP weakly correlated with other top SNPs (r2≤0.165)
# rs174546 r2 with rs174584: 1.0 CHS, 1.0 CHB, 0.94 JPT, 0.97 Chinese DAI (CDX), 0.97 KHV (vietnam)

###################
# plots Europeans#
###################
Traits1<-c("AA_to_DGLA")
Traits2<-c("GLA_to_LA")
Traits3<-c("DHA_to_DPA_n3")
# Traits4<-c("OA_to_SA")
Traits5<-c("POA_to_PA")
# Traits_list<-c(Traits1,Traits2,Traits3,Traits4,Traits5)
Traits_list<-c(Traits3)

for(k in 1:length(Traits_list)){
	print(k)
	Traits<-"POA_to_PA"
	print(Traits)
	Pos<-unlist(lapply(1:length(Traits),FUN=function(x)
	grep(Traits[x],fa.tab1$file)))
	Dat_fa.tab<-fa.tab1[Pos,]
	unique(Dat_fa.tab$file)
	# Dat_fa.tab<-Dat_fa.tab[Dat_fa.tab$SNP != "rs603424",]
	# grep("GLA_to_LA",unique(fa.tab1$file))
	# grep(Traits[x],unique(fa.tab1$file))
	

	# j<-3
	for(i in 1:length(gene_tab1$gene)){
		for(j in 1:length(Studies)){
			print(i)
			print(j)
			# plot_strategy="alltraits_1study_alltissues"
			regional_association_plot(gene=gene_tab1$gene[4],study=Studies[1],gtex_tissues=c("whole blood","adipose subcutaneous","liver","adipose visceral omentum"),fix_charge=FALSE,gene_tab=gene_tab1,ref=ref1,fa.tab=Dat_fa.tab,gtex_data=gtex_data1,eqtlgen_data=eqtlgen_data1,bbj_eqtl_data=NULL,gwis_file=gwis_file1,region=500000,Top.marker="rs603424")
			# rs3734398
			# ,Plot_height=500
			# rs6584379
			# 
			# eqtlgen_data1
			# gtex_data1
			# 
		}
	}
}
###################
# plots Cancer#
###################

# Dat_fa.tab[Dat_fa.tab$SNP=="rs2524299",]

Trait<-c("AA:DGLA")
# Trait<-c("GLA:LA")
genes<-c("FADS1","FADS2")
Traits<-c("AA:DGLA","GLA:LA")
gene<-genes[which(Traits==Trait)]
# gene<-"FADS2"
# Trait<-c("GLA:LA")

# Traits<-c("POA_to_PA")
# Pos<-unlist(lapply(1:length(Traits),FUN=function(x)
# 	grep(Traits[x],fa.tab1$file)))
# Dat_fa.tab<-fa.tab1[Pos,]
# unique(Dat_fa.tab$file)



# which(bbj_eqtl_data1$SNP == "rs61896141")
# ref1[ref1$V2=="rs174546",]
# for(study in c("CHARGE","Framingham","Shin")){
	# print(study)

Top.markers<-c("rs174546","rs2727271","rs2524299","rs968567","rs61896141","rs5792235","rs174575","rs603424","rs10509744","rs10883508","rs4244338","rs11190588","rs683854","rs3853519","rs55825790","rs1054411")	
Studies<-c("CHARGE","Framingham","Shin")

Tissues<-c("colon sigmoid","colon transverse")

Tissues<-Tissues[!Tissues %in% c("artery coronary","colon sigmoid","colon transverse","lung")]
Tissues2<-Tissues[!Tissues %in% c("artery coronary")]

# for(i in 1:length(Top.markers)){
# for(i in 1:5){
		# for(k in 8:18){
			# head(Dat_fa.tab)
			# plot_strategy="alltraits_1study_alltissues"
# crc_data.eas1<-data.frame(data_list3[9],stringsAsFactors=F)
# fa.tab.eas1<-data.frame(data_list3[3],stringsAsFactors=F)
# gwis_file.eas1<-unlist(data_list3[7])


# snp_info_fail_files<-c("N6meta2031.tbl.fixed.tab_snps_fail_info_score.txt","	N6meta2041.tbl.fixed.tab_snps_fail_info_score.txt")
# snps_info_fail1<-readLines(paste0("~/fatty-acids/colocalisation/data/",snp_info_fail_files[1]))
# snps_info_fail2<-readLines(paste0("~/fatty-acids/colocalisation/data/",snp_info_fail_files[2]))

# snps_info_fail1<-snps_info_fail1[!duplicated(snps_info_fail1)]
# snps_info_fail2<-snps_info_fail2[!duplicated(snps_info_fail2)]
# snps_info_fail<-c(snps_info_fail1,snps_info_fail2)
# snps_info_fail<-snps_info_fail[!duplicated(snps_info_fail)]
# fa.tab2<-fa.tab1[!fa.tab1$SNP %in% snps_info_fail,]
# fa.tab3<-fa.tab2[which(fa.tab2$eaf>0.05),]

Studies<-c("CHARGE","Framingham","Shin")
Tissues<-unique(gtex_data1$tissue)
Tissues<-Tissues[Tissues %in% c("adipose subcutaneous","adipose visceral omentum","liver","whole blood")]
Tissues<-Tissues[!Tissues %in%c("esophagus mucosa", "adipose visceral omentum"   ,"liver","whole blood"   ,"adipose subcutaneous","skin not sun exposed suprapubic" ,"skin sun exposed lower leg","artery coronary","colon transverse","esophagus muscularis")]
Tissues<-Tissues[!Tissues %in% c("adipose subcutaneous","adipose visceral omentum","liver","whole blood","artery coronary")]

Tissues<-c("adipose visceral omentum","liver","whole blood","adipose subcutaneous","skin not sun exposed suprapubic" ,"skin sun exposed lower leg")
# cancer_data2<-cancer_data1[cancer_data1$trait %in% c("Esophageal adenocarcinoma", "Lung cancer"),]
# cancer_data2<-cancer_data2[cancer_data2$study !="BJ",]
Tissues<-c("colon transverse","colon sigmoid")
cancer_data2<-cancer_data1[cancer_data1$trait %in% c("Malignant skin cancer","Malignant non-melanoma skin cancer","Basal cell carcinoma"),]

regional_association_plot(gene="FADS1",study=Studies[1],gtex_tissues=NULL,fix_charge=FALSE,gene_tab=gene_tab1,ref=ref1,fa.tab=fa.tab1[fa.tab1$trait %in% c("AA:DGLA","GLA:LA"),],gtex_data=NULL,eqtlgen_data=NULL,bbj_eqtl_data=NULL,gwis_file=gwis_file1,lun_data=NULL,crc_data=NULL,cancer_data=cancer_data2,region=500000,Top.marker="rs174546",fix_log10=FALSE,ld_eas_pop=NULL)
	# }

	# rs174570
	# rs2727263
		# "colon sigmoid","colon transverse"
		# lun_data1
		# ,"adipose subcutaneous","liver",
		# "adipose visceral omentum"
		# "whole blood",
		# eqtlgen_data1
	}

	regional_association_plot(gene=gene_tab1$gene[2],study=Studies[1],plot_strategy="alltraits_1study_alltissues",gtex_tissues=c("whole blood","adipose subcutaneous","liver","adipose visceral omentum"),fix_charge=FALSE,gene_tab=gene_tab1,ref=ref1,fa.tab=Dat_fa.tab,gtex_data=gtex_data1,eqtlgen_data=eqtlgen_data1,bbj_eqtl_data=NULL,gwis_file=gwis_file1,region=500000,Top.marker="rs174546")
	


# }

# RS_number	rs61896141	rs174546	rs968567	rs2727271	rs2524299
# rs61896141	1.0	0.38	0.981	0.004	0.005
# rs174546	0.38	1.0	0.381	0.265	0.235
# rs968567	0.981	0.381	1.0	0.004	0.005
# rs2727271	0.004	0.265	0.004	1.0	0.929
# rs2524299	0.005	0.235	0.005	0.929	1.0

)