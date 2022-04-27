ieugwasr::get_access_token()
library(TwoSampleMR)
library(ieugwasr)
ao<-available_outcomes()

Token<-ieugwasr::check_access_token()
 
# datasets from mrbase / OpenGWAS
Cancer_mrbase<-get_mrbase_cancer_data()
save(Cancer_mrbase,file="cancers_mrbase.Rdata")
write.table(Cancer_mrbase,"cancers_mrbase.txt",sep="\t",col.names=T,row.names=F,quote=F)



ao<-gwasinfo(access_token = Token)

get_mrbase_cancer_data<-function(){

	

	# Breast cancer - harmonised already via MR-Base
	# Glioma, gliomascan - harmonised already via MR-Base
	# Lung cancer - harmonised already via MR-Base
	# Ovarian cancer - harmonised already via MR-Base
	# Prostate cancer - harmonised already via MR-Base
	# Neuroblastoma - harmonised already via MR-Base (another version available obtained by correspondence)

	# "Neuroblastoma"
	Cancers<-c("Breast cancer","Lung cancer","Ovarian cancer")
# "Prostate cancer"
	
	Pos<-unlist(lapply(Cancers,FUN=function(i) grep(i,ao$trait)))
	ID.tricl<-ao$id[grep("TRICL",ao$consortium)]
	# ao[ao$id %in% ID.tricl,c("trait","ncase")]
	ID.ocac<-ao$id[grep("OCAC",ao$consortium)]
	IDS<-ao$id[Pos]
	ID.tricl<-ID.tricl[!ID.tricl %in% IDS]
	ID.ocac<-ID.ocac[!ID.ocac %in% IDS]
	
	ID.other<-c(ao$id[ao$trait=="Neuroblastoma"],ao$id[ao$trait=="Thyroid cancer"],ao$id[ao$trait=="Upper gastrointestinal cancers"],ao$id[ao$trait=="Gallbladder cancer"])

	# ID.other<-c(ao$id[ao$trait=="Glioma" & ao$filename=="for_Philip_data_delivery.txt.new.tab"],ao$id[ao$trait=="Neuroblastoma"],ao$id[ao$trait=="Upper gastrointestinal cancers"],ao$id[ao$trait=="Thyroid cancer"],ao$id[ao$trait=="Gallbladder cancer"])

	snplist<-read.table("~/fatty-acids-mr/instruments/snplist_Europeans.txt",sep="\t",head=T,stringsAsFactors=F)
	# snplist2<-read.table("~/fatty-acids-mr/instruments/snplist_East_Asians.txt",sep="\t",head=T,stringsAsFactors=F)
	snps<-unique(snplist$SNP)
	# IDS2<-ao$id[ao$subcategory=="Cancer"]
	# c(IDS,IDS2)

	IDS1<-IDS[1:10]
	IDS2<-IDS[11:20]
	IDS3<-IDS[21:25]
	IDS4<-IDS[26:35]
	IDS5<-IDS[36:44]
	length(IDS)
	cancer_dat <- extract_outcome_data(
	    snps = snps,
	    outcomes = IDS1)
	cancer_dat2 <- extract_outcome_data(
	    snps = snps,
	    outcomes = IDS2)
	cancer_dat3 <- extract_outcome_data(
	    snps = snps,
	    outcomes = IDS3)
	cancer_dat4 <- extract_outcome_data(
	    snps = snps,
	    outcomes = IDS4)
	cancer_dat5 <- extract_outcome_data(
	    snps = snps,
	    outcomes = IDS5)
	cancer_dat6 <- extract_outcome_data(
	    snps = snps,
	    outcomes = ID.tricl)
	cancer_dat7 <- extract_outcome_data(
	    snps = snps,
	    outcomes = ID.ocac)
	cancer_dat8 <- extract_outcome_data(
	    snps = snps,
	    outcomes = ID.other)

	IDS.all<-c(IDS,ID.tricl,ID.ocac,ID.other)

	# names(cancer_dat8)[names(cancer_dat8) == "SNP"]<-"rsid"
	# names(cancer_dat8)[names(cancer_dat8) == "beta.outcome"]<-"lnor"
	# names(cancer_dat8)[names(cancer_dat8) == "se.outcome"]<-"se"
	# names(cancer_dat8)[names(cancer_dat8) == "ncase.outcome"]<-"ncase"
	# names(cancer_dat8)[names(cancer_dat8) == "ncontrol.outcome"]<-"ncontrol"
	# names(cancer_dat8)[names(cancer_dat8) == "pval.outcome"]<-"p"
	# names(cancer_dat8)[names(cancer_dat8) == "eaf.outcome"]<-"eaf"
	# names(cancer_dat8)[names(cancer_dat8) == "effect_allele.outcome"]<-"Effect.Allele"
	# names(cancer_dat8)[names(cancer_dat8) == "other_allele.outcome"]<-"Other.Allele"
	# names(cancer_dat8)[names(cancer_dat8) == "consortium.outcome"]<-"study"
	# names(cancer_dat8)[names(cancer_dat8) == "pmid.outcome"]<-"pmid"
	# names(cancer_dat8)[names(cancer_dat8) == "originalname.outcome"]<-"outcome"
	# names(cancer_dat8)<-gsub(".outcome","",names(cancer_dat8))
	# cancer_dat8$effect_allele_confirmed<-TRUE
	# load("cancers_mrbase.Rdata")
	# Cancer_mrbase<-rbind(Cancer_mrbase,cancer_dat8)

	# names(cancer_dat8)[!names(cancer_dat8) %in% names(Cancer_mrbase)]
	# ao[ao$id %in% IDS.all,c("trait","ncase","ncontrol")]
	# ao[ao$id %in% IDS,c("trait","consortium","ncase","ncontrol","id")]

	Cancer_mrbase<-do.call(rbind,list(cancer_dat,cancer_dat2,cancer_dat3,cancer_dat4,cancer_dat5,cancer_dat6,cancer_dat7))

	names(Cancer_mrbase)[names(Cancer_mrbase) == "SNP"]<-"rsid"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "beta.outcome"]<-"lnor"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "se.outcome"]<-"se"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "ncase.outcome"]<-"ncase"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "ncontrol.outcome"]<-"ncontrol"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "pval.outcome"]<-"p"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "eaf.outcome"]<-"eaf"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "effect_allele.outcome"]<-"Effect.Allele"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "other_allele.outcome"]<-"Other.Allele"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "consortium.outcome"]<-"study"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "pmid.outcome"]<-"pmid"
	names(Cancer_mrbase)[names(Cancer_mrbase) == "originalname.outcome"]<-"outcome"
	names(Cancer_mrbase)<-gsub(".outcome","",names(Cancer_mrbase))
	Cancer_mrbase$effect_allele_confirmed<-TRUE
	# names(Cancer_mrbase)[names(Cancer_mrbase) == "proxy.outcome"]<-"proxy"
	# names(Cancer_mrbase)[names(Cancer_mrbase) == "target_snp.outcome"]<-"target_snp"
	return(Cancer_mrbase)
}





# all files stored here on rdsf if want to do direct lookups
# cd /projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc

# TRICL_ICR_MDACC_IARC_NIH_GWAMA_Adeno.rsid.txt.tab
# TRICL_ICR_MDACC_IARC_NIH_GWAMA_overall.rsid.txt.tab
# TRICL_ICR_MDACC_IARC_NIH_GWAMA_Squamous.rsid.txt.tab
# TRICL_Meta_032114_Adeno.csv.tab
# TRICL_Meta_032114_Overall.csv.tab
# TRICL_Meta_032114_Squamous.csv.tab

# Ovarian_clearcell.txt.tab
# Ovarian_endometrioid.txt.tab
# Ovarian_mucinous.txt.tab
# Ovarian_overall.txt.tab
# Ovarian_serous_hg.txt.tab
# Ovarian_serouslowgrade.txt.tab

# oncoarray_bcac_public_release_oct17_GWAS_ER_negative_disease.txt.tab
# oncoarray_bcac_public_release_oct17_GWAS_ER_positive_disease.txt.tab
# oncoarray_bcac_public_release_oct17_GWAS.txt.tab
# oncoarray_bcac_public_release_oct17_iCOGS_ER_negative_disease.txt.tab
# oncoarray_bcac_public_release_oct17_iCOGS_ER_positive_disease.txt.tab
# oncoarray_bcac_public_release_oct17_iCOGS.txt.tab
# oncoarray_bcac_public_release_oct17_meta-analysis_ER_negative_disease.txt.tab
# oncoarray_bcac_public_release_oct17_meta-analysis_ER_positive_disease.txt.tab
# oncoarray_bcac_public_release_oct17_meta-analysis.txt.tab
# oncoarray_bcac_public_release_oct17_oncoarray_ER_negative_disease.txt.tab
# oncoarray_bcac_public_release_oct17_oncoarray_ER_positive_disease.txt.tab
# oncoarray_bcac_public_release_oct17_oncoarray.txt.tab
# oncoarray_bcac_public_release_oct17_oncoarray.txt.tab.sig

# Onco_TRICL_032116_Adeno.csv.tab
# Onco_TRICL_032116_Ever.csv.tab
# Onco_TRICL_032116_Never.csv.tab
# Onco_TRICL_032116_Overall.csv.tab
# Onco_TRICL_032116_Small.csv.tab
# Onco_TRICL_032116_Squam.csv.tab

# results for prostate cancer here
# cd /projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/Practical_data_app_262/App262_MRbase/MetaOverall/
