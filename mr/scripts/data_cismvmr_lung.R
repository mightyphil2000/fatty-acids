source("~/fatty-acids/colocalisation/scripts/regional_association_plots_functions.R")
source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/outcome_data/scripts/functions_combine_and_format_outcomes.R")

library(plyr)
library(gassocplot)
library(biomaRt)

data_list<-load_data(
	cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata",
	cancer3="~/fatty-acids/colocalisation/data/cancer_data_colocalisation.RData",
	# bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",
	eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_imputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
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

tissues<-unique(gtex_data1$tissue)
tissues<-tissues[!tissues %in%c("artery coronary","lung","colon transverse")]

table(cancer_data1$trait)

fads1<-format_data3(lun_data=NULL,crc_data=crc_data1,cancer_data=NULL,gtex_data=gtex_data1,fa.tab=fa.tab1[fa.tab1$trait %in% c("GLA:LA","AA:DGLA"),],ref=ref1,eqtlgen_data=eqtlgen_data1,gene_tab=gene_tab1,ld_eas_pop=NULL,gene="FADS1",region=500000,gtex_tissues=tissues,gwis_file=gwis_file1,studies="CHARGE",fix_log10=FALSE,bbj_eqtl_data=NULL,pre_formatted_data=TRUE) #fix_log10 rescales the Z scores when -log10Pvalue is greater than 1000. Should only be needed for the regional association plots 
fads2<-format_data3(lun_data=NULL,crc_data=crc_data1,cancer_data=NULL,gtex_data=gtex_data1,fa.tab=fa.tab1[fa.tab1$trait %in% c("GLA:LA","AA:DGLA"),],ref=ref1,eqtlgen_data=eqtlgen_data1,gene_tab=gene_tab1,ld_eas_pop=NULL,gene="FADS2",region=500000,gtex_tissues=tissues,gwis_file=gwis_file1,studies="CHARGE",fix_log10=FALSE,bbj_eqtl_data=NULL,pre_formatted_data=TRUE) 

ref <- read.table(ref_dat,stringsAsFactors=F,head=F)

data_list<-prep_data(Dat=fads1)
Dat<-data.frame(data_list[1],stringsAsFactors=F)
ref_dat<-data.frame(data_list[2],stringsAsFactors=F)
Dat2<-harmonise_dat(Dat=Dat,ref_dat=ref_dat,marker="marker")
fads1_data<-format_results()

data_list<-prep_data(Dat=fads2)
Dat<-data.frame(data_list[1],stringsAsFactors=F)
ref_dat<-data.frame(data_list[2],stringsAsFactors=F)
Dat2<-harmonise_dat(Dat=Dat,ref_dat=ref_dat,marker="marker")
fads2_data<-format_results()

fads2_data<-fads2_data[grep("FADS2",fads2_data$trait),]
Dat<-rbind(fads1_data,fads2_data)
Dat[Dat$marker=="rs174546",]
table(Dat$trait)
dim(Dat)
head(Dat)

Dat<-Dat[order(Dat$trait2),]

save(Dat,file="~/fatty-acids/mr/data/data_for_cismvmr_temp.RData")
load("~/fatty-acids/mr/data/data_for_cismvmr.RData")



max(Dat$eaf[which(!is.na(Dat$eaf))])

# fads1_dat<-format_data(Dat=fads1)
# fads2_dat<-format_data(Dat=fads2)

# # fads1_dat[fads1_dat$SNP=="rs174546",]
# fads2_dat<-fads2_dat<-fads2_dat[fads2_dat$trait == "FADS2 expression",]
# fads_dat<-rbind(fads1_dat,fads2_dat)

# save(list=c("fads_dat","ld.matrix"),file="~/fatty-acids/mr/data/data_for_cismvmr.RData")
# load(file="~/fatty-acids/mr/data/data_for_cismvmr.RData")

format_results<-function(){
	Dat3<-rbind.fill(ref_dat,Dat2)
	Dat3$study[is.na(Dat3$study)]<-"GECCO/CORECT/CCFR"	
	Dat4<-Dat3[!is.na(Dat3$ma_count),]
	Dat5<-Dat3[is.na(Dat3$ma_count),]
	Dat4<-Dat4[Dat4$ma_count>=10,]
	Dat<-rbind(Dat5,Dat4)
	Dat<-Dat[,c("marker","effect_allele","other_allele","eaf","beta","se","study","trait","maf","trait2")]
	return(Dat)
}

prep_data<-function(Dat=NULL){		
	gwis<-data.frame(Dat[1],stringsAsFactors=F)
	gtex<-data.frame(Dat[2],stringsAsFactors=F)
	eqtl<-data.frame(Dat[3],stringsAsFactors=F)
	crc<-data.frame(Dat[4],stringsAsFactors=F)
	gtex2<-effect_allele_gtex_data(dat=gtex)
	ref_dat<-gwis[gwis$trait == "AA:DGLA / D5D",]
	gwis2<-gwis[gwis$trait != "AA:DGLA / D5D",]

	aa_info<-read.table("~/fatty-acids/colocalisation/data/aa_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	dgla_info<-read.table("~/fatty-acids/colocalisation/data/dgla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	gla_info<-read.table("~/fatty-acids/colocalisation/data/gla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	la_info<-read.table("~/fatty-acids/colocalisation/data/la_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	
	snps1<-snps_pass(Dat=aa_info)
	snps2<-snps_pass(Dat=dgla_info)

	Pos<-ref_dat$marker %in% snps1 & ref_dat$marker %in% snps2	
	ref_dat<-ref_dat[which(Pos),]		
	
	Dat<-do.call(rbind.fill,list(gwis2,gtex2,eqtl,crc))
	
	return(list(Dat,ref_dat))
}
	

 
snps_pass<-function(Dat=NULL){
	Dat1<-Dat[Dat$snp!=".",]
	# Dat1[Dat1$V2=="rs174546",]
	# ref[ref$V2=="rs174546",]
	Dat2<-Dat[Dat$snp==".",]
	chr<-unlist(strsplit(ref$V1,split="chr"))
	ref$chr<-chr[chr!=""]
	Dat3<-merge(Dat2,ref,by.x=c("chr","bp"),by.y=c("chr","V2"))
	snps1<-Dat1$snp
	snps2<-Dat3$V4	
	snps<-c(snps1,snps2)
	return(snps)
}
format_data<-function(Dat=NULL){
	Markers<-data.frame(Dat[3])	
	B.matrix<-data.frame(Dat[5])		
	SE.matrix<-data.frame(Dat[6])
	Trait_names<-unlist(Dat[4])
	Dat_list<-NULL	
	for(i in 1:ncol(B.matrix)){
		print(i)
		Dat_list[[i]]<-do.call(cbind,list(B.matrix[,i],SE.matrix[,i],Trait_names[i]))
	}
	Dat2<-data.frame(do.call(rbind,Dat_list),stringsAsFactors=F)
	rownames(Markers)<-NULL
	Dat<-cbind(Markers,Dat2)
	names(Dat)<-c("SNP","chr","pos_grch37","effect_allele","other_allele","eaf","b","se","trait_detail")
	Dat$effect_allele<-toupper(Dat$effect_allele)
	Dat$other_allele<-toupper(Dat$other_allele)
    Pos<-grep("FADS1",Dat$trait_detail)    
    Dat$trait<-NA
    Dat$trait[Pos]<-"FADS1 expression"
    Pos<-grep("FADS2",Dat$trait_detail)
    Dat$trait[Pos]<-"FADS2 expression"    
    Dat$trait[Dat$trait_detail =="Colorectal cancer GECCO/CORECT/CCFR" ]<-"Colorectal cancer"
    Dat$study<-NA
    Dat$study[Dat$trait_detail =="Colorectal cancer GECCO/CORECT/CCFR" ]<-"GECCO/CORECT/CCFR"
    Dat$study[grep("eQTLGen",Dat$trait_detail)]<-"eQTLGen"
	Dat$study[grep("GTEx",Dat$trait_detail)]<-"GTEx"
	Dat$study[grep("CHARGE",Dat$trait_detail)]<-"CHARGE"   
    Dat$trait[Dat$trait_detail =="AA:DGLA / D5D (CHARGE)"  ]<-"AA:DGLA"
    Dat$trait[Dat$trait_detail =="GLA:LA / D6D (CHARGE)"  ]<-"GLA:LA" 
    Dat1<-Dat[grep("expression",Dat$trait_detail),]   
    tissue<-trimws(unlist(strsplit(Dat1$trait_detail,split="in")))
	Dat1$tissue<-tissue[seq(2,length(tissue),by=4)]   	
	Dat2<-Dat[grep("expression",Dat$trait_detail,invert=T),]   
    Dat2$tissue<-NA
    Dat<-rbind(Dat1,Dat2)
    return(Dat)
}

