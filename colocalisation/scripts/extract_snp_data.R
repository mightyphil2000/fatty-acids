# cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/
# source("~/fatty-acids/colocalisation/scripts/extract_snps_function.R")
ref <- read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",head=F,stringsAsFactors=F)
Chrs<-c(11,6,6,10,2,16,20)
Starts<-c(61043499,10480992,52632196,101606881,27219709,14568448,12489627)
Ends<-c(62159523,11544547,53713947,102624591,28246554,15733196,13647411)
Snplist<-lapply(1:length(Chrs),FUN=function(x)
	find_snp_positions(Chr=Chrs[x],Position=c(Starts[x],Ends[x]))
	)
Snplist<-unlist(Snplist)
Snplist<-strsplit(Snplist,split="_")
# x<-29995
Snplist2<-unlist(lapply(1:length(Snplist),FUN=function(x)
	unlist(Snplist[x])[1])
	)
write.table(Snplist2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_ukb.txt",col.names=F,row.names=F,quote=F)

##################
# extract regions#
##################

# Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/ratios/"
Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/ratios/"

# Dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/"

setwd(Dir)
Files<-dir()[grep(".tab",dir())]

# #snplist_coloc_ukb.txt contains 59588 SNPs +/- 500kb of each of 7 index SNPs selected as instruments for fatty acid traits. The SNPs reside at 7 genomic regions: FADS, ELOVL2, ELOVL5, SCD, GCKR, PDXDC1, SPTLC3
# Charge
 
Charge<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0(Dir,Files[i])	
	Charge[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_ukb.txt",File=target_file,exact_match=TRUE,file_sep="\t") 
}
Charge1<-do.call(rbind,Charge)

#######
# Shin#
#######

# Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/ratios/"
Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/imputed/ratios/"

# Dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/shin/"
setwd(Dir)
Files<-dir()[grep(".tab",dir())]
Shin<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0(Dir,Files[i])
	Shin[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_ukb.txt",File=target_file,exact_match=TRUE,file_sep="\t")  
}
Shin1<-do.call(rbind,Shin)

# Tintle
# Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/ratios/"
Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/imputed/ratios/"

# Dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/tintle/"
setwd(Dir)
Files<-dir()[grep(".tab",dir())]
Tin<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0(Dir,Files[i])
	Tin[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_ukb.txt",File=target_file,exact_match=TRUE,file_sep="\t")  
}
Tin1<-do.call(rbind,Tin)


########################
# CHARGE East Asia####
########################

ChargeEA1<-extract_snp_data(Dir="/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge_asia/ratios/",all_files=

########################
# Dorajoo East Asia####
########################

Dor1<-extract_snp_data(Dir="/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_pc/ratios/")


Dor2<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_ukb.txt",File="/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_lnD5D_pooled_allchr_qc1.tab",exact_match=TRUE,file_sep="\t") 
Dor3<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_ukb.txt",File= "/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_lnD6D_pooled_allchr_qc1.tab",exact_match=TRUE,file_sep="\t") 

library(plyr)
Dor4<-do.call(rbind.fill,list(Dor1,Dor2,Dor3))

 
###########
#BBJ eQTLs#
###########

bbj_bcells<-extract_snp_data(Dir="/projects/MRC-IEU/users/ph14916/eqtl_bbj/B_cells/fattyacids/",all_files=TRUE)
bbj_cd8tcells<-extract_snp_data(Dir="/projects/MRC-IEU/users/ph14916/eqtl_bbj/CD8+T_cells/fattyacids/",all_files=TRUE)
bbj_nkcells<-extract_snp_data(Dir="/projects/MRC-IEU/users/ph14916/eqtl_bbj/NK_cells/fattyacids/",all_files=TRUE)
bbj_cd4tcells<-extract_snp_data(Dir="/projects/MRC-IEU/users/ph14916/eqtl_bbj/CD4+T_cells/fattyacids/",all_files=TRUE)
bbj_mono<-extract_snp_data(Dir="/projects/MRC-IEU/users/ph14916/eqtl_bbj/Monocytes/fattyacids/",all_files=TRUE)
bbj_blood<-extract_snp_data(Dir="/projects/MRC-IEU/users/ph14916/eqtl_bbj/Peripheral_blood/fattyacids/",all_files=TRUE)

bbj_list<-ls()[grep("bbj",ls())]

########################
# Dorajoo East Asia####
########################
############
# save data####
############

if(sum(grep("imputed",Dir))==1){
	save(list=c("Charge1","Shin1","Tin1"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_ukb_gwis_ratios_imputed.Rdata")	
}
if(sum(grep("imputed",Dir))==0){
	save(list=c("Charge1","Shin1","Tin1"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata")
}

save(list=c("ChargeEA1","Dor4"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata")

save(list=bbj_list,file="/projects/MRC-IEU/users/ph14916/eqtl_bbj/snplist_coloc_ukb_eqtlbbj.Rdata")


cd ~/fatty-acids/colocalisation/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_ukb_gwis_ratios_imputed.Rdata .

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/eqtl_bbj/snplist_coloc_ukb_eqtlbbj.Rdata .


########################
# extract ELOVL5 region#
########################


# Elov5_charge<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/DGLA_to_GLA.tab",exact_match=TRUE,file_sep="\t")
# Elov5_charge$study<-"charge"
# Elov5_fram<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/tintle/DGLA_to_GLA.tab",exact_match=TRUE,file_sep="\t")
# Elov5_fram$study<-"fram"
# Elov5_shin<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/shin/DGLA_to_GLA.tab",exact_match=TRUE,file_sep="\t")
# Elov5_shin$study<-"shin"
# elov5<-do.call(rbind,list(Elov5_charge,Elov5_fram,Elov5_shin))
# save(elov5,file="elovl5_dgla_gla.Rdata")

# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/elovl5_dgla_gla.Rdata .

extract_snp_data<-function(Dir=NULL,all_files=FALSE){	
	setwd(Dir)
	if(!all_files){
		Files<-dir()[grep(".tab",dir())]
	}else{
		Files<-dir()
	}
	Dat<-NULL 
	for(i in 1:length(Files)){
		print(i)
		print(Files[i])
		target_file<-paste0(Dir,Files[i])
		Dat[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_ukb.txt",File=target_file,exact_match=TRUE,file_sep="\t")  
	}
	Dat1<-do.call(rbind,Dat)
	return(Dat1)
}

b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}

format_dat<-function(dat=NULL,exposure=NULL,consortium="CHARGE",file=NULL,population="European",mechanism="FA biosynthesis",gene=NULL){
	z<-dat$beta/dat$se
	n<-dat$samplesize
	maf<-dat$eaf
	if(maf>0.5){
		maf<-1-maf
	}
	beta.sd<-b_sd(z,maf,n)
	dat$se.sd<-abs(beta.sd /z)
	dat$beta.sd<-beta.sd
	dat$r2<-2*beta.sd^2*maf*(1-maf)
	dat$exposure<-exposure
	dat$trait<-exposure
	dat$consortium<-consortium
	dat$type.fa<-"other"
	dat$file<-file
	dat$population<-population
	load("~/fatty-acids-mr/instruments/fatty_acid_instruments_ldc_v3.Rdata") 
	dat$Chr<-unique(Dat.ldc$Chr[Dat.ldc$SNP == dat$SNP])
	dat$SNP.GRCh38.p12 <-  unique(Dat.ldc$SNP.GRCh38.p12[Dat.ldc$SNP == dat$SNP])
	dat$mechanism<-mechanism
	dat$gene <- gene
	return(dat)
}


