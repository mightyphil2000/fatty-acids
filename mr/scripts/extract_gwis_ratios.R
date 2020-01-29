
# source("~/fatty-acids-mr/instruments/Extract_SNPs_function.R")
# ref <- read.table("~/fatty-acids-mr/instruments/data_maf0.01_rs.bim.gz")
# Snplist<-find_snp_positions(Chr=6,Position=c(53132196,53213977))
# /write.table(Snplist,"~/fatty-acids-mr/instruments/ELOVL5")

########################
# extract FADS region#
########################

setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge")
Files<-dir()[grep(".tab",dir())]

# Charge
Charge<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/",Files[i])
	Charge[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt",File=target_file,exact_match=TRUE,file_sep="\t") #snplist_coloc.txt contains 17381 SNPs +/- 500kb of each of 9 index SNPs selected as instruments for fatty acid traits. The SNPs reside at 6 genomic regions: FADS, ELOVL2, SCD, GCKR, PDXDC1, SPTLC3 
}
Charge1<-do.call(rbind,Charge)

# Shin
setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/shin")
Files<-dir()[grep(".tab",dir())]
Shin<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/shin/",Files[i])
	Shin[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt",File=target_file,exact_match=TRUE,file_sep="\t") #snplist_coloc.txt contains 17381 SNPs +/- 500kb of each of 9 index SNPs selected as instruments for fatty acid traits. The SNPs reside at 6 genomic regions: FADS, ELOVL2, SCD, GCKR, PDXDC1, SPTLC3 
}
Shin1<-do.call(rbind,Shin)

# Tintle
setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/tintle")
Files<-dir()[grep(".tab",dir())]
Tin<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/tintle/",Files[i])
	Tin[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt",File=target_file,exact_match=TRUE,file_sep="\t") #snplist_coloc.txt contains 17381 SNPs +/- 500kb of each of 9 index SNPs selected as instruments for fatty acid traits. The SNPs reside at 6 genomic regions: FADS, ELOVL2, SCD, GCKR, PDXDC1, SPTLC3 
}
Tin1<-do.call(rbind,Tin)

save(list=c("Charge1","Shin1","Tin1"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_gwis_ratios.Rdata")

cd ~/fatty-acids-mr/coloc
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/snplist_coloc_gwis_ratios.Rdata .

########################
# extract ELOVL5 region#
########################


Elov5_charge<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/DGLA_to_GLA.tab",exact_match=TRUE,file_sep="\t")
Elov5_charge$study<-"charge"
Elov5_fram<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/tintle/DGLA_to_GLA.tab",exact_match=TRUE,file_sep="\t")
Elov5_fram$study<-"fram"
Elov5_shin<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt",File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/shin/DGLA_to_GLA.tab",exact_match=TRUE,file_sep="\t")
Elov5_shin$study<-"shin"
elov5<-do.call(rbind,list(Elov5_charge,Elov5_fram,Elov5_shin))
save(elov5,file="elovl5_dgla_gla.Rdata")

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/elovl5_dgla_gla.Rdata .

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


