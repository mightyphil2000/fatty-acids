# rm(list=ls())
library(TwoSampleMR)
library(plyr)
# devtools::install_github('MRCIEU/TwoSampleMR')

# source(extrac_snps_functions.R)
instruments<-prep_instruments(load_ref=TRUE)

ref4<-data.frame(instruments[1])
bed37<-unlist(instruments[2])
bed37_2<-unlist(instruments[3])
snps<-ref4$V2

#######################
# extract fatty acids#
######################
# cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/
Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed"
Files<-dir(Dir)
Files<-Files[grep(".tab",Files)]
Files1<-paste(Dir,"/",Files,sep="")

Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_pc"
Files<-dir(Dir)
Files<-Files[grep(".tab",Files)]
Files2<-paste(Dir,"/",Files,sep="")

Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tanaka/hg19/imputed"
Files<-dir(Dir)
Files<-Files[grep(".tab",Files)]
Files3<-paste(Dir,"/",Files,sep="")
Files<-c(Files1,Files2,Files3)

Res<-NULL
for(i in 1:length(Files)){
	print(Files[i])
	Res[[i]]<-extract_snps(snplist=snps,File=Files[i],exact_match=TRUE,file_sep="\t")
}

Faa<-do.call(rbind.fill,Res)

write.table(Faa,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/pufas_sensitivityhp.txt",sep="\t",col.names=T,row.names=F,quote=F)
cd ~/fatty-acids/mr/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/pufas_sensitivityhp.txt .

pufas<-read.table("~/fatty-acids/mr/data/pufas_sensitivityhp.txt",head=T,stringsAsFactors=F,sep="\t")


fat<-read.table("~/fatty-acids/mr/data/fatty_acid_GWASanalysis_table2.txt",sep="\t",stringsAsFactors=F,head=T)

fat$filename<-gsub(".txt",".tab",fat$filename)
fat$filename<-gsub("pooled_allchr_qc","pooled_allchr_qc1",fat$filename)
fat$filename<-gsub("\xca","",fat$filename)
fat$filename<-trimws(fat$filename)

fat$filename[!fat$filename %in% pufas$file]
pufas$file[!pufas$file %in% fat$filename]

pufas.m<-merge(pufas,fat[,c("pmid","author","consortium","trait",	"chain","chain.length","SigSNPs","	sample_size.analysis","population","sex","year","filename")],by.x="file",by.y="filename")

names(pufas.m)[names(pufas.m)=="effect_allele_freq"]<-"eaf"
names(pufas.m)[names(pufas.m)=="trait"]<-"exposure"
names(pufas.m)[names(pufas.m)=="snp"]<-"SNP"

# ID<-paste(pufas.m$exposure,pufas.m$SNP)
# ID_dup<-unique(ID[duplicated(ID)])
# pufas.m[ID %in% ID_dup,]
pufas.m<-pufas.m[order(pufas.m$n,decreasing=T),]

ID<-paste(pufas.m$population,pufas.m$exposure,pufas.m$SNP)

pufas.m<-pufas.m[!duplicated(ID),]
pufas.m$exposure<-paste(pufas.m$exposure,pufas.m$population)

exposure_dat<-define_exposure()
define_exposure<-function(Dat=pufas.m){
	# PDXDC1/NTAN1 & arachidonic acid 
	Dat<-Dat[Dat$SNP %in% c("rs16966952","rs4985155","rs4985167","rs2280018")	,]
	Dat<-Dat[grep("arachidonic",Dat$exposure),]
	return(Dat)
}

write.table(exposure_dat,"~/fatty-acids/mr/data/exposure_dat.txt",sep="\t",col.names=T,row.names=F,quote=F)



###############
#Biobank Japan#
###############

# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/

setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/auto_rsq07.mac10")
Files<-"CoCa.auto.rsq07.mac10.txt" 
Res<-extract_snps(snplist=bed37,File=Files,exact_match=FALSE,file_sep=" ")

bbj<-merge(Res,ref4,by.x=c("CHR","POS"),by.y=c("V1","V4"),all.x=T)
bbj$outcome<-"Colorectal cancer"
names(bbj)[names(bbj) == "Allele2"]<-"Effect.Allele"
names(bbj)[names(bbj) == "Allele1"]<-"Other.Allele"
names(bbj)[names(bbj) == "AF_Allele2"]<-"eaf"
names(bbj)[names(bbj) == "BETA"]<-"lnor"
names(bbj)[names(bbj) == "SE"]<-"se"
names(bbj)[names(bbj) == "p.value"]<-"p"
names(bbj)[names(bbj) == "Rsq"]<-"info"
names(bbj)[names(bbj) == "V2"]<-"rsid"
bbj$study<-"Biobank Japan"
bbj$population<-"East Asian"
bbj$pmid<-"biorxiv"
bbj$effect_allele_confirmed<-TRUE
bbj$ncase[bbj$file == "CoCa.auto.rsq07.mac10.txt"]<-7062
bbj$ncontrol[bbj$file == "CoCa.auto.rsq07.mac10.txt"]<-195745 
bbj$outcome[bbj$file == "CoCa.auto.rsq07.mac10.txt"]<-"Colorectal cancer"
write.table(bbj,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/bbj_senshp.txt",sep="\t",col.names=T,row.names=F,quote=F)

cd ~/fatty-acids/mr/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/bbj_senshp.txt .

#############
# UK Biobank#
#############
"rs174546"
setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/")
setwd("/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/overall_colorectal_cancer")
setwd("/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results")
Files<-"/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/overall_colorectal_cancer/5c93555a-faa4-4053-9c19-1f5dfa4381ac/overall_colorectal_cancer_imputed.txt.gz"
Files2<-"/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/overall_pan_cancer/34499930-c7ef-4be8-8647-263dbeecd2fc/overall_pan_cancer_imputed.txt.gz"

Res<-extract_snps(snplist=snps,File=Files2,exact_match=TRUE,file_sep="\t")
Res$outcome<-"Overall cancer" 
Res1<-extract_snps(snplist=snps,File=Files,exact_match=TRUE,file_sep="\t")
Res1$outcome<-"Colorectal cancer"
ukb<-rbind(Res,Res1)

ukb$study<-"UK Biobank"
ukb$population<-"European"
ukb$effect_allele_confirmed<-TRUE
ukb$pmid<-"unpublished"

ukb$ncase<-NA
ukb$ncontrol<-NA
ukb$ncase[ukb$outcome=="Colorectal cancer"]<-5713
ukb$ncontrol[ukb$outcome=="Colorectal cancer"]<-372016
ukb$ncase[ukb$outcome=="Overall cancer"]<-52400
ukb$ncontrol[ukb$outcome=="Overall cancer"]<-372016


beta<-ukb$BETA
se<-ukb$SE
U<-ukb$ncase/(ukb$ncontrol+ukb$ncase)
ukb$lnor<-beta/(U*(1-U)) 
ukb$se<-se/(U*(1-U))

names(ukb)[names(ukb) == "SNP"]<-"rsid"
names(ukb)[names(ukb) == "ALLELE1"]<-"Effect.Allele" 
names(ukb)[names(ukb) == "ALLELE0"]<-"Other.Allele" 
names(ukb)[names(ukb) == "A1FREQ"]<-"eaf" 
names(ukb)[names(ukb) == "INFO"]<-"info" 
names(ukb)[names(ukb) == "P_BOLT_LMM_INF"]<-"p" #the P value calcualted from the beta and standard error match P_BOLT_LMM_INF

write.table(ukb,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/ukb_senshp.txt",sep="\t",col.names=T,row.names=F,quote=F)

# temp<-read.table("~/MR_FattyAcids/data/summary_data/Final/ukbiobank.txt",sep="\t",head=T,stringsAsFactors=F)
# temp[temp$rsid=="rs174546",]
# ukb[ukb$rsid=="rs174546",c("outcome","lnor","se")]

cd ~/fatty-acids/mr/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/ukb_senshp.txt .

###############
# CORECT/GECCO#
###############

instruments<-prep_instruments(load_ref=TRUE,local=TRUE)
ref4<-data.frame(instruments[1])
snps<-ref4$V2

load("~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata")
crg<-CRC[CRC$rs_number %in% ref4$V2,]
crg$outcome<-"Colorectal cancer"
crg$study<-"CORECT/GECCO"
crg$pmid<-30510241 
crg$ncase<-58221	
crg$ncontrol<-67694
crg$outcome<-"Colorectal cancer"
crg$population<-"European"
crg$UKbiobank<-TRUE

names(crg)[names(crg) == "rs_number"]<-"rsid"
names(crg)[names(crg) == "Effect"]<-"lnor"
names(crg)[names(crg) == "StdErr"]<-"se"
names(crg)[names(crg) == "P.value"]<-"p"
names(crg)[names(crg) == "HetISq"]<-"I2"
names(crg)[names(crg) == "HetChiSq"]<-"Q"
names(crg)[names(crg) == "HetPVal"]<-"Phet"
names(crg)[names(crg) == "Allele1"]<-"Effect.Allele"
names(crg)[names(crg) == "Allele2"]<-"Other.Allele"
names(crg)[names(crg) == "Freq1"]<-"eaf"
crg$effect_allele_confirmed<-TRUE

#######
# ACCC#
#######

accc<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt",sep="\t",head=T,stringsAsFactors=F)
SNP<-gregexpr(":",accc$snp)
Test<-unlist(lapply(1:length(SNP),FUN=function(x)
	length(unlist(SNP[x]))))
Pos<-which(Test==3)
accc2<-accc[Pos,]
SNP<-unlist(strsplit(accc2$snp,split=":"))
accc2$chr<-SNP[seq(1,length(SNP),by=4)]
accc2$bp<-SNP[seq(2,length(SNP),by=4)]

acc<-merge(accc2,ref4,by.x=c("chr","bp"),by.y=c("V1","V4"))
acc$outcome<-"Colorectal cancer"
names(acc)[names(acc) == "Allele1"]<-"Effect.Allele"
names(acc)[names(acc) == "Allele2"]<-"Other.Allele"
names(acc)[names(acc) == "Freq1"]<-"eaf"
names(acc)[names(acc) == "Effect"]<-"lnor"
names(acc)[names(acc) == "StdErr"]<-"se"
names(acc)[names(acc) == "P.value"]<-"p"
names(acc)[names(acc) == "HetPVal"]<-"Phet"
names(acc)[names(acc) == "V2"]<-"rsid"
acc$study<-"ACCC" #The Asian Colorectal Cancer Consortium?
acc$population<-"East Asian"
acc$pmid<-"31826910"
acc$effect_allele_confirmed<-TRUE
acc$UKbiobank<-FALSE
acc$ncase<-23572
acc$ncontrol<-48700 

ukb<-read.table("~/fatty-acids/mr/data/ukb_senshp.txt",sep="\t",head=T,stringsAsFactors=F)
bbj<-read.table("~/fatty-acids/mr/data/bbj_senshp.txt",sep="\t",head=T,stringsAsFactors=F)


crc<-do.call(rbind.fill,list(ukb,bbj,crg,acc))
crc.m<-merge(crc,ref4,by.x="rsid",by.y="V2")
crc.m<-crc.m[,names(crc.m) != "chr"]
names(crc.m)[names(crc.m) == "V1"]<-"chr"
names(crc.m)[names(crc.m) == "V4"]<-"bp_GRCh37"

# crc.m[,c("rsid","chr","bp_GRCh37","Effect.Allele","Other.Allele","eaf","p","file","outcome","study","population","effect_allele_confirmed","pmid","ncase","ncontrol","lnor","se")]


crc.m$outcome<-paste(crc.m$outcome,crc.m$study)
write.table(crc.m,"~/fatty-acids/mr/data/outcome_dat.txt",sep="\t",col.names=T,row.names=F,quote=F)



prep_instruments<-function(load_ref=FALSE,local=FALSE){
	ins<-c("rs174546", #FADS 
	"rs3734398", #ELOVL2
	"rs603424", #SCD (POA:PA)
	"rs6584379", #second POA:PA hit
	"rs4985155", #PDXDC1/NTAN1
	"rs4985167",#PDXDC1/NTAN1
	"rs780093", #GCKR
	"rs680379", #SPTLC3
	"rs12210577", #ELOVL5
	"rs16966952", #PDXDC1/NTAN1
	"rs2280018", #PDXDC1/NTAN1
	"rs10740118", # LA was associated with multiple SNPs on chromosome 10 in a region that included nuclear receptor binding factor 2 (NRBF2), jumonji domain containing 1C (JMJD1c), and receptor accessory protein 3 (REEP3) (Figure 3a). no associations with other n6 fatty acids in Guan paper
	"rs3134950") #For AGPAT1 on chromosome 6, we observed an association of the rs3134950 SNP and AdrA, following adjustment for its fatty acid precursor, AA. 


	eqtls<-c("rs78244690", #eqTL ELOVL5 in adipose 
	"rs209533", #eqTL ELOVL5 in blood
	"rs2562896", #sQTL ELOVL5 in colon sigmoid
	"rs2562896", #sQTL ELOVL5 in colon transverse
	"rs9382192", #sQTL ELOVL5 in liver 
	"rs2057024", #sQTL ELOVL5 in adipose subcutaneous
	"rs2057023", #sQTL ELOVL5 in adipose visceral
	"rs2057024", #sQTL ELOVL5 in blood
	"rs239585", #eQTL for ELOVL4 in adipose subcutaneous
	"rs240237", #eQTL for ELOVL4 in adipose visceral
	"rs72894911") #eQTL for ELOVL4 in blood


	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K_bed_hg38.txt.gz",  

	if(load_ref==TRUE){
		ref_dat="/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim" #GRCh37
		if(local){
			ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz"
		}
		ref <- read.table(ref_dat,stringsAsFactors=F,head=FALSE)
	}

	# all(ins %in% ref$V2  )
	# all(eqtls %in% ref$V2 )
	ref2<-ref[ref$V2 %in% ins,]
	ref3<-ref[ref$V2 %in% eqtls,]
	ref4<-rbind(ref2,ref3)
	bed37<-paste(ref4$V1,ref4$V4,sep="_")
	bed37_2<-gsub("_",":",bed37)
	return(list(ref4,bed37,bed37_2))
}

