# map rsids to ld matrix freq file
# restrict file to snps present in ukbiobank bim file
# remove positions mappign to >1 rsid
# remove SNPs with 0 MAF / NaN correlation
# retrieve SNPs from fatty acid ratio files (AA:DGLA and LA:GLA)
# restrict all files to common set of SNPs. 
# make sure effect allele in fatty acid results files are same as effect alleles in LD files

chs<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/CHS_ld_reference_panel/CHS_EA_HRC_LD_11_61043499_62159523.frq",head=T,stringsAsFactors=F)
# nrow(chs)-length(which(chs$MAF==0))

ukb<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",head=F,stringsAsFactors=F)

chs.ld<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/CHS_ld_reference_panel/CHS_EA_HRC_LD_11_61043499_62159523.ld",head=F,stringsAsFactors=F)

ukb1<-ukb[ukb$V1 == 11,]
ukb1$SNP<-paste(ukb1$V1,":",ukb1$V4,sep="")
names(ukb1)[names(ukb1) == "V2"]<-"rsid"
names(ukb1)[names(ukb1) == "V5"]<-"allele1_ukb"
names(ukb1)[names(ukb1) == "V6"]<-"allele2_ukb"
names(ukb1)[names(ukb1) == "V4"]<-"pos_grch37"
names(ukb1)[names(ukb1) == "V1"]<-"chr"
ukb1<-ukb1[,names(ukb1) != "V3"]

Pos<-ukb1$pos_grch37[duplicated(ukb1$pos_grch37)]
length(Pos) # 409 duplicate positions
ukb2<-ukb1[!ukb1$pos_grch37 %in% Pos,]
# ukb2<-ukb1[ukb1$pos_grch37 %in% Pos,]


Pos<-match(ukb2$SNP,chs$SNP)
# Match finds the position of the SNPs in first argument in the second argument. The positions correspond to the second argument. 
# Dat<-c("abc","gef","bcd","fhi")
# Pos<-match(c("bcd","abc","gef","mnb","xcv","AAA","bbbb"),Dat)

Pos<-Pos[!is.na(Pos)]
chs1<-chs[Pos,]
chs.ld2<-chs.ld[Pos,Pos]
chs.ld3<-cbind(chs1$SNP,chs.ld2)


ukb3<-ukb2[ukb2$SNP %in% chs1$SNP,]
# head(ukb3)
# head(chs1)
all(chs1$SNP == ukb3$SNP)
any(chs1$SNP != ukb3$SNP)
all(ukb3$SNP == chs1$SNP )
any(ukb3$SNP != chs1$SNP )

chs2<-cbind(chs1,ukb3$rsid)
names(chs.ld3)[names(chs.ld3) == "chs1$SNP"]<-"SNP"
all(chs2$SNP == chs.ld3$SNP)
any(chs2$SNP != chs.ld3$SNP)
all(chs.ld3$SNP == chs2$SNP)
any(chs.ld3$SNP != chs2$SNP)

chs.ld4<-chs.ld3[,names(chs.ld3) != "SNP"]

Pos<-which(chs2$MAF == 0)
chs2[Pos,]
# head(chs2)
all(chs.ld4[Pos,]=="NaN")
all(chs.ld4[,Pos]=="NaN")
Pos<-which(chs2$MAF != 0)
any(chs2[Pos,]=="NaN")
any(chs.ld4[Pos,Pos]=="NaN")
chs3<-chs2[Pos,]
chs.ld5<-chs.ld4[Pos,Pos]
any(chs.ld5=="NaN")
any(chs3=="NaN")
names(chs3)[names(chs3) =="ukb3$rsid"]<-"rsid"
chs3$rsid<-as.character(chs3$rsid)

Temp<-ukb2[ukb2$SNP %in% chs1$SNP,]
all(Temp$SNP == chs1$SNP)
any(Temp$SNP != chs1$SNP)
all(chs1$SNP == Temp$SNP)
any(chs1$SNP  != Temp$SNP)

######################################################################################
# Extract SNPs for AA:DGLA (FADS1 enzyme activity) and LA:GLA (FADS2 enzyme activity)#
######################################################################################

setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge")

d5d<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/AA_to_DGLA.tab",sep="\t",head=T,stringsAsFactors=F)
d6d<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/GLA_to_LA.tab",sep="\t",head=T,stringsAsFactors=F)
d6d[which(d6d$SNP == "rs11231053"),]
# d5d2<-d5d[d5d$SNP %in% snps,]

any(duplicated(d5d$SNP))
any(duplicated(d6d$SNP))

Pos<-match(chs3$rsid,d5d$SNP)
Pos<-Pos[!is.na(Pos)]
d5d2<-d5d[Pos,]

Pos<-match(chs3$rsid,d5d$SNP)
Pos<-Pos[!is.na(Pos)]
d5d2<-d5d[Pos,]
Pos<-match(chs3$rsid,d6d$SNP)
Pos<-Pos[!is.na(Pos)]
d6d2<-d6d[Pos,]

# d6d2 has 2 more SNPs than d5d2
# restrict d6d2 to same SNPs as d5d2. and make order of SNPs in d6d2 same as order of SNPs in d5d2
Pos<-match(d5d2$SNP,d6d2$SNP) #find position of SNPs in argument 1 in argument 2. Position relates to d6d2 but order of positions corresponds to argument 1
d6d3<-d6d2[Pos,]
all(d6d3$SNP == d5d2$SNP)
any(d6d3$SNP != d5d2$SNP)
all(d5d2$SNP== d6d3$SNP )
any(d5d2$SNP!= d6d3$SNP )

Pos<-which(chs3$rsid %in% d5d2$SNP)
# using match would get the positinos in ld matrix same order as fatty acid files
chs4<-chs3[Pos,]
chs.ld6<-chs.ld5[Pos,Pos]

all(chs4$rsid == d5d2$SNP)
any(chs4$rsid != d5d2$SNP)
all(d5d2$SNP== chs4$rsid )
any(d5d2$SNP!= chs4$rsid )

names(chs4)[names(chs4) == "A1"]<-"effect_allele"
names(chs4)[names(chs4) == "A2"]<-"other_allele"
names(chs4)[names(chs4) == "MAF"]<-"eaf"
names(chs4)[names(chs4) == "SNP"]<-"chr_bp"
names(chs4)[names(chs4) == "rsid"]<-"SNP"


d6d3_h<-harmonise_dat(Dat=d6d3,ref_dat=chs4)
d5d2_h<-harmonise_dat(Dat=d5d2,ref_dat=chs4)
dim(d6d3_h)
dim(d5d2_h)

all(d6d3_h$SNP == d5d2_h$SNP)
all(d5d2_h$effect_allele == d5d2_h$ref_ea)
all(d5d2_h$other_allele == d5d2_h$ref_oa)

dim(d5d2_h)

Pos<-which(chs4$SNP %in% d5d2_h$SNP)
chs5<-chs4[Pos,]
chs.ld7<-chs.ld6[Pos,Pos]

dim(chs4)
dim(chs5)
dim(chs.ld6)
dim(chs.ld7)

all(d5d2_h$SNP == d6d3_h$SNP)
all(d5d2_h$SNP == chs5$SNP)
any(d5d2_h$SNP != chs5$SNP)
all(chs5$SNP== d5d2_h$SNP )
any(chs5$SNP!= d5d2_h$SNP )

all(d6d3_h$SNP == chs5$SNP)
any(d6d3_h$SNP != chs5$SNP)
all(chs5$SNP== d6d3_h$SNP )
any(chs5$SNP!= d6d3_h$SNP )

dim(d5d2_h)
dim(d6d3_h)
dim(chs5)
dim(chs.ld7)

write.table(d5d2_h,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/AA_to_DGLA_fine_mapping.tab",sep="\t",col.names=T,row.names=F,quote=F)
write.table(d6d3_h,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/GLA_to_LA_fine_mapping.tab",sep="\t",col.names=T,row.names=F,quote=F)
write.table(chs5,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/CHS_EA_HRC_LD_11_61043499_62159523_fine_mapping.frq",sep="\t",col.names=T,row.names=F,quote=F)
write.table(chs.ld7,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/CHS_EA_HRC_LD_11_61043499_62159523_fine_mapping.ld",sep="\t",col.names=T,row.names=F,quote=F)
write.table(chs5$SNP,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/snplist_fine_mapping.txt",col.names=F,row.names=F,quote=F)


#####################################################################################################################
# Extract SNPs for lung cancer and colorectal cancer and harmonise with d5d, d6d and chs ld files for colocalisation#
#####################################################################################################################

snplist<-readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/snplist_fine_mapping.txt")
d5d<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/AA_to_DGLA_fine_mapping.tab",sep="\t",head=T,stringsAsFactors=F)
d6d<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/GLA_to_LA_fine_mapping.tab",sep="\t",head=T,stringsAsFactors=F)
chs<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/CHS_EA_HRC_LD_11_61043499_62159523_fine_mapping.frq",sep="\t",head=T,stringsAsFactors=F)
chs.ld<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/CHS_EA_HRC_LD_11_61043499_62159523_fine_mapping.ld",sep="\t",head=T,stringsAsFactors=F)

Lun<-read.table("/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc/Onco_TRICL_032116_Overall.csv.tab",sep="\t",head=T,stringsAsFactors=F)
Crc<-read.csv("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/colorectal_cancer/1255_MarginalResults_HRC125K_20191105.csv",head=T,stringsAsFactors=F)

Lun2<-Lun[Lun$snp %in% snplist,]

if(any(duplicated(Crc$SNP_Name))) warning("duplicate positions")
if(any(duplicated(Crc$rs_number))) warning("duplicate rsids")
Dups<-Crc$rs_number[duplicated(Crc$rs_number)]
Crc2<-Crc[!Crc$rs_number %in% Dups,]
Crc2<-Crc2[Crc2$rs_number %in% snplist,]
Crc2<-Crc2[Crc2$rs_number %in% Lun2$snp,]

all(Crc2$rs_number %in% Lun2$snp)

dim(Crc2)
Pos<-match(Crc2$rs_number,Lun2$snp)
Lun3<-Lun2[Pos,]
all(Lun3$snp == Crc2$rs_number)
any(Lun3$snp != Crc2$rs_number)
all(Crc2$rs_number == Lun3$snp)
any(Crc2$rs_number != Lun3$snp)

Pos<-match(Crc2$rs_number,d5d$SNP)
d5d2<-d5d[Pos,]
all(Crc2$rs_number == d5d2$SNP)
any(Crc2$rs_number != d5d2$SNP)

Pos<-match(Crc2$rs_number,d6d$SNP)
d6d2<-d6d[Pos,]
all(Crc2$rs_number == d6d2$SNP)
any(Crc2$rs_number != d6d2$SNP)

Pos<-match(Crc2$rs_number,chs$SNP)
chs2<-chs[Pos,]
chs.ld2<-chs.ld[Pos,Pos]
all(Crc2$rs_number == chs2$SNP)
any(Crc2$rs_number != chs2$SNP)



names(Lun3)[names(Lun3) == "snp"]<-"SNP"
names(Lun3)[names(Lun3) == "effect_allele_freq"]<-"eaf"
names(Lun3)[names(Lun3) == "p"]<-"pval"
names(Lun3)[names(Lun3) == "n"]<-"samplesize"
Lun3$ncases<-29266
Lun3$ncontrols<-56450
Lun3$trait<-"Lung cancer"

names(Crc2)[names(Crc2) == "rs_number"]<-"SNP"
names(Crc2)[names(Crc2) == "Allele1"]<-"effect_allele"
names(Crc2)[names(Crc2) == "Allele2"]<-"other_allele"
names(Crc2)[names(Crc2) == "Freq1"]<-"eaf"
names(Crc2)[names(Crc2) == "Effect"]<-"beta"
names(Crc2)[names(Crc2) == "StdErr"]<-"se"
SNP_Name Chromosome
names(Crc2)[names(Crc2) == "P.value"]<-"pval"
names(Crc2)[names(Crc2) == "Chromosome"]<-"CHR"
Crc2$trait<-"Colorectal cancer"
Crc2$ncases<-58221	
Crc2$ncontrols<-67694
Crc2$effect_allele<-toupper(Crc2$effect_allele)
Crc2$other_allele<-toupper(Crc2$other_allele)

d5d2$trait<-"ratio_AA_to_DGLA"
d6d2$trait<-"ratio_GLA_to_LA"

dim(d5d2)
dim(d6d2)
dim(chs2)
dim(chs.ld2)
dim(Lun3)
dim(Crc2)

d5d2[d5d2$eaf>0.50,]
chs2[chs2$eaf>0.50,]
chs2[chs2$SNP == "rs741887",]
d5d2[d5d2$SNP == "rs741887",]
d6d2[d6d2$SNP == "rs741887",]

all(d6d2$effect_allele == d5d2$effect_allele)
all(chs2$effect_allele == d5d2$effect_allele)


Crc_h<-harmonise_dat(Dat=Crc2,ref_dat=chs2)
Lun_h<-harmonise_dat(Dat=Lun3,ref_dat=chs2)

all(Crc_h$effect_allele == chs2$effect_allele)
all(Lun_h$effect_allele == chs2$effect_allele)

all(Crc_h$SNP == chs2$SNP)
any(Crc_h$SNP != chs2$SNP)
all(chs2$SNP== Crc_h$SNP )
any(chs2$SNP!= Crc_h$SNP )

all(Lun_h$SNP == chs2$SNP)
any(Lun_h$SNP != chs2$SNP)
all(chs2$SNP == Lun_h$SNP)
any(chs2$SNP != Lun_h$SNP)

write.table(Crc_h,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/colorectal_cancer_forcoloc.tab",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Lun_h,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/lung_cancer_forcoloc.tab",sep="\t",col.names=T,row.names=F,quote=F)
write.table(d6d2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/GLA_to_LA_forcoloc.tab",sep="\t",col.names=T,row.names=F,quote=F)
write.table(d5d2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/AA_to_DGLA_forcoloc.tab",sep="\t",col.names=T,row.names=F,quote=F)
write.table(chs2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/CHS_LD_forcoloc.frq",sep="\t",col.names=T,row.names=F,quote=F)
write.table(chs.ld2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/CHS_LD_forcoloc.ld",sep="\t",col.names=T,row.names=F,quote=F)

dim(Crc_h)
dim(Lun_h)
dim(d6d2)
dim(d5d2)
dim(chs2)
dim(chs.ld2)

# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/* .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/* .

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/snplist_fine_mapping.txt .

scp ph14916@bluecrystalp3.acrc.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/snplist_fine_mapping.txt .

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/fine_mapping/data/
# Dat<-d6d3
# ref_dat<-chs4
# Dat=Lun3
# ref_dat=chs2

harmonise_dat<-function(Dat=NULL,ref_dat=NULL){
	names(ref_dat)[names(ref_dat) == "effect_allele"]<-"ref_ea"
	names(ref_dat)[names(ref_dat) == "other_allele"]<-"ref_oa"
	names(ref_dat)[names(ref_dat) == "eaf"]<-"ref_eaf"
	Dat<-merge(Dat,ref_dat,by="SNP")
	Alleles<-paste(Dat$effect_allele,Dat$other_allele,sep="")
	Pos<-Alleles %in% c("GC","CG","TA","AT")
	Dat2<-Dat[Pos,] #palindromic SNPs
	Dat3<-Dat[!Pos,] #non palindromic SNPs
	
	# harmonise non palindromic SNPs
	ea<-Dat3$effect_allele
	oa<-Dat3$other_allele
	eaf<-Dat3$eaf
	beta<-Dat3$beta
	if(any(ea != Dat3$ref_ea & ea != Dat3$ref_oa)) stop("some SNPs are on different strands to the reference")
	Pos<-which(ea != Dat3$ref_ea) #positions where effect allele is different from effect allele in reference set
	Dat3$effect_allele[Pos]<-oa[Pos]
	Dat3$other_allele[Pos]<-ea[Pos]
	Dat3$eaf[Pos]<-1-eaf[Pos]
	Dat3$beta[Pos]<-beta[Pos]*-1
	all(Dat3$effect_allele == Dat3$ref_ea)
	
	# deal with palindromic SNPs
	# check that palindromic SNPs "appear to be on same strand". 
	ea<-Dat2$effect_allele
	if(any(ea != Dat2$ref_ea & ea != Dat2$ref_oa)) warning("some palindromic SNPs appear to be on different strands to the reference, indicating allele errors. these SNPs will be dropped")
	Pos<-ea != Dat2$ref_ea & ea != Dat2$ref_oa
	Dat2_1<-Dat2[!Pos,]# drop SNPs that are not palindromic in both the test and reference datasets. E.g. one rs11231053 is C/G in CHARGE GLA:LA but is A/G UK biobank bim file and dbSNP.  
	# ukb[which(ukb$V2 == "rs11231053"),]

	# code all palindromic SNPs so that effect allele is the minor allele. The effect allele is always ≤0.5 in the reference dataset (i.e. CHS LD file). Drop palindromic SNPs with MAF close to 0.5 defined as >0.42
	if(any(ref_dat$eaf>0.5)) stop("some effect alleles are not minor allele in reference dataset")
	if(any(Dat2_1$ref_eaf>0.5)) stop("some effect alleles of palindromic SNPs are not the minor allele")
	ea<-Dat2_1$effect_allele
	oa<-Dat2_1$other_allele
	eaf<-Dat2_1$eaf
	beta<-Dat2_1$beta

	Pos<-eaf>0.5
	Dat2_1$eaf[Pos]<-1-eaf[Pos]
	Dat2_1$beta[Pos]<-beta[Pos]*-1
	Dat2_1$effect_allele[Pos]<-oa[Pos]
	Dat2_1$other_allele[Pos]<-ea[Pos]
	Dat2_1[,c("eaf","ref_eaf")]
	# exclude palindromic SNPs with MAF >0.42
	Pos<-Dat2_1$ref_eaf>=0.43 | Dat2_1$eaf>=0.43
	# Dat2_1[which(Pos),]
	# Lun3[Lun3$SNP %in% c("rs7112985","rs7946441"),]
	# Crc2[Crc2$SNP %in% c("rs7112985","rs7946441"),]
	# chs2[chs2$SNP %in% c("rs7112985","rs7946441"),]
	# d5d2[d5d2$SNP %in% c("rs7112985","rs7946441"),]
	# d6d2[d6d2$SNP %in% c("rs7112985","rs7946441"),]
	Dat2_2<-Dat2_1[!Pos,]

	# join non palindromic and palindromic SNPs together
	Har<-rbind(Dat2_2,Dat3)
	# make sure order of SNPs in Har is same as in ref_dat
	Pos<-match(ref_dat$SNP,Har$SNP)
	Pos<-Pos[!is.na(Pos)]
	Har<-Har[Pos,]
	all(ref_dat$SNP[ref_dat$SNP %in% Har$SNP] == Har$SNP)
	return(Har)
}





