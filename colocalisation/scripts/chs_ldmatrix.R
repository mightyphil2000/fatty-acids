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

ukb$SNP<-paste(ukb$V1,":",ukb$V4,sep="")
ukb1<-ukb[ukb$V1 == 11,]
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

dim(chs.ld5)
dim(chs3)
Pos<-which(chs3$rsid %in% c("rs174537","rs174546"))
Pos<-which(chs3$rsid %in% c("rs102275","rs174546"))
chs3[Pos,]
chs.ld5[Pos,Pos]
Temp<-ukb2[ukb2$SNP %in% chs1$SNP,]
all(Temp$SNP == chs1$SNP)
any(Temp$SNP != chs1$SNP)
all(chs1$SNP == Temp$SNP)
any(chs1$SNP  != Temp$SNP)
