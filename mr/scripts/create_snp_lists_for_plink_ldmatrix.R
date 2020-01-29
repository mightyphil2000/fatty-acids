# UK biobank and 1000 genomes in build 37 format
# Create list of SNPs +/- 500kb of index SNP for each region 

# 			GRCh37							GRCh38		
# Region	Chr	start		stop			start		stop
# FADS		11	61043499	62,159,523		61,276,027	62,392,051
# ELOVL2	6	10,480,992	11,544,547		10,480,759	11,544,305
# ELOVL5	6	52,632,196	53,713,947		52,767,398	53,849,179
# SCD		10	101,606,881	102,624,591		99,847,233	100,864,826
# GCKR		2	27,219,709	28,246,554		26,996,839	28,023,684
# PDXDC1	16	14,568,448	15,733,196		14,474,591	15,639,339
# SPTLC3	20	12,489,627	13,647,411		12,508,972	13,669,103

# cd /projects/MRC-IEU/users/ph14916/
# mv /projects/MRC-IEU/scratch/for_Phil/* /projects/MRC-IEU/users/ph14916/ukb_bed

ukb<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",head=F,stringsAsFactors=F)
# kg1<-read.table("/projects/MRC-IEU/users/ph14916/1000genomes/data_maf0.01_rs.bim",head=F,stringsAsFactors=F)


###############################################
# snplist for FADS region (rs174528 index SNP)#
###############################################


Pos<-which(ukb$V1 == 11)
ukb2<-ukb[Pos,]
Pos<-which(ukb2$V4>61043499 & ukb2$V4<62159523)
ukb2<-ukb2[Pos,]
fads<-ukb2$V2
write.table(fads,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs174528_ukb.txt",col.names=F,row.names=F,quote=F)

# Pos<-which(kg1$V1 == 11)
# kg2<-kg1[Pos,]
# Pos<-which(kg2$V4>61043499 & kg2$V4<62159523)
# kg2<-kg2[Pos,]
# fads2<-kg2$V2

# length(fads)-length(fads2). 1333 more SNPs in bim file for UK biobank


##################################################
# snplist for ELOVL2 region (rs3734398 index SNP)#
##################################################

# ukb<-read.table("/projects/MRC-IEU/scratch/for_Phil/UKBB_10K.bim",head=F,stringsAsFactors=F)

# kg1<-read.table("/projects/MRC-IEU/users/ph14916/1000genomes/data_maf0.01_rs.bim",head=F,stringsAsFactors=F)

Pos<-which(ukb$V1 == 6)
ukb2<-ukb[Pos,]
Pos<-which(ukb2$V4>10480992 & ukb2$V4<11544547)
ukb2<-ukb2[Pos,]
elovl2<-ukb2$V2
write.table(elovl2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs3734398_ukb.txt",col.names=F,row.names=F,quote=F)


###################################################
# snplist for ELOVL5 region (rs12210577 index SNP)#
####################################################

Pos<-which(ukb$V1 == 6)
ukb2<-ukb[Pos,]
Pos<-which(ukb2$V4>52632196 & ukb2$V4<53713947)
ukb2<-ukb2[Pos,]
elovl5<-ukb2$V2
write.table(elovl5,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs12210577_ukb.txt",col.names=F,row.names=F,quote=F)


##############################################
# snplist for SCD region (rs603424 index SNP)#
##############################################

Pos<-which(ukb$V1 ==10)
ukb2<-ukb[Pos,]
Pos<-which(ukb2$V4>101606881 & ukb2$V4<102624591)
ukb2<-ukb2[Pos,]
scd<-ukb2$V2
write.table(scd,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs603424_ukb.txt",col.names=F,row.names=F,quote=F)

##############################################
# snplist for GCKR region (rs780093 index SNP)#
##############################################

Pos<-which(ukb$V1 ==2)
ukb2<-ukb[Pos,]
Pos<-which(ukb2$V4>27219709 & ukb2$V4<28246554)
ukb2<-ukb2[Pos,]
gckr<-ukb2$V2
write.table(gckr,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs780093_ukb.txt",col.names=F,row.names=F,quote=F)


# Pos<-which(kg1$V1 == 2)
# kg2<-kg1[Pos,]
# Pos<-which(kg2$V4>27219709 & kg2$V4<28246554)
# kg2<-kg2[Pos,]
# gckr2<-kg2$V2

# length(gckr)-length(gckr2). 1019 more SNPs in bim file for UK biobank


##############################################
# snplist for PDXDC1 region ( rs4985155 index SNP)#
##############################################

Pos<-which(ukb$V1 ==16)
ukb2<-ukb[Pos,]
Pos<-which(ukb2$V4>	14568448 & ukb2$V4<15733196)
ukb2<-ukb2[Pos,]
pdx<-ukb2$V2
write.table(pdx,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs4985155_ukb.txt",col.names=F,row.names=F,quote=F)



##############################################
# snplist for SPTLC3 region (rs680379 index SNP)#
##############################################

Pos<-which(ukb$V1 ==20)
ukb2<-ukb[Pos,]
Pos<-which(ukb2$V4>12489627 & ukb2$V4<13647411)
ukb2<-ukb2[Pos,]
spt<-ukb2$V2
write.table(spt,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs680379_ukb.txt",col.names=F,row.names=F,quote=F)

