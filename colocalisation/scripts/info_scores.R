
# Extract the info scores for summary-data-based imputations in the fatty acid GWAS summary data
# linoleic acid (18:2n6)
# gamma-linolenic acid (18:3n6)
# dihomo-gamma-linolenic acid (20:3n6)
# arachidonic acid (20:4n6)
# scp ph14916@bluecrystalp3.acrc.bris.ac.uk:/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/N6meta* . 
ref<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",sep="\t",head=F,stringsAsFactors=F)

setwd("/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/")
Files<-c("N6meta2041.tbl.fixed.tab","N6meta2031.tbl.fixed.tab","N6meta1821.tbl.fixed.tab", "N6meta1831.tbl.fixed.tab"  )
# File<-Files[1]

system("wc N6meta1831.tbl.fixed.tab")
snps_file1<-snps_pass_info(File=Files[1])
snps_file2<-snps_pass_info(File=Files[2])
snps_file3<-snps_pass_info(File=Files[3])
snps_file4<-snps_pass_info(File=Files[4])

outfile1<-paste0("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/",Files[1],"_snps_pass_info_score90.txt")
outfile2<-paste0("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/",Files[2],"_snps_pass_info_score90.txt")
outfile3<-paste0("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/",Files[3],"_snps_pass_info_score90.txt")
outfile4<-paste0("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/",Files[4],"_snps_pass_info_score90.txt")

write.table(snps_file1,outfile1,col.names=F,row.names=F,quote=F)
write.table(snps_file2,outfile2,col.names=F,row.names=F,quote=F)
write.table(snps_file3,outfile3,col.names=F,row.names=F,quote=F)
write.table(snps_file4,outfile4,col.names=F,row.names=F,quote=F)

snps_pass_info<-function(File=NULL){
	Dat<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
	Dat1<-Dat[Dat$snp!=".",]
	# Dat1[Dat1$V2=="rs174546",]
	# ref[ref$V2=="rs174546",]
	Dat2<-Dat[Dat$snp==".",]
	Dat3<-merge(Dat2,ref,by.x=c("chr","bp"),by.y=c("V1","V4"))
	snps1<-Dat1$snp[Dat1$info>=0.8]
	snps2<-Dat3$V2[Dat3$info>=0.8]
	snps<-c(snps1,snps2)
	return(snps)
}

cd ~/fatty-acids/colocalisation/data
scp ph14916@bluecrystalp3.acrc.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/*snps_pass_info_score* . 