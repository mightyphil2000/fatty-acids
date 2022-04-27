cox_snps<-c("rs174546","rs62575596","rs689470","rs115693689","rs71635405","rs113455966")
#cox SNPs plus rs174546 which is a FADS SNP
write.table(cox_snps,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox_snps.txt",col.names=FALSE,row.names=FALSE,quote=F)

