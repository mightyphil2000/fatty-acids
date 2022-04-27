# devtools::install_github("explodecomputer/genetics.binaRies")
# genetics.binaRies::get_plink_binary()
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/GSCAN_incUKBiobank_results/
# library(dplyr)
library(ieugwasr)
Cpd<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/GSCAN_incUKBiobank_results/cpd_no23andme_sig.txt",sep="\t",head=T,stringsAsFactors=F)

Cpd$id <-"cpd"

Clump_relaxed<-ld_clump( clump_r2 = 0.3,
    dplyr::tibble(rsid=Cpd$RSID, pval=Cpd$PVALUE, id=Cpd$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict01<-ld_clump( clump_r2 = 0.01,
    dplyr::tibble(rsid=Cpd$RSID, pval=Cpd$PVALUE, id=Cpd$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict<-ld_clump( clump_r2 = 0.001,
    dplyr::tibble(rsid=Cpd$RSID, pval=Cpd$PVALUE, id=Cpd$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

write.table(Clump_relaxed,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_relaxed.txt",sep="\t",col.names=T,row.names=F,quote=F)

write.table(Clump_strict01,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_strict01.txt",sep="\t",col.names=T,row.names=F,quote=F)

head /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_strict01.txt

wc /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_strict.txt

write.table(Clump_strict,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_strict.txt",sep="\t",col.names=T,row.names=F,quote=F)

snplist<-unique(c(Clump_relaxed$rsid,Clump_strict01$rsid,Clump_strict$rsid))
# Cpd2<-Cpd[Cpd$RSID %in% snplist,c("CHROM","POS","RSID")]

write.table(snplist,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/snplist_rsid_cpd_clumped.txt",col.names=F,row.names=F,quote=F)


Dat1<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_relaxed.txt",sep="\t",head=T,stringsAsFactors=F)
Dat1$clump_r2 <- 0.3

Dat2<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_strict01.txt",sep="\t",head=T,stringsAsFactors=F)
Dat2$clump_r2<-0.01

Dat3<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/gscan_cpd_incUKB_snps_clump_strict.txt",sep="\t",head=T,stringsAsFactors=F)
Dat3$clump_r2<-0.001

Clump<-do.call(rbind,list(Dat1,Dat2,Dat3))
Clump2<-merge(Clump,Cpd,by.x=c("rsid","id"),by.y=c("RSID","id"))

write.table(Clump2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_cpd_sig_clumped.txt",col.names=T,row.names=F,quote=F,sep="\t")

Clump3<-Clump2[Clump2$clump_r2!=0.3,]

write.table(Clump3,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_cpd_sig_clumped_strict.txt",col.names=T,row.names=F,quote=F,sep="\t")

Res<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_cpd_sig_clumped_strict.txt",sep="\t",stringsAsFactors=F,head=T)

write.table(Res[Res$clump_r2==0.001,],"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_cpd_sig_clumped_vstrict.txt",sep="\t",col.names=T,row.names=F,quote=F)


############################
# Clump CSI SNPs UK Biobank#
############################

Dat<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/CSI_results/csi_sig.txt",sep=" ",head=T,stringsAsFactors=F)
Dat$id<-"csi"


Clump_relaxed<-ld_clump( clump_r2 = 0.3,
    dplyr::tibble(rsid=Dat$SNP, pval=Dat$P, id=Dat$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict01<-ld_clump( clump_r2 = 0.01,
    dplyr::tibble(rsid=Dat$SNP, pval=Dat$P, id=Dat$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict<-ld_clump( clump_r2 = 0.001,
    dplyr::tibble(rsid=Dat$SNP, pval=Dat$P, id=Dat$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_relaxed$clump_r2<-0.3
Clump_strict$clump_r2<-0.001
Clump_strict01$clump_r2<-0.01
Clump<-do.call(rbind,list(Clump_relaxed,Clump_strict,Clump_strict01))

Csi_clump<-merge(Dat,Clump,by.x=c("SNP","id"),by.y=c("rsid","id"))

write.table(Csi_clump,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_clumped.txt",col.names=T,row.names=F,quote=F)

Temp<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_clumped.txt",sep=" ",head=TRUE,stringsAsFactors=FALSE)
snplist<-unique(Temp$SNP)
write.table(snplist,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_unique_rsidonly_clumped.txt",col.names=F,row.names=F,quote=F)

wc /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_unique_rsidonly_clumped.txt
head /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_clumped.txt


# Instruments for analysis of consortia summary data
Cpd<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_cpd_sig_clumped_strict.txt",sep="\t",head=T,stringsAsFactors=F)
Csi<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_clumped.txt",sep=" ",head=T,stringsAsFactors=F)
length(unique(Cpd$rsid))
Csi2<-Csi[Csi$clump_r2!=0.3,]
dim(Cpd)
snplist<-unique(c(Cpd$rsid,Csi2$SNP))
length(snplist)
head(Cpd)
Csi3<-Csi2[,c("SNP","CHR","BP","id","clump_r2")]
names(Csi3)<-c("rsid","CHROM","POS","id","clump_r2")
Cpd_csi<-rbind(Cpd[,c("rsid","CHROM","POS","id","clump_r2")],Csi3)
dim(Cpd_csi)
Cpd_csi2<-unique(Cpd_csi[,c("rsid","CHROM","POS")])

head(Cpd_csi2)

# snplist2<-unique(Cpd_csi$rsid)
# snplist2<-snplist2[!snplist2 %in% snplist]

# wtf<-unique(Cpd_csi[Cpd_csi$rsid %in% snplist2,c("rsid")])

# wtf %in% Cpd$RSID
length(snplist)
write.table(Cpd_csi,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cpd_csi.txt",col.names=T,row.names=F,quote=F)
write.table(Cpd_csi2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cpd_csi_coordinates.txt",col.names=T,row.names=F,quote=F)

write.table(snplist,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cpd_csi_snplist.txt",col.names=F,row.names=F,quote=F)

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cpd_csi* .

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/GSCAN_incUKBiobank_results/cpd_csi_coordinates.txt.txt .


scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/GSCAN_incUKBiobank_results/cpd_csi_snplist.txt .
