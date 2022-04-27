# devtools::install_github("explodecomputer/genetics.binaRies")
# genetics.binaRies::get_plink_binary()

# library(dplyr)

################
# clump eqtlgen#
################

library(ieugwasr)
File1<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox1_ptgs1_sig_eqtlgen.txt"
File2<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox2_ptgs2_sig_eqtlgen.txt"
Cox1<-read.table(File1,sep="\t",head=T,stringsAsFactors=F)
Cox1$id<-"Cox1"
Cox2<-read.table(File2,sep="\t",head=T,stringsAsFactors=F)
Cox2$id<-"Cox2"
Cox<-rbind(Cox1,Cox2)

Clump_relaxed<-ld_clump( clump_r2 = 0.3,
    dplyr::tibble(rsid=Cox$snp, pval=Cox$p, id=Cox$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict01<-ld_clump( clump_r2 = 0.01,
    dplyr::tibble(rsid=Cox$snp, pval=Cox$p, id=Cox$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict<-ld_clump( clump_r2 = 0.001,
    dplyr::tibble(rsid=Cox$snp, pval=Cox$p, id=Cox$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

write.table(Clump_relaxed,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox_eqtlgen_clump_relaxed.txt",sep="\t",col.names=T,row.names=F,quote=F)

write.table(Clump_strict01,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox_eqtlgen_clump_strict01.txt",sep="\t",col.names=T,row.names=F,quote=F)

write.table(Clump_strict,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox_eqtlgen_clump_strict.txt",sep="\t",col.names=T,row.names=F,quote=F)

Clump_relaxed$clump_r2<-0.3
Clump_strict$clump_r2<-0.001
Clump_strict01$clump_r2<-0.01
Clump<-do.call(rbind,list(Clump_relaxed,Clump_strict,Clump_strict01))

Cox_clump<-merge(Cox,Clump,by.x=c("snp","id"),by.y=c("rsid","id"))
snplist<-unique(Cox_clump$snp)
write.table(Cox_clump,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_eqtlgen_sig_clumped.txt",col.names=T,row.names=F,quote=F)

Temp<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_eqtlgen_sig_clumped.txt",sep=" ",head=T,stringsAsFactors=F)
snplist<-unique(Temp$snp)
write.table(snplist,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_eqtlgen_unique_rsidonly_sig_clumped.txt",col.names=F,row.names=F,quote=F)


#############
# clump GTEx#
#############

/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg38.txt

ref<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg38.txt",sep="\t",head=F,stringsAsFactors=F)
ref1<-ref[ref$V1 %in% c("chr1","chr9"),]

setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/")
Files<-dir()[grep("PTGS[[:digit:]]_COX[[:digit:]]",dir())]
i<-1
Dat_list<-NULL
for(i in 1:length(Files)){
	Dat<-read.table(Files[i],sep="\t",head=T,stringsAsFactors=F)	
	if(nrow(Dat)>0){
		Dat$ID<-Files[i]
		Dat_list[[i]]<-Dat
	}
}

Dat<-do.call(rbind,Dat_list)
head(Dat)
VarID<-unlist(strsplit(Dat$variant_id,split="_"))

# https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/README_eQTL_v8.txt
# gtex ALT allele is effect allele
# variant_id: variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b38

Dat$Chr<-VarID[seq(1,length(VarID),by=5)]
Dat$pos_b38<-VarID[seq(2,length(VarID),by=5)]
Dat$other_allele<-VarID[seq(3,length(VarID),by=5)]
Dat$effect_allele<-VarID[seq(4,length(VarID),by=5)]

Dat2<-merge(Dat,ref1,by.x=c("Chr","pos_b38"),by.y=c("V1","V2"),all.x=TRUE)


Clump_relaxed<-ld_clump( clump_r2 = 0.3,
    dplyr::tibble(rsid=Dat2$V4, pval=Dat2$pval_nominal , id=Dat2$ID),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict01<-ld_clump( clump_r2 = 0.01,
    dplyr::tibble(rsid=Dat2$V4, pval=Dat2$pval_nominal , id=Dat2$ID),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

Clump_strict<-ld_clump( clump_r2 = 0.001,
    dplyr::tibble(rsid=Dat2$V4, pval=Dat2$pval_nominal , id=Dat2$ID),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)



Clump_relaxed$clump_r2<-0.3
Clump_strict$clump_r2<-0.001
Clump_strict01$clump_r2<-0.01
Clump<-do.call(rbind,list(Clump_relaxed,Clump_strict,Clump_strict01))

Cox_clump<-merge(Dat2,Clump,by.x=c("V4","ID"),by.y=c("rsid","id"))
names(Cox_clump)[names(Cox_clump) == "V4"]<-"rsid"
names(Cox_clump)[names(Cox_clump) == "Allele1"]<-"other_allele"
names(Cox_clump)[names(Cox_clump) == "Allele2"]<-"effect_allele"
snplist<-unique(Cox_clump$rsid)

write.table(Cox_clump,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_gtex_sig_clumped.txt",col.names=T,row.names=F,quote=F)
write.table(snplist,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_gtex_unique_rsidonly_sig_clumped.txt",col.names=F,row.names=F,quote=F)



/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_gtex_sig_clumped.txt

