setwd("/projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/OncoArray_LungCancer/data")

snps<-c("rs7937840","rs174546","rs2524299","rs2072113")
snps<-c("rs2524299","rs2072113")
Res<-read.csv("Onco_TRICL_032116_Overall_Chr11.csv",head=T,stringsAsFactors=F)
Res1<-read.csv("Onco_TRICL_032116_Ever_Chr11.csv",head=T,stringsAsFactors=F)
Res2<-read.csv("Onco_TRICL_032116_Never_Chr11.csv",head=T,stringsAsFactors=F)

Res[Res$rs_number %in% snps,c("rs_number","Pvalue_fixed","N_study")]
Res1[Res1$rs_number %in% snps,c("rs_number","Pvalue_fixed","N_study")]
Res1[Res1$rs_number %in% snps,]
Res2[Res2$rs_number %in% snps,]

grep rs174546 Onco_TRICL_032116_Adeno_Chr11.csv

summary(Res$N_Cases)