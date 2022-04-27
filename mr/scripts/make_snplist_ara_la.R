# plan - test associations of arachidonic acid and linoleic acid with selected cancers in univariable and multivariable MR
devtools::install_github("mrcieu/ieugwasr")

la_nmr<-ieugwasr::tophits(id="met-d-LA",clump=0)
# la_nmr <- ieugwasr::associations(id="met-d-LA", variants=snplist,proxies=0)  
la_nmr$id <-"LA"
la_clump<-ieugwasr::ld_clump( clump_r2 = 0.001,
    dplyr::tibble(rsid=la_nmr$rsid, pval=la_nmr$p, id=la_nmr$id),
    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
)

la_nmr_clump<-la_nmr[la_nmr$rsid %in% la_clump$rsid,]
la_snplist<-la_nmr_clump$rsid

data.frame(la_nmr_clump[la_nmr_clump$chr==11,])

 # rs780093 
ara_snplist<-c("rs4985155","rs174546","rs3734398")
snplist<-unique(c(ara_snplist,la_snplist))
class(snplist)
# snplist<-ara_snplist

ara_snplist

write.table(snplist,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/ara_la_snplist_independent.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
