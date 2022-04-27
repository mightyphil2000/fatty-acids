# cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/
# setwd(Dir)
# Files<-dir()[grep(".tab",dir())]


Chargei<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"  ,File="/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/N6meta2041.tbl.fixed.tab",exact_match=TRUE,file_sep="\t") 

Chargeg<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"  ,File="/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/N6meta2041.tbl.fixed.tab",exact_match=TRUE,file_sep="\t") 

Dorajoo<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"  ,File="/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c203n6_pooled_allchr_qc1.tab",exact_match=TRUE,file_sep="\t") 
Dorajoo_aa<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt"  ,File="/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c204n6_pooled_allchr_qc1.tab",exact_match=TRUE,file_sep="\t") 


# Dat.m<-merge(Dorajoo,Dorajoo_aa,by="snp")
# plot(Dat.m$effect_allele_freq.y,Dat.m$effect_allele_freq.x)
# Dat.m[which(Dat.m$effect_allele.y != Dat.m$effect_allele.x),]
# Diff<-Dat.m$effect_allele_freq.y-Dat.m$effect_allele_freq.x
# head(Dat.m)
# Dat.m[which(Diff>.2),]

# Dat.m[Dat.m$snp=="rs174546",]

Charge<-Chargei
save(Charge,file="~/ref_dat_charge_imputed.RData")

SCHS<-Dorajoo_aa
save(SCHS,file="~/ref_dat_schs.RData")


# Dat<-merge(Chargeg,Chargei,by="snp")
# plot(Dat$effect_allele_freq.y,Dat$effect_allele_freq.x)
scp  ph14916@bluecrystalp3.acrc.bris.ac.uk:~/ref_dat_charge_imputed.RData  .
scp  ph14916@bluecrystalp3.acrc.bris.ac.uk:~/ref_dat_schs.RData  .

