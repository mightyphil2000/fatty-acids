library(devtools)
library(CheckSumStats)
library(ieugwasr)
load_all()
document()

# File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Hu28298293_Zhu26932504/MUFA_and_SFA/SFA_PMID_26932504/18_0.TBL"
File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.gz.tab"

File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.fixed.tab"

rs12285167
rs760306
rs259874
rs472031

rs3741259
rs10026364

Dat<-extract_sig_snps(path_to_target_file=File,p_val_col_number=7)
Dat$id <- ""
Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,dplyr::tibble(rsid=Dat$snp, pval=Dat$p, id=Dat$id),pop="EUR")

# gwas_hits<-c("rs174528","rs16829840","rs10026364","rs17077488","rs10488885", "rs10819512","rs10938476","rs10986188")
gc_list2<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid)
gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)")
gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,efo="omega-6 polyunsaturated fatty acid measurement")
gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,efo="arachidonic acid measurement")


length(gc_list$in_gc)
length(gc_list$not_in_gc)
ensembl[ensembl$refsnp_id %in%  unlist(gc_list$in_gc),] 
# trait=c("Saturated fatty acid levels (stearic acid C18:0)","Stearic acid (18:0) levels")

snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=TRUE,snplist_user=Dat2$rsid)

puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")

Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")

Plot1<-make_plot_maf(ref_1000G="EUR",target_dat=Dat,Title_size=10,Title_xaxis_size=10)
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat,Title_size=10,Title_xaxis_size=6)

Plot3<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",beta="beta",se="se",Title_size=10,Title_xaxis_size=10)

Plot3_2<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",force_all_trait_study_hits=TRUE,beta="beta",se="se",Title_size=10,Title_xaxis_size=10)

