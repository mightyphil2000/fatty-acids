# figure 2
setwd("~/mrQC")
library(devtools)
library(CheckSumStats)
load_all()
document()

File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.gz.tab"
Dat<-extract_sig_snps(path_to_target_file=File,p_val_col_number=7)
Dat$id <- ""
Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,dplyr::tibble(rsid=Dat$snp, pval=Dat$p, id=Dat$id),pop="EUR")

gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",efo=c("arachidonic acid measurement","omega-6 polyunsaturated fatty acid measurement"))
gc_list2<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,efo=c("arachidonic acid measurement","omega-6 polyunsaturated fatty acid measurement"))
gc_list3<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,efo="arachidonic acid measurement")

snps<-Dat2$rsid
snps1<-unique(c(gc_list$in_gc,gc_list2$in_gc,gc_list3$in_gc))
snps[!snps %in% snps1]
snps[snps %in% snps1]

snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=TRUE,snplist_user=Dat2$rsid)
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")

Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")

dim(Dat)

Plot1<-make_plot_maf(ref_1000G="EUR",target_dat=Dat,Title_size=10,Title_xaxis_size=10)
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat,Title_size=10,Title_xaxis_size=8)
Plot3<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",beta="beta",se="se",Title_size=10,Title_xaxis_size=10,force_all_trait_study_hits=TRUE)
Dat<-predict_beta_sd(dat=Dat)
Plot4<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE,Title_size=10,Title_xaxis_size=10)
Plot5<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=TRUE,Title_size=10,Title_xaxis_size=10)
Plot6<-zz_plot(dat=Dat,beta="beta",se="se",Title_size=10,Title_xaxis_size=10)
# Plot6

Plot_list2<-c("Plot1","Plot2","Plot3", "Plot4","Plot5","Plot6")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report2.png")

