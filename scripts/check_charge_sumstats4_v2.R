# install.packages("devtools")
library(devtools)
devtools::install_github("MRCIEU/CheckSumStats")
library(CheckSumStats)

#supplementary figure 6. Relationship between reported and expected effect sizes for SNPs associated with arachidonic acid in CHARGE, before and after filtering out low quality SNPs 
setwd("~/mrQC")
library(devtools)
library(CheckSumStats)
load_all()
document()

File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.gz.tab"

Dat<-extract_sig_snps(path_to_target_file=File,p_val_col_number=7)
Dat$id <- ""
Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,dplyr::tibble(rsid=Dat$snp, pval=Dat$p, id=Dat$id),pop="EUR")
snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=TRUE,snplist_user=Dat2$rsid)
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")

Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")
Dat<-predict_beta_sd(dat=Dat)
Plot4<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE,Title_size=10,Title_xaxis_size=10)

# file cleaned
File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.fixed.tab"
Dat<-extract_sig_snps(path_to_target_file=File,p_val_col_number=7)
Dat$id <- ""
Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,dplyr::tibble(rsid=Dat$snp, pval=Dat$p, id=Dat$id),pop="EUR")
snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=TRUE,snplist_user=Dat2$rsid)
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")

Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")
Dat<-predict_beta_sd(dat=Dat)
Plot4_cleaned<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE,Title_size=10,Title_xaxis_size=10)

Plot_list2<-c("Plot4","Plot4_cleaned")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report_charge_cleaned.png")

