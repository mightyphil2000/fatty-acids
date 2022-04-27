File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.gz.tab"
File_cleaned<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/N6meta2041.tbl.fixed.tab"


setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/")
library(CheckSumStats)
load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/sigSNPs_N6meta2041_clumped.Rdata")

snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=TRUE,snplist_user=charge_top_hits)

File<-system.file("extdata", "ara_test_dat.txt", package = "mrQC")
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")
Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")

Plot1<-make_plot_maf(ref_1000G="EUR",target_dat=Dat,Title_size=10,Title_xaxis_size=10)
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat,Title_size=10,Title_xaxis_size=6)
Plot3<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",beta="beta",se="se",Title_size=10,Title_xaxis_size=10)
snplist<-make_snplist(ref1000G_superpops=FALSE,snplist_user=charge_top_hits) 
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t") 
Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")
Dat<-predict_beta_sd(dat=Dat)
Plot4<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE,Title_size=10,Title_xaxis_size=10)

snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=FALSE) 
snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=FALSE,snplist_user=charge_top_hits) 
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t") 
Dat2<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")
Dat2<-predict_beta_sd(dat=Dat2)
Plot4_2<-make_plot_pred_effect(dat=Dat2,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE,Title_size=10,Title_xaxis_size=10)
# Plot4_2

Plot5<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=TRUE,Title_size=10,Title_xaxis_size=10)

snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",snplist_user=charge_top_hits)
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")
Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")
Plot3_2<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",force_all_trait_study_hits=TRUE,beta="beta",se="se",Title_size=10,Title_xaxis_size=10)
# Plot3_2


load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/guan24823311/sigSNPs_N6meta2041_fixed_clumped.Rdata")
File_cleaned<-system.file("extdata", "ara_test_cleaned_dat.txt", package = "mrQC")
data(charge_top_hits_cleaned) 
# charge_top_hits_cleaned<-Clump_relaxed$rsid
snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=FALSE,snplist_user=charge_top_hits_cleaned)
puf<-extract_snps(snplist=snplist,path_to_target_file=File_cleaned,path_to_target_file_sep="\t")
Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",pmid=24823311,study="CHARGE",ncontrol="n",rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")
Plot3_3<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",force_all_trait_study_hits=TRUE,beta="beta",se="se",Title_size=10,Title_xaxis_size=10)
# Plot3_3

Dat<-predict_beta_sd(dat=Dat)
Plot4_3<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",Title_size=10,Title_xaxis_size=10)
Plot5_2<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=TRUE,Title_size=10,Title_xaxis_size=10)
# Plot5_2

Plot6<-zz_plot(dat=Dat,beta="beta",se="se",Title_size=10,Title_xaxis_size=10)
# Plot6

Plot_list2<-c("Plot1","Plot2","Plot3", "Plot4_2","Plot5","Plot6")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report2.png")

Plot_list2<-c("Plot1","Plot2","Plot3_3", "Plot4_2","Plot5_2","Plot6")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report3.png")

png("~/qc_report4.png")
	Plot3_2
dev.off()

png("~/qc_report5.png")
	Plot4_2
dev.off()


png("~/qc_report6.png")
	Plot4_3
dev.off()

png("~/qc_report7.png")
	Plot3_3
dev.off()
