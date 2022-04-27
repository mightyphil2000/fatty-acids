library(CheckSumStats)
# install.packages("devtools")
# devtools::install_github("MRCIEU/CheckSumStats")

devtools::load_all()
devtools::document()

# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/Dorajoo26584805/

# load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Dorajoo26584805/sigSNPs_score_c204n6_pooled_allchr_qc1_clumped.Rdata")
# File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Dorajoo26584805/score_c204n6_pooled_allchr_qc1.tab"
File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Dorajoo26584805/beta_log/score_c204n6_pooled_allchr_qc1.tab"
# File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Dorajoo26584805/beta_log/score_c204n6_control_allchr_qc.tab"
# File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Dorajoo26584805/beta_log/score_c204n6_case_allchr_qc.tab"

Dat<-extract_sig_snps(path_to_target_file=File,p_val_col_number=7)
Dat$id<-""
Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,
    dplyr::tibble(rsid=Dat$snp, pval=Dat$p, id=Dat$id),pop="EAS")

# gc_list2<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid)

gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",efo="arachidonic acid measurement")

snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",efo="arachidonic acid measurement",ref1000G_superpops=TRUE,snplist_user=Dat2$rsid)
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")
Dat<-format_data(dat=puf,outcome="arachidonic acid",population="East Asian",study="SCHS",ncontrol=1361,UKbiobank=FALSE,rsid="snp",effect_allele="effect_allele",other_allele="other_allele",beta="beta",se="se",eaf="effect_allele_freq",p="p")

Plot1<-make_plot_maf(ref_1000G="EAS",target_dat=Dat,Title_size=10,Title_xaxis_size=10)
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat,Title_size=10,Title_xaxis_size=8)
Plot3<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",efo="arachidonic acid measurement",beta="beta",se="se",force_all_trait_study_hits=TRUE)
Dat<-predict_beta_sd(dat=Dat)
Plot4<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE)
Plot5<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=TRUE)
Plot6<-zz_plot(dat=Dat,beta="beta",se="se")

Plot_list2<-c("Plot1","Plot2","Plot3", "Plot4","Plot5","Plot6")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report_schs.png")

