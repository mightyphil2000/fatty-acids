# library(CheckSumStats)
devtools::load_all()
devtools::document()

load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Tintle25500335/sigSNPs_FHS_RBC_C20_4n6_AGE_SEX_COHORT_ADJUST_filtered_clumped.Rdata")
File<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Tintle25500335/Results/FHS_RBC_C20_4n6_AGE_SEX_COHORT_ADJUST_filtered.txt"

Dat<-extract_sig_snps(path_to_target_file=File,p_val_col_number=4)
Dat$id<-""
Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,
    dplyr::tibble(rsid=Dat$rsid, pval=Dat$pval_SNP, id=Dat$id),pop="EUR")

Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,pval=1e-8,
    dplyr::tibble(rsid=Dat$rsid, pval=Dat$pval_SNP, id=Dat$id),pop="EUR")

gc_list2<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid)
gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,efo="arachidonic acid measurement")


snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=TRUE,snplist_user=Clump_relaxed$rsid)
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t")
Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",study="Framingham",ncontrol=2555,UKbiobank=FALSE,rsid="rsid",effect_allele="Al1",other_allele="Al2",beta="beta_SNP",se="SE_SNP",eaf="EAF",p="pval_SNP")

Plot1<-make_plot_maf(ref_1000G="EUR",target_dat=Dat)
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat)
Plot3<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",beta="beta",se="se")

snplist<-make_snplist(trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",ref1000G_superpops=FALSE,snplist_user=Clump_relaxed$rsid)
puf<-extract_snps(snplist=snplist,path_to_target_file=File,path_to_target_file_sep="\t") 

Dat<-format_data(dat=puf,outcome="arachidonic acid",population="European",study="Framingham",ncontrol=2555,UKbiobank=FALSE,rsid="rsid",effect_allele="Al1",other_allele="Al2",beta="beta_SNP",se="SE_SNP",eaf="EAF",p="pval_SNP")

Dat<-predict_beta_sd(dat=Dat)

Plot4<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE)

Plot5<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=TRUE)
Plot3_2<-make_plot_gwas_catalog(dat=Dat,trait="Plasma omega-6 polyunsaturated fatty acid levels (arachidonic acid)",force_all_trait_study_hits=TRUE,beta="beta",se="se")
# Plot3_2
Plot6<-zz_plot(dat=Dat,beta="beta",se="se")


Plot_list2<-c("Plot1","Plot2","Plot3", "Plot4","Plot5","Plot6")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report_framingham.png")
