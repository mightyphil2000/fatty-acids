# install.packages("devtools")
# devtools::install_github("MRCIEU/CheckSumStats")
library(devtools)
library(CheckSumStats)
library(ieugwasr)
load_all()
document()

Temp<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Hu28298293_Zhu26932504/MUFA_and_SFA/SFA_PMID_26932504/18_0.TBL",sep="\t",head=TRUE)
Temp[Temp$MarkerName == "rs1488080",]
# grep -w rs1488080 head /projects/MRC-IEU/users/ph14916/fatty_acids_summary/Hu28298293_Zhu26932504/MUFA_and_SFA/SFA_PMID_26932504/18_0.TBL

Dat<-extract_sig_snps(path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Hu28298293_Zhu26932504/MUFA_and_SFA/SFA_PMID_26932504/18_0.TBL",p_val_col_number=10)
 
Dat$id<-""
Dat2<-ieugwasr::ld_clump( clump_r2 = 0.01,
    dplyr::tibble(rsid=Dat$MarkerName, pval=Dat$P.value, id=Dat$id),pop="EAS")
traits<-c("Saturated fatty acid levels (stearic acid C18:0)","Stearic acid (18:0) levels")

# gc_list2<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid)
gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,trait=traits)

# gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat2$rsid,efo="fatty acid measurement")

Dat<-Dat[Dat$MarkerName %in% Dat2$rsid,]
# trait<-c("Saturated fatty acid levels (stearic acid C18:0)"
snplist<-make_snplist(trait=traits,ref1000G_superpops=TRUE,snplist_user=Dat$MarkerName)
sfa<-extract_snps(snplist=snplist,path_to_target_file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/Hu28298293_Zhu26932504/MUFA_and_SFA/SFA_PMID_26932504/18_0.TBL",path_to_target_file_sep="\t")
Dat<-format_data(dat=sfa,outcome="stearic acid",population="East Asian",study="NHAPC/MESA",ncontrol=3521,UKbiobank=FALSE,rsid="MarkerName",effect_allele="Allele1",other_allele="Allele2",beta="Effect",se="StdErr",eaf="Freq1",p="P.value")

Plot1<-make_plot_maf(ref_1000G="EAS",target_dat=Dat)
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat)

Plot3<-make_plot_gwas_catalog(dat=Dat,trait=traits,beta="beta",se="se")

Dat<-predict_beta_sd(dat=data.frame(Dat))
Plot4<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE)

Plot5<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=TRUE)
Plot3_2<-make_plot_gwas_catalog(dat=Dat,trait=traits,force_all_trait_study_hits=TRUE,beta="beta",se="se")
# Plot3_2
Plot6<-zz_plot(dat=Dat,beta="beta",se="se")

gc_list<-gwas_hit_in_gwas_catalog(gwas_hits=Dat$rsid,trait=traits)

Plot_list2<-c("Plot1","Plot2","Plot3", "Plot4","Plot5","Plot6")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report_nhapc_mesa_chi.png")
