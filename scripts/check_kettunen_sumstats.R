# devtools::install_github("MRCIEU/CheckSumStats")
library(devtools)
library(CheckSumStats)
library(ieugwasr)
# load_all()
# document()
# check()

ao<-gwasinfo()
ID<-ao$id[ao$trait=="18:2, linoleic acid (LA)"]
# ao$pmid[ao$trait=="18:2, linoleic acid (LA)"]
# sort(unique(ao$trait[ao$pmid==27005778]))

traits<-c("Metabolite levels (lipid measures)",	"Metabolite levels (small molecules and protein measures)",	"Metabolite levels (lipoprotein measures)")
# traits<-"Other polyunsaturated fatty acids than 18:2"  
# ID<-ao$id[ao$trait==traits]
# efo_id<-c("EFO_0006807","EFO_0004529")

top_hits<-tophits(id=ID,r2=0.01,pval=2.3e-9,force_server = TRUE )
# top_hits<-tophits(id=ID,r2=0.01,pval=5e-8,force_server = TRUE )
length(top_hits)
gc_list<-find_hits_in_gwas_catalog(gwas_hits=top_hits$rsid,trait=traits)
# top_hits[,c("rsid","chr")]

# puf <- ieugwasr::associations(id=ID, variants=top_hits$rsid,proxies=0)  
# snplist<-make_snplist(efo_id="EFO_0006807",ref1000G_superpops=TRUE,snplist_user=top_hits$rsid)
snplist<-make_snplist(trait=traits,ref1000G_superpops=TRUE,snplist_user=top_hits$rsid)
puf <- ieugwasr::associations(id=ID, variants=snplist,proxies=0)  
Dat<-format_data(dat=puf,outcome="linoleic acid",population="European",study="Kettunen et al",ncontrol=13524,UKbiobank=FALSE,rsid="rsid",effect_allele="ea",other_allele="nea",beta="beta",se="se",eaf="eaf",p="p")

gc_list<-find_hits_in_gwas_catalog(gwas_hits=Dat$rsid,trait=traits)

Plot1<-make_plot_maf(ref_1000G="EUR",target_dat=Dat)
Plot2<-make_plot_maf(ref_1000G=c("AFR","AMR","EAS","EUR","SAS","ALL"),target_dat=Dat)
Plot3<-make_plot_gwas_catalog(dat=Dat,efo_id="EFO_0006807",beta="beta",se="se")

Dat<-predict_beta_sd(dat=data.frame(Dat))
Plot4<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=FALSE,sd_est=NULL)
Plot5<-make_plot_pred_effect(dat=Dat,pred_beta = "beta_sd",pred_beta_se="se_sd",beta="beta",se="se",bias=TRUE,sd_est=NULL)
# Plot3_2<-make_plot_gwas_catalog(dat=Dat,efo_id="EFO_0006807",force_all_trait_study_hits=TRUE,beta="beta",se="se")
Plot3_2<-make_plot_gwas_catalog(dat=Dat,trait="Metabolite levels (lipid measures)",force_all_trait_study_hits=TRUE,beta="beta",se="se")
# Plot3_2

Plot6<-zz_plot(dat=Dat,beta="beta",se="se")


Plot_list2<-c("Plot1","Plot2","Plot3", "Plot4","Plot5","Plot6")
Plot_list<-lapply(1:length(Plot_list2),FUN=function(x) eval(parse(text=Plot_list2[x])))
combine_plots(Plot_list=Plot_list,out_file="~/qc_report_kettunen.png")
