# BiocManager::install("multtest")
library(ggforestplot)
library(ggplot2)
library(TwoSampleMR)
library(devtools)
setwd("~/TwoSampleMR")
load_all()
document()

source("~/fatty-acids/mr/scripts/mr_functions.R")
# exposure
source("~/fatty-acids/mr/scripts/mr_functions.R")
load("~/fatty-acids/mr/data/instruments.Rdata")

exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")

exp_dat<-exp_dat[!exp_dat$FADS,]

#Confirmed one SNP per FADS region per exposure
# Temp<-exp_dat[exp_dat$FADS,c("chr","position","id.exposure","SNP","exposure")]
# Temp[order(Temp$id.exposure),]

# Outcome
# Dat2
load("~/fatty-acids/mr/data/outcome_dat_secondary_pufas.Rdata")

exp_dat[exp_dat$SNP=="rs3798713",]

# harmonise data
dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat =outcome_dat,tolerance=0.05)

# outcome_dat$study[which(outcome_dat$SNP=="rs1741")]

# outcome_dat[outcome_dat$id == 141,"ncase.outcome"]
# res_mre<-mr(dat,method_list=c("mr_ivw_mre"))
# res<-mr(dat,method_list="mr_ivw")
# res_fe<-mr(dat,method_list="mr_ivw_fe")
# if(!all(res_mre$se ==res$se )) warning("underdispersion in effects, fixed effects model recommended in this case")

# Pos<-!res_mre$se ==res$se 
# res_mre[Pos[1],]
# res[Pos[1],]
# res_fe[Pos[1],]
# dat[dat$id.exposure == "exp10" & dat$id.outcome==1,  ]

# # Pos<-which(res_mre$exposure=="Arachidonic acid (20:4n6)")
# res_mre[Pos,]
# res[Pos[1],]
# res_fe[Pos[1],]

# res1<-mr(dat,method_list=c("mr_ivw_mre"))
# res2<-mr(dat,method_list=c("mr_ivw_fe"))
# res3<-mr(dat,method_list=c("mr_wald_ratio"))

# res<-rbind(res1,res2)
# res<-rbind(res,res3)
# res$id.analysis<-paste0(res$id.outcome,res$id.exposure)
# res<-res[!duplicated(res$id.analysis),]
# mr_method_list()$obj

sens<-c("mr_egger_regression","mr_weighted_median","mr_weighted_mode","mr_weighted_mode_nome")

# sens2<-c("mr_raps","mr_ivw_radial")

# exposures with ≥5 instruments and therefore possible to run mr_raps and mr_ivw_radial
# "Docosahexaenoic acid (22:6n3)"
# sen_exp<-c("Omega-3 fatty acids","Linoleic acid (18:2n6)","Omega-6 fatty acids","Other polyunsaturated fatty acids than 18:2","Ratio of omega-6 fatty acids to omega-3 fatty acids","Gamma-linolenic acid (18:3n6)")

# res<-mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_ivw","mr_wald_ratio",sens))
res<-mr(dat,method_list=c("mr_ivw","mr_wald_ratio","mr_ivw_fe"))
# res<-mr(dat[dat$id.outcome == 60,])

res$id.analysis<-paste0(res$id.outcome,res$id.exposure)

# TwoSampleMR::mr_method_list()
pos1<-which(res$nsnp==1)
pos2<-which(res$nsnp==2)
pos3<-which(res$nsnp>2)
res1_1<-res[pos2,] #where nsnps == 2, use ivw fe
res1_1<-res1_1[res1_1$method == "Inverse variance weighted (fixed effects)",]
res1_2<-res[which(res$nsnp!=2),]

res1_2<-subset_on_method(mr_res=res1_2,multi_snp_method="Inverse variance weighted")

res1<-rbind(res1_1,res1_2)
load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
res2<-format_results(res=res1)
# res2all<-format_results(res=res)

# res2$outcome2[res2$id.outcome==32]<-"Overall cancer32" #ovreall cancer FinGen 9792 cases 
# res2$outcome2[res2$id.outcome==33]<-"Overall cancer33" #overall cancer FinGen 31217 cases
res2$id<-paste0(res2$outcome2,res2$study.abbreviation,res2$id.exposure)
# res2$id<-paste0(res2$outcome2,res2$study.abbreviation,res2$id.exposure,res2$method)

dups<-res2$id[duplicated(res2$id)]
# res2[res2$id == dups[1],]
# res2[res2$id %in% dups,]
res2<-res2[order(res2$ncase.outcome,decreasing=TRUE),]
res2<-res2[!duplicated(res2$id),]

res3<-res2[res2$population=="European",]
res3$outcome2[res3$outcome2=="Basal cell skin carcinoma"]<-"Basal cell carcinoma"
res3$id.outcome
res3[res3$exposure == "Alpha-linolenic acid (18:3n3)" & res3$id.outcome %in% c(1,135),]

dat.meta_list<-NULL
exposures<-unique(res3$exposure)
for(i in 1:length(exposures)){
	# i<-12
	print(exposures[i])
	input_dat<-res3[res3$exposure == exposures[i],]
	# dups<-unique(input_dat$outcome2[duplicated(input_dat$outcome2)])
	print(input_dat[order(input_dat$outcome2),c("outcome2","nsnp","study.abbreviation")])
	dat.meta_list[[i]]<-meta_analysis3(dat=input_dat,exclude_finngen=TRUE)
} 

dat.meta<-do.call(rbind,dat.meta_list)

res4<-format_results4(dat=res3,dat.meta=dat.meta,exclude_finngen=TRUE)

individual_studies_sensitvity_analysis(dat=res3,exclude_finngen=TRUE,excl_fads=TRUE)

individual_studies_sensitvity_analysis(dat=res3,exclude_finngen=FALSE,excl_fads=TRUE)

dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat =outcome_dat,tolerance=0.05)
het<-mr_heterogeneity(dat)

het_fisherp<-meta_p(dat=het)
res_single <- mr_singlesnp(dat)
res_single2<-format_results(res=res_single)
# res_single2[res_single2$outcome2=="Colorectal cancer",]
# res_single[res_single$outcome=="Colorectal cancer | 60",]
# res_wald<-mr(dat,method_list="mr_wald_ratio")

# unique(res4$outcome2)
# load("~/fatty-acids/mr/results/res_single_secondary_pufas_exclfads.Rdata")
# save(res2all,file="~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_sensitivityanalyses.Rdata")

# version 3 relaxes tolerance threshold for palindromic SNPs from 0.08 (MAF<0.42 acceptable) to 0.05 (MAF<0.45 acceptable). This allows inclusion of one additional SNP for EPA from outside the FADS region, which was otherwise excluded (in versions 1 and 2). This relaxed threshold also includes an additional SNP from outside FADS region for LA and total omega 6
save(res4,file="~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_v3.Rdata")
save(res_single2,file="~/fatty-acids/mr/results/res_single_secondary_pufas_exclfads_v3.Rdata")
save(het,file="~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity_v3.Rdata")
save(het_fisherp,file="~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity_fisher_p_v3.Rdata")

load("~/fatty-acids/mr/results/res_single_secondary_pufas_exclfads_v3.Rdata")
# write.table(het,"~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
# write.table(res2all,"~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_sensitivityanalyses.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# res2<-format_results()
# head(res2)
# res2$study[which(res2$id.outcome %in% c("32","33","46","54"))]

# res2[res2$outcome2=="Colorectal cancer",]


# # Pos<-grep("East Asian",disc.tab9$population,invert=TRUE)
# ID<-disc.tab9$ID[Pos]



