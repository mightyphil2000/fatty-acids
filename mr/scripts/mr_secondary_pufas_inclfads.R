# BiocManager::install("multtest")
library(ggforestplot)
library(ggplot2)
library(TwoSampleMR)
library(devtools)
setwd("~/TwoSampleMR")
load_all()
document()
# check()
# install()


# mv_extract_exposures_local
# mv_extract_exposures
source("~/fatty-acids/mr/scripts/mr_functions.R")
# exposure
source("~/fatty-acids/mr/scripts/mr_functions.R")
load("~/fatty-acids/mr/data/instruments.Rdata")
exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")


exp_dat[exp_dat$exposure=="Eicosapentaenoic acid (20:5n3)",c("id","exposure","SNP","id.exposure","chr","eaf.exposure")]

# rs3798713 is a palindromic SNP associated with EPA with MAF ~0.43 that by default gets dropped. We will relax the threshold to 0.44 to retain the SNP. this results in an additional SNP being included for LA and omega 6 (rs28601761). The MAF looks very similar across across the outcome datasets. 
check_palindromic_snp_threshold()
check_palindromic_snp_outcomedat(snp="rs3798713")

#Confirmed one SNP per FADS region per exposure
# Temp<-exp_dat[exp_dat$FADS,c("chr","position","id.exposure","SNP","exposure")]
# Temp[order(Temp$id.exposure),]

# Outcome
# Dat2
load("~/fatty-acids/mr/data/outcome_dat_secondary_pufas.Rdata")

# harmonise data
dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat =outcome_dat,tolerance = 0.05) #relax tolerance from 0.08 to 0.05
dat<-exclude_ukb135_ala() #not necessary to include this in mr_secondary_pufas_exclfads.R. 

# table(dat$id.exposure)
# table(dat$exposure,dat$id.outcome)
# table(dat$exposure,dat$outcome)
# length(which(dat$id.outcome == "60"))
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

sens<-c("mr_egger_regression","mr_weighted_median","mr_weighted_mode","mr_weighted_mode_nome")

# res<-mr(dat,method_list=c("mr_ivw_fe","mr_ivw_mre","mr_ivw","mr_wald_ratio",sens))
res<-mr(dat,method_list=c("mr_ivw","mr_wald_ratio","mr_ivw_fe"))

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
dups<-res2$id[duplicated(res2$id)]
# res2[res2$id == dups[1],]
# res2[res2$id %in% dups,]
res2<-res2[order(res2$ncase.outcome,decreasing=TRUE),]
res2<-res2[!duplicated(res2$id),]
res3<-res2[res2$population=="European",]
res3$outcome2[res3$outcome2=="Basal cell skin carcinoma"]<-"Basal cell carcinoma"

dat.meta_list<-NULL
exposures<-unique(res3$exposure)
# res3$ncase.outcome[res3$study.abbreviation=="FinnGen"]

# i<-2
for(i in 1:length(exposures)){
	print(exposures[i])
	input_dat<-res3[res3$exposure == exposures[i],]
	dups<-unique(input_dat$outcome2[duplicated(input_dat$outcome2)])	
	print(input_dat[order(input_dat$outcome2),c("outcome2","nsnp","study.abbreviation")])
	dat.meta_list[[i]]<-meta_analysis3(dat=input_dat,exclude_finngen=TRUE)
} 
dat.meta<-do.call(rbind,dat.meta_list)
dat.meta$id.outcome
dat.meta[dat.meta$exposure == "Alpha-linolenic acid (18:3n3)" & dat.meta$id.outcome %in% c("1; 135"),]
# dat.meta[dat.meta$outcome2=="Lung cancer" & dat.meta$exposure=="Linoleic acid (18:2n6)",]

# Temp<-res3[res3$outcome2=="Lung cancer" & res3$study.abbreviation!="FinnGen",c("study.abbreviation","exposure","nsnp")]
# Temp[order(Temp$exposure),]
# which(res3$study.abbreviation=="UKB")


individual_studies_sensitvity_analysis(dat=res3,exclude_finngen=TRUE)
individual_studies_sensitvity_analysis(dat=res3,exclude_finngen=FALSE)

# Temp<-res3[res3$id.outcome %in% c(1,135,70),c("exposure","id.outcome","pval")]
# Temp[order(Temp$exposure),]

res4<-format_results4(dat=res3,dat.meta=dat.meta,exclude_finngen=TRUE)

dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat =outcome_dat,tolerance=0.05)
dat<-exclude_ukb135_ala()
het<-TwoSampleMR::mr_heterogeneity(dat)
het[het$exposure=="Eicosapentaenoic acid (20:5n3)",]
het_fisherp<-meta_p(dat=het)
het_fisherp[het_fisherp$exposure=="Eicosapentaenoic acid (20:5n3)",]
# het_fisherp[het_fisherp$outcome2 == ]

res_single <- mr_singlesnp(dat)
res_single2<-format_results(res=res_single)
# res_single2[res_single2$outcome2=="Colorectal cancer",]
# res_single[res_single$outcome=="Colorectal cancer | 60",]
# res_wald<-mr(dat,method_list="mr_wald_ratio")


load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads_v3.Rdata")
load("~/fatty-acids/mr/results/res_single_secondary_pufas_inclfads_v3.Rdata")
# load("~/fatty-acids/mr/results/res_single_secondary_pufas_inclfads.Rdata")
# version 3 relaxes tolerance threshold for palindromic SNPs from 0.08 (MAF<0.42 allowed) to 0.05 (MAF<0.45 allowed)
save(res4,file="~/fatty-acids/mr/results/mr_secondary_pufas_inclfads_v3.Rdata")
save(res_single2,file="~/fatty-acids/mr/results/res_single_secondary_pufas_inclfads_v3.Rdata")
save(het,file="~/fatty-acids/mr/results/secondary_pufas_inclfads_heterogeneity_v3.Rdata")
save(het_fisherp,file="~/fatty-acids/mr/results/secondary_pufas_inclfads_heterogeneity_fisherp_v3.Rdata")

# res2<-format_results()
# head(res2)
# res2$study[which(res2$id.outcome %in% c("32","33","46","54"))]

# res2[res2$outcome2=="Colorectal cancer",]


# # Pos<-grep("East Asian",disc.tab9$population,invert=TRUE)
# ID<-disc.tab9$ID[Pos]





check_palindromic_snp_threshold<-function(){
	maf<-exp_dat$eaf
	Pos<-maf>0.5
	maf[Pos]<-1-maf[Pos]
	exp_dat$maf<-maf
	alleles<-paste0(exp_dat$effect_allele.exposure,exp_dat$other_allele.exposure)
	Pos<-which(alleles %in% c("GC","CG","AT","TA","cg","gc","at","ta"))

	exp1<-exp_dat[Pos,]
	Pos<-exp1$maf>=0.42 & exp1$maf<0.45
	exp1<-exp1[Pos,c("SNP","exposure","effect_allele.exposure","other_allele.exposure","maf","eaf.exposure")]
	return(exp1)
	# snps == "rs3798713"
	# Pos<-maf>0.41 & maf<0.45
	# snps[Pos]
	# exp_dat[exp_dat$SNP=="rs3798713",c("effect_allele.exposure","other_allele.exposure","eaf.exposure")]
}


check_palindromic_snp_outcomedat<-function(snp="rs28601761"){
	outcome_dat<-outcome_dat[outcome_dat$population=="European",]
	maf<-outcome_dat$eaf.outcome[outcome_dat$SNP==snp]
	maf<-1-maf[maf>0.5]
	sort(maf)
	return(maf)
}
