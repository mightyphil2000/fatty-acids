library(TwoSampleMR)
# mv_extract_exposures_local
# mv_extract_exposures
source("~/fatty-acids/mr/scripts/mr_functions.R")
# exposure
load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")

# Outcome
# Dat2
load("~/fatty-acids/mr/data/outcome_dat_ara_la.Rdata")



# harmonise data
dat <- harmonise_data(exposure_dat = exposure_dat,outcome_dat =outcome_dat)

res_mre<-mr(dat,method_list=c("mr_ivw_mre"))
res<-mr(dat,method_list="mr_ivw")
if(!all(res_mre$b ==res$b )) warning("underdispersion in effects, fixed effects model recommended in this case")

# res1<-mr(dat,method_list=c("mr_ivw_mre"))
# res2<-mr(dat,method_list=c("mr_ivw_fe"))
# res3<-mr(dat,method_list=c("mr_wald_ratio"))

# res<-rbind(res1,res2)
# res<-rbind(res,res3)
# res$id.analysis<-paste0(res$id.outcome,res$id.exposure)
# res<-res[!duplicated(res$id.analysis),]

res<-mr(dat,method_list=c("mr_ivw_mre","mr_wald_ratio"))
res<-subset_on_method(mr_res=res,    multi_snp_method="Inverse variance weighted (multiplicative random effects)")
res$id.analysis<-paste0(res$id.outcome,res$id.exposure)
load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
res2<-format_results(res=res)
res2$id<-paste0(res2$outcome2,res2$study.abbreviation,res2$id.exposure)

dups<-res2$id[duplicated(res2$id)]
res2[res2$id %in% dups,]
res2<-res2[order(res2$ncase.outcome,decreasing=TRUE),]
res2<-res2[!duplicated(res2$id),]
res3<-res2[res2$population=="European",]
res3$outcome2[res3$outcome2=="Basal cell skin carcinoma"]<-"Basal cell carcinoma"
ara<-res3[res3$exposure == "Arachidonic acid",]
dat.meta1<-meta_analysis3(dat=ara)
la<-res3[res3$exposure == "Linoleic acid",]
any(duplicated(la$outcome2))
la$outcome2[duplicated(la$outcome2)]
dat.meta2<-meta_analysis3(dat=la)
dat.meta<-rbind(dat.meta1,dat.meta2)
res4<-format_results4(dat=res3)

dat <- harmonise_data(exposure_dat = exposure_dat,outcome_dat =outcome_dat)
het<-mr_heterogeneity(dat)
res_single <- mr_singlesnp(dat)
res_single2<-format_results(res=res_single)
# res_single2[res_single2$outcome2=="Colorectal cancer",]
# res_single[res_single$outcome=="Colorectal cancer | 60",]
# res_wald<-mr(dat,method_list="mr_wald_ratio")
res4[res4$exposure=="Arachidonic acid",c("exposure","outcome2","b","se","pval","study")]
res4[res4$exposure=="Linoleic acid",c("exposure","outcome2","b","se","pval","study")]
# unique(res4$outcome2)
load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
save(res4,file="~/fatty-acids/mr/results/mr_ara_la.Rdata")
save(res_single2,file="~/fatty-acids/mr/results/res_single_ara_la.Rdata")


# res2<-format_results()
# head(res2)
# res2$study[which(res2$id.outcome %in% c("32","33","46","54"))]

# res2[res2$outcome2=="Colorectal cancer",]


# # Pos<-grep("East Asian",disc.tab9$population,invert=TRUE)
# ID<-disc.tab9$ID[Pos]



