library(ieugwasr)
library(TwoSampleMR)

load("~/fatty-acids/mr/data/exposure_dat_ara_la_mvmr.RData")  


# id_exposure <- c("ieu-a-299", "ieu-a-300", "ieu-a-302")
# id_outcome <- "ieu-a-7"
# exposure_dat <- mv_extract_exposures(id_exposure)
exposure_dat$exposure
snps<-unique(exposure_dat$SNP)

load("~/fatty-acids/mr/data/outcome_dat_ara_la.Rdata")

# outcome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome)
ids<-unique(outcome_dat$id.outcome)
res_list<-NULL
for(i in 1:length(ids)){
	print(ids[i])
	mvdat <-mv_harmonise_data(exposure_dat, outcome_dat[outcome_dat$id.outcome==ids[i],],harmonise_strictness=2)
	class(mvdat)
	res_list[[i]] <-data.frame(mv_multiple(mvdat,instrument_specific=FALSE,pval_threshold=1))
}
res<-do.call(rbind,res_list)
# res[res$result.pval<0.05,]
res2<-format_results_mvmr(res=res)
res2<-res2[order(res2$ncase.outcome,decreasing=TRUE),]
# res2<-res2[!duplicated(res2$outcome),]
res2$id<-paste(res2$outcome2,res2$study,res2$result.exposure)
res2<-res2[!duplicated(res2$id),]

names(res2)<-gsub("result.","",names(res2))
res2$outcome2[res2$outcome2=="Basal cell skin carcinoma"]<-"Basal cell carcinoma"
dat.meta2<-meta_analysis3(dat=res2[res2$exposure=="Arachidonic acid",])
dat.meta3<-meta_analysis3(dat=res2[res2$exposure=="Linoleic acid",])
dat.meta<-rbind(dat.meta2,dat.meta3)
res2<-res2[res2$population=="European",]
res4<-format_results4(dat=res2)

save(res2,file="~/fatty-acids/mr/results/results_ara_la_mvmr.Rdata")
save(res4,file="~/fatty-acids/mr/results/results_ara_la_mvmr_meta_analysis.Rdata")