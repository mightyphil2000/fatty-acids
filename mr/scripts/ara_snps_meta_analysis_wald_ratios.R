source("~/fatty-acids/mr/scripts/mr_functions.R")
library(ggforestplot)
library(ggplot2) 

# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")

load("~/fatty-acids/mr/data/outcome_dat_ara_la.Rdata")
load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")

res_single2$outcome2[res_single2$outcome2=="Basal cell skin carcinoma"]<-"Basal cell carcinoma"
outcomes<-unique(res_single2$outcome2)

res_single2$id_study<-paste(res_single2$outcome2,res_single2$study,res_single2$SNP)
res_single2<-res_single2[order(res_single2$ncase.outcome,decreasing=TRUE),]
res_single2<-res_single2[!duplicated(res_single2$id_study),]

res_single2<-res_single2[!res_single2$SNP %in%  "All - MR Egger"   ,]
res_single2<-res_single2[res_single2$exposure=="Arachidonic acid",]
res_single2<-res_single2[res_single2$population=="European",]
res_single2<-res_single2[res_single2$SNP != "rs3734398",]
res_single2<-res_single2[!is.na(res_single2$b),]
meta_list<-NULL
one_study_list<-NULL
for(i in 1:length(outcomes))
{
	test_dat<-res_single2[res_single2$outcome2==outcomes[i],]
	
	if(length(unique(test_dat$study))>1){
		meta_list[[i]]<-meta_analysis_wr(dat=test_dat)	
	}
	if(length(unique(test_dat$study))==1)
	{
		one_study_list[[i]]<-test_dat
	}
}

meta<-do.call(rbind,meta_list)
one_study<-do.call(rbind,one_study_list)
names(meta)[!names(meta) %in% names(one_study)]
names(one_study)[!names(one_study) %in% names(meta)]
meta$exposure<-"Arachidonic acid"
meta$outcome<-meta$outcome2
one_study$nstudies<-1
one_study$Q.p<-NA
one_study$studies<-one_study$study
# meta$study<-one_study$study
names(one_study)[names(one_study) == "p"]<-"pval"
meta_dat<-plyr::rbind.fill(meta,one_study)


# meta_dat2<-meta_dat[meta_dat$outcome2 %in% c("Overall cancer","Basal cell carcinoma"),]
# meta_dat<-meta_dat[!meta_dat$outcome2 %in% c("Overall cancer","Basal cell carcinoma"),]
res_single2$nstudies<-1
# res_single2[res_single2$outcome2=="Malignant skin cancer"    ,]
# res_single3<-res_single2[res_single2$outcome2 %in% c("Malignant skin cancer","Overall cancer","Basal cell carcinoma"),]

# temp<-plyr::rbind.fill(meta_dat2,res_single3)
# temp$id.temp<-paste(temp$SNP,temp$outcome2)
# temp<-temp[order(as.numeric(temp$ncase.outcome),decreasing=TRUE),]
# temp<-temp[!duplicated(temp$id.temp),]

names(meta)[!names(meta) %in% names(res_single2)]
names(res_single2)[!names(res_single2) %in% names(meta)]

names(res_single2)[names(res_single2) == "p"]<-"pval"
res_single2$studies<-res_single2$study
res_single2$Q.p<-NA

meta_dat<-plyr::rbind.fill(meta_dat,res_single2)
meta_dat$id.temp<-paste(meta_dat$SNP,meta_dat$outcome2)
meta_dat<-meta_dat[order(as.numeric(meta_dat$ncase.outcome),decreasing=TRUE),]
meta_dat<-meta_dat[!duplicated(meta_dat$id.temp),]


# meta[meta$outcome2=="Malignant skin cancer"    ,"SNP"]
save(meta_dat,file="~/fatty-acids/mr/results/results_ara_snps_meta_analysis_wald_ratios.Rdata") #meta analysis of o wald raiot results only plus outcomes where only one study available (e.g. UKB overall cancer excluded non-menlanoma skin cancer and UKB non melanoma skin cancer)
