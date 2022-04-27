source("~/fatty-acids/mr/scripts/mr_functions.R")
library(ggforestplot)
library(ggplot2)

load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
res_single2<-res_single2[res_single2$exposure == "Arachidonic acid",]
res_single2<-res_single2[res_single2$SNP != "All - Inverse variance weighted",]
res_single2<-res_single2[res_single2$SNP != "All - MR Egger",]

# res_single2<-res_single2[res_single2$SNP =="rs3734398",]

res_single2$outcome2[res_single2$outcome2=="Basal cell skin carcinoma"]<-"Basal cell carcinoma"
outcomes<-unique(res_single2$outcome2)

res_single2$id_study<-paste(res_single2$outcome2,res_single2$study,res_single2$SNP)
res_single2<-res_single2[order(res_single2$ncase.outcome,decreasing=TRUE),]
res_single2<-res_single2[!duplicated(res_single2$id_study),]

res_single2<-res_single2[res_single2$population == "European",]
meta_list<-NULL
one_study_list<-NULL
for(i in 1:length(outcomes))
{
	outcomes[i]
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

names(one_study)[names(one_study) == "beta.outcome"]<-"b"
names(one_study)[names(one_study) == "se.outcome"]<-"se"
names(one_study)[names(one_study) == "pval.outcome"]<-"pval"
names(one_study)[names(one_study) == "effect_allele.outcome"]<-"effect_allele"
names(one_study)[names(one_study) == "other_allele.outcome"]<-"other_allele"
names(one_study)[names(one_study) == "eaf.outcome"]<-"eaf"

one_study$nstudies<-1
one_study$studies<-one_study$study
one_study$Q.p<-NA

meta_dat<-plyr::rbind.fill(meta,one_study)

names(meta)[!names(meta) %in% names(res_single2)]
names(res_single2)[!names(res_single2) %in% names(meta)]

names(res_single2)[names(res_single2) == "beta.outcome"]<-"b"
names(res_single2)[names(res_single2) == "se.outcome"]<-"se"
names(res_single2)[names(res_single2) == "pval.outcome"]<-"pval"
names(res_single2)[names(res_single2) == "effect_allele.outcome"]<-"effect_allele"
names(res_single2)[names(res_single2) == "other_allele.outcome"]<-"other_allele"
names(res_single2)[names(res_single2) == "eaf.outcome"]<-"eaf"
res_single2$nstudies<-1
res_single2$studies<-res_single2$study
res_single2$Q.p<-NA
meta_dat2<-plyr::rbind.fill(meta,res_single2)

meta_dat2[meta_dat2$outcome2=="Overall cancer",c("SNP","study")]

save(meta_dat,file="~/fatty-acids/mr/results/results_wald_ratios_ara_meta_analysis_europeans.Rdata") #meta analysis result only plus outcomes where only one study available (e.g. UKB overall cancer excluded non-menlanoma skin cancer and UKB non melanoma skin cancer)
save(meta_dat2,file="~/fatty-acids/mr/results/results_wald_ratios_ara_meta_analysis_europeans2.Rdata") #meta analysis results plus individual study results
# crc<-res_single2[res_single2$outcome2 == "Colorectal cancer",]
# crc[,c("SNP","study","beta.outcome","se.outcome")]


