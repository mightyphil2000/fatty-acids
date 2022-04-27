library(ggforestplot)
library(ggplot2)

load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
load("~/fatty-acids/mr/data/outcome_dat_ara_la.Rdata")
load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")

ara_snps<-exposure_dat$SNP[exposure_dat$exposure=="Arachidonic acid"]
outcome_dat<-outcome_dat[outcome_dat$SNP %in% ara_snps,]

outcome_dat$outcome2[outcome_dat$outcome2=="Basal cell skin carcinoma"]<-"Basal cell carcinoma"
outcomes<-unique(outcome_dat$outcome2)

outcome_dat$id_study<-paste(outcome_dat$outcome2,outcome_dat$study,outcome_dat$SNP)
outcome_dat<-outcome_dat[order(outcome_dat$ncase.outcome,decreasing=TRUE),]
outcome_dat<-outcome_dat[!duplicated(outcome_dat$id_study),]

meta_list<-NULL
one_study_list<-NULL
for(i in 1:length(outcomes))
{
	outcomes[i]
	test_dat<-outcome_dat[outcome_dat$outcome2==outcomes[i],]
	if(length(unique(test_dat$study))>1){
		meta_list[[i]]<-meta_analysis_snp(dat=test_dat)	
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
which(one_study$outcome2=="Colorectal cancer")
head(meta_dat)

names(meta)[!names(meta) %in% names(outcome_dat)]
names(outcome_dat)[!names(outcome_dat) %in% names(meta)]

names(outcome_dat)[names(outcome_dat) == "beta.outcome"]<-"b"
names(outcome_dat)[names(outcome_dat) == "se.outcome"]<-"se"
names(outcome_dat)[names(outcome_dat) == "pval.outcome"]<-"pval"
names(outcome_dat)[names(outcome_dat) == "effect_allele.outcome"]<-"effect_allele"
names(outcome_dat)[names(outcome_dat) == "other_allele.outcome"]<-"other_allele"
names(outcome_dat)[names(outcome_dat) == "eaf.outcome"]<-"eaf"
outcome_dat$nstudies<-1
outcome_dat$studies<-outcome_dat$study
outcome_dat$Q.p<-NA
meta_dat2<-plyr::rbind.fill(meta,outcome_dat)

save(meta_dat,file="~/fatty-acids/mr/results/results_ara_snps_meta_analysis.Rdata") #meta analysis result only plus outcomes where only one study available (e.g. UKB overall cancer excluded non-menlanoma skin cancer and UKB non melanoma skin cancer)
save(meta_dat2,file="~/fatty-acids/mr/results/results_ara_snps_meta_analysis2.Rdata") #meta analysis results plus individual study results
# crc<-outcome_dat[outcome_dat$outcome2 == "Colorectal cancer",]
# crc[,c("SNP","study","beta.outcome","se.outcome")]


