library(purrr)
library(dplyr)
library(plyr)
library(LDlinkR)

load("~/fatty-acids/mr/results/fads_region_lookups_gwascatalog.Rdata")

Res<-dplyr::left_join(Variants@variants,Assoc@risk_alleles,
			by="variant_id") %>% 
	dplyr::left_join(Assoc@associations,
			by = 'association_id') %>% 
	dplyr::left_join(study2association,
			by="association_id") %>%
   dplyr::left_join(studies@studies,
			by = c('study_id')) %>%
   dplyr::left_join(studies@publications,
			by = c('study_id')) %>%
   dplyr::left_join(studies@ancestries,
			by="study_id") %>% 
   dplyr::left_join(studies@ancestral_groups,
		   	by=c("study_id","ancestry_id"))


Pos<-unlist(lapply(c("autoimmune","inflammation","inflammatory","inflammatory bowel disease", "Crohn’s disease" ,"Crohns disease" , "ulcerative colitis","asthma","esophagitis","eczema","atopic dermatitis","Chronic obstructive pulmonary disease","asbestosis","silicosis","bronchitis","cystitis","bladder inflammation","gingivitis","lichen planus","lichen sclerosus","pancreatitis","Barett’s esophagus","Baretts esophagus","sialadentis","Sjogren","thyroiditis","skin inflammation","psoriasis","hives"),FUN=function(x) grep(x,Res$reported_trait,ignore.case=TRUE)))
Pos<-unique(Pos)
Res<-Res[Pos,]
Res<-sample_size(Dat=Res)
# IDS<-Res$association_id
# DUPS<-unique(IDS[duplicated(IDS)])
# Res$study_id[Res$association_id %in% DUPS]
# data.frame(Res[Res$association_id %in% DUPS,c("association_id","ancestral_group","initial_sample_size","type")])
# Res[Res$association_id %in% DUPS,c("pvalue","association_id","type","ancestral_group")]


Res<-Res[!duplicated(Res$association_id),]

Res$rs174546_grch38<- 61802358
Res$dis_rs174546_grch38<-Res$chromosome_position - Res$rs174546_grch38

Var2<-unique(Res$variant_id)

Ld_eur<-lapply(Var2,FUN=function(x) LDpair(var1="rs174546",var=x,pop = "EUR",token = "a086cd05a12d",output = "table", file = FALSE)) 
Ld_eur<-do.call(rbind.fill,Ld_eur)[,c("r2","corr_alleles")] 

Ld_eas<-lapply(Var2,FUN=function(x) LDpair(var1="rs174546",var=x,pop = "EAS",token = "a086cd05a12d",output = "table", file = FALSE))
Ld_eas<-do.call(rbind.fill,Ld_eas)[,c("r2","corr_alleles")] 

r2_table<-data.frame(cbind(Ld_eur,Ld_eas) %>% cbind("variant_id"=unlist(Var2))) %>% 
	 dplyr::left_join(Res,by="variant_id")
names(r2_table)[names(r2_table) == "r2"]<-"r2_eur"
names(r2_table)[names(r2_table) == "r2.1"]<-"r2_eas"
names(r2_table)[names(r2_table) == "corr_alleles"]<-"corr_alleles_eur"
names(r2_table)[names(r2_table) == "corr_alleles.1"]<-"corr_alleles_eas"

Res1<-r2_table[,c("reported_trait","variant_id","risk_allele","risk_frequency","or_per_copy_number", "beta_number","standard_error","range","pvalue","chromosome_name","chromosome_position","dis_rs174546_grch38","ancestral_group","type","r2_eur","corr_alleles_eur","r2_eas","corr_alleles_eas","initial_sample_size","N","study_id","pubmed_id")]


write.table(Res1,"~/fatty-acids/mr/results/inflammation_lookups_gwascatalog.txt",sep="\t",col.names=T,row.names=F,quote=F)


sample_size<-function(Dat=NULL){
	Dat$initial_sample_size<-gsub(",","",Dat$initial_sample_size)
	N<-regmatches(Dat$initial_sample_size,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",Dat$initial_sample_size))
	Dat$N<-unlist(lapply(1:length(N),FUN=function(i) sum(as.numeric(unlist(N[i])))))
	Dat<-Dat[order(Dat$N,decreasing=T),]
	# Dat[,c("N","initial_sample_size")]
	return(Dat)
}


Res$pubmed_id[which(Res$reported_trait =="Chronic obstructive pulmonary disease-related biomarkers")]


