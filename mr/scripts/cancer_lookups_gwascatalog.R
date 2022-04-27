# devtools::install_github("CBIIT/LDlinkR")
library(LDlinkR)
library(gwasrapidd)
library(purrr)
library(dplyr)
# snps<-read.csv("~/fatty-acids/mr/data/regions_chr11_61276027-62392051-associations-2021-01-18.csv",head=T,stringsAsFactors=F,quote="")
# snps<-snps[,1]
# snps1<-snps[grep(" x ",snps)]
# snps1<-trimws(unlist(strsplit(snps1,split=" x ")))
# snps1<-gsub("\"","",snps1)
# snps2<-snps[grep(" x ",snps,invert=T)]
# snps<-c(snps1,snps2)
# snps<-unlist(strsplit(snps,split="-"))
# snps<-trimws(snps[seq(1,length(snps),by=2)])
# snps<-gsub("\"","",snps)
# snps<-snps[!duplicated(snps)]
Variants<-get_variants(genomic_range=list(chromosome="11",start=61276027,end=62392051))
snps<-Variants@variants$variant_id
Assoc<-get_associations(variant_id=snps)
assoc_id<-unique(Assoc@associations$association_id)
studies<-get_studies(association_id=assoc_id)

study_ids <- studies@studies$study_id
    names(study_ids) <- study_ids
    associations <-
      purrr::map(study_ids, ~ get_associations(study_id = .x))
    study2association <-
      purrr::imap_dfr(
        associations,
        ~ tibble::tibble(
          study_id = .y,
          association_id = .x@associations$association_id
        )
      )

IDS<-Assoc@associations$association_id
Dups<-unique(IDS[duplicated(Assoc@associations$association_id)])
unique(data.frame(Assoc@associations[Assoc@associations$association_id %in% Dups ,]))

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

Res<-unique(unlist(lapply(c("cancer","neoplasm","carcinoma","emia","malignancy","malignant"),FUN=function(x) grep(x,Res$reported_trait,ignore.case=T)))) %>% 
	Res[.,]

# Res1<-Res[Res$reported_trait == "Core binding factor acute myeloid leukemia",]
# Res2<-Res[Res$reported_trait != "Core binding factor acute myeloid leukemia",]
Res<-Res[!duplicated(Res$association_id),]
Res$rs174546_grch38<- 61802358
Res$dis_rs174546_grch38<-Res$chromosome_position - Res$rs174546_grch38
Res<-sample_size(Dat=Res)

Anc<-c("European","East Asian")
Pop<-c("EUR",'EAS')


Var2<-unique(Res$variant_id)

Ld_eur<-lapply(Var2,FUN=function(x) LDpair(var1="rs174546",var=x,pop = "EUR",token = "a086cd05a12d",output = "table", file = FALSE)$r2)
Ld_eas<-lapply(Var2,FUN=function(x) LDpair(var1="rs174546",var=x,pop = "EAS",token = "a086cd05a12d",output = "table", file = FALSE)$r2)

r2_table<-data.frame(cbind("r2_eur"=unlist(Ld_eur),"r2_eas"=unlist(Ld_eas)) %>% cbind("variant_id"=unlist(Var2))) %>% 
	 dplyr::left_join(Res,by="variant_id")


Res1<-r2_table[,c("variant_id","risk_allele","risk_frequency","or_per_copy_number", "beta_number","standard_error","range","pvalue","reported_trait","chromosome_name","chromosome_position","dis_rs174546_grch38","ancestral_group","type","r2_eas","r2_eur","initial_sample_size","N","study_id","pubmed_id")]

save.image(file="~/fatty-acids/mr/results/fads_region_lookups_gwascatalog.Rdata")

write.table(Res1,"~/fatty-acids/mr/results/cancer_lookups_gwascatalog.txt",sep="\t",col.names=T,row.names=F,quote=F)


# rs174549 hit for laryngeal cancer
# LDpair("rs174546", "rs174549", pop = "CHB", token = "a086cd05a12d", output = "table", file = FALSE)

LDpair("rs174546", var2, pop = "EUR", token = "a086cd05a12d", output = "table", file = FALSE)

LDpair("rs174546", "rs7937840", pop = "EAS", token = "a086cd05a12d", output = "table", file = FALSE)



sample_size<-function(Dat=NULL){
	Dat$initial_sample_size<-gsub(",","",Dat$initial_sample_size)
	N<-regmatches(Dat$initial_sample_size,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",Dat$initial_sample_size))
	Dat$N<-unlist(lapply(1:length(N),FUN=function(i) sum(as.numeric(unlist(N[i])))))
	Dat<-Dat[order(Dat$N,decreasing=T),]
	# Dat[,c("N","initial_sample_size")]
	return(Dat)
}

harmonise_cancer_names<-function(){		
	Dat$cancer<-NA
	Dat$cancer[unlist(lapply(c("colorectal","colon","rectal"),FUN=function(i) grep(i,Dat$reported_trait,ignore.case=TRUE)))]<-"Colorectal cancer"
	Dat$cancer[grep("breast",Dat$reported_trait,ignore.case=TRUE)]<-"Breast cancer"
	Dat$cancer[grep("lung",Dat$reported_trait,ignore.case=TRUE)]<-"Lung cancer"
	Dat$cancer[grep("prostate",Dat$reported_trait,ignore.case=TRUE)]<-"Prostate cancer"
	Dat$cancer[grep("esophageal",Dat$reported_trait,ignore.case=TRUE)]<-"Esophageal cancer"
	Dat$cancer[grep("pancreatic",Dat$reported_trait,ignore.case=TRUE)]<-"Pancreatic cancer"
	Dat$cancer[grep("diffuse large B-cell lymphoma",Dat$reported_trait,ignore.case=TRUE)]<-"Non-Hodgkin's lymphoma"
	Dat$cancer[grep("multiple myeloma",Dat$reported_trait,ignore.case=TRUE)]<-"Multiple myeloma"
	Dat$cancer[grep("osteosarcoma",Dat$reported_trait,ignore.case=TRUE)]<-"Osteosarcoma"
	Dat$cancer[grep("head and neck",Dat$reported_trait,ignore.case=TRUE)]<-"Head and neck cancer"
	Dat$cancer[grep("hepatocellular carcinoma",Dat$reported_trait,ignore.case=TRUE)]<-"Liver cancer"
	Dat$cancer[grep("ovarian",Dat$reported_trait,ignore.case=TRUE)]<-"Ovarian cancer"
	Dat$cancer[grep("Cholangiocarcinoma",Dat$reported_trait,ignore.case=TRUE)]<-"Bile duct cancer"
	Dat<-Dat[!is.na(Dat$cancer),]
	return(Dat)
}

sample_size<-function(){
	Dat$initial_sample_size<-gsub(",","",Dat$initial_sample_size)
	N<-regmatches(Dat$initial_sample_size,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",Dat$initial_sample_size))
	Dat$N<-unlist(lapply(1:length(N),FUN=function(i) sum(as.numeric(unlist(N[i])))))
	Dat<-Dat[order(Dat$N,decreasing=T),]
	return(Dat)
}


make_table<-function(){
	ancestries <-dplyr::left_join(gwas_studies@ancestries[gwas_studies@ancestries$type == "initial",],
	                       gwas_studies@ancestral_groups,
	                       by = c('study_id', 'ancestry_id'))

	Dat<-dplyr::left_join(gwas_studies@studies,ancestries,
	                       by = c('study_id')) %>% 
	   dplyr::left_join(gwas_studies@publications,
	                       by = c('study_id'))
	   return(Dat)
}

collapse_ancestry<-function(){
	Dat$ancestry<-NA
	Study_ids<-unique(Dat$study_id)
	for(i in 1:length(Study_ids)){
		Pos<-which(Dat$study_id == Study_ids[i])
		Dat$ancestry[Pos]<-paste(Dat$ancestral_group[Pos],collapse="; ")
	}
	Dat<-Dat[!duplicated(Dat$study_id),]
	Dat<-Dat[,!names(Dat) %in% c("ancestral_group","ancestry_id","type")]
	return(Dat)
}

collapse_cancer<-function(){
	Cancers<-unique(Dat1$cancer)
	Dat_list<-NULL
	for(i in 1:length(Cancers)){
		Dat2<-Dat1[Dat1$cancer == Cancers[i],]
		N_studies<-nrow(Dat2)
		Dat2<-Dat2[!duplicated(Dat2$cancer),]
		Dat2$N_studies<-N_studies
		Dat_list[[i]]<-	Dat2	
	}
	Dat3<-tibble(do.call(rbind,Dat_list))	
	return(Dat3)
}


