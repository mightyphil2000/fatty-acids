######################
# Combine all datasets#
#######################
library(plyr)
source("~/fatty-acids/outcome_data/scripts/functions_combine_and_format_outcomes.R")
# library(meta)
# setwd("~/MR_FattyAcids/data/summary_data")

# four batches: 
# Can_corr - data obtained by correspondence or from manually currated/bespoke GWAS analyses of cancer in UK Biobank
# mrbase - data downloaded from MR-Base in 2019
# ukb_dat1 - ukbiobank summary data downloaded from OpenGWAS in 2020 from MRC IEU automated pipeline
# fin_dat1 - finGen summary data downloaded from OpenGWAS in 2020 (SAIGE)
# merge with meta data

load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
cor.met<-format_meta_dat1()
unique(cor.met[cor.met$Cancer.Group=="Esophageal cancer",c("study.abbreviation")])
mrb.met<-format_meta_dat2()
fin.met<-format_meta_dat3()
ukb.met<-format_meta_dat4()

dat.met<-do.call(rbind.fill,list(cor.met,mrb.met,fin.met,ukb.met))
table(dat.met$set)

dat.met2<-dat.met[dat.met$rsid == "rs174546",]
dat.met3<-dat.met[!dat.met$ID %in% unique(dat.met2$ID),]
unique(dat.met3[,c("cancer","study.abbreviation","population")])
dim(dat.met2)
# dat.met2[which(dat.met2$Cancer.Group == "Esophageal cancer"),c("cancer","study.abbreviation")]
#######################################
# find proxies for missing target SNPS#
#######################################

dat.proxies<-format_proxies(Dat=dat.met3,snp="rs174546")
dat.met4<-rbind.fill(dat.met2,dat.proxies)
dim(dat.met4)
length(which(!is.na(dat.met4$ID)))

######################
# find missing allele#
######################

Dat1<-dat.met4[is.na(dat.met4$Other.Allele),]
Dat2<-dat.met4[!is.na(dat.met4$Other.Allele),]
Dat3<-get_missing_allele(Dat=Dat1,snp="rs174546")
dat.met5<-rbind(Dat2,Dat3)

############################
# find allele coding errors#
############################

Dat1<-dat.met5[!is.na(dat.met5$eaf),]
Dat2<-dat.met5[is.na(dat.met5$eaf),]

allele_errors<-find_allele_error(Dat1=Dat1)
allele_errors

# meta analysis 

dat.met6<-harmonise_ea(dat=dat.met5,ea="T",oa="C")
# unique(dat.met6[,c("Effect.Allele","Other.Allele","eaf","subpopulation")])

dat.meta<-meta_analysis(dat=dat.met6,IDS=disc.tab9$ID)
dat_outcomes_final<-format_final(dat=dat.meta,dat2=dat.met6)

dat_outcomes_final[which(dat_outcomes_final$Cancer.Group == "Esophageal cancer"),c("cancer","study.abbreviation")]

save(dat_outcomes_final,file="~/fatty-acids/mr/data/dat_outcomes_final.Rdata")
write.table(dat_outcomes_final,"~/fatty-acids/mr/data/dat_outcomes_final.txt",sep="\t",col.names=T,row.names=F,quote=F)


Dat<-load_plinkfrq(File<-"fatty_acid_snps_eur.frq",population<-"EUR")