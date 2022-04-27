source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/mr/scripts/mr_functions.R")
library(plyr)


ref<-read.table("~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",sep=" ",head=F,stringsAsFactors=F)

# 23andMe basal cell carcinoma
setwd("~/fatty-acids/colocalisation/data/cancer/23andme")
bcc<-read.table("23me_bcc_fads.txt",sep="\t",head=T,stringsAsFactors=F)
# scc<-read.table("23me_scc_fads.txt",sep="\t",head=T,stringsAsFactors=F)
all<-read.table("annotation.txt",sep="\t",head=T,stringsAsFactors=F)
gt<-read.table("gt_annotations.txt",sep="\t",head=T,stringsAsFactors=F)
im<-read.table("im_annotations.txt",sep="\t",head=T,stringsAsFactors=F)

bcc2<-format_23andme2(dat=bcc,ncase=12945, ncontrol=274252,outcome="Basal cell carcinoma",pmid=27539887,ID=1)

setwd("~/MR_FattyAcids/data/summary_data/BCC")

bcc<-read.table("chr11results.txt",sep="\t",head=T,stringsAsFactors=F)
bcc_hardvard<-format_data(Dat=bcc,outcome="Basal cell carcinoma",population="European",pmid=27539887,study="HNMSC",ncase=4242,ncontrol=12802,UKbiobank=FALSE,Effect.Allele="Allele1",Other.Allele="Allele2",lnor="Effect",se="StdErr",eaf="Freq1",p="Pvalue",effect_allele_confirmed=TRUE,ref="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",chr="Chr",pos="Pos",ID=70,Direction="Direction",phet="HetPVal",rsid="SNP")

# squamous cell carcinoma esophagous
escc1<-read.table("~/fatty-acids/colocalisation/data/cancer/bbj_escc_fads.txt",sep=" ",head=T,stringsAsFactors=F)
escc1<-format_data(Dat=escc1,outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=1300,ncontrol=195745,UKbiobank=FALSE,Effect.Allele="Allele2",Other.Allele="Allele1",lnor="BETA",se="SE",eaf="AF_Allele2",p="p.value",info="Rsq",effect_allele_confirmed=TRUE,ref="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",chr="CHR",pos="POS",ID=14)

escc2<-read.table("~/fatty-acids/colocalisation/data/cancer/nci_ugi_escc_fads.txt",sep="\t",head=T,stringsAsFactors=F)
escc2<-format_data(Dat=escc2,outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=25129146,study="N-UGC",
	ncase=2013,ncontrol=2701,UKbiobank=FALSE,rsid="rs",Effect.Allele="risk_allele",Other.Allele="reference_allele",lnor="beta",se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info="info",chr="chromosome",pos="position",ID=99)

escc<-rbind.fill(escc1,escc2)
escc<-escc[escc$info>=0.80,]
escc<-escc[escc$eaf>0.01,]

escc$alleles<-paste(escc$Effect.Allele,escc$Other.Allele,sep="")
escc<-escc[which(nchar(escc$alleles) == 2),]
Pos<-escc$alleles %in% c("GC","CG","TA","AT")
escc1<-escc[Pos,]
escc1<-escc1[which(escc1$eaf <0.40 | escc1$eaf > 0.60), ]
escc2<-escc[!Pos,]
escc<-rbind(escc1,escc2)

# escc$eaf<-c(0.05,0.55)
# escc$Effect.Allele<-c("C","A")
Dat2<-meta_analysis2(Dat=escc)
Dat2<-Dat2[Dat2$nstudies!=1,]
Dat3<-Dat2[!Dat2$alleles %in% c("GC","CG","AT","TA"),]
Dat3$study<-paste0(unique(Dat2$study)," excluding palindromic SNPs")
escc2<-rbind.fill(Dat2,escc)
escc3<-rbind.fill(escc2,Dat3)
escc<-escc3
save(escc,file="~/fatty-acids/colocalisation/data/escc.RData")

#bbj lung cancer
luca_bbj<-read.table("~/fatty-acids/colocalisation/data/cancer/bbj_luca_fads.txt",sep=" ",head=T,stringsAsFactors=F)
luca_bbj<-format_data(Dat=luca_bbj,outcome="Lung cancer",population="East Asian",pmid=32514122,study="BJ",ncase=4050,ncontrol=208403,UKbiobank=FALSE,Effect.Allele="Allele2",Other.Allele="Allele1",lnor="BETA",se="SE",eaf="AF_Allele2",p="p.value",info="Rsq",effect_allele_confirmed=TRUE,ref="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",chr="CHR",pos="POS",ID=17)

#ukb lung cancer
luca_ukb<-read.table("~/fatty-acids/colocalisation/data/cancer/ukb_luca_fads.txt",sep="\t",head=T,stringsAsFactors=F)
luca_ukb<-format_data(Dat=luca_ukb,outcome="Lung cancer",population="European",pmid="unpublished",study="UKB",ncase=2671,ncontrol=372016,UKbiobank=TRUE,Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",eaf="A1FREQ",p="P_BOLT_LMM_INF",info="INFO",effect_allele_confirmed=TRUE,rsid="SNP",chr="CHR",pos="BP",ID=149)

#ukb non melanoma skin cancer
nmskin_ukb<-read.table("~/fatty-acids/colocalisation/data/cancer/ukb_nmskin_fads.txt",sep="\t",head=T,stringsAsFactors=F)
nmskin_ukb<-format_data(Dat=nmskin_ukb,outcome="Malignant non-melanoma skin cancer",population="European",pmid="unpublished",study="UKB",ncase=23694,ncontrol=372016,UKbiobank=TRUE,Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",eaf="A1FREQ",p="P_BOLT_LMM_INF",info="INFO",effect_allele_confirmed=TRUE,rsid="SNP",chr="CHR",pos="BP",ID=156)

save(nmskin_ukb,file="~/fatty-acids/colocalisation/data/nmskin_ukb.RData")


####################################################################
# lung cancer, respiratory cancer and skin cancers from OpenGWAS####
#####################################################################
       
load("~/fatty-acids/colocalisation/data/cancer/cancer_dat.RData")

can_dat1<-cancer_dat[cancer_dat$id=="ukb-d-C3_SKIN",]
can_dat1<-format_data(Dat=can_dat1,outcome="Malignant skin cancer",population="European",pmid="ukb-d-C3_SKIN",study="UKB",
	ncase=16531,ncontrol=344663,UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position",ID=152)
can_dat2<-cancer_dat[cancer_dat$id=="ukb-d-C3_RESPIRATORY_INTRATHORACIC",]
can_dat2<-format_data(Dat=can_dat2,outcome="Respiratory and intrathoracic cancer",population="European",pmid="ukb-d-C3_RESPIRATORY_INTRATHORACIC",study="UKB",
	ncase=1944,ncontrol=359250,UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position",ID=160)

can_dat3<-cancer_dat[cancer_dat$id=="ukb-b-8837" ,]
can_dat3<-format_data(Dat=can_dat3,outcome="Basal cell carcinoma",population="European",pmid="ukb-b-8837",study="UKB",ID=135,ncase=4290,ncontrol=458643,UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat4<-cancer_dat[cancer_dat$id=="finn-a-C3_RESPIRATORY_INTRATHORACIC" ,]
can_dat4<-format_data(Dat=can_dat4,outcome="Respiratory and intrathoracic cancer",population="European",pmid="finn-a-C3_RESPIRATORY_INTRATHORACIC",study="FinnGen",ID=54,ncase=615,ncontrol=95884,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat5<-cancer_dat[cancer_dat$id== "finn-a-C3_SKIN"   ,]
can_dat5<-format_data(Dat=can_dat5,outcome="Malignant skin cancer",population="European",pmid="finn-a-C3_SKIN",study="FinnGen",ID=46,ncase=895,ncontrol=95604,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat6<-cancer_dat[cancer_dat$id== "ieu-a-984"     ,]
can_dat6<-format_data(Dat=can_dat6,outcome="Lung adenocarcinoma",population="European",pmid=28604730,study="ILCCO",ID=72,ncase=11245,ncontrol=54619,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat7<-cancer_dat[cancer_dat$id==  "ieu-a-988"      ,]
can_dat7<-format_data(Dat=can_dat7,outcome="Small cell lung carcinoma",population="European",pmid=28604730,study="ILCCO",ID=78,ncase=2791,ncontrol=20580,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")


can_dat8<-cancer_dat[cancer_dat$id==  "ieu-a-985"      ,]
can_dat8<-format_data(Dat=can_dat8,outcome="Lung cancer in ever smokers",population="European",pmid=28604730,study="ILCCO",ID=76,ncase=23848,ncontrol=16605,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat9<-cancer_dat[cancer_dat$id==  "ieu-a-986",]
can_dat9<-format_data(Dat=can_dat9,outcome="Lung cancer in never smokers",population="European",pmid=28604730,study="ILCCO",ID=77,ncase=2303,ncontrol=6995,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat10<-cancer_dat[cancer_dat$id==   "ieu-a-989" ,]
can_dat10<-format_data(Dat=can_dat10,outcome="Squamous cell lung cancer",population="European",pmid=28604730,study="ILCCO",ID=79,ncase=7704,ncontrol=54763,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat11<-cancer_dat[cancer_dat$id==   "ieu-a-987" ,]
can_dat11<-format_data(Dat=can_dat11,outcome="Lung cancer",population="European",pmid=28604730,study="ILCCO",ID=75,ncase=29863,ncontrol=55586,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,chr="chr",pos="position")

can_dat12<-cancer_dat[cancer_dat$id==   "finn-a-LUNG_CANCER" ,]
can_dat12<-format_data(Dat=can_dat12,outcome="Lung cancer",population="European",pmid="finn-a-LUNG_CANCER",study="FinnGen",ID=42,ncase=673,ncontrol=95826,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,ref="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",chr="chr",pos="position")


# esophageal adenocarcinoma
ea<-read.table("~/fatty-acids/colocalisation/data/cancer/ukb_ea_fads.txt",sep="\t",head=T,stringsAsFactors=F)
ea<-format_data(Dat=ea,outcome="Esophageal adenocarcinoma",population="European",pmid="unpublished",study="UKB",
	ncase=740,ncontrol=372016,UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",eaf="A1FREQ",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,chr="CHR",pos="BP",ID=144)

ea2<-read.table("~/fatty-acids/colocalisation/data/cancer/ea_eas_fads.txt",sep=" ",head=T,stringsAsFactors=F)
ea2<-format_data(Dat=ea2,outcome="Esophageal adenocarcinoma",population="European",pmid=27527254,study="EAS",
	ncase=4112,ncontrol=17159,UKbiobank=FALSE,rsid="MarkerName",Effect.Allele="effect_allele",Other.Allele="non.effect_allele",lnor="Effect",se="StdErr",eaf=NA,p="P.value",Direction="Direction",effect_allele_confirmed=TRUE,chr="CHR",pos="position",ID=24)

# esophageal adenocarcinoma
Dat=list(ea,ea2)
Dat1<-clean_and_combine_data(Dat=Dat)

Dat1<-harmonise_dat(Dat=Dat1,ref_dat=Dat1[Dat1$study == "UKB",],marker="rsid",effect_allele="Effect.Allele",other_allele="Other.Allele",EAF="eaf",effect="lnor",assume_same_strand=TRUE)
Dat1<-Dat1[Dat1$rsid !="rs4963283",]
ea_dat<-meta_analysis2(Dat=Dat1)
ea_dat2<-format_results(Dat=ea_dat)
eac<-do.call(rbind.fill,list(ea_dat2,ea,ea2))
save(eac,file="~/fatty-acids/colocalisation/data/eac.RData")



# lung cancer
can_dat11<-data.frame(can_dat11)
can_dat12<-data.frame(can_dat12)
# can_dat11[can_dat11$rsid == "rs1692122",]
Dat=list(luca_bbj,luca_ukb,can_dat11,can_dat12)
Dat1<-clean_and_combine_data(Dat=Dat)


# can_dat11[can_dat11$rsid == "rs1692122",]
# Dat1[which(Dat1$rsid=="rs1692122"), c("study","Effect.Allele","Other.Allele","eaf","lnor","se","p")]

Dat1<-harmonise_dat(Dat=Dat1,ref_dat=Dat1[Dat1$study == "UKB",],marker="rsid",effect_allele="Effect.Allele",other_allele="Other.Allele",EAF="eaf",effect="lnor",assume_same_strand=FALSE)

Dat2<-Dat1[Dat1$population!="East Asian",]
Dat3<-Dat1[Dat1$study!="FinnGen",]
Dat4<-Dat3[Dat3$population!="East Asian",]

luca_dat<-meta_analysis2(Dat=Dat1)
luca_dat_excleas<-meta_analysis2(Dat=Dat2)
luca_dat_exclFinGenn<-meta_analysis2(Dat=Dat3)
luca_dat_exclFinGenn_excleas<-meta_analysis2(Dat=Dat4)

luca_dat2<-format_results(Dat=luca_dat)
luca_dat_excleas2<-format_results(Dat=luca_dat_excleas)
luca_dat_exclFinGenn2<-format_results(Dat=luca_dat_exclFinGenn)
luca_dat_exclFinGenn_excleas2<-format_results(Dat=luca_dat_exclFinGenn_excleas)
luca<-do.call(rbind.fill,list(luca_dat2,luca_dat_excleas2,luca_dat_exclFinGenn2,luca_dat_exclFinGenn_excleas2,luca_bbj,luca_ukb,can_dat11,can_dat12))
# luca[which(luca$rsid=="rs1692122"), c("study","Effect.Allele","Other.Allele","eaf","lnor","se","p")]
save(luca,file="~/fatty-acids/colocalisation/data/luca.RData")

# Malignant skin cancer
can_dat1<-data.frame(can_dat1)
can_dat5<-data.frame(can_dat5)

Dat=list(can_dat1,can_dat5)
Dat1<-clean_and_combine_data(Dat=Dat)

Dat1<-harmonise_dat(Dat=Dat1,ref_dat=Dat1[Dat1$study == "UKB",],marker="rsid",effect_allele="Effect.Allele",other_allele="Other.Allele",EAF="eaf",effect="lnor",assume_same_strand=FALSE)
mskin_dat<-meta_analysis2(Dat=Dat1)
mskin_dat<-format_results(Dat=mskin_dat)
mskin<-do.call(rbind.fill,list(mskin_dat,can_dat1,can_dat5))
save(mskin,file="~/fatty-acids/colocalisation/data/mskin.RData")

# Respiratory and intrathoracic cancer

Dat=list(data.frame(can_dat2),data.frame(can_dat4))
Dat1<-clean_and_combine_data(Dat=Dat)

Dat1<-harmonise_dat(Dat=Dat1,ref_dat=Dat1[Dat1$study== "UKB",],marker="rsid",effect_allele="Effect.Allele",other_allele="Other.Allele",EAF="eaf",effect="lnor",assume_same_strand=FALSE)
resp_dat<-meta_analysis2(Dat=Dat1)
resp_dat<-format_results(Dat=resp_dat)
resp<-do.call(rbind.fill,list(resp_dat,can_dat2,can_dat4))
save(resp,file="~/fatty-acids/colocalisation/data/resp.RData")

# basal cell carcinoma
# bcc2$study
# can_dat3$study
# bcc_hardvard$study

can_dat3<-data.frame(can_dat3)
Dat=list(bcc2,can_dat3,bcc_hardvard)
Dat1<-clean_and_combine_data(Dat=Dat)

Dat1<-harmonise_dat(Dat=Dat1,ref_dat=Dat1[Dat1$study== "UKB",],marker="rsid",effect_allele="Effect.Allele",other_allele="Other.Allele",EAF="eaf",effect="lnor",assume_same_strand=FALSE)
bcc_dat<-meta_analysis2(Dat=Dat1)
bcc_dat<-format_results(Dat=bcc_dat)
bcc<-do.call(rbind.fill,list(bcc_dat,bcc2,can_dat3,bcc_hardvard))
save(bcc,file="~/fatty-acids/colocalisation/data/bcc.RData")



# ILCCO
can_dat6<-data.frame(can_dat6)
can_dat7<-data.frame(can_dat7)
can_dat8<-data.frame(can_dat8)
can_dat9<-data.frame(can_dat9)
can_dat10<-data.frame(can_dat10)
can_dat11<-data.frame(can_dat11)
Dat=list(can_dat6,can_dat7,can_dat8,can_dat9,can_dat10,can_dat11)
Dat1<-clean_and_combine_data2(Dat=Dat)
ilcco<-Dat1
save(ilcco,file="~/fatty-acids/colocalisation/data/ilcco.RData")



load("~/fatty-acids/colocalisation/data/luca.RData")
luca<-luca[luca$study %in%  c("ILCCO/UKB","BJ") ,]
load("~/fatty-acids/colocalisation/data/resp.RData")
resp<-resp[resp$study == "UKB",]
load("~/fatty-acids/colocalisation/data/bcc.RData")
bcc<-bcc[bcc$study == "UKB/HNMSC/23andMe",]

load("~/fatty-acids/colocalisation/data/mskin.RData")
# mskin$study[mskin$rsid=="rs174546"]
mskin<-mskin[mskin$study == "UKB",]

load("~/fatty-acids/colocalisation/data/escc.RData")
# escc[escc$rsid=="rs174546",c("ncase","study")]
escc<-escc[escc$study == "BJ/N-UGC",]

load("~/fatty-acids/colocalisation/data/nmskin_ukb.RData")

load("~/fatty-acids/colocalisation/data/eac.RData")
eac<-eac[eac$study == "EAS/UKB" ,]

load("~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata")

crc_accc<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt",sep="\t",head=T,stringsAsFactors=F)

cancer_data_colocalisation<-do.call(rbind.fill,list(bcc,escc,eac,luca,mskin,nmskin_ukb,resp,CRC,crc_accc))
save(cancer_data_colocalisation,file="~/fatty-acids/colocalisation/data/cancer_data_colocalisation.RData")

# ilcco_dat<-format_results(Dat=Dat1)
# Dat1<-harmonise_dat(Dat=Dat1,ref_dat=Dat1[Dat1$outcome== "Lung cancer",],marker="rsid",effect_allele="Effect.Allele",other_allele="Other.Allele",EAF="eaf",effect="lnor",assume_same_strand=FALSE)
# ilcco_dat<-meta_analysis2(Dat=Dat1)

load("~/fatty-acids/colocalisation/data/cancer_data_colocalisation.RData")
acc<-preformat_accc_4()
dat<-format_data(Dat=acc,outcome="Colorectal cancer",population="East Asian",pmid=31826910,study="ACCC",ncase=23572,ncontrol=48700,UKbiobank=FALSE,rsid="rsid",Effect.Allele="Allele1",Other.Allele="Allele2",eaf="Freq1",lnor="Effect",se="StdErr",p="P.value",effect_allele_confirmed=TRUE,ID=3,all_summary_stats=FALSE,summary_set="FAregions",Direction="Direction",phet="HetPVal",I2="HetISq",Q="HetChiSq")

cancer_data_colocalisation_v2<-rbind.fill(dat,cancer_data_colocalisation)
save(cancer_data_colocalisation_v2,file="~/fatty-acids/colocalisation/data/cancer_data_colocalisation_v2.RData")

preformat_accc_4<-function(){
	Crc_accc1<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt",sep="\t",stringsAsFactors=F,head=T)	
	Crc_accc2<-read.table("~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_batch2_regions.txt",sep="\t",stringsAsFactors=F,head=T)	
	accc<-rbind(Crc_accc1,Crc_accc2)
	accc<-accc[which(accc$region == "FADS"),]
	
	
	# Dir<-accc$Direction
	# Dir2<-gsub("\\?","",Dir)
	# Nstudies<-nchar(Dir2)
	# accc2<-accc[Nstudies>=median(Nstudies),]
	# accc2<-accc2[accc2$Freq1>0.01,]

	accc2<-accc[which(accc$N >= 69175),] #3rd quartile
	accc2<-accc2[which(accc2$Freq1>0.01),]	

	SNP<-gregexpr(":",accc2$snp)
	Test<-unlist(lapply(1:length(SNP),FUN=function(x)
		length(unlist(SNP[x]))))
	Pos<-which(Test==3)
	accc2<-accc2[Pos,]
	SNP<-unlist(strsplit(accc2$snp,split=":"))
	accc2$chr<-SNP[seq(1,length(SNP),by=4)]
	accc2$bp<-SNP[seq(2,length(SNP),by=4)]
	accc2$chr<-paste0("chr",accc2$chr)
	acc<-find_rsids(dat=accc2,ref_dat=TRUE)
	return(acc)
}


find_rsids<-function(dat=NULL,ref_dat=FALSE){	
	if(ref_dat){
		ref<-read.table("~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt",head=F,stringsAsFactors=F,sep=" ")
	}
	dat2<-merge(dat,ref,by.x=c("chr","bp"),by.y=c("V1","V2"))
	names(dat2)[names(dat2) == "V4"]<-"rsid"
	return(dat2)
}
