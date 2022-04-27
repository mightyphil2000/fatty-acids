
library(ggforestplot)
library(ggplot2)

library(plyr)
# source("~/fatty-acids-mr/instruments/Extract_SNPs_function.R")
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
snplist<-c("rs7937840","rs174546","rs2524299","rs7937840", "rs2277283","rs10897266")

##############
# lung cancer#
##############

dat <- ieugwasr::associations(id="ieu-a-987",variants=snplist,proxies=0)  #ILCCO / ID=75
dat<-format_data(Dat=data.frame(dat,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid=28604730,ncase=29266,ncontrol=56450,study="ILCCO",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=75,open_gwas=TRUE,efo = "lung carcinoma")
ilc<-dat

bbj<-ieugwasr::associations(id= "bbj-a-133", variants=snplist)  #BJ / ID = 17
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Lung cancer",population="East Asian",pmid=32514122,study="BJ",ncase=4050,ncontrol=208403,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=17,open_gwas=TRUE,efo="lung carcinoma")
bbj<-dat

fin <- ieugwasr::associations(id="finn-a-LUNG_CANCER_MESOT", variants=snplist,proxies=0)  # FinnGen / ID = 42
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Lung cancer",population="European",pmid="finn-a-LUNG_CANCER_MESOT",ncase=673,ncontrol=95826,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=42,open_gwas=TRUE,efo="lung carcinoma")
fin<-dat
ukb<-read.table("~/fatty-acids/mr/results/rs7937840_INCENP_luc149.txt",sep="\t",head=T,stringsAsFactors=F)
Lun<-do.call(rbind.fill,list(ilc,bbj,fin,ukb))

##########################################
# esophageal squamous cell carcinoma
##########################################

bbj<-ieugwasr::associations(id= "bbj-a-117", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=32514122,study="BJ",ncase=1300,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=14,open_gwas=TRUE,efo="esophageal squamous cell carcinoma")
bbj<-dat
ugi<-read.table("~/fatty-acids/mr/results/INCENP_eso99.txt",sep="\t",head=T,stringsAsFactors=F)
Eso<-rbind.fill(bbj,ugi)

##########################################
# respiratory and intrathoracic cancer#######
##########################################

fin <- ieugwasr::associations(id="finn-a-C3_RESPIRATORY_INTRATHORACIC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="finn-a-C3_RESPIRATORY_INTRATHORACIC",ncase=615,ncontrol=95884,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=54,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
fin<-dat
ukb <- ieugwasr::associations(id="ukb-d-C3_RESPIRATORY_INTRATHORACIC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Respiratory and intrathoracic cancer",population="European",pmid="ukb-d-C3_RESPIRATORY_INTRATHORACIC",ncase=1944,ncontrol=359250,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=160,open_gwas=TRUE,efo=c("mouth neoplasm","oral cavity cancer","laryngeal squamous cell carcinoma","lung carcinoma","malignant pleural mesothelioma","pharynx cancer","hypopharynx cancer"))
ukb<-transform_betas(dat=dat)

Rip<-rbind.fill(fin,ukb)


###########################
# non melanoma skin cancer###
##############################

Nmc<-read.table("~/fatty-acids/mr/results/INCENP_nmc156.txt",sep="\t",head=T,stringsAsFactors=F)



############################
# basal cell carcinoma#######
############################
Tme<-read.table("~/fatty-acids/mr/results/INCENP_bcc1.txt",sep="\t",head=T,stringsAsFactors=F)
Tme$outcome<-"Basal cell carcinoma"
ukb <- ieugwasr::associations(id="ukb-b-8837", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Basal cell carcinoma",population="European",pmid="ukb-b-8837",ncase=4290,ncontrol=458643,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=135,open_gwas=TRUE,efo = "basal cell carcinoma")
dat<-transform_betas(dat=dat)
ukb<-dat
Bcc<-rbind.fill(Tme,ukb)
Bcc$Effect.Allele[Bcc$Effect.Allele==TRUE]<-"T"


############################
# malignant skin cancer#######
############################

fin <- ieugwasr::associations(id="finn-a-C3_SKIN", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="finn-a-C3_SKIN",ncase=895,ncontrol=95604,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=46,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
fin<-dat
ukb <- ieugwasr::associations(id="ukb-d-C3_SKIN", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(ukb,stringsAsFactors=F),outcome="Malignant skin cancer",population="European",pmid="ukb-d-C3_SKIN",ncase=16531,ncontrol=344663,study="UKB",UKbiobank=TRUE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=152,open_gwas=TRUE,efo=c("non-melanoma skin carcinoma","cutaneous squamous cell carcinoma","basal cell carcinoma","cutaneous melanoma"))
dat<-transform_betas(dat=dat)
ukb<-dat
Msc<-rbind.fill(fin,ukb)


##################
#colorectal cancer#
####################

bbj<-ieugwasr::associations(id= "bbj-a-107", variants=snplist) 
dat<-format_data(Dat=data.frame(bbj,stringsAsFactors=F),outcome="Colorectal cancer",population="East Asian",pmid=32514122,study="BJ",ncase=7062,ncontrol=195745,UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",p="p",effect_allele_confirmed=TRUE,eaf="eaf",all_summary_stats=TRUE,ID=12,open_gwas=TRUE,efo="colorectal cancer")
bbj<-dat
fin <- ieugwasr::associations(id="finn-a-CUSTOM_COLORECTAL_CANCER_EXALLC", variants=snplist,proxies=0)  
dat<-format_data(Dat=data.frame(fin,stringsAsFactors=F),outcome="Colorectal cancer",population="European",pmid="finn-a-CUSTOM_COLORECTAL_CANCER_EXALLC",ncase=843,ncontrol=95656,study="FinnGen",UKbiobank=FALSE,rsid="rsid",Effect.Allele="ea",Other.Allele="nea",lnor="beta",se="se",eaf="eaf",p="p",effect_allele_confirmed=TRUE,all_summary_stats=TRUE,ID=36,open_gwas=TRUE,efo= "colorectal cancer")
fin<-dat
ukb<-read.table("~/fatty-acids/mr/results/INCENP_crc143.txt",sep="\t",head=T,stringsAsFactors=F)

Crc<-do.call(rbind.fill,list(bbj,fin,ukb))


Dat<-do.call(rbind.fill,list(Crc,Lun,Eso,Rip,Nmc,Bcc,Msc))
Dat<-Dat[order(Dat$rsid),]

Dat<-harmonise_effect_allele(Dat=Crc[Crc$rsid =="rs7937840",])

Can<-do.call(rbind.fill,list(Lun,Eso,Rip,Nmc,Bcc,Msc,Crc))
Pos<-Can$eaf>0.5
EA<-Can$Effect.Allele[Pos]
OA<-Can$Other.Allele[Pos]
lnor<-Can$lnor[Pos]
eaf<-Can$eaf[Pos]
Can$Effect.Allele[Pos]<-OA
Can$Other.Allele[Pos]<-EA
Can$lnor[Pos]<-lnor*-1
Can$eaf[Pos]<-1-eaf

Met<-do.call(rbind,lapply(unique(Can$outcome),FUN=function(x) meta_analysis2(Can[Can$outcome==x,])))

Met<-Met[Met$rsid %in% c("rs7937840","rs2277283","rs10897266"),c("outcome","rsid","lnor","se","p","study","Q.p","ncase","ncontrol","eaf","Effect.Allele","Other.Allele","nstudies")]
Met<-Met[as.numeric(Met$nstudies)>1,]


Can<-Can[,c("outcome","rsid","lnor","se","p","study","ncase","ncontrol","eaf","Effect.Allele","Other.Allele")]
Can.all<-rbind.fill(Met,Can)
Can.all<-Can.all[Can.all$rsid %in% c("rs7937840","rs2277283","rs10897266"),]


Can.all$lnor<-as.numeric(Can.all$lnor)
Can.all$se<-as.numeric(Can.all$se)
Can.all$eaf<-as.numeric(Can.all$eaf)
Can.all$p<-as.numeric(Can.all$p)

Can.all$name<-paste0(Can.all$rsid,"\nNo. cases=",Can.all$ncase)
Can.all<-Can.all[order(Can.all$ncase,decreasing=T),]
Can.all<-Can.all[order(Can.all$outcome),]

# png("~/fatty-acids/mr/results/plots/inflammation_lookups_opengwas.png", width = 600, height = 800)
	# forestplot(df = Can.all,logodds = TRUE,name=name,
	#   estimate=lnor, se=se,shape=NULL,
	# 	colour = rsid,xlab = "OR (95% CI) per copy of the minor allele")

# mr_res<-format_bystudy(cancer<-"Respiratory and intrathoracic cancer",overall_effect=TRUE)
# mr_res[,c("study.abbreviation","Q.p")]
Can.all$weight<-1/Can.all$se/10
Can.all$shape<-15
Can.all$shape[!is.na(Can.all$nstudies)]<-23
Can1<-Can.all[Can.all$shape ==23,]
Can2<-Can.all[Can.all$shape !=23,]
Can.all<-rbind(Can1,Can2)
Can.all$study[which(as.numeric(Can.all$nstudies)>1)]<-"Overall effect"

plotlist<-NULL
Plotlist2<-NULL 
snps<-unique(Can.all$rsid)
Outcomes<-unique(Can.all$outcome)

# j<-1
for(j in 1:length(Outcomes)){
	plot.dat2<-Can.all[Can.all$outcome == Outcomes[j],]
	for(i in 1:length(snps)){
		plot.dat<-plot.dat2[plot.dat2$rsid == snps[i],]
		plotlist[[i]]<-forestplot(df = plot.dat,logodds = TRUE,name=study,
					  estimate=lnor, se=se,shape=NULL,
						colour = NULL,xlab = "OR (95% CI) per minor allele") +
						geom_point(shape=plot.dat$shape,size=plot.dat$weight,fill=c("black"))+
						theme(plot.title = element_text(size = ""),text = element_text(size=10))+
						ggtitle(snps[i]) 
	}

	
	Plot<-cowplot::plot_grid(plotlist=plotlist,nrow=3,ncol=1)
	# Title<-"Lung cancer"
	Title<-unique(plot.dat$outcome)	
	Title<-gsub("Esophageal squamous cell carcinoma","Esophageal SCC",Title)
	Title_size<-20
	title <- cowplot::ggdraw() + 
					cowplot::draw_label(
						Title,
						# fontface = 'bold',
						fontface = 'plain',
						x = 0,
						hjust = 0,
						size=Title_size)  +
					ggplot2::theme(
					# add margin on the left of the drawing canvas,
					# so title is aligned with left edge of first plot
						plot.margin = ggplot2::margin(0, 0, 0, 7)
						)

	Plotlist2[[j]]<-cowplot::plot_grid(title, Plot,ncol = 1,rel_heights = c(0.05, 1))
}

Plot4<-cowplot::plot_grid(plotlist=Plotlist2[c(1,5,7)],nrow=1,ncol=3)
png("~/fatty-acids/mr/results/plots/incenp_lookups1.png", width = 1200, height = 1500)
	Plot4
dev.off()

Plot4<-cowplot::plot_grid(plotlist=Plotlist2[c(2,3,4,6)],nrow=1,ncol=4)
png("~/fatty-acids/mr/results/plots/incenp_lookups2.png", width = 1200, height = 1500)
	Plot4
dev.off()

Can.all<-rbind.fill(Met,Can)
Can.all<-Can.all[Can.all$rsid %in% c("rs7937840","rs2277283","rs10897266"),]
write.table(Can.all,"~/fatty-acids/mr/results/incenp_lookups.txt",sep="\t",col.names=T,row.names=F,quote=F)
########################
######### on server########
########################

################
# lung cancer########
################

ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_lung_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE) # UKB / ID=149
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Lung cancer",efo="lung carcinoma",ID=149)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/rs7937840_INCENP_luc149.txt",sep="\t",col.names=T,row.names=F,quote=F)

########################################
# esophageal squamous cell carcinoma########
################################################

ugi<-extract_snps(snplist=snplist,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/summary_chr_all.txt",exact_match=TRUE,file_sep="\t",Comment="")
dat<-format_data(Dat=ugi,outcome="Esophageal squamous cell carcinoma",population="East Asian",pmid=25129146,ncase=2013,ncontrol=2701,study="N-UGC",UKbiobank=FALSE,rsid="rs",Effect.Allele="risk_allele",Other.Allele="reference_allele",lnor="beta",se="standard_error_of_beta",eaf="risk_allele_frequency",p="pvalue",effect_allele_confirmed=TRUE,info="info",ID=99,all_summary_stats=TRUE,efo="esophageal squamous cell carcinoma") 
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_eso99.txt",sep="\t",col.names=T,row.names=F,quote=F)


###########################
# non melanoma skin cancer###
##############################

ukb<-extract_snps(File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_nm_skin_cancer_imputed.txt.gz",snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE)
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Non-melanoma skin cancer",efo="non-melanoma skin carcinoma",ID=156)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_nmc156.txt",sep="\t",col.names=T,row.names=F,quote=F)

############################
# basal cell carcinoma#######
############################
dat<-extract_snps_and_format_bcc_23andMe(snplist=snplist)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_bcc1.txt",sep="\t",col.names=T,row.names=F,quote=F)

########################
# Colorectal cancer####
########################

ukb<-extract_snps(snplist=snplist,exact_match=TRUE,file_sep="\t",Test.gz=TRUE,
	File="/projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/maf_overall_colorectal_cancer_imputed.txt.gz")
dat<-format_data(Dat=ukb,population="European",pmid="ukb-cancer",study="UKB",ncase="ncase",ncontrol="ncontrol",UKbiobank=TRUE,rsid="SNP",Effect.Allele="ALLELE1",Other.Allele="ALLELE0",lnor="beta",se="se",p="P_BOLT_LMM_INF",effect_allele_confirmed=TRUE,eaf="A1FREQ",info="INFO",all_summary_stats=TRUE,
	outcome="Colorectal cancer",efo="colorectal cancer",ID=143
	)
write.table(dat,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_crc143.txt",sep="\t",col.names=T,row.names=F,quote=F)

