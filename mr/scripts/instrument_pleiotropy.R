setwd("~/fatty-acids-mr/instruments")
source("functions_for_replication.R")
library(ggforestplot)
library(ggplot2)

load("~/fatty-acids-mr/instruments/fatty_acid_instruments_ldc_v3.Rdata") #CHARGE instrumentss
library(TwoSampleMR)
ao<-available_outcomes()
snps<-unique(Dat.ldc$SNP)

ao.rf<-ao[ao$category == "Risk factor",]
ao.ket<-ao[ao$author == "Kettunen",]

snps<-c("rs174546","rs1741","rs881803")
# ieugwasr::get_access_token() 

# ao[ao$trait == "Lung cancer","id"]
# lc_dat <- extract_outcome_data(
#     snps = "rs1051730",
#     outcomes =  "ieu-a-966"
# )

rf_dat <- extract_outcome_data(
    snps = snps,
    outcomes = ao.rf$id
)

# Pos<-grep("smoke",rf_dat$outcome)
# ao$consortium[ao$id %in% c("ieu-a-961","ieu-a-962","ieu-a-963")]
# rf_dat[Pos,c("SNP","beta.outcome","se.outcome","outcome","pval.outcome")]
# names(rf_dat)



ket_dat <- extract_outcome_data(
    snps = snps,
    outcomes = ao.ket$id
)

# i<-which(snps=="rs3734398")
rf_dat2<-merge(rf_dat,ao,by.x="id.outcome",by.y="id")
unique(rf_dat2$subcategory[rf_dat2$SNP=="rs881803"])

Plot1<-format_dat(Dat=rf_dat2[rf_dat2$SNP=="rs174546",])
Plot1<-Plot1[order(abs(Plot1$b_sd),decreasing=TRUE),]

write.table(Plot1[,c("trait","sig","b_sd","lci","uci","id.outcome","samplesize.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome","pval.outcome","population","Category")],paste0(out.file1,".txt"),sep="\t",col.names=TRUE,row.names=FALSE)

Plot2<-format_dat(Dat=rf_dat2[rf_dat2$SNP=="rs1741",])
write.table(Plot2[,c("trait","sig","b_sd","lci","uci","id.outcome","samplesize.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome","pval.outcome","population","Category")],paste0(out.file2,".txt"),sep="\t",col.names=TRUE,row.names=FALSE)

Plot3<-format_dat(Dat=rf_dat2[rf_dat2$SNP=="rs881803",])
write.table(Plot3[,c("trait","sig","b_sd","lci","uci","id.outcome","samplesize.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome","pval.outcome","population","Category")],paste0(out.file3,".txt"),sep="\t",col.names=TRUE,row.names=FALSE)

library(plyr)
Plot4<-do.call(rbind.fill,list(Plot1,Plot2,Plot3))
Plot4<-Plot4[!Plot4$Category %in% c("PUFA desaturase activity","Arachidonic acid","Dihomo-gamma-linolenic acid" ),]


# sort(unique(rf_dat$outcome))
# rf_dat[grep("smoked",rf_dat$outcome,ignore.case=T),]

# Plot1[Plot1$pval.outcome<0.05/36,"trait"]
# Plot2[Plot2$pval.outcome<0.05/36,"trait"]
# Plot3[Plot3$pval.outcome<0.05/36,"trait"]
dim(Plot2)
out.file1<-paste("~/fatty-acids/mr/results/rf_pleiotropy_",unique(Plot1$SNP),sep="")
Title.plot<-unique(Plot1$SNP)
table(Plot1[,"Pval<0.0014"])

Plot.pdf1<-forestplot(df = Plot1,
	name=trait,
	  logodds = FALSE,
	  estimate=b_sd,
	  se=se.sd,
	  shape=sig,
	  colour = Category,
	   xlab = "SD change (95% CI) in risk factor per copy of the major allele")+
# labs(title=Title.plot,size=2)+
theme(plot.title = element_text(size = 10))

out.file2<-paste("~/fatty-acids/mr/results/rf_pleiotropy_",unique(Plot2$SNP),sep="")
Title.plot<-unique(Plot2$SNP)
Plot.pdf2<-forestplot(df = Plot2,
	name=trait,
	  logodds = FALSE,
	  estimate=b_sd,
	  se=se.sd,
	  shape=sig,
	  colour = Category,
	   xlab = "SD change (95% CI) in risk factor per AA raising allele")+
# labs(title=Title.plot,size=2)+
theme(plot.title = element_text(size = 10))

out.file3<-paste("~/fatty-acids/mr/results/rf_pleiotropy_",unique(Plot3$SNP),sep="")
Title.plot<-unique(Plot3$SNP)
Plot.pdf3<-forestplot(df = Plot3,
	name=trait,
	  logodds = FALSE,
	  estimate=b_sd,
	  se=se.sd,
	  shape=sig,
	  colour = Category,
	   xlab = "SD change (95% CI) in risk factor per DGLA raising allele")+
# labs(title=Title.plot,size=2)+
theme(plot.title = element_text(size = 10))
	

out.file4<-paste("~/fatty-acids/mr/results/rf_pleiotropy_",paste(unique(Plot4$SNP),collapse="_"),sep="")
Title.plot<-unique(Plot4$SNP)
Plot.pdf4<-forestplot(df = Plot4,
	name=trait,
	  logodds = FALSE,
	  estimate=b_sd,
	  se=se.sd,
	  shape=SNP,
	  colour = Category,
	   xlab = "SD change (95% CI) in risk factor per PUFA raising allele")+
# labs(title=Title.plot,size=2)+
theme(plot.title = element_text(size = 10))

pdf(out.file1)
	Plot.pdf1
dev.off()
pdf(out.file2)
	Plot.pdf2
dev.off()
pdf(out.file3)
	Plot.pdf3
dev.off()
pdf(out.file4)
	Plot.pdf4
dev.off()

save.image("~/fatty-acids/mr/results/instrument_pleiotropy.RData")
# load("~/fatty-acids/mr/results/instrument_pleiotropy.RData")

# a<-Plot1[Plot1$trait %in% c("Heart rate","Platelet count"),]
# a[a$SNP=="rs174546",]
# z<-a$beta.outcome/a$se.outcome
# n<-a$sample_size
# maf<-a$eaf.outcome
#  sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)

format_dat<-function(Dat=NULL){
	Dat1<-Dat[order(Dat$samplesize.outcome,decreasing=T),]
	Dat1$originalname.outcome[Dat1$originalname.outcome=="Haemoglobin concentration"] <- "Hemoglobin concentration"
	Dat1$originalname.outcome[Dat1$originalname.outcome=="Mean corpuscular hemoglobin"] <- "Mean cell haemoglobin"
	Dat1$originalname.outcome[Dat1$originalname.outcome=="Mean corpuscular hemoglobin concentration"] <- "Mean cell haemoglobin concentration"
	Dat1$originalname.outcome[Dat1$originalname.outcome=="Mean corpuscular volume"] <- "Mean cell volume"

	Dat2<-Dat1[!duplicated(paste(Dat1$originalname.outcome,Dat1$SNP)),]

	rf_haem<-Dat2[Dat2$subcategory=="Haemotological",]
	rf_Xhaem<-Dat2[Dat2$subcategory!="Haemotological",]
	rf_haem1<-rf_haem[rf_haem$originalname.outcome %in% c("White blood cell count (basophil)","Eosinophil counts","Granulocyte count","Hematocrit","Hemoglobin concentration","Lymphocyte counts","Mean cell haemoglobin","Mean cell haemoglobin concentration","Mean cell volume","Mean platelet volume","Monocyte count","Myeloid white cell count","Neutrophil count","Packed cell volume","Platelet count","Plateletcrit","Red blood cell count","Reticulocyte count","White blood cell count"),]
	Dat3<-rbind(rf_haem1,rf_Xhaem)

	Dat4<-Dat3[which(!Dat3$consortium %in% c("CHARGE","ADIPOGen")),]
	Dat5<-Dat4[which(!Dat4$subcategory %in% c("Personality","Protein","Psychiatric / neurological","Metal")),]
	Dat6<-Dat5[!Dat5$originalname.outcome %in% c("Childhood intelligence", "College completion",#education
		"Obesity class 1", "Obesity class 2","Obesity class 3" ,"Birth length",#anthropometric
		"Extreme height", "Extreme body mass index","Childhood obesity","Difference in height between childhood and adulthood"  ,"Infant head circumference","Extreme waist-to-hip ratio","Difference in height between adolescence and adulthood",
		"Parents' age at death", #age at death
		"Father's age at death","Top 1 % survival",
		"Urate",
		 "Overweight", #binary risk factors
		 "Alcohol dependence",
		 "Oligoclonal band status",
		 "Microalbuminuria",
		 "Leptin", #hormone
		 "Plasma cortisol"
		 ),]

	Dat7<-Dat6[!is.na(Dat6$beta.outcome),]
	Dat8<-Dat7[which(Dat7$samplesize.outcome>5000),]
	Dat8<-Dat8[order(Dat8$subcategory),]

	Dat8$subcategory[Dat8$originalname.outcome %in% c("telomere length", "Years of schooling","Heart rate","C-reactive protein", "Percent emphysema")]<-"Other"
	Dat8$subcategory[Dat8$subcategory %in% c("Lipid","Glycemic")]<-"Lipid/glycemic"
	Dat8$subcategory[Dat8$subcategory %in% c("Kidney")]<-"Kidney function"
	Dat8$subcategory[Dat8$subcategory %in% c("Bone")]<-"Bone health"
	Dat8$subcategory[Dat8$subcategory %in% c("Reproductive aging", "Sleeping")]<-"Other"
	Sort.dat<-data.frame(matrix(c(as.factor(c("Anthropometric","Behavioural","Lipid/glycemic","Haemotological","Kidney function","Bone health","Other") ),c("Anthropometric","Behavioural","Lipid/glycemic","Haemotological","Kidney function","Bone health","Other")),nrow=7,ncol=2))
	names(Sort.dat)<-c("Sort","subcategory")

	Sort.dat$Sort<-1:nrow(Sort.dat)
	Sort.dat$subcategory<-as.character(Sort.dat$subcategory)
	Dat9<-merge(Dat8,Sort.dat,by="subcategory")
	Dat9<-Dat9[order(Dat9$Sort),]
	Dat9$population
	names(Dat9)[names(Dat9) == "subcategory"]<-"Category"
	if(unique(Dat9$SNP) == "rs174546" ){
		exposure<-aa_dgla_instrument()
		exposure$Category<-"PUFA desaturase activity"
	}
	if(unique(Dat9$SNP) == "rs1741" ){
		exposure<-aa_PDXDC1_instrument()
		exposure$Category<-"Arachidonic acid"
	}

	if(unique(Dat9$SNP) == "rs881803" ){
		exposure<-dgla_PDXDC1_instrument()
		exposure$Category<-"Dihomo-gamma-linolenic acid"
		exposure<-exposure[exposure$SNP =="rs881803", ]
		exposure$effect_allele
	}


	Dat9<-plyr::rbind.fill(exposure,Dat9)
	Z<-Dat9$beta.outcome/Dat9$se.outcome
	Pos<-which(Dat9$beta.outcome==0)
	Z2<-qnorm(Dat9$pval.outcome/2,lower.tail=FALSE)
	Dat9$b_sd<-b_sd(z=Z,maf=Dat9$eaf.outcome,n=Dat9$samplesize.outcome)
	Dat9$b_sd2<-b_sd(z=Z2,maf=Dat9$eaf.outcome,n=Dat9$samplesize.outcome)
	Dat9$se.sd<-abs(Dat9$b_sd /Z)
	Dat9$se.sd2<-abs(Dat9$b_sd2 /Z2)
	Pos<-grep("SD",Dat9$unit)
	Pos2<-which(is.na(Dat9$b_sd))
	Pos<-Pos2[Pos2 %in% Pos]
	Dat9$b_sd[Pos]<-Dat9$beta.outcome[Pos]
	Dat9$se.sd[Pos]<-Dat9$se.outcome[Pos]
	Dat9<-Dat9[!is.na(Dat9$b_sd),]
	Pos<-Dat9$b_sd==0
	Dat9$b_sd[Pos]<-Dat9$b_sd2[Pos]
	Dat9$se.sd[Pos]<-Dat9$se.sd2[Pos]
	Dat9$id<-paste(Dat9$originalname.outcome,Dat9$SNP)
	Dat9$effect_allele.outcome<-toupper(Dat9$effect_allele.outcome)
	Dat9$other_allele.outcome<-toupper(Dat9$other_allele.outcome)
	if(length(unique(Dat9$effect_allele.outcome))>1) stop("effect allele inconsistent across traits")

	if(exposure$beta.outcome<0){
		Dat9$b_sd<-Dat9$b_sd *-1	
		effect_allele<-unique(Dat9$effect_allele.outcome)
		other_allele<-unique(Dat9$other_allele.outcome)
		Dat9$effect_allele.outcome<-other_allele
		Dat9$other_allele.outcome<-effect_allele
		Dat9$eaf.outcome<-1-Dat9$eaf.outcome
	}

	Dat9$trait<-gsub("(^[[:alpha:]])", "\\U\\1", Dat9$trait, perl=TRUE)
	Dat9$sig<- FALSE
	Dat9$sig[Dat9$pval.outcome<0.05/36]<-TRUE
	Dat9$lci<-Dat9$b_sd-1.96*Dat9$se.sd
	Dat9$uci<-Dat9$b_sd+1.96*Dat9$se.sd
	return(Dat9)
}

aa_dgla_instrument<-function(){
	source("~/fatty-acids/mr/scripts/mr_functions.R")
	exp<-read.table("~/fatty-acids/mr/data/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",head=T,stringsAsFactors=F)
	exp[exp$population=="East Asian",]
	exposure_dat<-format_exposure2(dat=exp,standardise_beta=TRUE)

	exposure_dat<-exposure_dat[exposure_dat$SNP == "rs174546",]
	eur<-exposure_dat[exposure_dat$population == "European",]
	eur<-eur[eur$exposure == "AA:DGLA",]
	names(eur)[names(eur) == "exposure"]<-"trait"
	names(eur)<-gsub("exposure","outcome",names(eur))
	# names(eur)[names(eur) == "beta.exposure"]<-"beta.outcome"
	# names(eur)[names(eur) == "se.exposure"]<-"se.outcome"
	# names(eur)[names(eur) == "eaf.exposure"]<-"eaf.outcome"
	# names(eur)[names(eur) == "effect_allele.exposure"]<-"eaf.outcome"
	# names(eur)[names(eur) == "samplesize.exposure"]<-"samplesize.outcome"
	return(eur)
}


	

}

b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}

aa_PDXDC1_instrument<-function(){
	source("~/fatty-acids/mr/scripts/mr_functions.R")
	# exposure
	source("~/fatty-acids/mr/scripts/mr_functions.R")
	load("~/fatty-acids/mr/data/instruments.Rdata")

	exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")

	exp_dat<-exp_dat[!exp_dat$FADS,]
	exp_dat<-exp_dat[exp_dat$exposure == "Arachidonic acid (20:4n6)"   ,]
	names(exp_dat)<-gsub("exposure","outcome",names(exp_dat))
	names(exp_dat)[names(exp_dat)=="outcome"]<-"trait"
	return(exp_dat)
}


dgla_PDXDC1_instrument<-function(){
	source("~/fatty-acids/mr/scripts/mr_functions.R")
	# exposure
	source("~/fatty-acids/mr/scripts/mr_functions.R")
	load("~/fatty-acids/mr/data/instruments.Rdata")

	exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")

	exp_dat<-exp_dat[!exp_dat$FADS,]
	exp_dat<-exp_dat[exp_dat$exposure == "Dihomo-gamma-linolenic acid (20:3n6)"     ,]
	names(exp_dat)<-gsub("exposure","outcome",names(exp_dat))
	names(exp_dat)[names(exp_dat)=="outcome"]<-"trait"
	return(exp_dat)
}

