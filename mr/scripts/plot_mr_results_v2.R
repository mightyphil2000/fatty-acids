source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")

library(devtools)
library(TwoSampleMR)
library(plyr)
library(ggforestplot)
library(ggplot2)

library(tidyverse)

# overall cancer
# sort(unique(mr_res1$cancer))
mr_res<-format_bystudy(cancer<-"Cancer (all cause)",overall_effect=TRUE)
mr_res[,c("study.abbreviation","Q.p")]

mr_res$weight<-1/mr_res$se/10
# exclude the smaller fingen study
mr_res<-mr_res[mr_res$cases!="9792",]
P1<-forestplot(df = mr_res,logodds = TRUE,name=study.abbreviation,
				  estimate=b, se=se,shape=NULL,
					colour = NULL,xlab = "")+
					geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=20))
png("~/fatty-acids/mr/results/plots/ggforest_overall_cancer_bystudy.png", width = 600, height = 480)
	print(P1) 
dev.off()

# exclude non melanoma skin cancer
# sort(unique(mr_res1$cancer))
mr_res<-format_bystudy(cancer<-"Cancer (excluding non-melanoma skin cancer)",overall_effect=TRUE)
mr_res[,c("study.abbreviation","Q.p")]
mr_res$weight<-1/mr_res$se/10
# exclude the smaller fingen study
mr_res<-mr_res[mr_res$cases!="9792",]
P1<-forestplot(df = mr_res,logodds = TRUE,name=study.abbreviation,
				  estimate=b, se=se,shape=NULL,
					colour = NULL,xlab = "")+
					geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=20))
png("~/fatty-acids/mr/results/plots/ggforest_overall_cancer_bystudy.png", width = 600, height = 480)
	print(P1) 
dev.off()


# respiratory and intrathoracic cancer

mr_res<-format_bystudy(cancer<-"Respiratory and intrathoracic cancer",overall_effect=TRUE)
mr_res[,c("study.abbreviation","Q.p")]
mr_res$weight<-1/mr_res$se
P1<-forestplot(df = mr_res,logodds = TRUE,name=study.abbreviation,
				  estimate=b, se=se,shape=NULL,
					colour = NULL,xlab = "")+
					geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=20))

png("~/fatty-acids/mr/results/plots/ggforest_respiratory_intrathoracic_cancer_bystudy.png", width = 600, height = 480)
	print(P1) 
dev.off()


# Colorectal cancer

mr_res<-format_bystudy(cancer<-"Colorectal cancer",overall_effect=TRUE)
mr_res[,c("exposure","study.abbreviation","Q.p")]
mr_res$weight<-1/mr_res$se/10
P1<-forestplot(df = mr_res,logodds = TRUE,name=study.abbreviation,
				  estimate=b, se=se,shape=NULL,
					colour = NULL,xlab = "")+
					geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=15))

png("~/fatty-acids/mr/results/plots/ggforest_crc_bystudy_mr_results_rep_v3.png", width = 600, height = 480)
	print(P1) 
dev.off()


# lung cancer

mr_res<-format_bystudy(cancer<-"Lung cancer",overall_effect=TRUE)
mr_res[,c("exposure","study.abbreviation","Q.p")]
mr_res$weight<-1/mr_res$se/4
P1<-forestplot(df = mr_res,logodds = TRUE,name=study.abbreviation,
				  estimate=b, se=se,shape=NULL,
					colour = NULL,xlab = "")+
					geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=15))

png("~/fatty-acids/mr/results/plots/ggforest_lung_bystudy_mr_results_rep_v3.png", width = 600, height = 480)
	print(P1) 
dev.off()

# Esophageal squamous cell carcinoma
mr_res<-format_bystudy(cancer<-"Esophageal squamous cell carcinoma"     ,overall_effect=TRUE)
mr_res[,c("exposure","study.abbreviation","Q.p")]
mr_res$weight<-1/mr_res$se
P1<-forestplot(df = mr_res,logodds = TRUE,name=study.abbreviation,
				  estimate=b, se=se,shape=NULL,
					colour = NULL,xlab = "")+
					geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=15))

png("~/fatty-acids/mr/results/plots/ggforest_escc_bystudy_mr_results_rep_v3.png", width = 600, height = 480)
	print(P1) 
dev.off()
	
# Malignant skin cancer

mr_res<-format_bystudy(cancer<-"Malignant skin cancer"     ,overall_effect=TRUE)
mr_res[,c("exposure","study.abbreviation","Q.p")]
mr_res$weight<-1/mr_res$se/5
P1<-forestplot(df = mr_res,logodds = TRUE,name=study.abbreviation,
				  estimate=b, se=se,shape=NULL,
					colour = NULL,xlab = "")+
					geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=15))

png("~/fatty-acids/mr/results/plots/ggforest_msc_bystudy_mr_results_rep_v3.png", width = 600, height = 480)
	print(P1) 
dev.off()
	


mr_res<-format_dat2(dat=mr_res1,p_threshold=0.05/67,sort="system")

# mr_res<-format_dat2(dat=mr_res1,p_threshold=0.05/67,sort="system")
mr_res$OR<-round(exp(mr_res$b),2)
mr_res$LCI<-round(exp(mr_res$b-mr_res$se*1.96),2)
mr_res$UCI<-round(exp(mr_res$b+mr_res$se*1.96),2)
# mr_res$pval<-round(mr_res$pval,2)

mr_res[,c("outcome","OR","LCI","UCI","pval","id.outcome")]

mr_res<-format_dat2(dat=mr_res1,sort="system")
dim(mr_res)
# mr_res<-mr_res[which(mr_res$Cancer.Group == "Colorectal cancer"),]
# mr_res<-mr_res[which(mr_res$system == "Blood"),]
mr_res$OR<-round(exp(mr_res$b),2)
mr_res$LCI<-round(exp(mr_res$b-mr_res$se*1.96),2)
mr_res$UCI<-round(exp(mr_res$b+mr_res$se*1.96),2)
mr_res[,c("outcome","OR","LCI","UCI","pval")]


# Plot results for cancers that did not pass 0.05/67 and which had powr to detect odds ratios ≥1.05
mr_res<-format_nulls(power05=FALSE)
mr_res[,c("outcome","OR","LCI","UCI","pval","cases","controls")]

# mr_res$cancer

# dev.off()
P1<-forestplot(df = mr_res,logodds = TRUE,name=outcome,
				  estimate=b, se=se,
				  shape=NULL, colour = NULL,xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=14))


P1<-forest_plot_1_to_many(mr_res = mr_res,b = "b",se = "se",TraitM = "cancer",col1_width = 2.0,col1_title = "",
		exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = 0.90,up = 1.2, by = "Cancer.Group",xlab = "",addcols = "cases",
		addcol_widths = 0.5,addcol_titles = "", subheading_size = NULL,shape_points = 15,colour_scheme = "black",col_text_size = 3,weight = "weight" )	
# P1
####################################################################################
# Plot results for all cancers defined as discovery studies and passing Pbonf<0.05#
######################################################################################

# selectedd either AA:DGLA or GLA:LA for each cancer. AA:DGLA only vaid for European only analyses. GLA:LA only valid for East Asian or European+East Asian analyses

P1<-forestplot(df = mr_res,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = NULL,xlab = "")+
		geom_point(shape="square",size=1/mr_res$se/5,fill=c(mr_res$Colour),colour = mr_res$Colour)+
		theme(text = element_text(size=20))

mr_res$Colour<-NA
mr_res$Colour[mr_res$system == "Digestive"]<-"black"
mr_res$Colour[mr_res$system == "Respiratory"]<-"red"
mr_res$Colour[mr_res$system == "Integumentary"]<-"blue"

P1<-forestplot(df = mr_res,
			logodds = TRUE,
			name=cancer,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = NULL,
				   xlab = "")+
			theme(text = element_text(size=20))+
			geom_point(shape="square",size=1/mr_res$se/5,fill=c(mr_res$Colour),colour = mr_res$Colour)


			# theme(text = element_text(size=100))
# theme(plot.title = element_text(size = ""))+


png("~/fatty-acids/mr/results/plots/ggforest_cancerv2.png", width = 900, height = 1000)
	print(P1) 
dev.off()

###########################
# gastrointestinal cancers#
###########################

mr_res1$Cancer.Group[mr_res1$Cancer.Group %in% c("Liver & biliary tract cancer","Cancer of digestive organs","Pancreatic cancer","Digestive system cancer")]<-"Other"
mr_res<-format_dat2(dat=mr_res1,sort="Cancer.Group")

mr_res<-mr_res[mr_res$system == "Digestive",]


P1<-forestplot(df = mr_res,
			logodds = TRUE,
			name=cancer,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group,
				   xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=14))

png("~/fatty-acids/mr/results/plots/ggforest_digestive_cancer.png", width = 600, height = 480)
	print(P1) 
dev.off()

###########################
#Respiratory cancerss
###########################
mr_res1<-mr_res1[!is.na(mr_res1$system),]
mr_res1<-mr_res1[mr_res1$system =="Respiratory",]
mr_res1[,c("outcome","b","pval")]
mr_res<-format_dat2(dat=mr_res1,sort="Cancer.Group")


P1<-forestplot(df = mr_res,
			logodds = TRUE,
			name=cancer,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group,
				   xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=15))

png("~/fatty-acids/mr/results/plots/ggforest_respiratory_cancers.png", width = 900, height = 500)
	print(P1) 
dev.off()

###########################
# smoking associated cancers
###########################
mr_res1$Cancer.Group[mr_res1$Cancer.Group %in% c("Respiratory and intrathoracic cancer","Esophageal cancer")]<-"Other"
mr_res<-format_dat2(dat=mr_res1,sort="Cancer.Group")
Plot1<-mr_res[mr_res$system == "Respiratory",]
Plot2<-mr_res[mr_res$cancer == "Esophageal squamous cell carcinoma \ncases= 3313",]
mr_res<-rbind(Plot1,Plot2)


P1<-forestplot(df = mr_res,
			logodds = TRUE,
			name=cancer,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group,
				   xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=12))

png("~/fatty-acids/mr/results/plots/ggforest_smoking_cancer.png", width = 600, height = 600)
	print(P1) 
dev.off()


###########################
# skin cancers
###########################

mr_res<-format_dat2(dat=mr_res1)
mr_res<-mr_res[mr_res$system == "Integumentary",]

P1<-forestplot(df = mr_res,
			logodds = TRUE,
			name=cancer,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group,
				   xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=14))+
			theme(legend.position = "none")


png("~/fatty-acids/mr/results/plots/ggforest_skin_cancer.png", width = 400, height = 400)
	print(P1) 
dev.off()

