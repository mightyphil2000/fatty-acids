source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
# rm(list=ls())
# install.packages("devtools")

# devtools::install_github("MRCIEU/TwoSampleMR")

library(devtools)
library(TwoSampleMR)
library(plyr)
library(ggforestplot)
library(ggplot2)

load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
load("~/fatty-acids/mr/results/mr_results_discovery_v2.Rdata")

# load("~/fatty-acids/mr/results/mr_results_rep.Rdata")
# load("~/fatty-acids/mr/results/mr_results_discovery.Rdata")
table(mr_res_disc$Cancer.Group)
mr_res_disc[mr_res_disc$Cancer.Group == "Colorectal cancer",c("cancer","cases","controls","OR","LCI","UCI","study.abbreviation","pval")]
mr_res_rep[which(mr_res_rep$cancer == "Colorectal cancer"),c("cases","controls","OR","LCI","UCI","study.abbreviation","pval")]
mr_res_rep[which(mr_res_rep$Cancer.Group == "Colorectal cancer"),c("cancer","site","cases","controls","OR","LCI","UCI","study.abbreviation","pval")]

########################
# Colorectal cancer####
########################

Dat<-format_colorectal()

P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 2.9,col1_title = "\n",
		exponentiate = TRUE,trans = "identity",ao_slc = FALSE,lo =0.95,up = 1.15, by=NULL,xlab="",xlab_top = "OR (95% CI) per SD genetic\nincrease in FADS1 enzyme activity",addcols = "cases",addcol_widths = 1.0,addcol_titles = "No. of\ncases", subheading_size = 20,shape_points = 15,colour_scheme = "black",col_text_size = 6,weight ="weight",row_line_colour="white",x_axis_label_size=20)	

png("~/fatty-acids/mr/results/plots/colorectal_cancer_1tm.png", width = 960, height = 960)
	print(P1) 
dev.off()

Dat<-format_colorectal2()
# Dat[,c("cancer","OR","LCI","UCI")]

P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 2.9,col1_title = "\n",
		exponentiate = TRUE,trans = "identity",ao_slc = FALSE,lo =0.85,up = 1.20, by="crc_group",xlab="",xlab_top = "OR (95% CI) per SD genetic\nincrease in FADS1 enzyme activity",addcols = "cases",addcol_widths = 1.0,addcol_titles = "No. of\ncases", subheading_size = 20,shape_points = 15,colour_scheme = "black",col_text_size = 6,weight ="weight",row_line_colour="white",x_axis_label_size=20)	

png("~/fatty-acids/mr/results/plots/colorectal_cancer_1tm_v2.png", width = 960, height = 960)
	print(P1) 
dev.off()


##################################################
# Plot results for all cancers passing Pbonf<0.05#
###################################################

Dat<-format_dat(Dat=mr_res_disc)
names(mr_res_rep)[grep("population",names(mr_res_rep))]
exposure_dat[,c("exposure","id.exposure","population")]
format_dat2<-function(dat=NULL){
	unique(Dat[,c("exposure","id.exposure","population","beta.exposure")])
	head(Dat)
	Dat1<-Dat[Dat$population == "European; East Asian" & expusure == "GLA:LA", ]
	Dat2<-Dat[Dat$population == "East Asian" & expusure == "GLA:LA", ]
	Dat3<-Dat[Dat$population == "European" & expusure == "GLA:LA", ]
}
head(Dat)

head(exposure_dat)

P1<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = system,
				   xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=14))
png("~/fatty-acids/mr/results/plots/ggforest_cancer.png", width = 960, height = 960)
	print(P1) 
dev.off()

Dat2<-format_dat(Dat=mr_res_disc,p_threshold=0.05)
P2<-forestplot(df = Dat2,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = system,
				   xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=10))
pdf("~/fatty-acids/mr/results/plots/ggforest_cancer05.pdf")
	print(P2)
dev.off()

Dat3<-format_dat(Dat=mr_res_disc,p_threshold=0.05/67)
P3<-forestplot(df = Dat3,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = system,
				   xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=15))
png("~/fatty-acids/mr/results/plots/ggforest_cancerBon05.png", width = 500, height = 500)
	print(P3)
dev.off()

Dat3$cancer<-gsub("Esophageal squamous cell carcinoma","Esophageal squamous\ncell carcinoma",Dat3$cancer)
Dat3$cancer<-gsub("Malignant non-melanoma skin cancer","Malignant non-melanoma\nskin cancer",Dat3$cancer)
Dat3$cancer<-gsub("Respiratory and intrathoracic cancer","Respiratory & intrathoracic\ncancer",Dat3$cancer)
Dat3$system[Dat3$system == "Integumentary"]<-"Skin"
P4<-forest_plot_1_to_many(mr_res = Dat3,b = "b",se = "se",TraitM = "cancer",col1_width = 1.2,col1_title = c("",""),
		exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = 0.95,up = 1.25, by = "system",xlab = "",addcols = "cases",
		addcol_widths = 0.4,addcol_titles = "", subheading_size = 12,shape_points = 15,colour_scheme = "black",col_text_size = 4,weight = NULL )	
png("~/fatty-acids/mr/results/plots/1tmforest_cancerBon05.png",width = 700, height = 700)
	print(P4)
dev.off()

# mr_res[mr_res$pval<0.05/43,c("cancer","OR","LCI","UCI","power10")]
# Dat[Dat$system=="Blood",c("cancer","id.outcome")]
# Dat$cancer2<-Dat$cancer
# plot_dat(Dat=Dat)


# List.plots<-system_plots()
# List.plots[5]
# forest1m<-List.plots[1]
# ggforest<-List.plots[2]
# Systems<-List.plots[3]

# mr_res[mr_res$cancer =="Endometrioid ovarian cancer",c("cancer","Cancer.Group","system","site")]

################################################
# Plots for cancers by system or cancer group
##############################

####################
# skin cancers
####################

Dat<-format_skin()
P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "study.abbreviation",col1_width = 1,col1_title = "",
		exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up = NULL, by = "cancer",xlab = "",addcols = "cases",
		addcol_widths = 1,addcol_titles = "", subheading_size = 7,shape_points = 15,colour_scheme = "black",col_text_size = 2,weight = NULL )	
Dat<-format_skin3()
P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 1.7,col1_title = "",
		exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up = NULL, by = NULL,xlab = "",addcols = "cases",
		addcol_widths = 0.7,addcol_titles = "", subheading_size = NULL,shape_points = 15,colour_scheme = "black",col_text_size = 3,weight = NULL )	

File1<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$system),"1tm_v2.pdf")
pdf(File1)
	print(P1)
dev.off()

Dat<-format_skin3()
File2<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$system),"ggf.pdf")

P2<-forestplot(df = Dat,logodds = TRUE,name=cancer2,estimate=b,se=se, shape=NULL,colour = NULL,xlab = "")+theme(plot.title = element_text(size = 5))+theme(text = element_text(size=10))
pdf(File2)
	print(P2)
dev.off()


# ovarian
Dat<-format_site(Site="Ovary")
Dat<-format_ovarian()

dev.off()
P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 1,col1_title = "",		exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up = NULL, by = "ovarian_subtype", xlab = "", 
		addcols = "cases",addcol_widths = 0.3, addcol_titles = "", subheading_size = 7,  shape_points = 15,
		colour_scheme = "black", col_text_size = 2)

File1<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$site),"1tm.pdf")
pdf(File1)
	print(P1)
dev.off()

P2<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=Shape,
				  colour = ovarian_subtype,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))

File2<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$site),"ggf.pdf")
pdf(File2)
	print(P2)
dev.off()

# prostate
Dat<-format_prostate()
P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 1,col1_title = "",		exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up = NULL, by = NULL, xlab = "", 
		addcols = "cases",addcol_widths = 0.4, addcol_titles = "", subheading_size = 7,  shape_points = 15,
		colour_scheme = "black", col_text_size = 3)
File1<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$site),"1tm.pdf")
pdf(File1)
	print(P1)
dev.off()

P2<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = NULL,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))
File2<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$site),"ggf.pdf")
pdf(File2)
	print(P2)
dev.off()


###########################
# Digestive system cancers
##############################

# mr_res[mr_res$cancer == "Esophageal cancer",c("cancer","site","system","Cancer.Group")]
# mr_res[mr_res$system == "Digestive",c("cancer","site","system","Cancer.Group")]
# dat1[order(dat1$Cancer.Group),c("cancer","Cancer.Group","study.abbreviation")]
Dat<-format_digestive()
Dat[Dat$Cancer.Group=="Esophageal cancer",c("cancer","id.outcome","study.abbreviation","Cancer.Group")]
# mr_res_disc[mr_res_disc$Cancer.Group == "Esophageal cancer",c("cancer","id.outcome","study.abbreviation","Cancer.Group")]
# mr_res_rep[mr_res_rep$Cancer.Group == "Esophageal cancer",c("cancer","id.outcome","study.abbreviation","Cancer.Group")]
# mr_res_disc[mr_res_disc$cancer == "Esophageal squamous cell carcinoma",c("study.abbreviation","cases","b","se","pval")]
# mr_res_rep[mr_res_rep$Cancer.Group == "Esophageal cancer",c("study.abbreviation","cases","b","se","pval")]


P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 2.5 ,col1_title = "",exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up =NULL, by = "Cancer.Group", xlab = "", addcols = "cases",addcol_widths = 0.5, addcol_titles = "", subheading_size = 10,  shape_points = 15,colour_scheme = "black", col_text_size = 3)

File1<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$system),"1tm.pdf")
pdf(File1)
	print(P1)
dev.off()

P2<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))
File2<-paste0("~/fatty-acids/mr/results/plots/",unique(Dat$system),"ggf.pdf")
pdf(File2)
	print(P2)
dev.off()

############
# All cause#
#############

Dat<-format_allcause()
P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "allcause1",col1_width = 1.2 ,col1_title = "",exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up =NULL, by = NULL, xlab = "", addcols = "cases",addcol_widths = 0.5, addcol_titles = "", subheading_size = 10,  shape_points = 15,colour_scheme = "black", col_text_size = 3)
File1<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$Cancer.Group)),"1tm.pdf")
pdf(File1)
	print(P1)
dev.off()

P2<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))

File2<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$Cancer.Group)),"ggf.pdf")
pdf(File2)
	print(P2)
dev.off()

########################
#Respiratory cancers###
########################

Dat<-format_respiratory()
Dat<-format_respiratory2()

P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 1.5 ,col1_title = "",exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up =NULL, by = "Cancer.Group", xlab = "", addcols = "cases",addcol_widths = 0.4, addcol_titles = "", subheading_size = 10,  shape_points = 15,colour_scheme = "black", col_text_size = 3)
File1<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$system)),"1tm.pdf")
pdf(File1)
	print(P1)
dev.off()

P2<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))

File2<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$system)),"ggf.pdf")
pdf(File2)
	print(P2)
dev.off()

# Dat<-format_blood()
# P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 1.5 ,col1_title = "",exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up =NULL, by = "Cancer.Group", xlab = "", addcols = "cases",addcol_widths = 0.4, addcol_titles = "", subheading_size = 10,  shape_points = 15,colour_scheme = "black", col_text_size = 3)
# File1<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$system)),"1tm.pdf")
# pdf(File1)
# 	print(P1)
# dev.off()


##########
#Blood cancer#
###########

Dat<-format_blood()
Dat2<-format_blood2()

P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 2.2 ,col1_title = "",exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up =NULL, by = "Cancer.Group2", xlab = "", addcols = "cases",addcol_widths = 0.4, addcol_titles = "", subheading_size = 10,  shape_points = 15,colour_scheme = "black", col_text_size = 3)
P2<-forest_plot_1_to_many(mr_res = Dat2,b = "b",se = "se",TraitM = "cancer",col1_width = 2.2 ,col1_title = "",exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up =NULL, by = "Cancer.Group2", xlab = "", addcols = "cases",addcol_widths = 0.4, addcol_titles = "", subheading_size = 10,  shape_points = 15,colour_scheme = "black", col_text_size = 3)
File2<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$system)),"1tm.pdf")
pdf(File2)
	print(P2)
dev.off()

File1<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$system)),"1tm_incl_cases<1000.pdf")
pdf(File1)
	print(P1)
dev.off()


P2<-forestplot(df = Dat2,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group2,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))

File2<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$system)),"ggf.pdf")
pdf(File2)
	print(P2)
dev.off()

P1<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = Cancer.Group2,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))

File1<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$system)),"ggf_incl_cases<1000.pdf")
pdf(File1)
	print(P1)
dev.off()


#breast
Dat<-format_breast()

P1<-forest_plot_1_to_many(mr_res = Dat,b = "b",se = "se",TraitM = "cancer",col1_width = 2.2 ,col1_title = "",exponentiate = TRUE,trans = "log2",ao_slc = FALSE,lo = NULL,up =NULL, by = NULL, xlab = "", addcols = "cases",addcol_widths = 0.8, addcol_titles = "", subheading_size = 10,  shape_points = 15,colour_scheme = "black", col_text_size = 4)
File1<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$Cancer.Group)),"1tm.pdf")
pdf(File1)
	print(P1)
dev.off()

P2<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = NULL,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = 5))+
			theme(text = element_text(size=10))

File2<-paste0("~/fatty-acids/mr/results/plots/",gsub(" ","_",unique(Dat$Cancer.Group)),"ggf.pdf")
pdf(File2)
	print(P2)
dev.off()



################################################################################################################
#Plots for cancers separately by each available study and meta analysis#
####################################################################################

p<-NULL
Files<-NULL
for(i in 1:length(cancers)){
	File<-paste0("~/fatty-acids/mr/results/plots/rep/",cancers[i],".pdf")
	Files[[i]]<-gsub(" ","_",File)
	p[[i]]<-replication_plots(cancer=cancers[i])
	dev.off()
}

i<-12
pdf(Files[i])
p[i]
dev.off()	


head(meta.tab9)


replication_plots(cancer="Oral cavity and pharyngeal cancer")
