source("~/fatty-acids/mr/scripts/mr_functions.R")
library(ggforestplot)
library(ggplot2)

# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
# load("~/fatty-acids/mr/results/results_ara_snps_meta_analysis_wald_ratios.Rdata")
# source("~/fatty-acids/mr/scripts/wald_ratios_rs3734398_meta_analysis_europeans.R")

load(file="~/fatty-acids/mr/results/results_wald_ratios_ara_meta_analysis_europeans.Rdata") #meta analysis result only plus outcomes where only one study available (e.g. UKB overall cancer excluded non-menlanoma skin cancer and UKB non melanoma skin cancer)
load(file="~/fatty-acids/mr/results/results_wald_ratios_ara_meta_analysis_europeans2.Rdata") #meta analysis results plus individual study results
meta_dat3<-rbind(meta_dat,meta_dat2)

Plot_dat<-format_plot(Plot_dat=meta_dat3)

P1<-forestplot(df = Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,colour = SNP,xlab = "")+
		scale_x_continuous(breaks=c(0.70,1.0,1.5,2.0),labels=c(0.70,1.0,1.5,2.0),trans="log")
png("~/fatty-acids/mr/results/plots/ggforest_rs3734398_rs4985155_rs174546_wald_ratios_ara.png", width = 600, height = 480)
	print(P1) 
dev.off()


# meta analysis across overall, lung and crc
Plot_dat<-format_meta_analysis_rs3734398_wald_ratio2(meta_dat=meta_dat3)

unique(Plot_dat2$outcome2)
Plot_dat[,c("SNP","study","p","or","lci","uci")]

P1<-forestplot(df = Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,colour = SNP,xlab = "")

png("~/fatty-acids/mr/results/plots/ggforest_rs3734398_rs4985155_rs174546_wald_ratios_maxpower_ara.png", width = 600, height = 480)
	print(P1) 
dev.off()

# combine two plots
Plot_dat2<-format_plot(Plot_dat=meta_dat3)
Plot_dat2<-Plot_dat2[Plot_dat2$outcome2 %in% c("Non-melanoma skin cancer","Basal cell carcinoma","Malignant skin cancer"),]
Plot_dat<-format_meta_analysis_rs3734398_wald_ratio2(meta_dat=meta_dat3)
Plot_dat<-plyr::rbind.fill(Plot_dat,Plot_dat2)
unique(Plot_dat2$outcome2)

Plot_dat[,c("SNP","outcome2","study","p","or","lci","uci")]

Plot_dat$outcome_name<-gsub("Overall effect","Overall cancer2",Plot_dat$outcome_name)

P1<-forestplot(df = Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,colour = SNP,xlab = "")+
	scale_x_continuous(breaks=c(0.70,1.0,1.5,2.0),labels=c(0.70,1.0,1.5,2.0),trans="log")
png("~/fatty-acids/mr/results/plots/ggforest_rs3734398_rs4985155_rs174546_wald_ratios_maxpower_ara2.png", width = 600, height = 480)
	print(P1) 
dev.off()

format_plot<-function(Plot_dat=NULL){
	ID<-paste(Plot_dat$outcome2,Plot_dat$SNP)
	Plot_dat<-Plot_dat[!duplicated(ID),]
	Plot_dat$b<-as.numeric(Plot_dat$b)
	Plot_dat$se<-as.numeric(Plot_dat$se)
	Plot_dat$ncase.outcome<-as.numeric(Plot_dat$ncase.outcome)
	Plot_dat$shape<-15
	# Plot_dat$shape[Plot_dat$SNP == "Overall effect"]<-23		
	Plot_dat$weight<-1/Plot_dat$se/5
	
	Plot_dat<-Plot_dat[Plot_dat$outcome2 !="Overall cancer (excluding non-melanoma skin cancer)",]
	Plot_dat<-Plot_dat[Plot_dat$outcome2 !="Respiratory and intrathoracic cancer",]
	# Plot_dat$outcome3<-gsub("\\(excluding non-melanoma skin cancer)","\n(excl. nm-skin cancer)",Plot_dat$outcome2)
	# Plot_dat$outcome3<-gsub("Respiratory and","Respiratory &\n",Plot_dat$outcome3)		
	Plot_dat$ncase.outcome2<- unlist(lapply(Plot_dat$outcome2,FUN=function(x)
		median(Plot_dat$ncase.outcome[Plot_dat$outcome2==x])))
	Plot_dat$outcome_name<-paste0(Plot_dat$outcome2,"\nN.cases=",Plot_dat$ncase.outcome2)		
	Plot_dat$outcome_name[Plot_dat$outcome2=="Overall cancer"]<-"Overall cancer\nN.cases=70223-101440"
	Plot_dat<-Plot_dat[order(Plot_dat$ncase.outcome2,decreasing=TRUE),]
	return(Plot_dat)
}


 