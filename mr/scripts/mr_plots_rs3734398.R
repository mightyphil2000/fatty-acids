source("~/fatty-acids/mr/scripts/mr_functions.R")
library(ggforestplot)
library(ggplot2)

# load("~/fatty-acids/mr/results/results_ara_snps_meta_analysis.Rdata")#meta analysis result only plus outcomes where only one study available (e.g. UKB overall cancer excluded non-menlanoma skin cancer and UKB non melanoma skin cancer)
load("~/fatty-acids/mr/results/results_ara_snps_meta_analysis2.Rdata") #meta analysis results plus individual study results


# plot meta analysis of rs3734398
Plot_dat<-format_meta_analysis_rs3734398()

Plot_dat[,c("b","se","study","pval")]
P1<-forestplot(df = Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,color =NULL,xlab = "", title = NULL)+geom_point(shape=Plot_dat$shape,size=Plot_dat$weight,fill=c("black"))
png("~/fatty-acids/mr/results/plots/ggforest_rs3734398_unweightmr.png", width = 600, height = 480)
	print(P1) 
dev.off()

# ELOVL2 activity SD change per copy of allele T EAF=0.57	T	C/T
# source("raw_rs3734398_meta_analysis_europeans")
load("~/fatty-acids/mr/results/results_rs3734398_meta_analysis_europeans.Rdata")
Plot_dat<-format_meta_analysis_rs3734398_wald_ratio()
Plot_dat$weight<-1/Plot_dat$se.elovl2/2
# Plot of overall cancer, crc and lung cancer Wald ratios in Europeans plus combined effect 
P1<-forestplot(df = Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b.elovl2, se=se.elovl2,shape=NULL,color =NULL,xlab = "", title = NULL)+geom_point(shape=Plot_dat$shape,size=Plot_dat$weight,fill=c("black"))

png("~/fatty-acids/mr/results/plots/ggforest_rs3734398_wald_ratios_maxpower.png", width = 600, height = 480)
	print(P1) 
dev.off()

# Plot of individual cancer wald ratios overall cancer, crc and lung cancer Wald ratios in Europeans plus combined effect 
load("~/fatty-acids/mr/results/results_rs3734398_meta_analysis_europeans.Rdata")
meta_dat<-elovl2_function(Plot_dat=meta_dat)
meta_dat<-format_plot(meta_dat=meta_dat)
meta_dat$weight<-1/meta_dat$se.elovl2/2

P1<-forestplot(df = meta_dat,logodds = TRUE,name=outcome_name,
	  estimate=b.elovl2, se=se.elovl2,shape=NULL,color =NULL,xlab = "", title = NULL)+geom_point(shape=meta_dat$shape,size=meta_dat$weight,fill=c("black"))
png("~/fatty-acids/mr/results/plots/ggforest_rs3734398_wald_ratios.png", width = 600, height = 480)
	print(P1) 
dev.off()

# Plot of individual cancer wald ratios plus combined effect across overall cancer, crc and lung cancer
load("~/fatty-acids/mr/results/results_rs3734398_meta_analysis_europeans.Rdata")
meta_dat1<-format_meta_analysis_rs3734398_wald_ratio()
meta_dat1<-meta_dat1[meta_dat1$study=="Overall effect" ,]
meta_dat1$outcome2<-"Overall cancer2"
load("~/fatty-acids/mr/results/results_rs3734398_meta_analysis_europeans.Rdata")
meta_dat2<-elovl2_function(Plot_dat=meta_dat)
meta_dat<-plyr::rbind.fill(meta_dat1,meta_dat2)
meta_dat<-format_plot(meta_dat=meta_dat)
meta_dat$weight<-1/meta_dat$se.elovl2/4
meta_dat$shape<-15
P1<-forestplot(df = meta_dat,logodds = TRUE,name=outcome_name,
	  estimate=b.elovl2, se=se.elovl2,shape=NULL,color =NULL,xlab = "", title = NULL)+geom_point(shape=meta_dat$shape,size=meta_dat$weight,fill=c("black"))

png("~/fatty-acids/mr/results/plots/ggforest_rs3734398_wald_ratios2.png", width = 600, height = 480)
	print(P1) 
dev.off()


format_plot<-function(meta_dat=NULL){
	meta_dat$shape<-15
	meta_dat$shape[meta_dat$study == "Overall effect"]<-23		
	meta_dat$weight<-1/meta_dat$se.elovl2/20	
	meta_dat$outcome2<-gsub(" and ","&\n",meta_dat$outcome2)
	meta_dat<-meta_dat[meta_dat$outcome2 != "Overall cancer (excluding non-melanoma skin cancer)",]
	meta_dat<-meta_dat[order(as.numeric(meta_dat$ncase.outcome),decreasing=TRUE),]	
	meta_dat$outcome_name<-paste0(meta_dat$outcome2,"\nN. cases=",meta_dat$ncase.outcome)
	return(meta_dat)
}