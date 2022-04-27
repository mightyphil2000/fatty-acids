source("~/fatty-acids/mr/scripts/mr_functions.R")
# these scripts created forest plots for the individual SNP outcome effects unweighted by the ara effect size (ie are just the original sNP effect sizes. effect allele is the allele associated with higher ARA)
# this analysis is equivalent to an unweighted MR analysis
# only interested in this analysis in East Asians or meta analysis Europeans plus East Asians
library(ggforestplot)
library(ggplot2)


load("~/fatty-acids/mr/results/results_ara_snps_meta_analysis.Rdata")#meta analysis result only plus outcomes where only one study available (e.g. UKB overall cancer excluded non-menlanoma skin cancer and UKB non melanoma skin cancer)
load("~/fatty-acids/mr/results/results_ara_snps_meta_analysis2.Rdata") #meta analysis results plus individual study results



meta_dat2$lci<-round(exp(as.numeric(meta_dat2$b)-1.96*as.numeric(meta_dat2$se)),3)
meta_dat2$uci<-round(exp(as.numeric(meta_dat2$b)+1.96*as.numeric(meta_dat2$se)),3)
meta_dat2$or<-round(exp(as.numeric(meta_dat2$b)),3)

# Esophageal sCC

Plot_dat<-meta_dat2
Plot_dat<-Plot_dat[Plot_dat$outcome2=="Esophageal squamous cell carcinoma"    ,]
Plot_dat<-format_plot(rm_rs3734398=TRUE)
Plot_dat$lci<-round(exp(as.numeric(Plot_dat$b.x)-1.96*as.numeric(Plot_dat$se)),3)
Plot_dat$uci<-round(exp(as.numeric(Plot_dat$b.x)+1.96*as.numeric(Plot_dat$se)),3)
Plot_dat$or<-round(exp(as.numeric(Plot_dat$b.x)),3)


Min<-min(Plot_dat$lci)
Max<-max(Plot_dat$uci)
title_size<-10
text_size<-10
Plot_dat[,c("study","SNP","effect_allele.x","or","lci","uci","pval","outcome2")]



P3<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,
	  estimate=b.x, se=se,shape=study,color = NULL,xlab = "", title = Plot_dat$outcome_name)+
		coord_cartesian(xlim = c(Min, Max))+
		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))
		# geom_point(size=Plot_dat$weight)
png("~/fatty-acids/mr/results/plots/ggforest_escc_rs4985155_rs174546.png", width = 600, height = 480)
	print(P3) 
dev.off()


# Colorectal cancer
Plot_dat<-meta_dat2
Plot_dat<-Plot_dat[Plot_dat$outcome2=="Colorectal cancer",]
Plot_dat<-format_plot(rm_rs3734398=TRUE)

Plot_dat$lci<-round(exp(as.numeric(Plot_dat$b.x)-1.96*as.numeric(Plot_dat$se)),3)
Plot_dat$uci<-round(exp(as.numeric(Plot_dat$b.x)+1.96*as.numeric(Plot_dat$se)),3)
Plot_dat$or<-round(exp(as.numeric(Plot_dat$b.x)),3)


Min<-min(Plot_dat$lci)
Max<-max(Plot_dat$uci)
title_size<-10
text_size<-10
Plot_dat[,c("study","SNP","effect_allele.x","or","lci","uci","pval","outcome2")]


P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,
	  estimate=b.x, se=se,shape=study,color = NULL,xlab = "", title = Plot_dat$outcome_name)+
	coord_cartesian(xlim = c(Min, Max))+
	theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))
		# geom_point(size=Plot_dat$weight)

png("~/fatty-acids/mr/results/plots/ggforest_crc_rs4985155_rs174546.png", width = 600, height = 480)
	print(P1) 
dev.off()


# Lung cancer
Plot_dat<-meta_dat2
Plot_dat<-Plot_dat[Plot_dat$outcome2=="Lung cancer",]
Plot_dat<-format_plot(rm_rs3734398=TRUE)
Plot_dat$lci<-round(exp(as.numeric(Plot_dat$b.x)-1.96*as.numeric(Plot_dat$se)),3)
Plot_dat$uci<-round(exp(as.numeric(Plot_dat$b.x)+1.96*as.numeric(Plot_dat$se)),3)
Plot_dat$or<-round(exp(as.numeric(Plot_dat$b.x)),3)


Min<-min(Plot_dat$lci)
Max<-max(Plot_dat$uci)
title_size<-10
text_size<-10
Plot_dat[,c("study","SNP","effect_allele.x","or","lci","uci","pval","outcome2")]




P2<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,
	  estimate=b.x, se=se,shape=study,color = NULL,xlab = "", title = Plot_dat$outcome_name)+
	coord_cartesian(xlim = c(Min, Max))+
	theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))
png("~/fatty-acids/mr/results/plots/ggforest_luc_rs4985155_rs174546.png", width = 600, height = 480)
	print(P2) 
dev.off()


Plotlist<-list(P1,P2,P3)
Plot4<-cowplot::plot_grid(plotlist=Plotlist,nrow=3,ncol=1)
png("~/fatty-acids/mr/results/plots/cow_plot_single_snp_ea_eur_unweighted_mr.png", width = 600, height = 1000)
	Plot4
dev.off()



format_plot<-function(rm_rs3734398=FALSE){
	if(rm_rs3734398)
	{
		Plot_dat<-Plot_dat[!Plot_dat$SNP %in% "rs3734398",]
	}
	Plot_dat1<-Plot_dat[Plot_dat$study == "Overall fixed effect",]
	Plot_dat2<-Plot_dat[Plot_dat$study != "Overall fixed effect",]
	Plot_dat<-rbind(Plot_dat1,Plot_dat2)
	Plot_dat<-Plot_dat[order(Plot_dat$SNP),]
	Plot_dat$b<-as.numeric(Plot_dat$b)
	Plot_dat$se<-as.numeric(Plot_dat$se)
	Plot_dat$eaf<-as.numeric(Plot_dat$eaf)
	Plot_dat<-ara_increasing_effect_allele(Plot_dat=Plot_dat)	
	Plot_dat$shape<-15
	Plot_dat$shape[Plot_dat$study == "Overall fixed effect"]<-23		
	Plot_dat$weight<-1/Plot_dat$se/5

	Plot_dat$ncase_total<-NA
	Plot_dat$ncase_total<- unique(as.numeric(unlist(lapply(Plot_dat$outcome2,FUN=function(x)	Plot_dat$ncase.outcome[Plot_dat$outcome2==x & Plot_dat$study=="Overall fixed effect"]))))
	Plot_dat$outcome_name<-paste(Plot_dat$outcome2," N. cases=",Plot_dat$ncase_total)
	Plot_dat$study<-gsub("GECCO/CORECT/CCFR","GECCO/CORECT/\nCCFR",Plot_dat$study)
	return(Plot_dat)
}





# # Basal cell carcinoma
# Plot_dat<-meta_dat2
# Plot_dat<-Plot_dat[Plot_dat$outcome2=="Basal cell carcinoma",]
# Plot_dat<-format_plot(rm_rs3734398=TRUE)
# P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,
# 	  estimate=b.x, se=se,shape=study,color = NULL,xlab = "")
# 		# geom_point(size=Plot_dat$weight)
# png("~/fatty-acids/mr/results/plots/ggforest_bcc_rs4985155_rs174546.png", width = 600, height = 480)
# 	print(P1) 
# dev.off()

# # "Overall cancer"     
# Plot_dat<-meta_dat2
# Plot_dat<-Plot_dat[Plot_dat$outcome2=="Overall cancer" ,]
# Plot_dat<-format_plot(rm_rs3734398=TRUE)
# P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b.x, se=se,shape=study,color = NULL,xlab = "")

# png("~/fatty-acids/mr/results/plots/ggforest_oc_rs4985155_rs174546.png", width = 600, height = 480)
# 	print(P1) 
# dev.off()


# # "Malignant skin cancer"             
# Plot_dat<-meta_dat2
# Plot_dat<-Plot_dat[Plot_dat$outcome2=="Malignant skin cancer"        ,]
# Plot_dat<-format_plot(rm_rs3734398=TRUE)
# P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b.x, se=se,shape=study,color = NULL,xlab = "")

# png("~/fatty-acids/mr/results/plots/ggforest_msc_rs4985155_rs174546.png", width = 600, height = 480)
# 	print(P1) 
# dev.off()


# # "Respiratory and intrathoracic cancer"      
# Plot_dat<-meta_dat2
# Plot_dat<-Plot_dat[Plot_dat$outcome2=="Respiratory and intrathoracic cancer",]
# Plot_dat<-format_plot(rm_rs3734398=TRUE)
# P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b.x, se=se,shape=study,color = NULL,xlab = "")

# png("~/fatty-acids/mr/results/plots/ggforest_ric_rs4985155_rs174546.png", width = 600, height = 480)
# 	print(P1) 
# dev.off()


# # "Overall cancer (excluding non-melanoma skin cancer)"
# Plot_dat<-meta_dat2
# Plot_dat<-Plot_dat[Plot_dat$outcome2=="Overall cancer (excluding non-melanoma skin cancer)",]
# Plot_dat<-format_plot(rm_rs3734398=TRUE)
# P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b.x, se=se,shape=study,color = NULL,xlab = "")

# png("~/fatty-acids/mr/results/plots/ggforest_oc_excl_nmsc_rs4985155_rs174546.png", width = 600, height = 480)
# 	print(P1) 
# dev.off()

# #  "Non-melanoma skin cancer"  
# Plot_dat<-meta_dat2
# Plot_dat<-Plot_dat[Plot_dat$outcome2== "Non-melanoma skin cancer"  ,]
# Plot_dat<-format_plot(rm_rs3734398=TRUE)
# P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b.x, se=se,shape=study,color = NULL,xlab = "")

# png("~/fatty-acids/mr/results/plots/ggforest_nmsc_rs4985155_rs174546.png", width = 600, height = 480)
# 	print(P1) 
# dev.off()