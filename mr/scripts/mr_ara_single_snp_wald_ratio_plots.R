library(ggforestplot)
library(ggplot2)

# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
load("~/fatty-acids/mr/results/results_ara_snps_meta_analysis_wald_ratios.Rdata")

# load("~/fatty-acids/mr/data/outcome_dat_ara_la.Rdata")
# load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")
# crc<-outcome_dat[outcome_dat$outcome2 == "Colorectal cancer",]
# crc[,c("SNP","study","beta.outcome","se.outcome")]
# ara_snps<-exposure_dat$SNP[exposure_dat$exposure=="Arachidonic acid"]
# crc<-crc[crc$SNP %in% ara_snps,]
# crc1<-meta_analysis_snp(dat=crc)
# res4<-res4[res4$population=="European",]

# to diamond meta analysis shape
 # ggplot2::scale_shape_manual(
 #         values = c(23L, 21L, 21L, 21L, 21L),
 #         labels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS")
 #       )
     

# to make x axes same across different plots
# ggplot2::coord_cartesian(xlim = c(-0.3, 0.4))


# Colorectal cancer
dev.off()
meta_dat$lci<-round(exp(as.numeric(meta_dat$b)-1.96*as.numeric(meta_dat$se)),3)
meta_dat$uci<-round(exp(as.numeric(meta_dat$b)+1.96*as.numeric(meta_dat$se)),3)
meta_dat$or<-round(exp(as.numeric(meta_dat$b)),3)
Min<-min(meta_dat$lci)
Max<-max(meta_dat$uci)

meta_dat[meta_dat$SNP=="rs4985155" & meta_dat$outcome2 %in% c("Overall cancer","Colorectal cancer","Lung cancer","Non-melanoma skin cancer"," Malignant skin cancer","Basal cell carcinoma"),c("outcome2","SNP","or","lci","uci","pval")]

meta_dat[meta_dat$SNP=="rs174546" & meta_dat$outcome2 %in% c("Overall cancer","Colorectal cancer","Lung cancer","Non-melanoma skin cancer"," Malignant skin cancer","Basal cell carcinoma"),c("outcome2","SNP","or","lci","uci","pval")]


text_size<-15
title_size<-15
Plot_dat<-meta_dat
Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Colorectal cancer")

Plot_dat[,c("SNP","or","lci","uci","pval")]

P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",  title = Plot_dat$outcome_name3)+
	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
	coord_cartesian(xlim = c(Min, Max))+
		theme(text = element_text(size=text_size))+		
		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+
		scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")

		# theme(axis.text.x = element_text(face="bold", color="#993333",size=14, angle=45))+
		
png("~/fatty-acids/mr/results/plots/ggforest_crc_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
	print(P1) 
dev.off()


# Lung cancer
Plot_dat<-meta_dat
Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Lung cancer")
Plot_dat[,c("SNP","or","lci","uci","pval")]

P2<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",title = Plot_dat$outcome_name3)+
	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
	coord_cartesian(xlim = c(Min, Max))+
		theme(text = element_text(size=text_size))+
		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+
		scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")

		geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
		theme(text = element_text(size=20))


png("~/fatty-acids/mr/results/plots/ggforest_luc_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
print(P2) 
dev.off()

# "Basal cell carcinoma"
# rs174546 missing from UKB for basal cell carcinoma
# exclude BCC for now
unique(meta_dat$outcome2)
Plot_dat<-meta_dat
Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Basal cell carcinoma")
Plot_dat[,c("SNP","or","lci","uci","pval")]

P3<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",title = Plot_dat$outcome_name3)+
	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
	coord_cartesian(xlim = c(Min, Max))+
		theme(text = element_text(size=text_size))+
		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+
		scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")


png("~/fatty-acids/mr/results/plots/ggforest_bcc_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
	print(P3) 
dev.off()


# "Overall cancer"      
unique(meta_dat$outcome2)
Plot_dat<-meta_dat
Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Overall cancer")
Plot_dat
P4<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",title = Plot_dat$outcome_name3)+
	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
	coord_cartesian(xlim = c(Min, Max))+
		theme(text = element_text(size=text_size))+
		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+
		scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")

png("~/fatty-acids/mr/results/plots/ggforest_oc_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
	print(P4) 
dev.off()


# "Malignant skin cancer"  
unique(meta_dat$outcome2)
Plot_dat<-meta_dat
Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Malignant skin cancer"        )
P5<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",title = Plot_dat$outcome_name3)+
	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
	coord_cartesian(xlim = c(Min, Max))+
		theme(text = element_text(size=text_size))+
		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+
		scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")

png("~/fatty-acids/mr/results/plots/ggforest_msc_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
	print(P5) 
dev.off()


# # "Overall cancer (excluding non-melanoma skin cancer)"
# unique(meta_dat$outcome2)
# Plot_dat<-meta_dat
# Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Overall cancer (excluding non-melanoma skin cancer)"        )
# P6<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",title = Plot_dat$outcome_name3)+
# 	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
# 		theme(text = element_text(size=text_size))+
# 		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+
# 		scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")

# png("~/fatty-acids/mr/results/plots/ggforest_oc_excl_nmsc_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
# 	print(P6) 
# dev.off()

#"Non-melanoma skin cancer"  
unique(meta_dat$outcome2)
Plot_dat<-meta_dat
Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Non-melanoma skin cancer")
P7<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",title = Plot_dat$outcome_name3)+
	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
	coord_cartesian(xlim = c(Min, Max))+
		theme(text = element_text(size=text_size))+
		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")


png("~/fatty-acids/mr/results/plots/ggforest_nmsc_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
	print(P7) 
dev.off()

#"Respiratory and intrathoracic cancer" 
# unique(meta_dat$outcome2)
# Plot_dat<-meta_dat
# Plot_dat<-format_plot(rm_rs3734398=TRUE,cancer="Respiratory and intrathoracic cancer" )
# P8<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,estimate=b, se=se,shape=NULL,colour = NULL,xlab = "",title = Plot_dat$outcome_name3)+
# 	geom_point(shape=Plot_dat$shape,size=4,fill=c("black"))+
# 	coord_cartesian(xlim = c(Min, Max))+
# 		theme(text = element_text(size=text_size))+
# 		theme(plot.title = element_text(size = title_size),text = element_text(size=text_size))+
# 		scale_x_continuous(breaks=c(0.65,1.0,1.5,2.0),labels=c(0.65,1.0,1.5,2.0),trans="log")


# png("~/fatty-acids/mr/results/plots/ggforest_ric_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
# 	print(P8) 
# dev.off()

# cow plot / all plots 
# Plot_dat<-meta_dat
# Plot_dat<-format_plot(rm_rs3734398=TRUE)

# P_all<-forestplot(df = Plot_dat,logodds = TRUE,name=outcome_name2,estimate=b, se=se,shape=NULL,colour = SNP,xlab = "")

# png("~/fatty-acids/mr/results/plots/ggforest_all_rs4985155_rs174546_wald_ratios.png", width = 600, height = 480)
# 	print(P_all) 
# dev.off()


Plotlist<-list(P4,P1,P2,P7,P3,P5)
Plot4<-cowplot::plot_grid(plotlist=Plotlist,nrow=3,ncol=2)
png("~/fatty-acids/mr/results/plots/cow_plot_wald_ratio_ara.png", width = 1200, height = 1500)
	Plot4
dev.off()


format_plot<-function(rm_rs3734398=FALSE,cancer=NULL){
	Plot_dat<-Plot_dat[Plot_dat$SNP != "All - MR Egger",]
	Plot_dat$SNP[Plot_dat$SNP == "All - Inverse variance weighted"]<-"Overall effect"
	if(!is.null(cancer))
	{
		Plot_dat<-Plot_dat[Plot_dat$outcome2==cancer,]
	}
	
	Plot_dat<-Plot_dat[Plot_dat$population =="European",]
	Plot_dat<-Plot_dat[Plot_dat$exposure=="Arachidonic acid",]
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
	Plot_dat$ncase.outcome<-as.numeric(Plot_dat$ncase.outcome)
	Plot_dat$shape<-15
	Plot_dat$shape[Plot_dat$SNP == "Overall effect"]<-23		
	Plot_dat$weight<-1/Plot_dat$se/5
	
	Plot_dat<-Plot_dat[Plot_dat$outcome2 !="Overall cancer (excluding non-melanoma skin cancer)",]
	Plot_dat<-Plot_dat[Plot_dat$outcome2 !="Respiratory and intrathoracic cancer",]
	Plot_dat$outcome3<-gsub("\\(excluding non-melanoma skin cancer)","\n(excl. nm-skin cancer)",Plot_dat$outcome2)
	Plot_dat$outcome3<-gsub("Respiratory and","Respiratory &\n",Plot_dat$outcome3)		
	Plot_dat$ncase.outcome2<- unlist(lapply(Plot_dat$outcome3,FUN=function(x)
		median(Plot_dat$ncase.outcome[Plot_dat$outcome3==x])))
	Plot_dat$outcome_name2<-paste(Plot_dat$outcome3,"\nN. cases=",Plot_dat$ncase.outcome3)	
	Plot_dat$outcome_name3<-paste(Plot_dat$outcome2," (N. cases=",Plot_dat$ncase.outcome2,")")	
	Plot_dat<-Plot_dat[order(Plot_dat$ncase.outcome2,decreasing=TRUE),]
	return(Plot_dat)
}


 