library(ggforestplot)
library(ggplot2)

# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
load("~/fatty-acids/mr/results/mr_ara_la.Rdata")

# load("~/fatty-acids/mr/data/outcome_dat_ara_la.Rdata")
# load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")
# crc<-outcome_dat[outcome_dat$outcome2 == "Colorectal cancer",]
# crc[,c("SNP","study","beta.outcome","se.outcome")]
# ara_snps<-exposure_dat$SNP[exposure_dat$exposure=="Arachidonic acid"]
# crc<-crc[crc$SNP %in% ara_snps,]
# crc1<-meta_analysis_snp(dat=crc)
# res4<-res4[res4$population=="European",]


Plot_dat<-res4
Plot_dat<-format_plot()
P1<-forestplot(df = Plot_dat,logodds = TRUE,name=outcome_name,estimate=b, se=se,shape=NULL,colour = exposure,xlab = "",weight=Plot_dat$weight)

png("~/fatty-acids/mr/results/plots/ggforest_ara_la_univariableMR.png", width = 600, height = 480)
	print(P1) 
dev.off()

format_plot<-function(cancer=NULL){
	if(!is.null(cancer))
	{
		Plot_dat<-Plot_dat[Plot_dat$outcome2==cancer,]
	}
	# Plot_dat<-Plot_dat[Plot_dat$population =="European",]
	# Plot_dat<-Plot_dat[Plot_dat$exposure=="Arachidonic acid",]
	Plot_dat$b<-as.numeric(Plot_dat$b)
	Plot_dat$se<-as.numeric(Plot_dat$se)
	Plot_dat$ncase.outcome<-as.numeric(Plot_dat$ncase.outcome)
	Plot_dat<-Plot_dat[order(Plot_dat$ncase.outcome,decreasing=TRUE),]
	# Plot_dat$ncase.outcome<-as.numeric(Plot_dat$ncase.outcome)
	Plot_dat$shape<-15
	Plot_dat$weight<-1/Plot_dat$se/5
	# Plot_dat$name <- paste(Plot_dat$outcome2,"\n",Plot_dat$ncase.outcome)

	Plot_dat$outcome_name<-paste(Plot_dat$outcome2,"\nN. cases=",Plot_dat$ncase.outcome)
	Plot_dat<-Plot_dat[Plot_dat$outcome2 !="Overall cancer (excluding non-melanoma skin cancer)",]
	Plot_dat$outcome_name<-gsub("\\(excluding non-melanoma skin cancer)","\n(excl. nm-skin cancer)",Plot_dat$outcome_name)
	Plot_dat$outcome_name<-gsub("Respiratory and","Respiratory &\n",Plot_dat$outcome_name)		
	return(Plot_dat)
}
