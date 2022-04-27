library(ggforestplot)
library(ggplot2)

load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
load("~/fatty-acids/mr/data/outcome_dat_ara_la.Rdata")
load("~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata")
crc<-outcome_dat[outcome_dat$outcome2 == "Colorectal cancer",]
# crc[,c("SNP","study","beta.outcome","se.outcome")]
ara_snps<-exposure_dat$SNP[exposure_dat$exposure=="Arachidonic acid"]
crc<-crc[crc$SNP %in% ara_snps,]
meta_analysis_snp(dat=crc)
# res4<-res4[res4$population=="European",]
Plot_dat<-res_single2
Plot_dat$b<-as.numeric(res4$b)
Plot_dat$se<-as.numeric(res4$se)
Plot_dat<-Plot_dat[Plot_dat$exposure=="Arachidonic acid",]
Plot_dat<-Plot_dat[Plot_dat$outcome2=="Colorectal cancer",]
Plot_dat<-Plot_dat[Plot_dat$method!="Colorectal cancer",]
head(Plot_dat)
P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,
	  estimate=b, se=se,shape=NULL,
		colour = study.abbreviation,xlab = "")

head(Plot_dat)
Plot_dat<-res_single2[res_single2$outcome2=="Colorectal cancer" & res_single2$population == "European",]
Plot_dat<-res_single2[res_single2$outcome2=="Colorectal cancer" & res_single2$population == "East Asian",]
Plot_dat<-Plot_dat[order(Plot_dat$exposure),]
P1<-forestplot(df = Plot_dat,logodds = TRUE,name=SNP,
	  estimate=b, se=se,shape=study.abbreviation,
		colour = exposure,xlab = "")
		
		# geom_point(shape="square",size=1/res2$se/5,fill=c(res2$Colour),colour = mr_res$Colour)+
		# theme(text = element_text(size=20))

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