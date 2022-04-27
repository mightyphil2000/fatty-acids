library(ggforestplot)
source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")	
mr_res<-format_mrresults()
Cancers<-disc.tab9$cancer[disc.tab9$power05>0.8]
Cancers2<-disc.tab9$cancer[disc.tab9$power05>0.8 & !disc.tab9$cancer %in% c("Colorectal cancer","Lung cancer","Malignant non-melanoma skin cancer","Colon cancer","Basal cell carcinoma","Malignant skin cancer","Rectal cancer","Distal colorectal cancer")]

mr_res<-mr_res[ mr_res$cancer %in% Cancers ,c("cancer","exposure","b","se","OR","LCI","UCI","pval","cases","population","plot_name","system","Cancer.Group")]
mr_res<-mr_res[!duplicated(mr_res$cancer),]
mr_res[!duplicated(mr_res$Cancer.Group),c("cancer","OR","LCI","UCI","pval")]
mr_res$Colour<-NA 
mr_res$Colour[mr_res$cancer %in% Cancers2]<-"black"
mr_res$Colour[!mr_res$cancer %in% Cancers2]<-"red"
mr_res$system[mr_res$system %in% c("Respiratory","Urinary")]<-"Other"
mr_res$system[mr_res$system == "Multiple"]<-"All cause"

# forestplot(df = mr_res,logodds = TRUE,name=plot_name,
# 	  estimate=b, se=se,shape=system,
# 		colour = Colour,xlab = "")

P1<-forestplot(df = mr_res,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = NULL,xlab = "")+
		geom_point(shape="square",size=1/mr_res$se/5,fill=c(mr_res$Colour),colour = mr_res$Colour)+
		theme(text = element_text(size=20))

png("~/fatty-acids/mr/results/plots/ggforest_cancers_power05.png", width = 1200, height = 1600)
	print(P1) 
dev.off()


# P1<-forestplot(df = mr_res2,logodds = TRUE,name=plot_name,
# 	  estimate=b, se=se,shape=NULL,
# 		colour = NULL,xlab = "")+
# geom_point(shape=mr_res2$Shape,size=1/mr_res2$se/10,fill=c(mr_res2$Colour),colour = mr_res2$Colour)+
# theme(legend.position = "none")+
# theme(text = element_text(size=20))

# unique(disc.tab9$cancer)
# disc.tab9$cancer[disc.tab9$power05>0.8 & !disc.tab9$cancer %in% c("Colorectal cancer","Lung cancer","Malignant non-melanoma skin cancer","Colon cancer","Basal cell carcinoma","Malignant skin cancer","Rectal cancer","Distal colorectal cancer")]


