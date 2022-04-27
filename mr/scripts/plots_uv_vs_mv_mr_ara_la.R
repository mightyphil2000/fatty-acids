library(ggforestplot)
library(ggplot2)


# load("~/fatty-acids/mr/results/results_ara_la_mvmr.Rdata")
load("~/fatty-acids/mr/results/results_ara_la_mvmr_meta_analysis.Rdata")

res4$design<-"multivariable MR"
res4$model<-"adjusted for linoleic acid"
res4$model[res4$exposure=="Linoleic acid"]<-"adjusted for arachidonic acid"
res1<-res4
load("~/fatty-acids/mr/results/mr_ara_la.Rdata")

res4$design<-"univariable MR"
res4$model<-"unadjusted"
res4$model[res4$exposure=="Linoleic acid"]<-"unadjusted"
res2<-res4
res<-plyr::rbind.fill(res1,res2)


# Arachidonic acid 
Plot_dat<-res # results with independent outcomes combined by fixed effects meta analysis of MR result
Plot_dat<-format_plot(fa="Arachidonic acid")
P1<-forestplot(df =Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,colour =model ,xlab = "")

	 #   +
		# scale_size_manual(values=c(10,10,10,10,10,10,10,10,10,10,10,10,10,10))
		# scale_shape_manual(
		# 	size=5,
  #        values = c(23L),
  #        labels = c("Meta-analysis"))+



  #        , "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS")
  #      )
     

png("~/fatty-acids/mr/results/plots/ggforest_ara_uni_plus_mvmr.png", width = 900, height = 1000)
	print(P1) 
dev.off()

# Linoleic acid 
Plot_dat<-res # results with independent outcomes combined by fixed effects meta analysis of MR result
Plot_dat<-format_plot(fa="Linoleic acid")
P1<-forestplot(df =Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,colour =model ,xlab = "")

png("~/fatty-acids/mr/results/plots/ggforest_la_uni_plus_mvmr.png", width = 900, height = 1000)
	print(P1) 
dev.off()
# +
# 		geom_point(shape="square",size=10)+
# 		theme(text = element_text(size=20))
		# /1/Plot_data$se/5

		
		# fill=c(res2$Colour)
		# colour = mr_res$Colour

# ara and la
Plot_dat<-res # results with independent outcomes combined by fixed effects meta analysis of MR result
Plot_dat<-format_plot()

P1<-forestplot(df =Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=model,colour =Fatty_acid ,xlab = "")

png("~/fatty-acids/mr/results/plots/ggforest_ara_la_uni_plus_mvmr.png", width = 450, height = 500)
	print(P1) 
dev.off()

format_plot<-function(cancer=NULL,fa=NULL){
	
	if(!is.null(cancer))
	{
		Plot_dat<-Plot_dat[Plot_dat$outcome2==cancer,]
	}

	if(!is.null(fa))
	{
		Plot_dat<-Plot_dat[Plot_dat$exposure==fa,]
	}

	if(is.null(fa))
	{
		Plot_dat$model[grep("for",Plot_dat$model)]<-"adjusted"
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
	Plot_dat$Fatty_acid<-Plot_dat$exposure
	return(Plot_dat)
}
