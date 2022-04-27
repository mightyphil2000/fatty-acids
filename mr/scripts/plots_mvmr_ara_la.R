library(ggforestplot)
library(ggplot2)


# load("~/fatty-acids/mr/results/results_ara_la_mvmr.Rdata")
load("~/fatty-acids/mr/results/results_ara_la_mvmr_meta_analysis.Rdata")

head(res4)
Plot_dat<-res4 # results with independent outcomes combined by fixed effects meta analysis of MR result
# Plot_data<-res2 #no meta analysis
Plot_dat<-format_plot(fa="Arachidonic acid")
Plot_dat<-format_plot(fa="Linoleic acid")

P1<-forestplot(df =Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,colour =NULL ,xlab = "")+
		geom_point(shape="square",size=10)+
		theme(text = element_text(size=20))
		# /1/Plot_data$se/5

		
		# fill=c(res2$Colour)
		# colour = mr_res$Colour

png("~/fatty-acids/mr/results/plots/ggforest_ara_la_mvmr.png", width = 900, height = 1000)
	print(P1) 
dev.off()

Plot_dat<-res4
Plot_dat<-format_plot()
P1<-forestplot(df =Plot_dat,logodds = TRUE,name=outcome_name,
	  estimate=b, se=se,shape=NULL,colour =exposure ,xlab = "")
png("~/fatty-acids/mr/results/plots/ggforest_ara_la_mvMR.png", width = 600, height = 480)
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
