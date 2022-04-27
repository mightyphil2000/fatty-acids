library(ggforestplot)
library(ggplot2)

# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
# load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads.Rdata")
load("~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_sensitivityanalyses.Rdata")
fads1<-res4
fads1$fads<-TRUE

fads2<-res2all
fads2$fads<-FALSE
# Dat<-rbind(fads1,fads2)
Dat<-fads2

Plot_dat<-format_plot(dat=Dat,cancer=c("Colorectal cancer"))
Plot_dat<-Plot_dat[order(Plot_dat$PUFA,Plot_dat$chain.length),]
Plot_dat<-Plot_dat[Plot_dat$nsnp>=5,]
# Plot_dat<-Plot_dat[order(Plot_dat$chain.length),]
# Plot_dat<-Plot_dat[Plot_dat]

# Dups<-unique(Plot_dat$exposure[duplicated(Plot_dat$exposure)])
# Plot_dat<-Plot_dat[Plot_dat$exposure %in% Dups,]
Plot_dat<-Plot_dat[!Plot_dat$exposure %in% unique(Plot_dat$exposure)[5:6],]
Plot_dat<-Plot_dat[!Plot_dat$method %in% unique(Plot_dat$method)[c(1,5)],]
Plot_dat$FADS<-Plot_dat$fads
Plot_dat$FADS<-""
Plot_dat$FADS[Plot_dat$fads]<-"included"
Plot_dat$FADS[!Plot_dat$fads]<-"excluded"
P1<-forestplot(df = Plot_dat,logodds = TRUE,name=exposure,estimate=b, se=se,shape=study.abbreviation,colour = method,xlab = "",weight=Plot_dat$weight)

png("~/fatty-acids/mr/results/plots/ggforest_ric_secondary_pufas_univariableMR.png", width = 600, height = 600)
	print(P1) 
dev.off()

format_plot<-function(dat=NULL,cancer=NULL){
	if(!is.null(cancer))
	{
		dat<-dat[dat$outcome2 %in%cancer,]
	}
	dat$PUFA<-NA
	dat$PUFA[grep("n3",dat$exposure)]<-"Omega 3"
	dat$PUFA[dat$exposure=="Omega-3 fatty acids"]<-"Omega 3"
	dat$PUFA[grep("n6",dat$exposure)]<-"Omega 6"
	dat$PUFA[dat$exposure=="Omega-6 fatty acids"]<-"Omega 6"
	dat$PUFA[grep("n3 or n6",dat$exposure)]<-"other"
	dat$PUFA[is.na(dat$PUFA)] <- "other"

	dat1<-dat[grep(":",dat$exposure),]
	dat1<-dat1[dat1$exposure!= "Other polyunsaturated fatty acids than 18:2",]
	dat2<-dat[grep(":",dat$exposure,invert=TRUE),]
	dat2<-rbind(dat2,dat[dat$exposure== "Other polyunsaturated fatty acids than 18:2",])
	Pos<-gregexpr(":",dat1$exposure)
	dat1$chain.length<-as.numeric(substr(dat1$exposure,start=unlist(Pos)-2,stop=unlist(Pos)-1))
	dat2$chain.length<-NA
	dat<-rbind(dat1,dat2)
	# dat<-dat[dat$population =="European",]
	# dat<-dat[dat$exposure=="Arachidonic acid",]
	dat$b<-as.numeric(dat$b)
	dat$se<-as.numeric(dat$se)
	dat$ncase.outcome<-as.numeric(dat$ncase.outcome)
	dat<-dat[order(dat$ncase.outcome,decreasing=TRUE),]
	# dat$ncase.outcome<-as.numeric(dat$ncase.outcome)
	dat$shape<-15
	dat$weight<-1/dat$se/5
	# dat$name <- paste(dat$outcome2,"\n",dat$ncase.outcome)

	dat$outcome_name<-paste(dat$outcome2,"\nN. cases=",dat$ncase.outcome)
	# dat<-dat[dat$outcome2 !="Overall cancer (excluding non-melanoma skin cancer)",]
	dat$outcome_name<-gsub("\\(excluding non-melanoma skin cancer)","\n(excl. nm-skin cancer)",dat$outcome_name)
	dat$outcome_name<-gsub("Respiratory and","Respiratory &\n",dat$outcome_name)		
	return(dat)
}
