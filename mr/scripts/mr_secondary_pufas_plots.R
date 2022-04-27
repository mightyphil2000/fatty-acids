library(ggforestplot)
library(ggplot2)
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads_v2.Rdata")
fads1<-res4

fads1$fads<-TRUE
load("~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_v2.Rdata")
fads2<-res4
fads2$fads<-FALSE
Dat<-rbind(fads1,fads2)

# [8] 
unique(Dat$outcome2)
Plot_dat<-format_plot(dat=Dat,cancer=c("Colorectal cancer"   ))
Plot_dat<-Plot_dat[order(Plot_dat$PUFA,Plot_dat$chain.length),]
# Plot_dat<-Plot_dat[order(Plot_dat$chain.length),]
# Plot_dat<-Plot_dat[Plot_dat]

Dups<-unique(Plot_dat$exposure[duplicated(Plot_dat$exposure)])
Plot_dat<-Plot_dat[Plot_dat$exposure %in% Dups,]
Plot_dat$FADS<-Plot_dat$fads
Plot_dat$FADS<-""
Plot_dat$FADS[Plot_dat$fads]<-"included"
Plot_dat$FADS[!Plot_dat$fads]<-"excluded"
# Plot_dat<-Plot_dat[!Plot_dat$exposure %in% c("Omega-3 fatty acids","Omega-6 fatty acids"),]
Plot_dat<-Plot_dat[Plot_dat$PUFA !="other",]
P1<-forestplot(df = Plot_dat,logodds = TRUE,name=exposure,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=20))
P1<-P1+ggforce::facet_col(
    facets = ~PUFA,
    scales = "free_y",
    space = "free" )

png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3.png", width = 600, height = 600)
	print(P1) 
dev.off()


# png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3.png", width = 1000, height = 600)
# 	print(P1) 
# dev.off()

