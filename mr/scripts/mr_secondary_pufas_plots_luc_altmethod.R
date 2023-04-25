library(metafor)
library(ggforestplot)
library(ggplot2)
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
source("~/fatty-acids/mr/scripts/mr_functions.R")
# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads_altmethod.Rdata")
fads1<-res1
fads1$fads<-TRUE
load("~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_altmethod.Rdata")
fads2<-res1
fads2$fads<-FALSE
Dat<-rbind(fads1,fads2)
Plot_dat1<-format_plot2(dat=Dat,cancer="Lung cancer")
Plot_dat1$col<-"black"
Plot_dat1$Order<-1:nrow(Plot_dat1)
Plot_dat1$col[seq(2,nrow(Plot_dat1),2)]<-"red"
Plot_dat1$exposure_study[seq(2,nrow(Plot_dat1),2)]<-""
Plot_dat1$method <- "MR-pooled-study"


load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads_v3.Rdata")
fads1<-res4
#fads1[fads1$outcome2 == "Lung cancer",]
fads1$fads<-TRUE
load("~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_v3.Rdata")
fads2<-res4
fads2$fads<-FALSE
Dat<-rbind(fads1,fads2)
Plot_dat2<-format_plot(dat=Dat,cancer=c("Lung cancer"   ))
Plot_dat2$col<-"gray"
Plot_dat2$Order<-1:nrow(Plot_dat2)
Plot_dat2$col[seq(2,nrow(Plot_dat2),2)]<-"orange"
Plot_dat2$exposure_study[seq(2,nrow(Plot_dat2),2)]<-""
Plot_dat2$method <- "MR-by-study"
Plot_dat<-plyr::rbind.fill(Plot_dat1,Plot_dat2)
Min<-min(Plot_dat$lci)
Max<-max(Plot_dat$uci)
Plot_dat$FADS
# P<-forestplot(df = Plot_dat,logodds = TRUE,name=exposure,estimate=b,se=se,xlab = "",colour=FADS,shape=method)

png("~/fatty-acids/mr/results/plots/metaforest_luc_secondary_pufas_univariableMR_altmethod_for_reviewer.png", width = 1000, height = 700)
  forestplot(df = Plot_dat,logodds = TRUE,name=exposure,estimate=b,se=se,xlab = "",colour=FADS,shape=method)  
dev.off()

