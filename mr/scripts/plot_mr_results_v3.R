library(ggforestplot)
library(ggplot2)
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
plot_dat<-format_all_discovery_v2()
plot_dat<-format_plot_dat()
table(plot_dat$system)
Min<-min(plot_dat$LCI)
Max<-max(plot_dat$UCI)
Max<-Max+0.01
Min<-Min-0.01

plot_dat$Colour<-rep("black",nrow(plot_dat))
plot_dat$Colour[which(plot_dat$pval<0.05/67)]<-"red"
plot_dat$name2<-paste0(plot_dat$OR," (",plot_dat$LCI,"-",plot_dat$UCI,")")

plot_dat2[plot_dat2$cancer == "Colorectal",]

plot_dat1<-plot_dat[plot_dat$system %in% c("AReproductive","DRespiratory","EIntegumentary"),]
Colour1<-plot_dat1$Colour

plot_dat2<-plot_dat[!plot_dat$system %in% c("AReproductive","DRespiratory","EIntegumentary"),]
plot_dat2$outcome[grep("hodg",plot_dat2$outcome,ignore.case=TRUE)]
Colour2<-plot_dat2$Colour

# plot_dat$cancer[which(plot_dat$pval<0.05/67)]


P1<-forestplot(df = plot_dat1,logodds = TRUE,name=outcome,estimate=b,se=se,xlab = "", psignif = 0.05/nrow(plot_dat1))+theme(legend.position = "none")+
	geom_point(shape=15,size=1/plot_dat1$se/14,fill=Colour1,colour=Colour1) +scale_x_continuous(limits=c(Min, Max),trans='log10',breaks=c(0.85,1.0,1.25),labels=c("0.8","1.0","1.25"))+theme(plot.title = element_text(size = ""),text = element_text(size=15))

# P1
  # scale_x_continuous(trans='log10')
P1<-P1+ggforce::facet_col(
    facets = ~system,
    scales = "free_y",
    space = "free" )
P1

P2<-forestplot(df = plot_dat2,logodds = TRUE,name=outcome,estimate=b,se=se,xlab = "", psignif = 0.05/nrow(plot_dat2))+theme(legend.position = "none")+
	geom_point(shape=15,size=1/plot_dat2$se/14,fill=Colour2,colour=Colour2)  +scale_x_continuous(limits=c(Min, Max),trans='log10',breaks=c(0.85,1.0,1.25),labels=c("0.8","1.0","1.25"))+theme(plot.title = element_text(size = ""),text = element_text(size=15))
P2<-P2+ggforce::facet_col(
    facets = ~system,
    scales = "free_y",
    space = "free" )
P2

# geom_point(shape=mr_res$Shape,size=1/mr_res$se/10,fill=c(mr_res$Colour),colour = mr_res$Colour)+



png("~/fatty-acids/mr/results/plots/mr_results_all1_mr_results_rep_v4.png",width = 500, height = 1000)
	P1
dev.off()

png("~/fatty-acids/mr/results/plots/mr_results_all2_mr_results_rep_v4.png",width = 500, height = 1000)
	P2
dev.off()


###################################
# add in columns for or and 95% ci#
###################################

P1<-forestplot(df = plot_dat1,logodds = TRUE,name=name2,estimate=b,se=se,xlab = "", psignif = 0.05/nrow(plot_dat1))+theme(legend.position = "none")+
	geom_point(shape=15,size=1/plot_dat1$se/14,fill=Colour1,colour=Colour1) +scale_x_continuous(limits=c(Min, Max),trans='log10',breaks=c(0.85,1.0,1.25),labels=c("0.8","1.0","1.25"))+theme(plot.title = element_text(size = ""),text = element_text(size=15))
# P1
  # scale_x_continuous(trans='log10')
P1<-P1+ggforce::facet_col(
    facets = ~system,
    scales = "free_y",
    space = "free" )
P1

P2<-forestplot(df = plot_dat2,logodds = TRUE,name=name2,estimate=b,se=se,xlab = "", psignif = 0.05/nrow(plot_dat2))+theme(legend.position = "none")+
	geom_point(shape=15,size=1/plot_dat2$se/14,fill=Colour2,colour=Colour2)  +scale_x_continuous(limits=c(Min, Max),trans='log10',breaks=c(0.85,1.0,1.25),labels=c("0.8","1.0","1.25"))+theme(plot.title = element_text(size = ""),text = element_text(size=15))
P2<-P2+ggforce::facet_col(
    facets = ~system,
    scales = "free_y",
    space = "free" )
P2

png("~/fatty-acids/mr/results/plots/mr_results_all1_or95ci_mr_results_rep_v3.png",width = 500, height = 1000)
	P1
dev.off()

png("~/fatty-acids/mr/results/plots/mr_results_all2_or95ci_mr_results_rep_v3.png",width = 500, height = 1000)
	P2
dev.off()


# load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
# load("~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")
# mr_res1[mr_res1$cancer == "Colorectal cancer",c("cases",)]
# load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")

# library(devtools)
# library(TwoSampleMR)
# library(plyr)

# library(tidyverse)

# plot_dat<-format_all_discovery(Power=NULL)

# plot_dat<-plot_dat[order(plot_dat$Category),]
# plot_dat1<-plot_dat[1:41,]
# plot_dat2<-plot_dat[42:81,]
# plot_dat1[plot_dat1$outcome=="Esophageal squamous cell carcinoma\ncases=2013",c("or","lci","uci" )]
# table(plot_dat$Category)




