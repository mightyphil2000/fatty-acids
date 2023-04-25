library(metafor)
library(ggforestplot)
library(ggplot2)
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
source("~/fatty-acids/mr/scripts/mr_functions.R")
# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads_v3.Rdata")
fads1<-res4
#fads1[fads1$outcome2 == "Lung cancer",]
fads1$fads<-TRUE
load("~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_v3.Rdata")
fads2<-res4
fads2$fads<-FALSE
Dat<-rbind(fads1,fads2)
Plot_dat<-format_plot(dat=Dat,cancer=c("Lung cancer"   ))
Plot_dat[Plot_dat$exposure =="Omega-6 fatty acids",]
dim(Plot_dat)
Plot_dat$study.abbreviation
Plot_dat$ncase
Plot_dat$ncontrol

Min<-min(Plot_dat$lci)
Max<-max(Plot_dat$uci)
Plot_dat$col<-"black"
Plot_dat[,c("exposure","nsnp","col")]
Plot_dat$Order<-1:nrow(Plot_dat)
Plot_dat$col[seq(2,nrow(Plot_dat),2)]<-"red"

Plot_dat$exposure_study[seq(2,nrow(Plot_dat),2)]<-""


png("~/fatty-acids/mr/results/plots/metaforest_luc_secondary_pufas_univariableMR_v4.png", width = 1000, height = 700)
  forest(x=Plot_dat$b,sei=Plot_dat$se,
    xlim=c(-2.0,1.0),at=log(c(Min-0.01,1,Max+0.01)),atransf=exp,
    slab=Plot_dat$exposure_study,
    ilab=do.call(cbind,list(Plot_dat$nsnp,
    Plot_dat$pval,Plot_dat$Q_pval)),ilab.xpos=c(-1.15,-0.9,-0.6), ylim=c(2, 33),col=Plot_dat$col,rows=c(30:29,27:26,24:23,21:20,18:17,15:14,12:11,9:8,6:5,3:2),header=c("","Odds ratio (95% confidence interval) per SD increase in genetically proxied PUFA"),top=2,xlab="")+text(c(-1.15,-0.9,-0.6), 33.1, c("No.\nSNPs", "P\nvalue", "Q\nP value"))
dev.off()


# Plot_dat[Plot_dat$exposure=="Gamma-linolenic acid (18:3n6)"  ,c("outcome2","FADS","or","lci","uci")]
# Plot_dat[Plot_dat$exposure=="Dihomo-gamma-linolenic acid (20:3n6)" ,c("outcome2","FADS","or","lci","uci")]



# Plot_dat$Q_pval_name
# # P1<-forestplot(df = Plot_dat,logodds = TRUE,name=exposure,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=20))
# # P1<-P1+ggforce::facet_col(
# #     facets = ~PUFA,
# #     scales = "free_y",
# #     space = "free" )

# # P1<-forestplot(df = Plot_dat,logodds = TRUE,name=exposure2,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=15),legend.position="bottom")
# # P1<-P1+ggforce::facet_col(
# #     facets = ~PUFA,
# #     scales = "free_y",
# #     space = "free" )

# Plot_dat<-tweak_plot_dat()

# P5<-forestplot(df = Plot_dat,logodds = TRUE,name=exposure3,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=15),legend.position="bottom")
# P5<-P5+ggforce::facet_col(
#     facets = ~PUFA,
#     scales = "free_y",
#     space = "free" )


# P4<-forestplot(df = Plot_dat,logodds = TRUE,name=nsnps_name,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=15),legend.position="bottom")
# P4<-P4+ggforce::facet_col(
#     facets = ~PUFA,
#     scales = "free_y",
#     space = "free" )

# P2<-forestplot(df = Plot_dat,logodds = TRUE,name=pval,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=15),legend.position="bottom")
# P2<-P2+ggforce::facet_col(
#     facets = ~PUFA,
#     scales = "free_y",
#     space = "free" )

# P3<-forestplot(df = Plot_dat,logodds = TRUE,name=Q_pval_name,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=15),legend.position="bottom")
# P3<-P3+ggforce::facet_col(
#     facets = ~PUFA,
#     scales = "free_y",
#     space = "free" )



# png("~/fatty-acids/mr/results/plots/ggforest_luc_secondary_pufas_univariableMR_v3.png", width = 600, height = 700)
# 	print(P5) 
# dev.off()

# png("~/fatty-acids/mr/results/plots/ggforest_luc_secondary_pufas_univariableMR_v3_nsnps.png", width = 600, height = 700)
# 	print(P4) 
# dev.off()

# png("~/fatty-acids/mr/results/plots/ggforest_luc_secondary_pufas_univariableMR_v3_pval.png", width = 600, height = 700)
# 	print(P2) 
# dev.off()

# png("~/fatty-acids/mr/results/plots/ggforest_luc_secondary_pufas_univariableMR_v3_Q_pval.png", width = 600, height = 700)
# 	print(P3) 
# dev.off()

# message(unique(Plot_dat$exposure2))
# # png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3.png", width = 1000, height = 600)
# # 	print(P1) 
# # dev.off()

# tweak_plot_dat<-function()
# {
#     Plot_dat$nsnps_name[Plot_dat$exposure =="Omega-6 fatty acids"]<-gsub("44","43",Plot_dat$nsnps_name[Plot_dat$exposure =="Omega-6 fatty acids"])

#     # The issue here is that the effect is based on a meta-analysis of the MR result across UKB and ILCCO. For omega-6 fatty acids, the number of SNPs was 40 and 47 in ILCCO and UKB respectively (FADS included) and was 39 and 46 in ILCCO and UKB (FADS excluded). The number of SNPs is the average number of SNPs across the two, which gives 43.5 and 42. 43.5 is rounded to 44 by default by the round function, so becomes 44 which doesn't makes sense since it implies two SNPs were excluded when only one was excluded. 
#     # rounded 43.5 down to 43 rather than 44 (default in R is to round to nearest even number, hence 23.5 gets rounded down to 22 and not 24)
    
#    Plot_dat$nsnps_name[Plot_dat$exposure =="Docosahexaenoic acid (22:6n3)"]<-gsub("24","23",Plot_dat$nsnps_name[Plot_dat$exposure =="Docosahexaenoic acid (22:6n3)"])
#      # similar issue applied to DHA in UKB and ILCCO; rounding 23.5 down to 23 rather than 24
#     return(Plot_dat)
# }
