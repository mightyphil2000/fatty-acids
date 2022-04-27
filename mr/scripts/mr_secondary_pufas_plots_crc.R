# install.packages("forestplot")
# install.packages("metafor")
library(metafor)
library(forestplot)
# library(ggforestplot)
library(ggplot2)
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
source("~/fatty-acids/mr/scripts/mr_functions.R")
# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads_v3.Rdata")
fads1<-res4


fads1$fads<-TRUE
load("~/fatty-acids/mr/results/mr_secondary_pufas_exclfads_v3.Rdata")

res4[res4$exposure =="Eicosapentaenoic acid (20:5n3)", ]

res4[which(res4$exposure =="Omega-3 fatty acids" & res4$outcome=="Colorectal cancer | 60"), c("exposure","outcome","b","se","pval","nsnp","method")]

fads2<-res4
fads2$fads<-FALSE
Dat<-rbind(fads1,fads2)

Dat[Dat$exposure =="Eicosapentaenoic acid (20:5n3)", c("fads","outcome2")]
table(Dat$exposure[Dat$id.outcome == "60"])

Plot_dat<-format_plot(dat=Dat,cancer=c("Colorectal cancer"   ))
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


png("~/fatty-acids/mr/results/plots/metaforest_crc_secondary_pufas_univariableMR_v4.png", width = 1000, height = 700)
  forest(x=Plot_dat$b,sei=Plot_dat$se,
    xlim=c(-2.0,1.0),at=log(c(Min,1,Max+0.01)),atransf=exp,
    slab=Plot_dat$exposure_study,
    ilab=do.call(cbind,list(Plot_dat$nsnp,
    Plot_dat$pval,Plot_dat$Q_pval)),ilab.xpos=c(-1.15,-0.9,-0.6), ylim=c(2, 33),col=Plot_dat$col,rows=c(30:29,27:26,24:23,21:20,18:17,15:14,12:11,9:8,6:5,3:2),header=c("","Odds ratio (95% confidence interval) per SD increase in genetically proxied PUFA"),top=2,xlab="")+text(c(-1.15,-0.9,-0.6), 33.1, c("No.\nSNPs", "P\nvalue", "Q\nP value"))
dev.off()

# ,ilab.pos=4
# forestplot(df = Plot_dat,logodds = TRUE,name=exposure3,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=15),legend.position="bottom")

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

# Plot_dat$Q_pval
# Plot_dat$Q_pval_name
# P3<-forestplot(df = Plot_dat,logodds = TRUE,name=Q_pval_name,estimate=b, se=se,shape=FADS,colour = PUFA,xlab = "",weight=Plot_dat$weight)+theme(plot.title = element_text(size = ""),text = element_text(size=15),legend.position="bottom")
# P3<-P3+ggforce::facet_col(
#     facets = ~PUFA,
#     scales = "free_y",
#     space = "free" )


# png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3.png", width = 600, height = 700)
# 	print(P5) 
# dev.off()

# png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3_nsnps.png", width = 600, height = 700)
# 	print(P4) 
# dev.off()

# png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3_pval.png", width = 600, height = 700)
# 	print(P2) 
# dev.off()

# png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3_Q_pval.png", width = 600, height = 700)
# 	print(P3) 
# dev.off()

# message(unique(Plot_dat$exposure2))
# # png("~/fatty-acids/mr/results/plots/ggforest_crc_secondary_pufas_univariableMR_v3.png", width = 1000, height = 600)
# # 	print(P1) 
# # dev.off()

# format(as.numeric(Plot_dat$pval), scientific = TRUE)

# library(
# 	)
# library(dplyr)
# sci <- data.frame(sci_note = 2.8 * 10^(1:10))
# sci <- sci %>% mutate(translation = format(sci_note, scientific = FALSE, big.mark = ","))
# kable(sci, col.names = c("Scientific Notation", "Full digit equivalent"))

# numb <- c(0.05, 0.05671, 0.000000027)



# P1<-forest(x=Plot_dat$b,sei=Plot_dat$se,
#     xlim=c(-16, 5),  atransf=exp,
#     slab=Plot_dat$exposure,
#     ilab=do.call(cbind,list(Plot_dat$nsnp,
#     Plot_dat$pval,Plot_dat$Q_pval)),ilab.xpos=c(-7.7,-6.5,-4.2),ilab.pos=4,ylim=c(2, 33),col=Plot_dat$col,rows=c(30:29,27:26,24:23,21:20,18:17,15:14,12:11,9:8,6:5,3:2),header=c("","OR (95% CI) per SD increase\nin genetically proxied PUFA"),top=2,xlab="")

# # cex=0.75, 
# rm(exp)
# P1
#     # ,order=Plot_dat$Order)
# ?forest.default