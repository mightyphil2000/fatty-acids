# source("~/fatty-acids/mr/scripts/mr_functions.R")
setwd("~/fatty-acids/mr/results")
library(ggforestplot)
library(ggplot2)

Dat<-read.table("cismvMR_all_results_eas_ukb.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)
Dat<-Dat[Dat$panel == "UKB",]
Dat1<-rbind(Dat,Dat)

Dat1$b[1:22]<-Dat1$beta1_prune_FIVW[1:22]
Dat1$b[23:44]<-Dat1$beta2_prune_FIVW[23:44]

Dat1$se[1:22]<-Dat1$se1_prune_FIVW[1:22]
Dat1$se[23:44]<-Dat1$se2_prune_FIVW[23:44]

Dat1$exposure[1:22]<-Dat1$exposure1[1:22]
Dat1$exposure[23:44]<-Dat1$exposure2[23:44]
Dat1$tissue<-Dat1$exposure
Dat1$tissue<-gsub("FADS1 expression in ","",Dat1$tissue)
Dat1$tissue<-gsub("FADS2 expression in ","",Dat1$tissue)
Dat1$exposure[grep("FADS1",Dat1$exposure)]<-"FADS1 expression"
Dat1$exposure[grep("FADS2",Dat1$exposure)]<-"FADS2 expression"
Dat1$exposure<-gsub(" / D5D","",Dat1$exposure)
Dat1$exposure<-gsub(" / D6D","",Dat1$exposure)

Dat1$study<-c("GTEx")
Dat1$study[Dat1$exposure %in% c("AA:DGLA","GLA:LA")]<-"CHARGE"
Dat1$study[Dat1$tissue == "blood"]<-"eQTLGen"

Dat1$tissue[Dat1$study == "CHARGE"]<-"blood"
Dat1$tissue[Dat1$tissue=="whole blood"]<-"blood"
Dat2<-Dat1[Dat1$exposure %in% c("AA:DGLA","GLA:LA"),]
Dat2<-Dat2[order(Dat2$exposure),]
Dat3<-Dat1[!Dat1$exposure %in% c("AA:DGLA","GLA:LA"),]
Dat1<-rbind(Dat2,Dat3)
crc<-Dat1[Dat1$outcome == "Colorectal cancer",]
lc<-Dat1[Dat1$outcome == "Lung cancer",]
msc<-Dat1[Dat1$outcome == "Malignant skin cancer",]




# +geom_point(shape=Plot_dat$shape,size=Plot_dat$weight,fill=c("black"))
     
     # Forestplot
P1<-forestplot(
       df = crc,
       name=tissue,
       estimate = b,
       logodds = TRUE,
       colour = exposure,
       shape=study,
       xlab = "OR per SD increase in fatty acid ratio or gene expression"
     )

P2<-forestplot(
       df = lc,
       name=tissue,
       estimate = b,
       logodds = TRUE,
       colour = exposure,
       shape=study,
       xlab = "OR per SD increase in fatty acid ratio or gene expression"
     )

P3<-forestplot(
       df = msc,
       name=tissue,
       estimate = b,
       logodds = TRUE,
       colour = exposure,
       shape=study,
       xlab = "OR per SD increase in fatty acid ratio or gene expression"
     )     
     

png("~/fatty-acids/mr/results/plots/ggforest_cismvMR_crc.png", width = 600, height = 480)
	print(P1) 
dev.off()

png("~/fatty-acids/mr/results/plots/ggforest_cismvMR_lc.png", width = 600, height = 480)
	print(P2) 
dev.off()

png("~/fatty-acids/mr/results/plots/ggforest_cismvMR_msc.png", width = 600, height = 480)
	print(P3) 
dev.off()
