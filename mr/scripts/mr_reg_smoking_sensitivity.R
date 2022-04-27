library(plyr)
library(ggforestplot)
library(ggplot2)
library(meta)
library(metafor)

source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")


mr_res<-format_metareg_v2()

M<-create_correlation_matrix(study=mr_res)
M<-M$matrix
rownames(M)<-mr_res$study.abbreviation
colnames(M)<-mr_res$study.abbreviation
# corr_results_list2<-res$corr_results_list2
mr_res$se_decoupled<-decoupling(s=mr_res$se,C=M)
Weights<-1/mr_res$se^2
Weights_decoupled<-1/mr_res$se_decoupled^2

# test<-mr_res[c(1,9),]

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$infl,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$infl,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

table(mr_res$system_num2)

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$system2_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$system2_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$external_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$external_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))




# Pos<-grep("InterLymph",mr_res$study.abbreviation)
# mr_res[Pos,c("study.abbreviation","se","se_decoupled")]
# M[Pos,Pos]

# plot(mr_res$se_decoupled,mr_res$se)

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

write.table(M,"~/fatty-acids/mr/results/study_correlation_matrix.txt",sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)


M<-create_correlation_matrix(study=mr_res)
M

