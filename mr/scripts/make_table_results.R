source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
# load("~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")
load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")

length(unique(mr_res1$cancer))
mr_res<-format_dat2(dat=mr_res1,sort="system")
dim(mr_res)
cancer<-unlist(strsplit(mr_res$cancer,split="\n"))
mr_res$cancer<-trimws(cancer[seq(1,length(cancer),by=2)])

dim(meta.tab9)

excl<-c("cell","id.exposure","method","nsnp","have_gwas","total","Note","overlap","include","MAC_rs603424","MAC_rs3734398","MAC_rs174546","MAC100rs174546")
mr_res<-mr_res[,!names(mr_res) %in% excl]
mr_res$OR<-exp(mr_res$b)
mr_res$LCI<-exp(mr_res$b - 1.96*mr_res$se)
mr_res$UCI<-exp(mr_res$b + 1.96*mr_res$se)
length(which(!duplicated(mr_res$cancer)))
mr_res$id[grep("86",mr_res$id)]
write.table(mr_res,"~/fatty-acids/mr/results/mr_results_table.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# library(devtools)
# library(TwoSampleMR)
# library(plyr)
# library(ggforestplot)
# library(ggplot2)


# library(tidyverse)


