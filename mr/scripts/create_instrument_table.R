library(ggforestplot)
library(ggplot2)
library(TwoSampleMR)
# mv_extract_exposures_local
# mv_extract_exposures
source("~/fatty-acids/mr/scripts/mr_functions.R")
# exposure
source("~/fatty-acids/mr/scripts/mr_functions.R")
load("~/fatty-acids/mr/data/instruments.Rdata")
exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")

length(unique(exp_dat$SNP))
length(unique(exp_dat$exposure))

excl<-c("Stearidonic acid (18:4n3)","Eicosadienoic acid (20:2n6)","Adrenic acid (22:4n6)","Docosapentaenoic acid (22:5n6)","Tetradecadienoic acid (14:2n9)","Dihomo-linolenic acid (20:3n3 or n6)","Ratio of omega-6 fatty acids to omega-3 fatty acids","Other polyunsaturated fatty acids than 18:2")

exp_dat<-exp_dat[!exp_dat$exposure %in% excl,]
length(unique(exp_dat$exposure))
exp_dat<-exp_dat[,c("exposure","SNP","beta.exposure","se.exposure","eaf.exposure","consortium","FADS","chr","position","samplesize.exposure","median_n","pval.exposure","id","pmid")]
exp_dat1<-exp_dat[exp_dat$FADS,]
length(unique(exp_dat$SNP))
exp_dat1[,c("exposure","SNP","beta.exposure","se.exposure","eaf.exposure","consortium","FADS","chr","position")]
write.table(exp_dat,"~/fatty-acids/mr/results/instruments_table.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)






head(exp_dat)


Temp<-exp_dat[,c("SNP","chr","FADS")]
Temp1<-Temp[Temp$FADS,]
Temp2<-Temp[!Temp$FADS,]
length(unique(Temp2$SNP))

Temp<-Temp[!duplicated(Temp$SNP),]
Temp<-Temp[!duplicated(Temp$chr),]
dim(Temp)
