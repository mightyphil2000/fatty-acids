load("~/fatty-acids/mr/data/instruments.Rdata")
Dat<-unique(Dat[,c("snp","chr","position")])

bed_38<-read.table("/Users/ph14916/fatty-acids/colocalisation/data/UKBB_10K_bed_hg38.txt.gz")

Dat2<-merge(Dat,bed_38,by.x="snp",by.y="V4",all.x=TRUE)

names(Dat2)<-c("rsid","chr","bp_grch37","V1","bp_grch38","V3")
Dat2<-Dat2[order(as.numeric(Dat2$chr),Dat2$bp_grch37),c("rsid","chr","bp_grch37","bp_grch38")]

write.table(Dat2,"~/fatty-acids/mr/data/instruments_secondary_pufas_gecco.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)




readLines("~/fatty-acids/mr/data/instruments_secondary_pufas_gecco.txt")