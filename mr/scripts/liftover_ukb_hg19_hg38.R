options("scipen" = 10)
Res<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",head=F,stringsAsFactors=F)
 # chrN:start-end formats):
 
# Bed<-paste("chr",Res$V1,":",Res$V4,"-",Res$V4,sep="")
Res1<-Res[,c("V1","V4","V4","V2")]
Res1$V1<-paste("chr",Res1$V1,sep="")
Res1$V4.1 <-Res1$V4.1+1
write.table(Res1,"/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt",col.names=F,row.names=F,quote=F)

sys.cmd<-paste("~/liftover/liftOver /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19.txt ~/liftover/hg19ToHg38.over.chain /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg38.txt /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg19_unMapped.txt")

system(sys.cmd)


