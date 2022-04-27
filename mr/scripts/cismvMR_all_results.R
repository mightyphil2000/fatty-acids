setwd("~/fatty-acids/mr/results")
Files<-dir()
Files<-Files[grep("cismvMR",Files)]
Files<-Files[grep("all_results",Files,invert=TRUE)]

# x<-Files[15]
Dat<-lapply(Files,FUN=function(x) 
	read.table(x,sep="\t",head=TRUE,stringsAsFactors=FALSE))

Dat<-do.call(rbind,Dat)
Dat1<-Dat[Dat$panel %in% c("EAS","UKB"),]
Dat2<-Dat[!Dat$panel %in% c("EAS","UKB"),]

write.table(Dat1,"cismvMR_all_results_eas_ukb.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(Dat2,"cismvMR_all_results_chb_jpt.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)