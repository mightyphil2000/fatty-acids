ao<-ieugwasr::gwasinfo()

IDS<-ao$id[which(ao$consortium == "GSCAN")]
b<-ieugwasr::associations(id=IDS, variants=c("rs174546","rs2524299") ,proxies=0) 

# flip results to reflect major allele
b$or<-round(exp(b$beta*-1),3)
b$lci<-round(exp(b$beta*-1-1.96*b$se),3)
b$uci<-round(exp(b$beta*-1+1.96*b$se),3)

b$beta2<-round(b$beta*-1,3)
b$lci2<-round(b$beta*-1-1.96*b$se,3)
b$uci2<-round(b$beta*-1+1.96*b$se,3)

b<-data.frame(b)

data.frame(b[b$rsid == "rs2524299",])
data.frame(b[b$rsid == "rs174546",])

ao1<-ao[ao$id %in% IDS ,]

write.table(b,"~/fatty-acids/mr/results/snplookups_rs174546_rs2524299_gscan.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)