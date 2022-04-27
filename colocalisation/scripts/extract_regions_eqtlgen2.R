# first run extract_regions_eqtlgen.sh
setwd("~/fatty-acids/colocalisation/data")
Files<-dir()
Files<-Files[grep("eqtlgen_",Files)]
df1<-lapply(1:length(Files), FUN=function(i) read.table(Files[i],sep="\t",head=T,stringsAsFactors=F))
df2<-do.call(rbind,df1)
save("df2",file="~/fatty-acids/colocalisation/data/eqtlgen.Rdata")

