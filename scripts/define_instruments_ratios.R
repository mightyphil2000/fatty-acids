# install.packages("BiocManager")
# install.packages("bit")
# install.packages("bit64")
# install.packages("RSQLite")
# BiocManager::install("biomaRt")
# module add languages/R-3.5.1-ATLAS-gcc-6.1
library(biomaRt)
# install.packages("BiocFileCache")
# library(devtools)
# install_github("MRCIEU/TwoSampleMR")
# install.packages("curl")
# install.packages("plyr")
library(TwoSampleMR)


# Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
# Mart = useMart(host="grch37.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
# Mart <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_SNP",path="/biomart/martservice",dataset="hsapiens_snp")
# grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

# Mart<-useEnsembl(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="grch37.ensembl.org", version = NULL, GRCh = 37, mirror = "uswest.ensembl.org", verbose = FALSE)

Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
Attr<-listAttributes(Mart)

setwd()

sys.cmd<-dir("/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge_asia/ratios")
sys.cmd<-paste("cp /newprojects/mrcieu/research/data/MR/fatty_acids/data/charge_asia/ratios/",sys.cmd," /projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge_asia/",sep="")
system(sys.cmd[2])

sys.cmd<-dir("/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/ratios")
sys.cmd<-sys.cmd[grep("pooled",sys.cmd)]
sys.cmd<-paste("cp /newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/ratios/",sys.cmd," /projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/dorajoo/",sep="")
system(sys.cmd[3])
system(sys.cmd[4])





# cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/
# source("~/fatty-acids-mr/instruments/Extract_SNPs_function.R")
# ref <- read.table("~/fatty-acids-mr/instruments/data_maf0.01_rs.bim.gz")
# 
# ref <- read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/data_maf0.01_rs.bim")
# res<-read.table("AA_to_DGLA.tab" ,sep="\t",head=T,stringsAsFactors=F)
# res1<-res[res$pval<5e-8,]

Files<-paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/",dir("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge"),sep="")
Files1<-paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge_asia/",dir("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge_asia/"),sep="")
Files2<-paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/dorajoo/",dir("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/dorajoo/"),sep="")
Res<-lapply(1:length(Files),FUN=function(i) extract_data2(File=Files[i],file_sep="\t",extract_pvalue=c(5e-8,7)))
Res1<-lapply(1:length(Files1),FUN=function(i) extract_data2(File=Files1[i],file_sep="\t",extract_pvalue=c(5e-8,7)))
Res2<-lapply(1:length(Files2),FUN=function(i) extract_data2(File=Files2[i],file_sep="\t",extract_pvalue=c(5e-8,7)))


names(Res)<-gsub(".tab","",Files)
names(Res1)<-gsub(".tab","",Files1)
names(Res2)<-gsub(".tab","",Files2)

Dat<-do.call(rbind,Res)
Dat1<-do.call(rbind,Res1)
Dat2<-do.call(rbind,Res2)
Dat$study<-"CHARGE"
Dat$population<-"Europeans"
Dat1$study<-"Zhu"
Dat1$population<-"East Asians"

Dat2$study<-"Dorajoo"
Dat2$population<-"East Asians"
Dat3<-do.call(rbind,list(Dat,Dat1,Dat2))

String<-unlist(strsplit(row.names(Dat3),split="\\."))
Dat3$exposure<-String[seq(1,length(String),by=2)]
String<-gregexpr("/",Dat3$exposure)
Pos<-unlist(lapply(1:length(String),FUN=function(x) 
	unlist(String[x])[length(unlist(String[x]))]))
Dat3$exposure<-substring(Dat3$exposure,Pos+1,nchar(Dat3$exposure))
row.names(Dat3)<-NULL
Dat<-Dat3[Dat3$study=="CHARGE",]

ref<-extract_data2(snplist=unique(Dat$SNP),file_sep="\t",exact_match=TRUE,File="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/data_maf0.01_rs.bim",Head=FALSE)

Dat<-merge(Dat,ref,by.x="SNP",by.y="V2",all.x=T)
names(Dat)[names(Dat)=="V4"]<-"pos"
names(Dat)[names(Dat)=="V1"]<-"chr"

# snps<-unique(Res2$SNP[is.na(Res2$pos)])
# snp.pos<-getBM(attributes=c("refsnp_id","chr_name","chrom_start"),filters="snp_filter",values=unique(Res1$SNP),mart=Mart)

# Res$SNP[!Res$SNP %in% ref$V2]
Dat$chr[Dat$SNP == "rs4985167"]<-16
Dat$pos[Dat$SNP == "rs4985167"]<-15082865
Dat$chr[Dat$SNP == "rs7192552"]<-16
Dat$pos[Dat$SNP == "rs7192552"]<-15122195
Dat$z<-Dat$beta/Dat$se
Dat$b_sd<-b_sd(Dat$z,maf=Dat$eaf,n=Dat$samplesize)
Dat$se_sd<-Dat$b_sd/Dat$z
Dat$r2<-2*Dat$b_sd^2*Dat$eaf*(1-Dat$eaf)
Dat$r2inv<-1/(Dat$r2 * 1000)

# Res1$exposure<-"AA_to_DGLA"
FA1<-format_data(Dat,type="exposure",phenotype_col="exposure",snp_col="SNP", samplesize_col = "samplesize",chr_col="chr",pos_col="pos",pval_col="r2inv",beta_col = "b_sd",se_col = "se_sd", eaf_col = "eaf")
FA2<-clump_data(FA1,clump_r2=0.001,clump_kb=10000)

save(FA2,file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/clumped_snps_gwis_charge.Rdata")

# write.table(FA2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/clumped_snps_charge.txt",sep="\t",col.names=T,row.names=F,quote=F)

# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/clumped_snps_charge.Rdata .

# FA3<-clump_data(FA1,clump_r2=0.0001,clump_kb=10000)
Res2[which(Res2$SNP=="rs968567"),]
Temp<-Res2[Res2$exposure=="GLA_to_LA",c("SNP","r2inv","pval","z","b_sd","eaf")]
Temp<-Temp[order(Temp$r2inv),]
which(Temp$SNP=="rs968567")
Temp$N<-1:nrow(Temp)

unique(FA2[order(FA2$chr.exposure),c("SNP","exposure","chr.exposure")])

b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}