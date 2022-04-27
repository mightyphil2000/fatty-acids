# install.packages("BiocManager")
# BiocManager::install("biomaRt")
# install.packages("arrow")
# library(biomaRt)
ukb<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",head=F,stringsAsFactors=F)
ukb$SNP<-paste(ukb$V1,":",ukb$V4,sep="")

Chrs<-c(11,6,6,10,2,16,20)
# Files<-c(
# 	"CHS_EA_HRC_LD_11_61043499_62159523.ld", 
# 	"CHS_EA_HRC_LD_6_10480992_11544547.ld",
# 	"CHS_EA_HRC_LD_6_52632196_53713947.ld",
# 	"CHS_EA_HRC_LD_10_101606881_102624591.ld"
# 	"CHS_EA_HRC_LD_2_27219709_28246554.ld",
# 	"CHS_EA_HRC_LD_16_14568448_15733196.ld"
# 	"CHS_EA_HRC_LD_20_12489627_13647411.ld"
Starts<-c(61043499,10480992,52632196,101606881,27219709,14568448,12489627)
Ends<-c(62159523,11544547,53713947,102624591,28246554,15733196,13647411)
Regions<-c("FADS","ELOVL2","ELOVL5","SCD","GCKR","PDXDC1","SPTLC3")

# Chr<-Chrs[1]
# Start<-Starts[1]
# End<-Ends[1]
# Region<-Regions[1]

# Res<-lapply(1:length(Chrs),FUN=function(x)
	# format_dat(Chr=Chrs[x],Start=Starts[x],End=Ends[x],Region=Regions[x]))

list.res<-lapply(1:length(Chrs),FUN=function(x)
	format_dat(Chr=Chrs[x],Start=Starts[x],End=Ends[x],Region=Regions[x])
	)

# head(list.res[[2]])

# df1<-do.call(rbind,list.rs)

# for(i in 1:length(Chrs)){
# 	print(i)
# 	chs<-read.table(paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/CHS_ld_reference_panel/CHS_EA_HRC_LD_",Chr,"_",Start,"_",End,".frq",sep=""),head=T,stringsAsFactors=F)
# 	chs.ld<-read.table(paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/CHS_ld_reference_panel/CHS_EA_HRC_LD_",Chr,"_",Start,"_",End,".ld",sep=""),head=F,stringsAsFactors=F)
# 	Pos<-which(chs$MAF != 0)
# 	chs2<-chs[Pos,]
# 	chs.ld2<-chs.ld[Pos,Pos]
# 	return(chs2)
# }
format_dat<-function(Chr=NULL,Start=NULL,End=NULL,Region=NULL){
	print(Chr)
	print(which(Chrs==Chr))
	chs<-read.table(paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/CHS_ld_reference_panel/CHS_EA_HRC_LD_",Chr,"_",Start,"_",End,".frq",sep=""),head=T,stringsAsFactors=F)
	
	chs.ld<-read.table(paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/CHS_ld_reference_panel/CHS_EA_HRC_LD_",Chr,"_",Start,"_",End,".ld",sep=""),head=F,stringsAsFactors=F)

	Pos<-which(chs$MAF != 0)
	chs2<-chs[Pos,]
	chs.ld2<-chs.ld[Pos,Pos]
	ukb2<-ukb[!duplicated(ukb$SNP),]
	ukb2<-ukb2[ukb2$V1 == Chr,]
	ukb3<-ukb2[ukb2$V4>=Start & ukb2$V4<=End,]
	chs3<-merge(chs2,ukb3[,c("SNP","V2","V4")],by="SNP",all.x=T)
	names(chs3)[names(chs3)=="V2"]<-"rsid"

	summary(chs3$MAF[which(is.na(chs3$rsid))]) # 75% of missing SNPs have MAF<=0.00083
	summary(chs3$MAF[which(!is.na(chs3$rsid))]) #75% of non missing SNPs have MAF>=0.007

	if(!all(chs2$SNP == chs3$SNP)) stop("mismatches in rsid")
	if(any(chs2$SNP != chs3$SNP)) stop("mismatches in rsid")
	if(!all( chs3$SNP== chs2$SNP)) stop("mismatches in rsid")
	if(any(chs3$SNP!= chs2$SNP )) stop("mismatches in rsid")

	names(chs3)[names(chs3)=="SNP"]<-"chr_bp"
	names(chs3)[names(chs3)=="rsid"]<-"SNP"
	names(chs3)[names(chs3)=="V4"]<-"bp_hg19"

	Pos<-which(is.na(chs3$bp_hg19))
	bp_hg19<-chs3$chr_bp[Pos]
	bp_hg19<-unlist(strsplit(bp_hg19,split=":"))
	bp_hg19<-bp_hg19[seq(2,length(bp_hg19),by=2)]
	chs3$bp_hg19[Pos]<-bp_hg19
	chs3$region <- Region
	return(list(chs3,chs.ld2,Region))
}



# getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="external_gene_name",values=gene_tab$gene,mart=Mart) #unique(Fa.tab$gene_id)


# library("biomaRt")
# Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")
# Attr<-listAttributes(Mart)
# ensembl<-getBM(attributes=c("refsnp_id",filters="snp_filter",values=b1$SNP,mart=Mart)