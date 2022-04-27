# header for gtex. need to add chr and pos to end of file 
# gene_id	variant_id	tss_distance	ma_samples	ma_count	maf	pval_nominal	slope	slope_se

# library(arrow)
# df <- read_parquet("/projects/MRC-IEU/users/ph14916/gtex2/EUR/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_EUR_eQTL_all_associations_Whole_Blood.v8.EUR.allpairs.chr11.parquet")
# df1<-data.frame(df)
# String<-strsplit(df1$variant_id,split="_")
# String<-unlist(String)
# df1$Chr<-String[seq(1,length(String),by=5)]
# df1$bp<-String[seq(2,length(String),by=5)]
# df2<-df1[df1$bp>=61276027 & df1$bp<= 62392051,]
	
# ukb<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg38.txt",sep="\t",head=F,stringsAsFactors=F)
# chs<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/CHS_ld_reference_panel/CHS_EA_HRC_LD_11_61043499_62159523.frq",head=T,stringsAsFactors=F)
# df3<-merge(df2,ukb[,c("V1","V2","V4")],by.x=c("Chr","bp"),by.y=c("V1","V2"),all.x=T)

# names(df3)[names(df3) == "V4" ]<-"SNP"

# save(df3,file="/projects/MRC-IEU/users/ph14916/gtex2/EUR/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_EUR_eQTL_all_associations_Whole_Blood.v8.EUR.allpairs_chr11_fads.Rdata")

setwd("/projects/MRC-IEU/users/ph14916/gtex2/all_fatty_acid_regions")
Files<-dir()
Files<-Files[grep("GTEx",Files,invert=T)]
Files<-Files[grep("FADS",Files)]
Files<-Files[grep("\r",Files,invert=T)]

# dat<-df1
# File<-Files[1]
ref<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K_bed_hg38.txt",sep="\t",head=F,stringsAsFactors=F)
ref1<-ref[ref$V1 %in% "chr11",]
# ref19<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",sep="\t",head=F,stringsAsFactors=F)



df3<-lapply(1:length(Files),FUN=function(i)
	format_dat(File=Files[i],ref_dat=ref1)
	)

df5<-do.call(rbind,df3)

load("/projects/MRC-IEU/users/ph14916/gtex2/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata")

df6<-df4[df4$region == "FADS",]
df7<-df6[df6$tissue != "lung",]
df8<-df6[df6$tissue == "lung",]
df9<-rbind(df5,df7)


# df10<-df5[df5$tissue == "lung",]
# Dups<-paste(df10$SNP,df10$gene_id)[duplicated(paste(df10$SNP,df10$gene_id))]
# df10<-df10[!paste(df10$SNP,df10$gene_id) %in% Dups,]
# Temp<-merge(df8,df10,by=c("SNP","gene_id"))
# plot(Temp$slope.x,Temp$slope.y)
df4<-df9
save(df4,file="/projects/MRC-IEU/users/ph14916/gtex2/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_FADS_ukbsnps_v2.Rdata")


# save(df4,file="/projects/MRC-IEU/users/ph14916/gtex2/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata")

cd ~/fatty-acids/colocalisation/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/gtex2/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_FADS_ukbsnps_v2.Rdata .

# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/gtex2/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata .


format_dat<-function(File=NULL,ref_dat=NULL){
	print(File)
	print(which(Files == File))
	dat<-read.table(File,sep="\t",head=T,stringsAsFactors=F)
	print(dim(dat))
	if(!is.null(ref_dat)){
		ref1<-ref_dat
	}else{
		ref1<-ref[ref$V1 %in% unique(dat$chr),]
	}	
	df2<-dat[dat$bp_hg38 %in% ref1$V2 ,]
	tissue_gene<-unlist(strsplit(File,split="GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_"))
	tissue_gene<-tissue_gene[tissue_gene!=""]
	tissue_gene<-unlist(strsplit(tissue_gene,split=".allpairs.txt.gz_"))
	Tissue<-tissue_gene[1]
	Tissue<-gsub("_"," ",Tissue)
	Tissue<-tolower(Tissue)
	Region<-tissue_gene[2]
	Region<-gsub(".txt","",Region)
	df2$tissue<-Tissue
	df2$region<-Region
	df2<-merge(df2,ref1[,c("V1","V2","V4")],by.x=c("chr","bp_hg38"),by.y=c("V1","V2"))
	names(df2)[names(df2) == "V4"]<-"SNP"
	unique(df2[df2$SNP=="rs174546",c("chr","bp_hg38")])
	return(df2)
}


# duplicate positions seem to correspond to CNVs overlapping with SNPs in UKB bim file
# ref1[ref1$V4 == "rs386701190",]
# ref1[ref1$V4 == "rs77143770",]

# ref19[ref19$V2 == "rs386701190",]
# ref19[ref19$V2 == "rs77143770",]

# dbSNP
# rs386701190
# build 37 52902188 52902193 
# build 38 53037390 53037395 

# rs77143770
# build 37 52902193
# build 38 53037395 


# IDs<-paste(df3$chr,df3$bp_hg38,df3$tissue,df3$region,df3$gene_id)
# df4<-df3[!IDs %in% dups,]
# dups<-unique(IDs[duplicated(IDs)])
# any(duplicated(IDs))

