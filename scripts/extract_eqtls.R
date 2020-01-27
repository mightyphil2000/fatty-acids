# source("~/fatty-acids/scripts/extract_snps_functions.R")
##########
# eQTLGen#
##########

cd /projects/MRC-IEU/scratch/
# cp /projects/MRC-IEU/scratch/eQTL-gen-results.zip /projects/MRC-IEU/users/ph14916/eQTLGen 
# cd /projects/MRC-IEU/users/ph14916/eQTLGen 

#updated file from Jie Zheng combined-eQTLs_full.EAF.beta.se.chr.pos.txt 

cp  /projects/MRC-IEU/scratch/combined-eQTLs_full.EAF.beta.se.chr.pos.txt /projects/MRC-IEU/users/ph14916/eQTLGen  
cd /projects/MRC-IEU/users/ph14916/eQTLGen  
# gunzip ENSG00000149485.eQTL.tab.gz
# gunzip ENSG00000134824.eQTL.tab.gz
# head ENSG00000149485.eQTL.tab
# head ENSG00000134824.eQTL.tab
# rs174546

# grep rs968567 ENSG00000134824.eQTL.tab

# grep -w ENSG00000149485 combined-eQTLs_full.EAF.beta.se.chr.pos.txt | grep -w rs968567 
# grep -w rs968567 ENSG00000149485.eQTL.tab

# head -1 combined-eQTLs_full.EAF.beta.se.chr.pos.txt 
# 	EAF 0.148312947324573	beta 0.560779864142844	se 0.0107584286171977

# head -1 ENSG00000149485.eQTL.tab
# 	EAF 0.148312947324573	beta 0.799339649931238 se	0.0153351414958981

# grep -w ENSG00000134824 combined-eQTLs_full.EAF.beta.se.chr.pos.txt | grep -w rs968567 
# 	EAF 0.148312947324573	beta 0.999342359176102 se	0.0097290180677342	

# grep -w rs968567 ENSG00000134824.eQTL.tab
# 	EAF 0.148312947324573	beta 1.30122387159944 se 0.0126679615256129

#########
#ENSEMBL#
#chr_name & chrom_start are hg38/CRCh38 format from ENSEMBL
#minor allele frequency in EMSEMBL is calculated using all 1000 genomes samples (global MAF)
# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")

# head(Dat.ldc)
# unique(Dat.ldc[,c("SNP","Chr","SNP.GRCh38.p12") ])

# library("biomaRt")
# Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
# Attr<-listAttributes(Mart)
# gene_names<-c("FADS1","FADS2","ELOVL2")
# Test<-lapply(1:length(gene_names),FUN=function(i) getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","wikigene_name"),filters="external_gene_name",values=gene_names[i],mart=Mart))
# Ens<-do.call(rbind,Test)
# Ens$ensembl_gene_id
setwd("/projects/MRC-IEU/users/ph14916/eQTLGen")
Ens<-c("ENSG00000197977", "ENSG00000230314",#ELOVL2
	"ENSG00000149485" ,#FADS1 (FADS region)
	"ENSG00000134824", #FADS2 (FADS region)
	"ENSG00000221968", #FADS3 (FADS region)
	"ENSG00000168496", #FEN1 (FADS region)
	"ENSG00000134825",#TMEM258 (FADS region)
	"ENSG00000124920", #MYRF (FADS region)
	"ENSG00000134780", #DAGLA (FADS region)
	"ENSG00000167994", #RAB3IL1  (FADS region)
	"ENSG00000167995", #BEST1  (FADS region)
	"ENSG00000167996", #FTH1 (FADS region)
	"ENSG00000099194") #SCD
Pos<-unlist(lapply(1:length(Ens),FUN=function(i) grep(Ens[i],dir())))
Files<-dir()[Pos]
lapply(1:length(Ens),FUN=function(i) grep(Ens[i],Files)) #ELOVL2 missing (ELOVL2)
	
# # dir()[grep("ENSG00000168496",dir())]
# # Res<-read.table("ENSG00000149485.eQTL.tab",sep="\t",head=T,stringsAsFactors=F)
# eQTL_gen<-extract_data(file_dir = "/projects/MRC-IEU/users/ph14916/eQTLGen",snplist = "/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt",type="not_fatty_acids",wk_dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/eQTLGen",file_list=Files)
# dim(eQTL_gen)
eqtl_gen<-extract_data3(gene="ENSG00000134824",snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt",File="/projects/MRC-IEU/users/ph14916/eQTLGen/combined-eQTLs_full.EAF.beta.se.chr.pos.txt")

# grep -w ENSG00000134824 /projects/MRC-IEU/users/ph14916/eQTLGen/combined-eQTLs_full.EAF.beta.se.chr.pos.txt | grep -wf /projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_temp2.txt > /projects/MRC-IEU/users/ph14916/fatty_acids_summary/output_temp.txt

# eqtl_gen<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0("/projects/MRC-IEU/users/ph14916/eQTLGen/",Files[i])
	eqtl_gen[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt",File=target_file,exact_match=TRUE,file_sep="\t") #snplist_coloc.txt contains 17381 SNPs +/- 500kb of each of 9 index SNPs selected as instruments for fatty acid traits. The SNPs reside at 6 genomic regions: FADS, ELOVL2, SCD, GCKR, PDXDC1, SPTLC3 
}
eqtl_gen1<-do.call(rbind,eqtl_gen)


# ELOVL2 gene not present in eQTLGen
save("eqtl_gen1",file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/eqtlgen_snplist_coloc.Rdata")

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/eqtlgen_snplist_coloc.Rdata .
# grep("ENSG00000197977",dir())
######
#GTEx#
######
# /projects/MRC-IEU/research/data/broad/public/gtex/released/2018-10-05/data/GTEx_Analysis_v7_eQTL_all_associations.

cp GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz /projects/MRC-IEU/users/ph14916/gtex

cd /projects/MRC-IEU/users/ph14916/gtex
# cat "/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_bedformat.txt"
setwd("/projects/MRC-IEU/users/ph14916/gtex")
Files<-dir()

gtex<-NULL 
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	target_file<-paste0("/projects/MRC-IEU/users/ph14916/gtex/",Files[i])
	gtex[[i]]<-extract_data2(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_bedformat.txt",File=target_file,exact_match=FALSE,file_sep="\t") #snplist_coloc.txt contains 17381 SNPs +/- 500kb of each of 9 index SNPs selected as instruments for fatty acid traits. The SNPs reside at 6 genomic regions: FADS, ELOVL2, SCD, GCKR, PDXDC1, SPTLC3 
}
gtex1<-do.call(rbind,gtex)

save("gtex1",file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc.Rdata")

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc.Rdata .

grep -w ENSG00000149485.12 Whole_Blood.allpairs.txt  | grep -f /projects/MRC-IEU/users/ph14916/temp.txt 
 
cat /projects/MRC-IEU/users/ph14916/temp.txt 
 
61595564
61600342

# gtex files too large for github. I though maybe save the separate tissues in separate files would help but even separated by tissue the files are too large 
# load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc.Rdata")
# Temp1<-split(gtex1,f=gtex1$file)
# names(Temp1)<-gsub(".txt","",unique(gtex1$file))
# # paste("trait",1:length(Temp1),sep="")
# list2env(Temp1,envir=.GlobalEnv)

# # save(list=names(Temp1),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc_v2.Rdata")

# # load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc_v2.Rdata")

# rm(list=c("Temp1","gtex1"))
# Objects<-ls()
# for(i in 1:length(Objects)){
# 	print(Objects[i])
# 	File.save<-paste("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc_",Objects[i],".Rdata",sep="")
# 	Obj<-Objects[i]
# 	save("Adipose_Subcutaneous.allpairs",file=File.save)
# }


# scp -T username@ip.of.server.copyfrom:"file1.txt file2.txt" "~/yourpathtocopy"




# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc_Whole_Blood.allpairs.Rdata .


# gtex_snplist_coloc_Adipose_Subcutaneous.allpairs.Rdata      
# gtex_snplist_coloc_Heart_Atrial_Appendage.allpairs.Rdata
# gtex_snplist_coloc_Adipose_Visceral_Omentum.allpairs.Rdata  
# gtex_snplist_coloc_Heart_Left_Ventricle.allpairs.Rdata
# gtex_snplist_coloc_Artery_Aorta.allpairs.Rdata              
# gtex_snplist_coloc_Liver.allpairs.Rdata
# gtex_snplist_coloc_Artery_Coronary.allpairs.Rdata           
# gtex_snplist_coloc_Lung.allpairs.Rdata
# gtex_snplist_coloc_Artery_Tibial.allpairs.Rdata             
# gtex_snplist_coloc_Colon_Sigmoid.allpairs.Rdata             
# gtex_snplist_coloc_Small_Intestine_Terminal_Ileum.allpairs.Rdata
# gtex_snplist_coloc_Colon_Transverse.allpairs.Rdata          
# gtex_snplist_coloc_Whole_Blood.allpairs.Rdata

# load("gtex_snplist_coloc_Whole_Blood.allpairs.Rdata")
# for i in hello 1 * 2 goodbye 
# do
#   echo "Looping ... i is set to $i"
# done

# while read gtexfiles.txt
# do
# 	cat $f
# done

# for f in gtexfiles*.txt;  do echo ${f}; done;

# for i in $(cat $1); do
#     echo "tester: $i"
# done

# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snplist_coloc_Adipose_Subcutaneous.allpairs.Rdata  .

# # Files<-Files[1]
# # Gtex<-extract_data(file_dir = "/projects/MRC-IEU/users/ph14916/gtex",snplist = "/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc_bedformat.txt",type="not_fatty_acids",wk_dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex",file_list=Files,exact_match=FALSE)

# # save("Gtex",file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gtex_snps.Rdata")