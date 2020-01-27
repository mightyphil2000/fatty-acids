# source("~/fatty-acids/scripts/extract_snps_functions.R")

##########
# TRICL##
##########


cd /projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc

Files<-c("Onco_TRICL_032116_Adeno.csv.tab","Onco_TRICL_032116_Ever.csv.tab","Onco_TRICL_032116_Never.csv.tab","Onco_TRICL_032116_Overall.csv.tab","Onco_TRICL_032116_Small.csv.tab","Onco_TRICL_032116_Squam.csv.tab")

Tricl<-extract_data(file_dir = "/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc",snplist = "/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt",type="not_fatty_acids",wk_dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/tricl",file_list=Files)

save("Tricl",file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/tricl_snps.Rdata")


# GECCO/CORECt

    
# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt ~/MR_FattyAcids/data/summary\ data/

# summary data provided for SNPs from all 6 genomic regions of interset. Should work directly with colocalisation pipeline
CRC<-read.csv("~/MR_FattyAcids/data/summary data/colorectal_cancer/061119/1255_MarginalResults_HRC125K_20191105.csv",head=T,stringsAsFactors=F)

save("CRC",file="~/MR_FattyAcids/data/summary data/colorectal_cancer/061119/crc_snps.Rdata")

snplist<-readLines("~/MR_FattyAcids/data/summary data/snplist_coloc.txt")

table(CRC$Chromosome)