# region coordinates can be found here fatty-acids/mr/data/candidate_regions_coordinates.xlsx
cd /projects/MRC-IEU/users/ph14916
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/

# plink_dir/plink --bfile 1000genomes/data_maf0.01_rs --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt  --make-just-bim --out elolv5_bim
# plink_dir/plink --bfile 1000genomes/data_maf0.01_rs --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt --out plink_dir/snp_r_stats --r square

# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/plink_dir/snp_r_stats.ld .


######################################################################################
# Genetate ld correlation matrix using UK Biobank random 10K from Tom Richardson
#######################################################################################
# uk biobank positions in bim file are buld 37

# /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bed  UKBB_10K.bim  UKBB_10K.fam

################################################################
#FADS region on chromosome 11 defined by rs174528 index SNP
################################################################


/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 61043499 \
	--to-bp 62159523 \
	--chr 11 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/fads_r_matrix_ukb \
	--r square 
	 
/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 61043499 \
	--to-bp 62159523 \
	--chr 11 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/fads_r_matrix_ukb \
	--freq 

###############################################################
# ELOVL2 region on chromosome 6 defined by rs3734398 index SNP#
###############################################################

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 10480992	 \
	--to-bp 11544547 \
	--chr 6 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/elovl2_r_matrix_ukb \
	--r square 
	 
/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 10480992	 \
	--to-bp 11544547 \
	--chr 6 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/elovl2_r_matrix_ukb \
	--freq 

#################################################################
# SCD region on chromosome 10 defined by rs603424 index SNP#
#################################################################

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 101606881 \
	--to-bp 102624591 \
	--chr 10 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/scd_r_matrix_ukb \
	--r square 

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 101606881 \
	--to-bp 102624591 \
	--chr 10 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/scd_r_matrix_ukb \
	--freq

###############################################################
# GCKR region on chromosome 2 defined by rs780093 index SNP#########
########################################################################

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 27219709 \
	--to-bp 28246554 \
	--chr 2	 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gckr_r_matrix_ukb \
	--r square 

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 27219709 \
	--to-bp 28246554 \
	--chr 2	 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/gckr_r_matrix_ukb \
	--freq


################################################################
# PDXDC1 region on chromosome 16 defined by rs4985155 index SNP#
################################################################

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 14568448 \
	--to-bp 15733196 \
	--chr 16 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/pdxdc1_r_matrix_ukb \
	--r square 

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 14568448 \
	--to-bp 15733196 \
	--chr 16 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/pdxdc1_r_matrix_ukb \
	--freq

##################################################################
# SPTLC3 region on chromosome 20 defined by rs680379 index SNP######
########################################################################

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 12489627 \
	--to-bp 13647411 \
	--chr 20 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/sptlc3_r_matrix_ukb \
	--r square 

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 12489627 \
	--to-bp 13647411 \
	--chr 20 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/sptlc3_r_matrix_ukb \
	--freq

 
###################################################################### 
# ELOVL5 region on chromosome 6 defined by rs12210577 index SNP#####
######################################################################

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 52632196 \
	--to-bp 53713947 \
	--chr 6 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/elovl5_r_matrix_ukb \
	--r square 

/projects/MRC-IEU/users/ph14916/plink_dir/plink \
	--bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  \
	--from-bp 52632196 \
	--to-bp 53713947 \
	--chr 6 \
	--out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/elovl5_r_matrix_ukb \
	--freq


# wc data_maf0.01_rs.bim
# grep rs174528 /projects/MRC-IEU/scratch/for_Phil/UKBB_10K.bim