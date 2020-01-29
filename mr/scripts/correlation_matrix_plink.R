cd /projects/MRC-IEU/users/ph14916


# plink_dir/plink --bfile 1000genomes/data_maf0.01_rs --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt  --make-just-bim --out elolv5_bim
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/
plink_dir/plink --bfile 1000genomes/data_maf0.01_rs --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/elolv5_snplist.txt --out plink_dir/snp_r_stats --r square

cd /Users/ph14916/fatty-acids-mr/coloc
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/plink_dir/snp_r_stats.ld .


######################################################################################
# Genetate ld correlation matrix using UK Biobank random 10K from Tom Richardson
#######################################################################################
# uk biobank positions in bim file are buld 37

# /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bed  UKBB_10K.bim  UKBB_10K.fam

#FADS region on chromosome 11 defined by rs174528 index SNP
/projects/MRC-IEU/users/ph14916/plink_dir/plink --bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs174528_ukb.txt  --out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs174528_r_matrix_ukb --r square

# ELOVL2 region on chromosome 6 defined by rs3734398 index SNP
/projects/MRC-IEU/users/ph14916/plink_dir/plink --bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs3734398_ukb.txt  --out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs3734398_r_matrix_ukb --r square

# SCD region on chromosome 10 defined by rs603424 index SNP
/projects/MRC-IEU/users/ph14916/plink_dir/plink --bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs603424_ukb.txt  --out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs603424_r_matrix_ukb --r square

# GCKR region on chromosome 2 defined by rs780093 index SNP
/projects/MRC-IEU/users/ph14916/plink_dir/plink --bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs780093_ukb.txt  --out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs780093_r_matrix_ukb --r square

# PDXDC1 region on chromosome 16 defined by rs4985155 index SNP
/projects/MRC-IEU/users/ph14916/plink_dir/plink --bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs4985155_ukb.txt  --out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs4985155_r_matrix_ukb --r square

# SPTLC3 region on chromosome 20 defined by rs680379 index SNP
/projects/MRC-IEU/users/ph14916/plink_dir/plink --bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs680379_ukb.txt  --out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs680379_r_matrix_ukb --r square
 
# ELOVL5 region on chromosome 6 defined by rs12210577 index SNP
/projects/MRC-IEU/users/ph14916/plink_dir/plink --bfile /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K  --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs12210577_ukb.txt  --out /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/rs12210577_r_matrix_ukb --r square






cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation/


cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary

wc data_maf0.01_rs.bim
grep rs174528 /projects/MRC-IEU/scratch/for_Phil/UKBB_10K.bim