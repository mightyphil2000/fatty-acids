ssh -X ph14916@bc4login.acrc.bris.ac.uk 
ssh -Y ph14916@bluecrystalp3.acrc.bris.ac.uk -X

# "~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata"
# "~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata"
# "~/MR_FattyAcids/data/summary_data/colorectal_cancer_ACCC/ACCC_summary_data_GRCh37_6_regions.txt"
# /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC

#coordinates chr11 FADS hg19 / grch37 61043499	62159523

##########################################
# squamous cell esophageal cancer ######
################################################

# NCI-UGC
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC

zcat summary_chr11.txt.gz  | awk -F "\t" '{ if(($3 == 11) && ($4 >= 61043499 && $4 <= 62159523)) { print } }' > temp.txt
zcat summary_chr11.txt.gz | head -1 > head.txt
cat head.txt temp.txt > nci_ugi_escc_fads.txt
rm temp.txt
cd ~/fatty-acids/colocalisation/data/cancer
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC/nci_ugi_escc_fads.txt . 

# bbj GRCh37
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/esophageal_cancer
tar -xvf bbj.EsCa.rsq07.mac10.tar  
zcat EsCa.auto.rsq07.mac10.txt.gz  | awk -F " " '{ if(($1 == 11) && ($2 >= 61043499 && $2 <= 62159523)) { print } }' > temp.txt
zcat EsCa.auto.rsq07.mac10.txt.gz | head -1 > head.txt
cat head.txt temp.txt > bbj_escc_fads.txt
rm temp.txt
cd ~/fatty-acids/colocalisation/data/cancer
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/esophageal_cancer/bbj_escc_fads.txt . 

############################
# esophageal adenocarcinoma#
############################

#######
# EAS#
######

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/esophageal_adenocarcinoma

# grep -w rs174546 Oesophageal_adenocarcinoma_Lancet_Oncol2016_Gharahkhani_et_al.txt
awk -F " " '{ if(($2 == 11) && ($3 >= 61043499 && $3 <= 62159523)) { print } }' Oesophageal_adenocarcinoma_Lancet_Oncol2016_Gharahkhani_et_al.txt > temp.txt
head -1 Oesophageal_adenocarcinoma_Lancet_Oncol2016_Gharahkhani_et_al.txt > head.txt 
cat head.txt temp.txt > ea_eas_fads.txt
rm temp.txt
rm head.txt
cd ~/fatty-acids/colocalisation/data/cancer
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/esophageal_adenocarcinoma/ea_eas_fads.txt . 

#############
#UK Biobank#
############

cd /projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/
zcat maf_overall_oesoph_cancer_imputed.txt.gz  | awk -F "\t" '{ if(($2 == 11) && ($3 >= 61043499 && $3 <= 62159523)) { print } }' > temp.txt
zcat maf_overall_oesoph_cancer_imputed.txt.gz | head -1 > head.txt
cat head.txt temp.txt > ukb_ea_fads.txt
mv ukb_ea_fads.txt /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/esophageal_adenocarcinoma/
rm temp.txt
rm head.txt
cd ~/fatty-acids/colocalisation/data/cancer
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/esophageal_adenocarcinoma/ukb_ea_fads.txt .

##############
# lung cancer#
##############

# biobank japan GRCh37
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/lung_cancer
tar -xvf bbj.LuCa.rsq07.mac10.tar
zcat LuCa.auto.rsq07.mac10.txt.gz  | awk -F " " '{ if(($1 == 11) && ($2 >= 61043499 && $2 <= 62159523)) { print } }' > temp.txt
zcat LuCa.auto.rsq07.mac10.txt.gz | head -1 > head.txt
cat head.txt temp.txt > bbj_luca_fads.txt
rm temp.txt
tar -zcvf LuCa.auto.rsq07.mac10.txt.gz  .
cd ~/fatty-acids/colocalisation/data/cancer
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank/lung_cancer/bbj_luca_fads.txt . 

# UK Biobank
cd /projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/
zcat maf_overall_lung_cancer_imputed.txt.gz  | awk -F "\t" '{ if(($2 == 11) && ($3 >= 61043499 && $3 <= 62159523)) { print } }' > temp.txt
zcat maf_overall_lung_cancer_imputed.txt.gz | head -1 > head.txt
cat head.txt temp.txt > ukb_luca_fads.txt
mv ukb_luca_fads.txt /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/lung_cancer
rm temp.txt
cd ~/fatty-acids/colocalisation/data/cancer
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/lung_cancer/ukb_luca_fads.txt . 

# # ILCCO
# cd /mnt/storage/private/mrcieu/
# research/scratch/IGD/data/public/

# zcat maf_overall_lung_cancer_imputed.txt.gz | grep -w rs174546
# maf_overall_lung_cancer_unadj_imputed.txt.gz

# # ILCCO
# cp /projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc/Onc* /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/lung_cancer/
# cd /projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc
# ls Onc*
# wc Onco_TRICL_032116_Adeno.csv.tab 

##############
# skin cancer#
##############

# ukb
cd /projects/MRC-IEU/research/projects/icep1/wp1/030/working/data/results/
zcat maf_overall_nm_skin_cancer_imputed.txt.gz  | awk -F "\t" '{ if(($2 == 11) && ($3 >= 61043499 && $3 <= 62159523)) { print } }' > temp.txt
zcat maf_overall_nm_skin_cancer_imputed.txt.gz | head -1 > head.txt
cat head.txt temp.txt > /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/ukb_nmskin_fads.txt
rm temp.txt
cd ~/fatty-acids/colocalisation/data/cancer
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/ukb_nmskin_fads.txt . 


##########
# 23andme#
##########

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation
awk -F "\t" '{ if(($5 == "chr11") && ($6 >= 61043499 && $6 <= 62159523)) { print } }' all_snp_info-4.1.txt > temp.txt
head -1 all_snp_info-4.1.txt > head.txt 
cat head.txt temp.txt > annotation.txt
awk '{print $4}' annotation.txt > temp.txt
echo "$(tail -n +2 temp.txt)"  > temp2.txt ; mv temp2.txt temp.txt
awk '{print $1}' annotation.txt > all.data.id.txt
echo "$(tail -n +2 all.data.id.txt)"  > all.data.id2.txt ; mv all.data.id2.txt all.data.id.txt

fgrep -wf temp.txt gt_snp_stat-4.1.txt > temp2.txt
head -1 gt_snp_stat-4.1.txt > head.txt
cat head.txt temp2.txt > gt_annotations.txt

fgrep -wf temp.txt im_snp_stat-4.1.txt > temp2.txt
head -1 im_snp_stat-4.1.txt > head.txt
cat head.txt temp2.txt > im_annotations.txt

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_bcc/Chahal_2016_basal_cell_carcinoma-4.1/
fgrep -wf /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/all.data.id.txt basal_cell_carcinoma-4.1.dat > temp2.txt
head basal_cell_carcinoma-4.1.dat > head.txt
cat head.txt temp2.txt > 23me_bcc_fads.txt 

 
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_scc/
fgrep -wf /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/all.data.id.txt squamous_cell_carcinoma-4.1.dat  > temp2.txt
head squamous_cell_carcinoma-4.1.dat > head.txt
cat head.txt temp2.txt > 23me_scc_fads.txt 

cd ~/fatty-acids/colocalisation/data/cancer/23andme
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_scc/23me_scc_fads.txt  . 
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/chahal_bcc/Chahal_2016_basal_cell_carcinoma-4.1/23me_bcc_fads.txt   . 
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/im_annotations.txt . 
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/gt_annotations.txt . 
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/non_melanoma_skin_cancer_23andMe/annotation/annotation.txt . 


 

