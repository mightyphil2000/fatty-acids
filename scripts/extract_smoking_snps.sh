##########################################
# Select GWAS sig SNPs for clumping#######
##########################################

cd /projects/MRC-IEU/scratch/for_Philip
cp -r /projects/MRC-IEU/scratch/for_Philip/GSCAN* /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/
cp -r CSI_results /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/GSCAN_incUKBiobank_results

head -1  cpd_no23andme.txt  > head_temp.txt
awk '{if ($7<5e-8) print }' cpd_no23andme.txt > cpd_no23andme_sig_temp.txt
cat head_temp.txt cpd_no23andme_sig_temp.txt > cpd_no23andme_sig.txt
rm *temp.txt

# Perform clumping on the above SNPs. Then extract the SNP list from the GSCAN results including and excluding UK Biobank
Dir=/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/GSCAN_noUKBiobank_results

###############
# including UKB#
################
head -1 $Dir/GSCAN_incUKBiobank_results/cpd_no23andme_sig.txt > head_temp

awk '{print $1}' $Dir/gscan_cpd_incUKB_snps_clump_relaxed.txt > temp 
sed -i '1d' temp 
grep -wf temp $Dir/GSCAN_incUKBiobank_results/cpd_no23andme_sig.txt > dat_temp
cat head_temp dat_temp  > $Dir/cpd_incUKB_sig_clump_relaxed.txt

awk '{print $1}' $Dir/gscan_cpd_incUKB_snps_clump_strict01.txt > temp #clump r2 0.01
sed -i '1d' temp 
grep -wf temp $Dir/GSCAN_incUKBiobank_results/cpd_no23andme_sig.txt > dat_temp
cat head_temp dat_temp  > $Dir/cpd_incUKB_sig_clump_strict01.txt

awk '{print $1}' $Dir/gscan_cpd_incUKB_snps_clump_strict.txt > temp #clump r2=0.001
sed -i '1d' temp 
grep -wf temp $Dir/GSCAN_incUKBiobank_results/cpd_no23andme_sig.txt > dat_temp
cat head_temp dat_temp  > $Dir/cpd_incUKB_sig_clump_strict.txt

########################################################
# excluding UKB from effect sizes estimate calculations#
########################################################

head -1 $Dir/GSCAN_noUKBiobank_results/cpd_noUKBiobank.txt > head_temp

awk '{print $1}' $Dir/gscan_cpd_incUKB_snps_clump_relaxed.txt > temp 
sed -i '1d' temp 
grep -wf temp $Dir/GSCAN_noUKBiobank_results/cpd_noUKBiobank.txt > dat_temp #SNP list comes from discovery study that included UKB. But betas/SEs are estimated excluded UKB
cat head_temp dat_temp  > $Dir/cpd_exclUKB_sig_clump_relaxed.txt

awk '{print $1}' $Dir/gscan_cpd_incUKB_snps_clump_strict01.txt > temp #clump r2 0.01
sed -i '1d' temp 
grep -wf temp $Dir/GSCAN_noUKBiobank_results/cpd_noUKBiobank.txt > dat_temp #SNP list comes from discovery study that included UKB. But betas/SEs are estimated excluded UKB
cat head_temp dat_temp  > $Dir/cpd_exclUKB_sig_clump_strict01.txt

awk '{print $1}' $Dir/gscan_cpd_incUKB_snps_clump_strict.txt > temp #clump r2=0.001
sed -i '1d' temp 
grep -wf temp $Dir/GSCAN_noUKBiobank_results/cpd_noUKBiobank.txt > dat_temp #SNP list comes from discovery study that included UKB. But betas/SEs are estimated excluded UKB
cat head_temp dat_temp  > $Dir/cpd_exclUKB_sig_clump_strict.txt

#########################
# CSI results UK Biobank#
##########################
Dir=/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/CSI_results
ls 

$Dir/2019.10.02\ Lifetime\ Smoking\ GWAS\ Data\ Sheet\ 1.txt 

head -1 $Dir/2019.10.02\ Lifetime\ Smoking\ GWAS\ Data\ Sheet\ 1.txt  >  head_temp
awk '{if ($10<5e-8) print }' $Dir/2019.10.02\ Lifetime\ Smoking\ GWAS\ Data\ Sheet\ 1.txt  > dat_temp
# awk '{print $10 }' $Dir/2019.10.02\ Lifetime\ Smoking\ GWAS\ Data\ Sheet\ 1.txt  | head
cat head_temp dat_temp > $Dir/csi_sig.txt
rm *temp



