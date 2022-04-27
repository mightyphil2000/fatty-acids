
# bladder cancer
scp -r ~/MR_FattyAcids/data/summary_data/bladder_cancer/ ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/
gunzip disease/bladder_cancer/Haycock/*.gz
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/bladder_cancer/Haycock/
cat chr1_GWAS_NBCS_risico_jan2017_methodscore_4mds_info.out  > temp.out
for i in {2..22} 
do
	echo "chr${i}_GWAS_NBCS_risico_jan2017_methodscore_4mds_info.out"		 	
    sed 1,1d "chr${i}_GWAS_NBCS_risico_jan2017_methodscore_4mds_info.out"  >> temp.out
done

mv temp.out chrALL_GWAS_NBCS_risico_jan2017_methodscore_4mds_info.out
wc chrALL_GWAS_NBCS_risico_jan2017_methodscore_4mds_info.out
gzip chr*


########
# melanoma
########

scp -r ~/MR_FattyAcids/data/summary_data/melanoma/ ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/melanoma/Melanoma_meta_single_files/
gunzip *.gz


cat Melanoma_metaanalysis_CHR1_15052019_EAF_RSQ.txt   > temp.out
for i in {2..22} 
do
	echo "Melanoma_metaanalysis_CHR${i}_15052019_EAF_RSQ.txt"		 	
    sed 1,1d "Melanoma_metaanalysis_CHR${i}_15052019_EAF_RSQ.txt"  >> temp.out
done

sed 1,1d Melanoma_metaanalysis_CHRX_15052019_EAF_RSQ.txt  >> temp.out
mv temp.out Melanoma_metaanalysis_chrALL_15052019_EAF_RSQ.txt
gzip chr*


###################################
# Upper gastrointestinal cancers#####
########################################

scp -r Upper_gastrointestinal_cancers/ ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/


# CC 
# Gastric cardia adenocarcinoma
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_CC

cat summary_chr1.txt > temp.out
for i in {2..22} 
do
	echo "summary_chr$i.txt"		 
	sed 1,1d "summary_chr$i.txt"  >> temp.out
done

gzip summary*
mv temp.out summary_chr_all.txt

# ESCC
# Esophageal squamous cell carcinoma

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_ESCC

cat summary_chr1.txt > temp.out
for i in {2..22} 
do
	echo "summary_chr$i.txt"		 
	sed 1,1d "summary_chr$i.txt"  >> temp.out
done

gzip summary*
mv temp.out summary_chr_all.txt



# gastric
# Gastric adenocarcinoma

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_gastric

cat summary_chr1.txt > temp.out
for i in {2..22} 
do
	echo "summary_chr$i.txt"		 
	sed 1,1d "summary_chr$i.txt"  >> temp.out
done

gzip summary*
mv temp.out summary_chr_all.txt

# NC
# Noncardia gastric adenocarcinoma
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Upper_gastrointestinal_cancers/Asian_UGI_Study_NC

# grep -w rs140081212 summary_chr_all.txt #top hit for gastric adenocarcinoma in gwas catalog
# head -1 summary_chr_all.txt

# eaf=0.83202
# G 0.32921
# 1.38987

cat summary_chr1.txt > temp.out
for i in {2..22} 
do
	echo "summary_chr$i.txt"		 
	sed 1,1d "summary_chr$i.txt"  >> temp.out
done

gzip summary*
mv temp.out summary_chr_all.txt

######################
# Japanese Biobank####
########################
# now uploaded to Open GWAS 
# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/Japanese_Biobank
# # gunzip *
# gunzip auto_rsq07.mac10.tar.gz 
# tar -xvf auto_rsq07.mac10.tar 
# cd auto_rsq07.mac10
# gunzip *


# Files<-dir()[grep("tar",dir())]
# for(i in 1:length(Files)){
#     print(i)
#     print(Files[i])
#     system(paste0("tar -xvf ",Files[i]))    
# }


############################################
# Chronic lymphocytic leukemia Interlymph###
############################################
# source("~/fatty-acids/colocalisation/scripts/Extract_SNPs_function.R")
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/CLL
gunzip * 

# tail -n +1 result* > cll.out
# cat result*  > cll2.out

sed 1,1d result_chr1.out  > cll.out
for i in {2..22} 
	do sed 1,1d "result_chr$i.out"  >> cll.out
done

head -1 result_chr1.out > head.txt
cat head.txt cll.out > cll2.out; mv cll2.out cll.out

# sed '$ d' cll.out > foo.txt #delete last row



#####################################
# Follicular lymphoma Interlymph ####
#####################################

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/FL
gunzip * 

sed 1,1d result_chr1.out  > fl.out
for i in {2..22} 
	do sed 1,1d "result_chr$i.out"  >> fl.out
done

head -1 result_chr1.out > head.txt
cat head.txt fl.out > fl2.out; mv fl2.out fl.out

# sed '$ d' cll.out > foo.txt #delete last row


###########################################
# Diffuse large B cell lymphoma interlymph####
###########################################

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/DLBCL
gunzip * 

sed 1,1d result_chr1.out  > dlbcl.out
for i in {2..22} 
	do sed 1,1d "result_chr$i.out"  >> dlbcl.out
done

head -1 result_chr1.out > head.txt
cat head.txt dlbcl.out > dlbcl2.out; mv dlbcl2.out dlbcl.out

# sed '$ d' cll.out > foo.txt #delete last row


#####################################
# Marginal zone lymphoma        ####
#####################################


cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/interlymph/MZL
gunzip * 

sed 1,1d result_chr1.txt.out  > mzl.out
for i in {2..22} 
	do sed 1,1d "result_chr$i.txt.out"  >> mzl.out
done

head -1 result_chr1.txt.out > head.txt
cat head.txt mzl.out > mzl2.out; mv mzl2.out mzl.out

# sed '$ d' cll.out > foo.txt #delete last row


#################################################### 
# Ovarian_cancer_EastAsians_30898391 OCAC
####################################################

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Ovarian_cancer_EastAsians_30898391
cat SummaryResults_Asian_chr1.txt > temp.out
for i in {2..23} 
do
	echo "SummaryResults_Asian_chr$i.txt"
	 # echo "meta_v3_onco_euro_caseonly_chr${i}_1.txt"	 
	 sed 1,1d "SummaryResults_Asian_chr$i.txt"  >> temp.out
done

mv temp.out SummaryResults_Asian_chr_all.txt

# lung cancer Patel 27488534cd 
cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/
cp -r /projects/MRC-IEU/research/data/evidencehub/summary/gwas/raw/TRICL_LungCancer /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/
ls 

