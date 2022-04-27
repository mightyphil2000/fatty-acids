# https://github.com/MRCIEU/mr-base/wiki/Resource-list---Data-infrastructure#generate-1000-genomes-plink-files
cd ~/1000genomes/EAS
Out=/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cismvMR


grep CHS ~/1000genomes/1000GP_Phase3.sample > ~/1000genomes/EAS/1000GP_Phase3_chs.sample
wc ~/1000genomes/EAS/1000GP_Phase3_chs.sample
# awk '{print $1 " " $1 " " 0 " " 0}' ~/1000genomes/EAS/1000GP_Phase3_jpt.sample > ~/1000genomes/EAS/1000GP_Phase3_jpt.sample_plink 

awk '{print $1 " " $1}' ~/1000genomes/EAS/1000GP_Phase3_chs.sample > ~/1000genomes/EAS/1000GP_Phase3_chs.sample_plink 

wc ~/1000genomes/EAS/1000GP_Phase3_chs.sample_plink 

# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/

~/plink_dir/plink \
	--bfile ~/1000genomes/EAS/EAS  \
	--from-bp 61043499 \
	--to-bp 62159523 \
	--keep ~/1000genomes/EAS/1000GP_Phase3_chs.sample_plink  \
	--chr 11 \
	--out $Out/fads_r_matrix_chs \
	--r square 

~/plink_dir/plink \
	--bfile ~/1000genomes/EAS/EAS  \
	--from-bp 61043499 \
	--to-bp 62159523 \
	--keep ~/1000genomes/EAS/1000GP_Phase3_chs.sample_plink \
	--chr 11 \
	--out $Out/fads_r_matrix_chs \
	--freq
 
cd $Out
