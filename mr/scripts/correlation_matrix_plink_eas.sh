# https://github.com/MRCIEU/mr-base/wiki/Resource-list---Data-infrastructure#generate-1000-genomes-plink-files
cd /data/ph14916/1000genomes/EAS
Out=/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cismvMR

# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/

~/plink_dir/plink \
	--bfile EAS  \
	--from-bp 61043499 \
	--to-bp 62159523 \
	--chr 11 \
	--out $Out/fads_r_matrix_eas \
	--r square 
 
~/plink_dir/plink \
	--bfile EAS \
	--from-bp 61043499 \
	--to-bp 62159523 \
	--chr 11 \
	--out $Out/fads_r_matrix_eas \
	--freq 

cd $Out

