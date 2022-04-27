# source("~/fatty-acids/colocalisation/scripts/extract_snps_function.R")

cd ~/plink_dir
./plink --bfile ~/1000genomes/EUR/EUR --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt    --freq  --out ~/1000genomes/EUR/fatty_acid_snps_eur

./plink --bfile ~/1000genomes/EAS/EAS --extract /projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_East_Asians_rsidsonly2_nodups.txt --freq  --out ~/1000genomes/EAS/fatty_acid_snps_eas

scp ph14916@epi-franklin.epi.bris.ac.uk:~/1000genomes/EAS/fatty_acid_snps_eas.frq .
scp ph14916@epi-franklin.epi.bris.ac.uk:~/1000genomes/EUR/fatty_acid_snps_eur.frq .



Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/"
# head N6meta2031.tbl.fixed.txt
# head N6meta2041.tbl.fixed.txt   
# Dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/charge/"
