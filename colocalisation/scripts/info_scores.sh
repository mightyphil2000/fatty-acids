cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed/
awk -F "\t" '{ if(($2 == 11) && ($12 >=0.8)) { print } }' N6meta2041.tbl.fixed.tab > ~/temp.txt
head -1 N6meta2041.tbl.fixed.tab > ~/head.txt
cat ~/head.txt ~/temp.txt  > ~/aa_info_snps80.txt


awk -F "\t" '{ if(($2 == 11) && ($12 >=0.8)) { print } }' N6meta2031.tbl.fixed.tab  > ~/temp.txt
head -1 N6meta2031.tbl.fixed.tab   > ~/head.txt
cat ~/head.txt ~/temp.txt  > ~/dgla_info_snps80.txt

awk -F "\t" '{ if(($2 == 11) && ($12 >=0.8)) { print } }' N6meta1821.tbl.fixed.tab  > ~/temp.txt
head -1 N6meta1821.tbl.fixed.tab > ~/head.txt
cat ~/head.txt ~/temp.txt  > ~/la_info_snps80.txt

awk -F "\t" '{ if(($2 == 11) && ($12 >=0.8)) { print } }' N6meta1831.tbl.fixed.tab  > ~/temp.txt
head -1 N6meta1831.tbl.fixed.tab > ~/head.txt
cat ~/head.txt ~/temp.txt  > ~/gla_info_snps80.txt

cd ~/fatty-acids/colocalisation/data
scp ph14916@bluecrystalp3.acrc.bris.ac.uk:~/*info_snps80.txt . 


