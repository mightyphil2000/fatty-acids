cd /Users/ph14916/fatty-acids/colocalisation/data
awk -F " " '{ if(($1 == "chr11") && ($2 >= 61043499 && $2 <= 62159523)) { print } }' ~/fatty-acids/colocalisation/data/UKBB_10K_bed_hg19.txt > temp.txt
mv temp.txt ukb_fads_ref_hg19.txt


