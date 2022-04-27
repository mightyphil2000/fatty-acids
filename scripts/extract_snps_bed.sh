#!/usr/bin/env bash
#to run: extract_regions.sh bash_input_gtex2.txt

Input='Ovarian_cancer_EastAsians_30898391.txt'
echo "Processing "$Input

while read line;
do

  #Split out the line.
 
  arr=($(IFS="\t" echo "$line"))
   
  chrom=${arr[0]}
  start=${arr[1]}
  
  echo "Extracting chr"$chrom" and "$start""
   
  awk -v x=$chrom -v y=$start -F ',' -v OFS='\t' '{if($3==x&&$4==y)  print}' $Input >> /projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Ovarian_cancer_EastAsians_30898391/header.txt

done<$1
