#!/usr/bin/env bash
#to run: extract_regions.sh bash_input_gtex2.txt

Input='/projects/MRC-IEU/users/ph14916/eQTLGen/combined-eQTLs_full.EAF.beta.se.chr.pos.txt'
OutputDir='/projects/MRC-IEU/users/ph14916/eQTLGen/'
header='/projects/MRC-IEU/users/ph14916/eQTLGen/combined-eQTLs_full.EAF.beta.se.chr.pos_header.txt'
echo "Processing "$Input

while read line;
do

  #Split out the line.
 
  arr=($(IFS="\t" echo "$line"))
   
  chrom=${arr[0]}
  start=${arr[1]}
  end=${arr[2]}
  region=${arr[3]}
  

  echo "Extracting chr"$chrom" from "$start" to "$end" bp ..."
  echo $region 
  
  output=$OutputDir"eqtlgen_"$region".txt"
  awk -v x=$chrom -v y=$start -v z=$end -F '\t' -v OFS='\t' '{if($11==x&&$12>=y&&$12<=z)  print}' $Input | cat $header - > $output

done<$1
