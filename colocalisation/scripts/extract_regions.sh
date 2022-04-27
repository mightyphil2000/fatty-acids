#!/usr/bin/env bash
#to run: extract_regions.sh bash_input_gtex2.txt

InputDir='/projects/MRC-IEU/users/ph14916/gtex2/all/'
OutputDir='/projects/MRC-IEU/users/ph14916/gtex2/all_fatty_acid_regions/'
header='/projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt'

while read line;
do

  #Split out the line.
 
  arr=($(IFS="\t" echo "$line"))
 
  filename=${arr[0]}
  chrom=${arr[1]}
  start=${arr[2]}
  end=${arr[3]}
  region=${arr[4]}
  
  echo "Processing "$filename
  echo "Extracting chr"$chrom" from "$start" to "$end" bp ..."
  echo $region 
  
  input=$InputDir$filename
  output=$OutputDir$filename"_"$region".txt"
  gunzip -cd $input | awk -v x=$chrom -v y=$start -v z=$end -F '\t' -v OFS='\t' '{split($2,a,"_"); if(a[1]==x&&a[2]>=y&&a[2]<=z)  print $0, a[1], a[2]}' - | cat $header - > $output

done<$1
