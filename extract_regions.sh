#!/usr/bin/env bash

InputDir='/projects/MRC-IEU/users/ph14916/gtex2/all/'
OutputDir='/projects/MRC-IEU/users/ph14916/gtex2/all_fatty_acid_regions'
header='/projects/MRC-IEU/users/ph14916/gtex2/all/GTex.header.txt 
'

#eQTL-Gen

while read line;
do

  #Split out the line.
  arr=( $(IFS="\t" echo "$line") )


  filename=${arr[0]}
  chrom=${arr[1]}
  start=${arr[2]}
  end=${arr[3]}
  region=${arr[4]}

  echo "Processing "$filename"and"$chrom
  echo "Extracting chr"$chrom" from "$start" to "$end" bp ..."

  input=$InputDir$filename
  output=$OutputDir$filename"_"$region".txt"

  #Extract region using coloc.
  gunzip -cd $input | awk -v x=$chrom -v y=$start -v z=$end -F '\t' -v OFS='\t' '{split($2,a,"_"); if(a[1]==x &&a[2]>=y&&a[2]<=z)  print $0, a[1], a[2]}' - | cat $header - > $output

done<$1

#gene_list - unique list of ENSG ids to lookup in each GTEx file.

gunzip -cd GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Adipose_Subcutaneous.allpairs.txt.gz | cut -f1 > gene_list_temp.txt

sort gene_list_temp.txt -T '.' | uniq > uniq_gene_list.temp.txt
grep -f gene_list.txt gtex2/all/uniq_gene_list.temp.txt > genenames.versionid.txt

gunzip -cd /projects/MRC-IEU/users/ph14916/gtex2/all/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Adipose_Visceral_Omentum.allpairs.txt.gz | awk 'FNR==NR{a[$1]=1;next}($1 in a)' /projects/MRC-IEU/users/ph14916/genenames.versionid.txt - > /projects/MRC-IEU/users/ph14916/gtex2/all/GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Adipose_Visceral_Omentum.allpairs_mygenes.txt

 grep ENSG00000149485  gene_list_temp.txt > gene_temp2

  chr1_13550_G_A_b38  -16003  18  18  0.0154905 0.227454  0.243701  0.201673
ENSG00000227232.5 chr1_14671_G_C_b38  -14882  12  12  0.010320.0907104  0.419844  0.247715
ENSG00000227232.5 chr1_14677_G_A_b38  -14876  60  60  0.0516351 0.290965  -0.120766 0.114243
ENSG00000227232.5 chr1_16841_G_T_b38  -12712  50  50  0.0430293 0.523798  -0.0784818  0.123024
ENSG00000227232.5 chr1_16856_A_G_b38  -12697  13  13  0.0111876 0.338762  -0.224404 0.234362
ENSG00000227232.5 chr1_17005_A_G_b38  -12548  18  18  0.0154905 0.151864  -0.284061 0.197935
ENSG00000227232.5 chr1_17147_G_A_b38  -12406  18  18  0.0172414 0.219269  -0.242539 0.197189
ENSG00000227232.5 chr1_17407_G_A_b38