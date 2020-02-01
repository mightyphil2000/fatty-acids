#!/usr/bin/env bash

#$1 - outcomes to be extracted.

WorkingDir='/Users/epdab/Documents/Jay_eQTL_project/coloc/'

coloc_dir=$WorkingDir"exposure_1Mb/"
output_dir1=$WorkingDir"outcome_1Mb/"
output_dir2=$WorkingDir"harmon_1Mb/"

#eQTL-Gen
#coloc_dir='/Users/epdab/Documents/biib_genes_projectMarch2019/coloc/eQTLGen_20190322/'
#output_dir=$coloc_dir"harmonised/"

while read line;
do
  echo "Processing "$line

  #Split out the line.
  arr=( $(IFS=" " echo "$line") )

  region=${arr[1]}
  gene=${arr[0]}
  outcome=${arr[2]}

  input_file=$coloc_dir$gene"."$region".txt"
  output_file1=$output_dir1$gene"."$outcome"."$region".txt"
  output_file2=$output_dir2$gene"."$outcome"."$region".txt"

  echo "Input GWAS file: "$input_file
  echo "Outcome file: "$output_file1
  echo "Harmon file: "$output_file2

  Rscript extract_outcomes.coloc.R $gene $outcome $input_file $output_file1 $output_file2

done < $1
