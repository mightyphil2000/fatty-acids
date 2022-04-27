# ENSG00000073756 = PTGS2 / COX2
# ENSG00000095303 = PTGS1 / COX1

#############################
# Extract eQTLs from eQTLGen#
#############################

cd /projects/MRC-IEU/users/ph14916/eQTLGen
gunzip ENSG00000073756.eQTL.tab.gz
gunzip ENSG00000095303.eQTL.tab.gz

# head ENSG00000073756.eQTL.tab #PTGS2 / COX2
# head ENSG00000095303.eQTL.tab # PTGS1 / COX1

# awk '{print $7}' ENSG00000073756.eQTL.tab | head

# ENSG00000073756.eQTL.tab
head -1 ENSG00000073756.eQTL.tab > head_temp
awk '{if ($7<5e-8) print }' ENSG00000073756.eQTL.tab > dat_temp
cat head_temp dat_temp > /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox2_ptgs2_sig_eqtlgen.txt

# ENSG00000095303.eQTL.tab
head -1 ENSG00000095303.eQTL.tab > head_temp
awk '{if ($7<5e-8) print }' ENSG00000095303.eQTL.tab > dat_temp
cat head_temp dat_temp > /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox1_ptgs1_sig_eqtlgen.txt

rm *temp*

wc /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox2_ptgs2_sig_eqtlgen.txt
wc /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/cox1_ptgs1_sig_eqtlgen.txt

gzip ENSG00000073756.eQTL.tab
gzip ENSG00000095303.eQTL.tab


####################################
# Extract eQTLs from GTEx version 8#
####################################
indir=/projects/MRC-IEU/research/data/broad/public/gtex/released/2020-03-09/data/v8_eQTL_all_associations
outdir=/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking
# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking

# ENSG00000073756 = PTGS2 / COX2
for infile in {Whole_Blood.allpairs.txt.gz,Lung.allpairs.txt.gz,Colon_Sigmoid.allpairs.txt.gz,Colon_Transverse.allpairs.txt.gz,Esophagus_Gastroesophageal_Junction.allpairs.txt.gz,Esophagus_Mucosa.allpairs.txt.gz,Esophagus_Muscularis.allpairs.txt.gz,Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz,Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz};do
	echo "infile="$infile
	outfile=$outdir"/"$infile"_PTGS2_COX2.txt"	
	echo "outfile="$outfile
	gunzip -cd $indir"/"$infile | head -1 > $outdir/head_temp
	gunzip -cd $indir"/"$infile | grep ENSG00000073756 | awk '{if ($7<5e-4) print }' > $outdir/dat_temp
	cat $outdir/head_temp $outdir/dat_temp  > $outfile
	rm $outdir/*temp*
done

# ENSG00000095303 = PTGS1 / COX1
for infile in {Whole_Blood.allpairs.txt.gz,Lung.allpairs.txt.gz,Colon_Sigmoid.allpairs.txt.gz,Colon_Transverse.allpairs.txt.gz,Esophagus_Gastroesophageal_Junction.allpairs.txt.gz,Esophagus_Mucosa.allpairs.txt.gz,Esophagus_Muscularis.allpairs.txt.gz,Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz,Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz};do
	echo "infile="$infile
	outfile=$outdir"/"$infile"_PTGS1_COX1.txt"	
	echo "outfile="$outfile
	gunzip -cd $indir"/"$infile | head -1 > $outdir/head_temp
	gunzip -cd $indir"/"$infile | grep ENSG00000095303 | awk '{if ($7<5e-4) print }' > $outdir/dat_temp
	cat $outdir/head_temp $outdir/dat_temp  > $outfile
	rm $outdir/*temp*
done


# Whole_Blood.allpairs.txt.gz
# Lung.allpairs.txt.gz

# Colon_Sigmoid.allpairs.txt.gz
# Colon_Transverse.allpairs.txt.gz

# Esophagus_Gastroesophageal_Junction.allpairs.txt.gz
# Esophagus_Mucosa.allpairs.txt.gz
# Esophagus_Muscularis.allpairs.txt.gz

# Skin_Not_Sun_Exposed_Suprapubic.allpairs.txt.gz
# Skin_Sun_Exposed_Lower_leg.allpairs.txt.gz

# gunzip -cd $input | grep ENSG00000073756 | head 
# gunzip -cd $input | awk '{print $1}' 