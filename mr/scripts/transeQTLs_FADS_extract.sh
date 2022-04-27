cd /projects/MRC-IEU/users/ph14916/eQTLGen/trans

gunzip 2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt.gz

grep FADS 2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt > fads_sig.txt 

grep ELOVL2 2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt > elovl2_sig.txt 


head -1 2018-09-04-trans-eQTLsFDR0.05-CohortInfoRemoved-BonferroniAdded.txt > head.txt

cat head.txt fads_sig.txt > fads2; mv fads2 fads_sig.txt

grep FADS 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > fads_all.txt 
head -1 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > head.txt
cat head.txt fads_all.txt > fads2; mv fads2 fads_all.txt

scp ~/fatty-acids-mr/instruments/ALLSNPs_Dec18_europeans.txt ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/ 

scp ~/fatty-acids-mr/instruments/targetSNPs_table_Dec18_europeans.txt ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/

scp ~/fatty-acids-mr/instruments/targetSNPs_table_Dec18.txt ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/




# awk -F, '$1 ==rs735665 && $8==FADS1 {print NR, $0}' 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > rs735665_FADS1_all.txt #didn't seem to work

# awk -F, '$1 ==rs735665 && $8==FADS2 {print NR, $0}' 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > rs735665_FADS2_all.txt

# awk -F, '$1 ==rs964184 && $8==FADS1 {print NR, $0}' 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > rs964184_FADS1_all.txt


# awk -F, '$1 ==rs964184 && $8==FADS2 {print NR, $0}' 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > rs964184_FADS2_all.txt

# awk -F, '$1 ==rs7412 && $8==FADS1 {print NR, $0}' 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > rs7412_FADS1_all.txt

# awk -F, '$1 ==rs7412 && $8==FADS2 {print NR, $0}' 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt > rs7412_FADS2_all.txt

# awk -F, '$75==100 && $76==200 && $79==300{print NR, $0}' file


# awk -v x=$chrom -v y=$start -v z=$end -F '\t' -v OFS='\t' '{if($11==x&&$12>=y&&$12<=z)  print}' $Input | cat $header - > $output

# Res
# rs735665
# rs964184
# rs7412

# FADS1
# FADS2

