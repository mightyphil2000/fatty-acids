# AA and DGLA in Dorajoo
cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo
#allele b is effect allele
# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/Dorajoo26584805/

head -1 score_c204n6_pooled_allchr_qc1.tab 

# AA
grep -w rs174546 beta_log/score_c204n6_pooled_allchr_qc1.tab 
snp	effect_allele	other_allele	effect_allele_freq	beta	se	p	n
rs174546	T	C	0.6619	-0.6049	0.0399	5.521e-52	1361

# DGLA
grep -w rs174546 beta_log/score_c203n6_pooled_allchr_qc1.tab
snp	effect_allele	other_allele	effect_allele_freq	beta	se	p	n
rs174546	T	C	0.6618	-0.3923	0.0389	5.933e-24	1361

head -1 beta_log/score_c183n6_pooled_allchr_qc1.tab
# GLA
grep -w beta_log/score_c183n6_pooled_allchr_qc1.tab
rsid	chr	pos	allele_a	allele_b	info	hwe	p_value	beta	sbaf
rs174546	11	61569830	C	T	1	.93281001	0	-.83222002	.053908002	.66713482

# LA
grep -w rs174546 score_c182n6_case_allchr_qc.txt
rsid	chr	pos	allele_a	allele_b	info	hwe	p_value	beta	sbaf
rs174546	11	61569830	C	T	1	.93281001	7.7050997e-09	.32032001	.055468	.66713482

head -1 beta_pc/score_c204n6_pooled_allchr_qc1.tab
grep -w rs174546 beta_pc/score_c204n6_pooled_allchr_qc1.tab
snp	effect_allele	other_allele	effect_allele_freq	beta	se	p	n
rs174546	T	C	0.6621	-1.06	0.069	2.959e-53	1361


grep -w rs174546 beta_pc/score_c203n6_pooled_allchr_qc1.tab
snp	effect_allele	other_allele	effect_allele_freq	beta	se	p	n
rs174546	T	C	0.6616	-0.1074	0.0131	2.464e-16	1361

# ratios
cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_pc/ratios
zcat AA_to_DGLA_pooled.tab.gz | grep -w rs174546
zcat AA_to_DGLA_pooled.tab.gz | head -1

SNP	effect_allele	other_allele	eaf	beta	se	pval	samplesize
rs174546	T	C	0.6621	-0.232751815858763	0.14566243729192	0.1100685848019	1361


# Arachidonic acid and DGLA in Charge 
cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/

head -1 N6meta2041.tbl.fixed.tab
grep -w rs174546 N6meta2041.tbl.fixed.tab
snp	effect_allele	other_allele	effect_allele_freq	beta	se	p	n
rs174546	t	c	0.3207	-1.6871	0.0253	0	8631

grep -w rs174546 N6meta2031.tbl.fixed.tab

snp	effect_allele	other_allele	effect_allele_freq	beta	se	p	n
rs174546	t	c	0.3262	0.3544	0.0135	1.8e-151	8631

#AA and DGLA framingham
cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle
head -1 FHS_RBC_C20_3n6_AGE_SEX_COHORT_ADJUST_filtered.tab

# DGLA
grep -w rs174546 FHS_RBC_C20_3n6_AGE_SEX_COHORT_ADJUST_filtered.tab
snp	beta	se	p	effect_allele_freq	effect_allele	other_allele	chr_name	chrom_start	N
rs174546	-0.0762403	0.0021468	5.87211e-225	0.667374	C	T11	61569830	2555

# AA
grep -w  rs174546 FHS_RBC_C20_4n6_AGE_SEX_COHORT_ADJUST_filtered.tab

snp	beta	se	p	effect_allele_freq	effect_allele	other_allele	chr_name	chrom_start	N
rs174546	0.0277407	0.00146332	4.23033e-75	0.667374	C	T	11	61569830	2555



#AA:DGLA and GLA:LA in CHARGE

cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/ratios
head -1 AA_to_DGLA.tab   
grep -w rs968567 AA_to_DGLA.tab   

b<-0.0277407
se<-0.00146332
z<-0.0277407/0.00146332
maf<-1-0.667374
n<-2555



b_sd1<-b_sd(z=z,maf=maf,n=n)
se<-b_sd1/z

b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}





