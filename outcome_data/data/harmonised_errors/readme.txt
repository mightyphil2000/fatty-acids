These files will be deleted entirely and not replaced:


1. bnh5.txt

EAF does not match reference sample and there appears to be evidence for effect allele conflict in GWAS catalog. Will be very difficult to correct this information through correspondence because file was obtained via NAture Genetics editor. 

2. gli67

Effect allele frequency incorrect but effect allele correct. It will be easier to replace this with gli967 than to correct eaf

3. neu106 

effect allele looks wrong. The simplest solution is to replace this with neu107 which comes from the same study and which looks like it has the correct information. I merge neu106 and neu107 and effect direction looked conflicting for some SNPs, consistent with conflict observed with GWAS catalog, indicating effect allele wrong in neu16 file. 

these files will be replaced with corrected files

4. gli967

eaf and effect allele wrong. The simple solution is to swap effect and non effect alleles. 

5. glioma133

eaf column wrong. MAF was incorrectly assumed to be eaf. eaf set to NA 

6.ncc132 eaf wrong. set eaf to NA. 

7. All UK biobank datasets obtained via ukb batch in Open GWAS
convert beta to log odds ratios 
bcc135
bla136
cdo142
kic145
lym151
msc152
ric160
sbc161
scc162
utc164