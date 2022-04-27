sc 

scp ~/Downloads/UKB_telomere_gwas_summarystats.tsv.gz ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/OpenGWAS/ukb_telomere_length/ 


scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cismvMR/fads_r_matrix_*.ld . 

ls /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cismvMR/fads_r_matrix_*.ld

scp ~/Downloads/1000GP_Phase3.sample  ph14916@epi-franklin.epi.bris.ac.uk:~/1000genomes

scp ~/MR_FattyAcids/data/crc_test_dat.txt ph14916@epi-franklin.epi.bris.ac.uk:~/mrQC/inst/extdata/

scp ph14916@epi-franklin.epi.bris.ac.uk:~/mrQC/man/figures/README-example5.png . 

scp ph14916@epi-franklin.epi.bris.ac.uk:~/qc_report2.png . 

scp ~/fatty-acids/mr/data/exposure_dat_ara_la.Rdata ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/exposure_dat2_ara_la_mvmr_preharmonise.RData  . 

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_unique_rsidonly_clumped.txt . 

scp /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cancer_casecontrol_plus_genetic_data.RData /mnt/storage/scratch/ph14916 ph14916@bc4login.acrc.bris.ac.uk:/mnt/storage/home/ph14916/

scp ph14916@bc4login.acrc.bris.ac.uk:/mnt/storage/home/ph14916/sibling_results.Rdata . 

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/mr_crc_lc_msc_d5d.Rdata . 

scp coding38.tsv ph14916@epi-franklin.epi.bris.ac.uk:/data/ph14916/UkbCancerMortality/data/




scp /projects/MRC-IEU/users/ph14916/cancer_mortality/UKBB_cancer_mortality_incl_relateds_non_white_british.Rdata

cd ~/fatty-acids/outcome_data/data/harmonised
scp -r /Users/ph14916/ieugwasr_oauth ph14916@epi-franklin.epi.bris.ac.uk:~/mrQC/
scp -r /Users/ph14916/ieugwasr_oauth ph14916@epi-franklin.epi.bris.ac.uk:~/OpenGWAS/
scp ph14916@epi-franklin.epi.bris.ac.uk:~/qc_report*.png .
scp ph14916@epi-franklin.epi.bris.ac.uk:~/qc_report4.png .
scp ph14916@epi-franklin.epi.bris.ac.uk:~/qc_report5.png .
scp ph14916@epi-franklin.epi.bris.ac.uk:~/qc_report6.png .
cd /Users/ph14916/fatty-acids/exposure
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/cancer_mortality/*png .

scp ph14916@epi-franklin.epi.bris.ac.uk:~/qc_report_AA_GWAS.png . 
scp ph14916@epi-franklin.epi.bris.ac.uk:~/qc_report_UKB_GWAS.png . 

scp ~/ICEP_caseonly_WG/cancer_mortality_ukb/gwas/submit_spacox.sh ph14916@bc4login.acrc.bris.ac.uk:/mnt/storage/home/ph14916/ukb_gwas_mortality/spacox/


scp ~/ICEP_caseonly_WG/cancer_mortality_ukb/gwas/test.sh ph14916@bc4login.acrc.bris.ac.uk:/mnt/storage/home/ph14916/ukb_gwas_mortality/spacox/

scp ~/ICEP_caseonly_WG/cancer_mortality_ukb/gwas/run_SPAcox.R ph14916@bc4login.acrc.bris.ac.uk:/mnt/storage/home/ph14916/ukb_gwas_mortality/spacox/

ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/cancer_mortality/

scp ~/ICEP_caseonly_WG/cancer_mortality_ukb/gwas/ExtractsnpDose_UKB_MARK_all_overall_cancer.sh ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/cancer_mortality/

scp /projects/MRC-IEU/users/ph14916/cancer_mortality/UKBB_cancer_mortality_incl_relateds_non_white_british.Rdata ph14916@bc4login.acrc.bris.ac.uk:/mnt/storage/scratch/ph14916

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/cancer_mortality/overall_exclc44_genetic_ids.txt ph14916@bc4login.acrc.bris.ac.uk:~/ukb_gwas_mortality

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/cancer_mortality/multiple_diagnoses.txt .


scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gliomascan/glioma_test_dat.txt . 
scp 'ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/{
	glioma/gli66.txt,
	endometrial_cancer/enc25.txt,
	acute_lymphoblastic_leukemia_22076464/All21.txt,
	esophageal_adenocarcinoma/esa24.txt,
	bladder_cancer/Haycock/blc105.txt,
	uveal_melanoma/uvm165.txt,
	melanoma/mel95.txt,gwas_catalog/Cervical_cancer_28806749/cec92.txt,
	Bcell_nonHodgkinlymphoma_Chinese_23749188/bnh5.txt,UKbiobank/brc138.txt,
	UKbiobank/brc139.txt,
	UKbiobank/crc143.txt,
	UKbiobank/blc137.txt,
	UKbiobank/opc157.txt,
	UKbiobank/leu146.txt,
	UKbiobank/lbc147.txt,
	UKbiobank/lic148.txt,
	UKbiobank/luc149.txt,
	UKbiobank/luc1499.txt,
	UKbiobank/lle150.txt,
	UKbiobank/mum154.txt,
	UKbiobank/myl155.txt,
	UKbiobank/nmc156.txt,
	UKbiobank/esa144.txt,
	UKbiobank/ovc158.txt,
	UKbiobank/oac141.txt,
	UKbiobank/oac140.txt,
	UKbiobank/pro159.txt,
	UKbiobank/mel153.txt,
	Upper_gastrointestinal_cancers/gca102.txt,
	Upper_gastrointestinal_cancers/esc99.txt,
	Upper_gastrointestinal_cancers/gac101.txt,
	Upper_gastrointestinal_cancers/nga103.txt,
	gwas_catalog/Cervical_cancer_28806749/cec92.txt,
	Bcell_nonHodgkinlymphoma_Chinese_23749188/bnh5.txt,
	non_melanoma_skin_cancer_23andMe/scc2.txt,
	non_melanoma_skin_cancer_23andMe/bcc1.txt,
	ewings_sarcoma/ews27.txt,
	interlymph/cll83.txt,
	interlymph/fll85.txt,
	interlymph/dlb84.txt,
	interlymph/mzl86.txt,
	gwas_catalog/BRCA1_2_negative_high_risk_breast_cancer_30323354/hrb88.txt,
	gwas_catalog/Thyroid_cancer_30104761/thc163.txt,
	gwas_catalog/Kidney_cancer_31231134/kif90,gwas_catalog/Kidney_cancer_31231134/kim91,
	gwas_catalog/Ovarian_cancer_EastAsians_30898391/ovc120.txt,
	gliomascan/gli67.txt,
	practical/pro128.txt
	}' . 



# scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/gwas_catalog/Ovarian_cancer_EastAsians_30898391/ovc120.txt .
# cd ~/OpenGWAS/data/HNC_summStats_Oncoarray/
scp ph14916@epi-franklin.epi.bris.ac.uk:~/dat2.Rdata .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/rs7937840_INCENP_luc149.txt .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_eso99.txt .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_nmc156.txt .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_bcc1.txt .
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/disease/UKbiobank/INCENP_crc143.txt .

scp ~/fatty-acids/outcome_data/data/input_predlnor_sh_impinfo.Rdata ph14916@epi-franklin.epi.bris.ac.uk:~

scp ph14916@epi-franklin.epi.bris.ac.uk:~/Dat*_info_lnorsh.Rdata .