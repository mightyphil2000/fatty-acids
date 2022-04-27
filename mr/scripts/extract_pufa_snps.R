# cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/
# head N6meta2041.tbl.fixed.tab
# grep rs2581624 N6meta2041.tbl.fixed.tab
# grep rs174548 N6meta2041.tbl.fixed.tab

# install.packages("devtools")
# devtools::install_github("MRCIEU/CheckSumStats")
library(ieugwasr)
# GWAS hits missing. 
Files_exclude<-c(
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C20_5n3_AGE_SEX_COHORT_ADJUST_filtered.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_5n6_AGE_SEX_COHORT_ADJUST_filtered.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_6n3_AGE_SEX_COHORT_ADJUST_filtered.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c226n3_pooled_allchr_qc1.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M01105.metal.pos.txt.gz.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M17805.metal.pos.txt.gz.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M34035.metal.pos.txt.gz.tab"
)
Files1<-Files[!Files %in% Files_exclude]
# Files1<-Files1[24:length(Files1)]
# Files1<-Files1[7:length(Files1)]

List<-NULL
for(i in 1:length(Files1)){
	print(Files1[i])
	# system.cmd<-paste("head -1",Files[i])
	# print(i)
	# system(system.cmd)	
	Dat<-extract_sig_snps(path_to_target_file=Files1[i],p_val_col_number=7,p_threshold=5e-8)
	Dat$File<-Files1[i]
	List[[i]]<-Dat
}

# length(List)
# Files1
Dat<-do.call(rbind,List)
table(Dat$File)

# try to extract sig SNPs in files_exclude through read.table
# didn't make a difference, no SNPs found
# List2<-NULL
# # Files_exclude<-Files_exclude[5:7]
# for(i in 1:length(Files_exclude)){
# 	print(Files_exclude[i])
# 	Res<-read.table(Files_exclude[i],sep="\t",head=TRUE,stringsAsFactors=FALSE,quote="")
# 	head(Res)
# 	Pos<-which(Res$p<5e-8)	
# 	if(length(Pos)>0){
# 		Res<-Res[Pos,]
# 		Res$File<-Files_exclude[i]
# 		List2[[i]]<-Res
# 	}
# }



# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/
# cd /newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/ratios/



# Open GWAS
# Kettunen GWAS
ids1<-find_ids_kettunen()
ids1<-c("met-c-856","met-c-855","met-c-917", "met-c-893", "met-c-852")



# omega-3 fatty acids
# omega-6 fatty acids
# polyunsaturated fatty acids other than 18:2
# linoleic acid (18:2n6)
# docosahexaenoic acid (22:6n3)

ids2<-finds_ids_ukb()
ids2<-c("met-d-Omega_6_pct","met-d-Omega_6_by_Omega_3","met-d-Omega_3_pct","met-d-DHA_pct","met-d-LA_pct")

Dat2<-tophits(id=c(ids1,ids2))
Dat2<-fomat_dat_og()

save(Dat,file="~/extract_sig_snps_pufas.RData")
save(Dat2,file="~/extract_sig_snps_pufas_opengwas.RData")

scp ~/extract_sig_snps* ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_instruments/

# Circulating omega-3 (i.e. docosahexaenoic acid (DHA) and total omega-3) and omega-6 (i.e. linoleic acid (LA) and total omega-6) fatty acids concentration 


fomat_dat_og<-function(){
	Dat2<-data.frame(Dat2)	
	names(Dat2)[names(Dat2) == "ea"]<-"effect_allele"
	names(Dat2)[names(Dat2) == "nea"]<-"other_allele"
	names(Dat2)[names(Dat2) == "eaf"]<-"effect_allele_freq"
	names(Dat2)[names(Dat2) == "rsid"]<-"snp"
	return(Dat2)
}


extract_sig_snps<-function(path_to_target_file=NULL,p_val_col_number=NULL,p_threshold=5e-8){
    path_to_outfile<-file.path(tempdir(), "output.txt")
    path_to_filehead<-file.path(tempdir(), "filehead.txt")
    path_to_outfile_plus_filehead<-file.path(tempdir(), "output_head.txt")

    system(paste0("head -1  ",path_to_target_file," > ",path_to_filehead))
    system(paste0("awk '{if ($",p_val_col_number,"<",p_threshold,") print }' ",path_to_target_file," > ",path_to_outfile))
    system(paste0("cat ",path_to_filehead," ",path_to_outfile," > ",path_to_outfile_plus_filehead))
    sig_dat<-utils::read.table(path_to_outfile_plus_filehead,head=TRUE,stringsAsFactors=FALSE,fill=TRUE)
    return(sig_dat)
}


find_ids_kettunen<-function(){
	ao<-gwasinfo()
	Pos<-which(ao$author == "Kettunen")
	ao1<-ao[Pos,]
	Pos1<-grep("omega-6",ao1$trait,ignore.case=TRUE)
	Pos2<-grep("omega-3",ao1$trait,ignore.case=TRUE)
	Pos3<-grep("polyunsaturated",ao1$trait,ignore.case=TRUE)
	Pos4<-grep("linoleic",ao1$trait,ignore.case=TRUE)
	Pos5<-grep("docosahexaenoic",ao1$trait,ignore.case=TRUE)
	ids<-ao1$id[c(Pos1,Pos2,Pos3,Pos4,Pos5)]
	return(ids)
}


finds_ids_ukb<-function(){
	ao<-gwasinfo()
	Pos<-grep("met-d",ao$id)
	ao1<-ao[Pos,]
	ao1$trait
	Pos1<-grep("omega",ao1$trait)
	Pos2<-grep("docosahexaenoic",ao1$trait)
	Pos3<-grep("linoleic",ao1$trait)
	ids<-ao1$id[c(Pos1,Pos2,Pos3)]
	return(ids)
}


Files<-c("/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/CHARGE_N3_ALA_fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/CHARGE_N3_EPA_fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/CHARGE_N3_DPA_fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/CHARGE_N3_DHA_fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/N6meta1821.tbl.fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/N6meta1831.tbl.fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/N6meta2031.tbl.fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/N6meta2041.tbl.fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/N6meta2241.tbl.fixed.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C18_2n6_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C18_3n3_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C18_3n6_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C20_2n6_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C20_3n6_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C20_4n6_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C20_5n3_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_4n6_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_5n3_AGE_SEX_COHORT_ADJUST_filtered.tab","/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_5n6_AGE_SEX_COHORT_ADJUST_filtered.tab"
	,"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_6n3_AGE_SEX_COHORT_ADJUST_filtered.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c182n6_pooled_allchr_qc1.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c183n3_pooled_allchr_qc1.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c183n6_pooled_allchr_qc1.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c203n6_pooled_allchr_qc1.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c204n6_pooled_allchr_qc1.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c205n3_pooled_allchr_qc1.tab",
	"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_log/score_c226n3_pooled_allchr_qc1.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M01105.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M01110.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M17805.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M18467.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M19323.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M32504.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M32980.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M33969.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M34035.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M34674.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M35718.metal.pos.txt.gz.tab",
"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/M33884.metal.pos.txt.gz.tab")



# awk '{if ($7<5e-08) print }' /newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_5n6_AGE_SEX_COHORT_ADJUST_filtered.tab

# head /newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/FHS_RBC_C22_5n6_AGE_SEX_COHORT_ADJUST_filtered.tab

#  > /tmp/Rtmp5HVeHJ/output.txt

# [1] 