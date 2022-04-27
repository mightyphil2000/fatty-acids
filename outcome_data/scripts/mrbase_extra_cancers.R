library(TwoSampleMR)
ieugwasr::get_access_token()

snplist<-read.table("~/fatty-acids-mr/instruments/snplist_Europeans.txt",sep="\t",head=T,stringsAsFactors=F)
snps<-unique(snplist$SNP)


Cancers_selected_ukb<-c("Cancer code, self-reported: basal cell carcinoma","Diagnoses - main ICD10: C67 Malignant neoplasm of bladder","Malignant neoplasm of digestive organs","Type of cancer: ICD10: C64 Malignant neoplasm of kidney, except renal pelvis", "Lymphomas","Malignant neoplasm of skin","Malignant neoplasm of respiratory system and intrathoracic organs","Cancer code  self-reported: small intestine/small bowel cancer","Cancer code  self-reported: squamous cell carcinoma","Malignant neoplasm of urinary organs")
 
 Cancers_selected_fin<-c("Malignant neoplasm of bladder","Primary_lymphoid and hematopoietic malignant neoplasms","Malignant neoplasm of brain","Malignant neoplasm of breast","Any event in cancer register","Neoplasms","Malignant neoplasm of digestive organs","Malignant neoplasm of eye, brain and central nervous system","Colorectal cancer","Malignant neoplasm of endocrine gland","Malignant neoplasm of corpus uteri","malignant neoplasm of female genital organs","Follicular lymphoma","Malignant neoplasm of kidney, except renal pelvis","Lung cancer and mesothelioma","Malignant neoplasm of bronchus and lung","Lymphoid leukaemia","malignant neoplasm of male genital organs","Malignant neoplasm of skin","Multiple myeloma and malignant plasma cell neoplasms","Non-follicular lymphoma","Other and unspecified types of non-Hodgkin lymphoma","Malignant neoplasm of lip, oral cavity and pharynx","Malignant neoplasm of ovary","Malignant neoplasm of pancreas","Malignant neoplasm of prostate","Malignant neoplasm of respiratory system and intrathoracic organs","Malignant neoplasm of thyroid gland","Malignant neoplasm of urinary organs")

# summary(ao$ncontrol[grep("finn",ao$id)])
ao<-available_outcomes()
# ID<-ao$id[grep("multiple myeloma",ao$trait,ignore.case=T)]

# ID<-ID[grep("finn",ID)]
# Temp<-ao[ao$id %in% ID,]
# Temp[,c("ncase","ncontrol","trait","id")]

# # ID<-c("ukb-d-C_DIGESTIVE_ORGANS","finn-a-C3_DIGESTIVE_ORGANS_EXALLC","finn-a-C3_DIGESTIVE_ORGANS" ) 
# # Temp$id[Temp$trait=="Malignant neoplasm of breast"  ]
# # ID<-"ukb-d-C3_BREAST_3"

# ukb_dat <- extract_outcome_data(snps = "rs139402",outcomes = ID)
# ukb_dat[,c("originalname.outcome","beta.outcome","se.outcome","pval.outcome","eaf.outcome","effect_allele.outcome","other_allele.outcome","id.outcome")]
# beta<--0.00153779
# se<-0.000293362
# controls<-355504
# cases<-5690
# 	u<-cases/(cases+controls)
# 	lnor <- beta / (u * (1 - u))
# 	lnor.se<-se / (u * (1 - u)) 	

# # fingen A allele increases risk 
# # ukb A allele increases risk 
# A allele increases risk

# data.frame(ao[ao$trait == "malignant neoplasm of male genital organs (ICD C excluded)",])
ukb.id<-ao$id[which(ao$trait %in% Cancers_selected_ukb)]
ukb.id<-ukb.id[grep("ukb",ukb.id)]

ukb_dat <- extract_outcome_data(snps = snps,outcomes = ukb.id)

save(ukb_dat,file="~/MR_FattyAcids/data/summary_data/mrbase_extra/ukb_mrbase.Rdata")
library(TwoSampleMR)
ao<-available_outcomes()
# "Malignant neoplasm of breast"
# ID<-ao$id[ao$trait == "Lung cancer and mesothelioma"]
# ID<- c("finn-a-LUNG_CANCER","finn-a-LUNG_CANCER_MESOT")
# ao$ncase[ao$id %in% ID]
# ID<-"finn-a-CUSTOM_COLORECTAL_CANCER_EXALLC"
# snp <- "rs174546"
# fin_dat2 <- extract_outcome_data(snps = snp,outcomes = ID)



fin.id<-ao$id[which(ao$trait %in% Cancers_selected_fin)]
fin.id<-fin.id[grep("finn",fin.id)]
# fin.id<-fin.id[grep("Colorectal",fin.id,ignore.case=T)]
# ao$ncase[ao$id == fin.id]
# fin_dat <- extract_outcome_data(snps = "rs4939827",outcomes = fin.id)
# exp(0.1569)
fin_dat <- extract_outcome_data(snps = "rs174546",outcomes = fin.id)

fin_dat.m<-merge(fin_dat,ao[,c("id","ncase","ncontrol")],by.x="id.outcome",by.y="id")
Dups<-unique(fin_dat.m$originalname.outcome[duplicated(fin_dat.m$originalname.outcome)])
fin_dat.m[fin_dat.m$originalname.outcome %in% Dups,]
fin_dat.m<-fin_dat.m[order(fin_dat.m$ncase,fin_dat.m$ncontrol,decreasing=T),]
fin_dat.m<-fin_dat.m[!duplicated(fin_dat.m[,c("SNP","originalname.outcome")]),]
fin_dat<-fin_dat.m
save(fin_dat ,file="~/MR_FattyAcids/data/summary_data/mrbase_extra/fin_mrbase.Rdata")





# define  Cancers_selected_ukb  and Cancers_selected_fin

# 	ao1<-ao[ao$subcategory=="Cancer",]
# 	Terms<-c("cancer","malignancy","carcinoma","adenocarcinoma","glioma","leukaemia","leukemia","malignancies","lymphoma","sarcoma","mesothelioma","melanoma","neoplasm","meningioma","myeloma","neuroblastoma","blastoma")
# 	Terms2<-c("oma","emia")

# 	Pos<-unlist(lapply(Terms,FUN=function(x) grep(x,ao$trait,ignore.case=T)))
# 	# Pos2<-unlist(lapply(Terms2,FUN=function(x) grep(x,ao$trait,ignore.case=T)))
# 	# Pos2<-unique(Pos2)
# 	# Pos2<-Pos2[!Pos2 %in% Pos]
# 	# ao$trait[Pos2]

# 	Pos<-unique(Pos)
	
# 	ao2<-ao[Pos,]
# 	ao2<-ao2[!ao2$id %in% ao1$id, ]
# 	IDS<-unlist(strsplit(ao2$id,split="-"))
# 	IDS<-unique(IDS[seq(1,length(IDS),by=3)])
	
# 	ukb<-data.frame(ao2[grep("ukb",ao2$id),])

# 	sort(ukb$trait)
# 	excl<-c("Non-cancer","Illnesses of father","Illnesses of siblings","Illnesses of mother","Interpolated Age","Age at cancer diagnosis" ,	"Cancer year/age first occurred" ,	"Interpolated Year","screening","Number of self-reported cancers"   ,	"Operation code",	"Operative procedures","Chemotherapy session","Secondary malignant neoplasm","benign","Follow-up examination","Family history","Number of self-reported non-cancer illnesses","Intraductal carcinoma in situ")

# excl2<-c("Neoplasms","skin","colon","prostate","breast","ovary","bowel","basal cell carcinoma","squamous cell carcinoma","melanoma","lung","rectum","Cancer diagnosed by doctor","Reported occurrences of cancer","respiratory intrathoracic")

# # excl<-c(excl,excl2)

# 	Pos.excl<-unique(unlist(lapply(excl,FUN=function(x) grep(x,ukb$trait,ignore.case=T))))
# 	Pos<-1:nrow(ukb)
# 	Pos<-Pos[!Pos %in% Pos.excl]
# 	sort(ukb$trait[Pos])
# 	Cancers_selected_ukb<-c("Cancer code, self-reported: basal cell carcinoma","Diagnoses - main ICD10: C67 Malignant neoplasm of bladder","Malignant neoplasm of digestive organs","Type of cancer: ICD10: C64 Malignant neoplasm of kidney, except renal pelvis", "Lymphomas","Malignant neoplasm of skin","Malignant neoplasm of respiratory system and intrathoracic organs","Cancer code  self-reported: small intestine/small bowel cancer","Cancer code  self-reported: squamous cell carcinoma","Malignant neoplasm of urinary organs")
 
 
# 	  all(Cancers_selected_ukb %in% ukb$trait)
# 	ukb$trait[grep("Cancer code",ukb$trait)]
	
# 	Temp<-ukb[ukb$trait %in% Cancers_selected_ukb,c("ncase","ncontrol","trait")]
# 	write.table(ukb[Pos,c("trait","ncase")],"~/MR_FattyAcids/data/summary_data/ukbiobank/cancer_ukb.txt",sep="\t",col.names=T,row.names=F,quote=F)
# 	dim(ukb[ukb$ncase>1000,c("trait","ncase")])

# 	fin<-data.frame(ao2[grep("fin",ao2$id),])
# 	fin<-fin[order(fin$ncase,decreasing=T),]
# 	fin<-fin[grep("Benign",fin$trait,invert=T,ignore.case=T),]
# 	sort(fin$trait)
# 	write.table(fin[,c("trait","ncase")],"~/MR_FattyAcids/data/summary_data/FinnGen/cancer_fingen.txt",sep="\t",col.names=T,row.names=F,quote=F)



# Cancers_selected_fin<-c("Malignant neoplasm of bladder","Primary_lymphoid and hematopoietic malignant neoplasms","Malignant neoplasm of brain","Malignant neoplasm of breast","Any event in cancer register","Neoplasms","Malignant neoplasm of digestive organs","Malignant neoplasm of eye, brain and central nervous system","Colorectal cancer","Malignant neoplasm of endocrine gland","Malignant neoplasm of corpus uteri","malignant neoplasm of female genital organs","Follicular lymphoma","Malignant neoplasm of kidney, except renal pelvis","Lung cancer and mesothelioma","Malignant neoplasm of bronchus and lung","Lymphoid leukaemia","malignant neoplasm of male genital organs","Malignant neoplasm of skin","Multiple myeloma and malignant plasma cell neoplasms","Non-follicular lymphoma","Other and unspecified types of non-Hodgkin lymphoma","Malignant neoplasm of lip, oral cavity and pharynx","Malignant neoplasm of ovary","Malignant neoplasm of pancreas","Malignant neoplasm of prostate","Malignant neoplasm of respiratory system and intrathoracic organs","Malignant neoplasm of thyroid gland","Malignant neoplasm of urinary organs")
# 	Temp<-fin[fin$trait %in% Cancers_selected_fin ,c("ncase","ncontrol","trait")]
# 	Temp<-Temp[order(Temp$trait),]
# 	Temp<-Temp[!duplicated(Temp$trait),]
# 	write.table(Temp,"~/MR_FattyAcids/data/summary_data/FinnGen/cancer_fingen2.txt",sep="\t",col.names=T,row.names=F,quote=F)





# 	prot<-data.frame(ao2[grep("prot",ao2$id),])


# 	ebi<-data.frame(ao2[grep("ebi",ao2$id),])

# 	bbj<-data.frame(ao2[grep("bbj",ao2$id),])
# 	nrow(bbj)
# 	sort(bbj$trait)












