# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/

setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/")
Gen1<-read.table("ukb_snps.txt",sep=" ",head=TRUE, stringsAsFactors=FALSE,fill=T) 
Gen2<-read.table("ukb_cox_csi_snps.txt",sep=" ",head=TRUE, stringsAsFactors=FALSE,fill=T)

Gen<-merge(Gen1,Gen2,by="ID")
Gen<-Gen[,grep("y",names(Gen),invert=T)]
names(Gen)<-gsub(".x","",names(Gen))

Phen<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/UKBB_cancer_aspirin_smoking.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)
# Phen<-dat
Dat<-merge(Phen,Gen,by.x="geneticID",by.y="ID")

# Dat2<-Dat
# Dat<-Dat2
# table(Dat$sig_cancers)
Dat<-sig_can_function()
Dat<-lc_crc_function()
Dat<-lc_crc_nm_cancers_function()
Dat<-smoking_cancers_function()

Dat<-format_cancer_coding() #change to 0/1 coding
Dat<-format_smoking()
Dat<-format_aspirin()

# Dat$GS_d5d<-Dat$rs174546*-0.86
Dat<-d5d_function()
Dat<-csi_genetic_score()
Dat<-cpd_genetic_score()
Dat<-cox_gtex_genetic_score()
Dat<-cox_eqtlgen_genetic_score()


Datsmk1<-Dat[Dat$smoking2 == 1,]
Datsmk0<-Dat[Dat$smoking2 == 0,]

Datasp1<-Dat[Dat$aspirin2 == 1,]
Datasp0<-Dat[Dat$aspirin2 == 0,]

save(list=c("Dat","Datsmk1","Datsmk0","Datasp1","Datasp0"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/UKBB_cancer_aspirin_smoking_plusgenetic_data.RData")


format_aspirin<-function(){
	Dat$aspirin2<-NA
	Pos_yes<-which(Dat$aspirin_baseline_f_20003 == "yes" & Dat$aspirin_baseline_f_6154 == "yes")
	Pos_no<-which(Dat$aspirin_baseline_f_20003 == "no" & Dat$aspirin_baseline_f_6154 == "no")
	Dat$aspirin2[Pos_yes]<-1
	Dat$aspirin2[Pos_no]<-0
	return(Dat)
}

format_cancer_coding<-function(){
	Cancers1<-names(Dat)[grep("cancer",names(Dat))]
	Cancers2<-names(Dat)[grep("overall",names(Dat))]
	Cancers3<-names(Dat)[grep("incident",names(Dat))]
	Cancers<-unique(c(Cancers1,Cancers2,Cancers3))
	for(i in Cancers){
		print(i)
		Dat[,i]<-Dat[,i]-1
	}
	return(Dat)
}


d5d_function<-function(){
	# change coding so that effect allele is allele c (the major allele and D5D raising allele )
	Dat$d5d<-NA
	Dat$d5d[Dat$rs174546==2]<-0
	Dat$d5d[Dat$rs174546==1]<-1
	Dat$d5d[Dat$rs174546==0]<-2
	Dat$d5d<-Dat$d5d*0.86
	return(Dat)
	# # effect of C allele on D5D activity is 0.86 SD units
	# Dat$rs174546_3<-Dat$rs174546_2 * 0.86
}

format_smoking<-function(){
	Dat$smoking<-as.character(Dat$smoking)
	Dat$smoking[Dat$smoking %in% c("Current","Previous")]<-"Ever"	
	Dat$smoking[Dat$smoking == "Prefer not to answer"]<-NA
	Dat$smoking2<-NA
	Dat$smoking2[which(Dat$smoking=="Ever")]<-1
	Dat$smoking2[which(Dat$smoking=="Never")]<-0
	return(Dat)
}


cpd_genetic_score<-function(){
# https://conservancy.umn.edu/bitstream/handle/11299/201564/README.txt?sequence=29&isAllowed=y
	cpd<-read.table("cpd_no23andme_sig_clump_relaxed.txt",sep="\t",head=TRUE,stringsAsFactors=F)
	snptest<-read.table("ukb_cpd_snp_stats_all.txt",sep=" ",head=T,stringsAsFactors=F)
	cpd_snptest<-merge(snptest,cpd,by.x="rsid",by.y="RSID")	
	cpd_snptest<-harmonise_effect_allele(Dat_harmonise=cpd_snptest,ea.x="alleleB",ea.y="ALT",oa.x="alleleA",oa.y="REF",eaf.x="alleleB_frequency",eaf.y=NULL,b.y="BETA")            
	Dat$GS_cpd<-create_gs(gs_dat=cpd_snptest,BETA="BETA")
	return(Dat)
}

csi_genetic_score<-function(){
	csi<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Csi_sig_clumped.txt",sep=" ",head=TRUE,stringsAsFactors=F)
	csi2<-csi[csi$clump_r2 == 0.3,]
	snptest<-read.table("ukb_snp_stats_all.txt",sep=" ",head=T,stringsAsFactors=F)
	# snptest2<-snptest[,c("rsid","alleleA","alleleB","alleleB_frequency")]
	# names(snptest)[names(snptest) %in% c("alleleA","alleleB","alleleB_frequency")]<-c("Other.allele")
	csi_snptest<-merge(snptest,csi2,by.y="SNP",by.x="rsid")	
	csi_snptest<-harmonise_effect_allele(Dat_harmonise=csi_snptest,ea.x="alleleB",ea.y="EFFECT_ALLELE",oa.x="alleleA",oa.y="OTHER_ALLELE",eaf.x="alleleB_frequency",eaf.y="EAF",b.y="BETA")
	# now effect allele for csi beta is same as coded allele in UKB 
	# Pal<-paste0(Test$EFFECT_ALLELE,Test$OTHER_ALLELE)	
	# Test[which(Pal %in% c("GC","CG","TA","AT")),c("EAF","alleleB_frequency")]
	Dat$GS_csi<-create_gs(gs_dat=csi_snptest,BETA="BETA")		
	return(Dat)
}


cox_gtex_genetic_score<-function(){
	gtex<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_gtex_sig_clumped.txt",sep=" ",head=T,stringsAsFactors=F)
	gtex<-gtex[gtex$clump_r2==0.3,]	
	gtex$gene<-NA
	gtex$gene[gtex$gene_id=="ENSG00000095303.14"]<-"COX1"
	gtex$gene[gtex$gene_id=="ENSG00000073756.11"]<-"COX2"
	# cox1<-gtex[gtex$gene == "COX1",]
	# cox2<-gtex[gtex$gene == "COX2",]
	snptest<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/ukb_snp_stats_all.txt",sep=" ",head=T,stringsAsFactors=F)
	gtex_snptest<-merge(snptest,gtex,by="rsid")
	# cox2_snptest<-merge(snptest,cox2,by="rsid")

	# exp_dat$maf = ifelse(exp_dat$ref_factor == -1, 1 - exp_dat$maf, exp_dat$maf)
    # ref_factor:               '1', when the minor allele is the alt base, '-1' when the minor allele is the reference base

	gtex_snptest<-harmonise_effect_allele(Dat_harmonise=gtex_snptest,ea.x="alleleB",ea.y="effect_allele",oa.x="alleleA",oa.y="other_allele",eaf.x="alleleB_frequency",eaf.y=NULL,b.y="slope")
	# cox2_snptest<-harmonise_effect_allele(Dat=cox2_snptest,ea.x="alleleB",ea.y="effect_allele",oa.x="alleleA",oa.y="other_allele",eaf.x="alleleB_frequency",eaf.y=NULL,b.y="slope")
	# alleleB and effect allele already the same prior to harmonisation. alleleA always equal to other_allele. Assume on same strand and include palindromic SNPs
	# Pal<-paste0(cox1_snptest$effect_allele,cox1_snptest$other_allele)
	# length(which(!Pal %in% c("GC","CG","AT","TA")))
	
	IDS<-unique(gtex_snptest$ID)
	String<-unlist(strsplit(IDS,split=".allpairs.txt.gz_"))
	IDS_nice<-String[seq(1,length(String),by=2)]
	Cox<-String[seq(2,length(String),by=2)]
	IDS_nice[grep("COX1",Cox)]<-paste0("COX1_",IDS_nice[grep("COX1",Cox)])
	IDS_nice[grep("COX2",Cox)]<-paste0("COX2_",IDS_nice[grep("COX2",Cox)])

	for(i in 1:length(IDS)){
		print(IDS[i])
		Dat[,paste0("GS_",IDS_nice[i])]<-NA
		Dat[,paste0("GS_",IDS_nice[i])]<-create_gs(gs_dat=gtex_snptest[gtex_snptest$ID == IDS[i],],BETA="slope")		
	}

	IDS_cox1<-IDS[grep("COX1",IDS)]
	IDS_cox2<-IDS[grep("COX2",IDS)]
	Dat$GS_gtex_cox1<-create_gs_unweighted(gs_dat=gtex_snptest[gtex_snptest$ID %in% IDS_cox1,],BETA="slope")		
	Dat$GS_gtex_cox2<-create_gs_unweighted(gs_dat=gtex_snptest[gtex_snptest$ID %in% IDS_cox2,],BETA="slope")
	return(Dat)
}


# mylogit <- glm(cancer ~ d5d*smoking2 +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat3, family = "binomial")


cox_eqtlgen_genetic_score<-function(){
	eqtl<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/Cox_eqtlgen_sig_clumped.txt",sep=" ",head=T,stringsAsFactors=F)
	eqtl<-eqtl[eqtl$beta!="Inf",]
	eqtl<-eqtl[eqtl$clump_r2==0.3,]	
	eqtl$gene<-toupper(eqtl$id)
	snptest<-read.table("ukb_snp_stats_all.txt",sep=" ",head=T,stringsAsFactors=F)
	eqtl_snptest<-merge(snptest,eqtl,by.x="rsid",by.y="snp")

	eqtl_snptest<-harmonise_effect_allele(Dat_harmonise=eqtl_snptest,ea.x="alleleB",ea.y="effect_allele",oa.x="alleleA",oa.y="other_allele",eaf.x="alleleB_frequency",eaf.y="effect_allele_freq",b.y="beta")
	
	IDS<-unique(eqtl_snptest$gene)
	for(i in 1:length(IDS)){
		print(IDS[i])
		Dat[,paste0("GS_eqtlgen_",IDS[i])]<-NA
		Dat[,paste0("GS_eqtlgen_",IDS[i])]<-create_gs(gs_dat=eqtl_snptest[eqtl_snptest$gene == IDS[i],],BETA="beta")		
	}
	return(Dat)
}


create_gs<-function(gs_dat=NULL,BETA=NULL){
	Dat2<-Dat
	# Dat3<-Dat
	# gs_dat<-gs_dat[!duplicated(gs_dat$rsid),]
	snps<-gs_dat$rsid
	# j<-1
	if(any(duplicated(snps))) stop("duplicate SNPs present")
	for(j in 1:length(snps)){
		print(j)
		# Beta<-gs_dat[,BETA][gs_dat$rsid == snps[j]]
		# Dat2[1,snps[4]]*Beta
		# (2-Dat3[1,snps[4]])*Beta*-1
		# sum(Dat2[,snps[j]][1:10]*Beta)
		# sum((2-Dat3[,snps[j]][1:10])*Beta*-1)
		# sum(Dat2[,snps[j]][1:10])
		Beta<-gs_dat[,BETA][gs_dat$rsid == snps[j]]
		Test<-sign(Beta)
		if(Test<0){
			Beta<-Beta*-1
			# Dat2[,snps[j]]<-Dat2[,snps[j]]*gs_dat[,BETA][gs_dat$rsid == snps[j]]
			Dat2[,snps[j]]<-(2-Dat2[,snps[j]])*Beta		
			# Beta<-gs_dat[,BETA][gs_dat$rsid == snps[j]]*-1
			# Dat3[,snps[j]]<-(2-Dat3[,snps[j]])*Beta	
		}
		if(Test>0){
			Dat2[,snps[j]]<-Dat2[,snps[j]]*Beta		
		}

		# check allele frequencies same in snptest and data files
		# sum(Dat2[,snps[4]])/(nrow(Dat2)*2)
		# snptest[snptest$rsid == snps[4],"alleleB_frequency"]
		# print(sum(Dat2[,snps[j]]))
		# head(Dat2)
	}
	if(length(snps)>1){
		GS<-rowSums(Dat2[,names(Dat2) %in% snps])			
		# GS2<-rowSums(Dat3[,names(Dat3) %in% snps])			
		# plot(GS1,GS)
	}
	if(length(snps)==1){
		GS<-Dat2[,names(Dat2) %in% snps]
	}
	return(GS)	
}

create_gs_unweighted<-function(gs_dat=NULL,BETA=NULL){
	Dat2<-Dat
	gs_dat<-gs_dat[!duplicated(gs_dat$rsid),]
	snps<-gs_dat$rsid
	# j<-4
	if(any(duplicated(snps))) stop("duplicate SNPs present")
	for(j in 1:length(snps)){
		print(j)
		Test<-sign(gs_dat[,BETA][gs_dat$rsid == snps[j]])
		# make effect alleles reflect COX increasing allele
		if(Test<0){
			Dat2[,snps[j]]<-2-Dat2[,snps[j]]	
		}
		# Dat2[,snps[j]]<-Dat2[,snps[j]]*sign(gs_dat[,BETA][gs_dat$rsid == snps[j]])
		# print(sum(Dat2[,snps[j]]))
		# head(Dat2)
	}
	if(length(snps)>1){
		GS<-rowSums(Dat2[,names(Dat2) %in% snps]) #number of COX increasing alleles			
	}
	if(length(snps)==1){
		GS<-Dat2[,names(Dat2) %in% snps]
	}
	return(GS)	
}


flip_strand<-function(Dat=NULL){
	Pos<-Dat$Effect.Allele.x!=Dat$Effect.Allele.y	
	strand1<-c("A","T","G","C")
	strand2<-c("T","A","C","G")
	# lnor.y<-Dat$lnor.y[Pos]*-1
	# Dat$lnor.y[Pos]<-lnor.y
	ea<-Dat$Effect.Allele.y[Pos]
	oa<-Dat$Other.Allele[Pos]
	Dat$Effect.Allele.y[Pos]<-strand2[match(ea,strand1)]
	Dat$Other.Allele[Pos]<-strand2[match(oa,strand1)]				
	return(Dat)
}
		
harmonise_effect_allele<-function(Dat_harmonise=NULL,ea.x=NULL,ea.y=NULL,oa.x=NULL,oa.y=NULL,eaf.x=NULL,eaf.y=NULL,b.y=NULL,drop_palindromics=FALSE){
	Pos<-Dat_harmonise[,ea.x]!=Dat_harmonise[,ea.y]
	if(sum(Pos)==0){
		all(Dat_harmonise[,ea.x]==Dat_harmonise[,ea.y] & Dat_harmonise[,oa.x]==Dat_harmonise[,oa.y])
		print("effect and other alleles in datasets x and y are already the same")
		return(Dat_harmonise)
	}
	if(sum(Pos)!=0){
	# Dat_harmonise[,c(ea.x,ea.y,oa.x,oa.y)]
		# all(Dat_harmonise[,ea.x]!=Dat_harmonise[,ea.y] & Dat_harmonise[,oa.x]!=Dat_harmonise[,oa.y])
		# print("effect and other alleles in datasets x and y are already the same")
		beta.y<-Dat_harmonise[,b.y][Pos]*-1
		Dat_harmonise[,b.y][Pos]<-beta.y
		oa<-Dat_harmonise[,ea.y][Pos]
		ea<-Dat_harmonise[,oa.y][Pos]
		Dat_harmonise[,ea.y][Pos]<-ea
		Dat_harmonise[,oa.y][Pos]<-oa
		eaf<-1-Dat_harmonise[,eaf.y][Pos]
		Dat_harmonise[,eaf.y][Pos]<-eaf	
		Test<-all(Dat_harmonise[,ea.x]==Dat_harmonise[,ea.y] & Dat_harmonise[,oa.x]==Dat_harmonise[,oa.y])
		if(Test) print("effect and other alleles in datasets x and y are now the same")
		if(!Test) print("attempted harmonisation failed probably because the studies are on different strands")
		return(Dat_harmonise)
	}
}




# Pos<-unlist(lapply(c("x","y"),FUN=function(i) grep(i,names(Gen))))
# names(Gen)[Pos]
# all(Gen$rs115693689.x == Gen$rs115693689.y)
# SNPS<-c("rs689470","rs115693689", "rs62575596","rs3025383","rs2007153","rs143127868","rs141607049", "rs117832812","rs28661610","rs139648962", "rs45568238")
# length(SNPS)
# # i<-"rs117832812"
# for(i in SNPS){
# 	print(i)
# 	print(all(Gen[,paste0(i,".x")]==Gen[,paste0(i,".y")]))
# }


# Temp<-readLines("rs174546_snps.txt")
# Temp2<-strsplit(Temp,split=" ")
# Temp3<-data.frame(do.call(rbind,lapply(1:length(Temp2),FUN=function(x) unlist(Temp2[x])[1:2])))
# names(Temp3)<-paste(Temp3[1,]) 
# Gen1<-Temp3[2:nrow(Temp3),]

# snps1<-names(Gen1)
# snps2<-names(Gen2)
# snps1<-snps1[grep("rs",snps1)]
# snps2<-snps2[grep("rs",snps2)]


# From layngeal cancer GWAS but not MR results
# "overall_larynx_cancer"

# MR results with P<0.05 &b>0
#  [1] "Basal cell carcinoma"                
#  [2] "Early-onset prostate cancer"         
#  [3] "Cancer (all cause)"                  
#  [4] "Malignant skin cancer"               
#  [5] "Malignant non-melanoma skin cancer"  
#  [6] "Liver cancer"                        
#  [7] "Respiratory and intrathoracic cancer"
#  [8] "Acute lymphoblastic leukaemia"       
#  [9] "Colon cancer"                        
# [10] "Colorectal cancer"                   
# [11] "Distal colorectal cancer"            
# [12] "Proximal colorectal cancer"          
# [13] "Rectal cancer"                       
# [14] "ER- breast cancer"                   
# [15] "Lung cancer"                         
# [16] "Oral cavity and pharyngeal cancer"   
# [17] "Oropharyngeal cancer"                
# [18] "Esophageal squamous cell carcinoma" 

# MR results with P<0.05 &b<0
# [1] "Endometrioid ovarian cancer" 
# "Multiple myeloma"           

# sig_can_function1<-function(bd=NULL){
# 	cancers<-c("overall_headneck_cancer","overall_kidney_cancer","overall_colorectal_cancer","overall_lung_cancer")
# 	Pos_case<-unique(unlist(lapply(cancers, FUN=function(x) which(bd[,x]==1))))
# 	Pos_controls<-which(bd$overall_pan_inclc44_cancer==0)
# 	bd$sig_cancers1<-NA
# 	bd$sig_cancers1[Pos_case]<-1
# 	bd$sig_cancers1[Pos_controls]<-0
# 	return(bd)	
# }

# table(Dat$overall_mult_myel)

sig_can_function<-function(){
	cancers<-c("overall_headneck_cancer","overall_colorectal_cancer","overall_lung_cancer","overall_nm_skin_cancer","overall_larynx_cancer","overall_liver_cell_cancer","overall_acute_lymph_leuk_cancer")
	Pos_case<-unique(unlist(lapply(cancers, FUN=function(x) which(Dat[,x]==2))))
	Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
	Dat$sig_cancers<-NA
	Dat$sig_cancers[Pos_case]<-2
	Dat$sig_cancers[Pos_controls]<-1

	return(Dat)	
}


lc_crc_function<-function(){
	lc_crc_cancers<-c("overall_colorectal_cancer","overall_lung_cancer")
	Pos_case<-unique(unlist(lapply(lc_crc_cancers, FUN=function(x) which(Dat[,x]==2))))
	Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
	Dat$lc_crc_cancers<-NA
	Dat$lc_crc_cancers[Pos_case]<-2
	Dat$lc_crc_cancers[Pos_controls]<-1	
	return(Dat)	
}

lc_crc_nm_cancers_function<-function(){
	# names(Dat)[grep("overall",names(Dat))]
	  # "overall_breast_cancer","overall_prostate_cancer"   
	  #           "overall_endometrial_cancer" 
	  #           "overall_pan_inclc44_cancer"
	  #            "overall_pan_exclc44_cancer" 
	  #            "overall_melanoma"
	  #            "overall_ovarian_cancer"
	  #            "overall_nm_skin_cancer"      
	  #            "overall_myel_leuk_cancer","overall_lymph_leuk_cancer"
	  #            "overall_mult_myel" 
	  #            "overall_leuk_cancer"
	  #            "overall_haem_cancer"
	  #            "overall_brain_cancer"
	lc_crc_nm_cancers<-c("overall_colorectal_cancer","overall_lung_cancer","overall_nm_skin_cancer")
	Pos_case<-unique(unlist(lapply(lc_crc_nm_cancers, FUN=function(x) which(Dat[,x]==2))))
	# Dat$overall_pan_inclc44_cancer
	Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
	Dat$lc_crc_nm_cancers<-NA
	Dat$lc_crc_nm_cancers[Pos_case]<-2
	Dat$lc_crc_nm_cancers[Pos_controls]<-1	
	return(Dat)	
}

smoking_cancers_function<-function(){
	# names(Dat)[grep("overall",names(Dat))]
	  # "overall_breast_cancer","overall_prostate_cancer"   
	  #           "overall_endometrial_cancer" 
	  #           "overall_pan_inclc44_cancer"
	  #            "overall_pan_exclc44_cancer" 
	  #            "overall_melanoma"
	  #            "overall_ovarian_cancer"
	  #            "overall_nm_skin_cancer"      
	  #            "overall_myel_leuk_cancer","overall_lymph_leuk_cancer"
	  #            "overall_mult_myel" 
	  #            "overall_leuk_cancer"
	  #            "overall_haem_cancer"
	  #            "overall_brain_cancer"
	smoking_cancers<-c("overall_cervical_cancer", "overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer")

	smoking_cancers[!smoking_cancers %in% names(Dat)]
	Pos_case<-unique(unlist(lapply(smoking_cancers, FUN=function(x) which(Dat[,x]==2))))
	Dat$overall_pan_inclc44_cancer
	Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
	Dat$smoking_cancers<-NA
	Dat$smoking_cancers[Pos_case]<-2
	Dat$smoking_cancers[Pos_controls]<-1
	Dat$non_smoking_cancers_inclc44<-Dat$overall_pan_inclc44_cancer
	Dat$non_smoking_cancers_inclc44[Pos_case]<-NA
	Dat$non_smoking_cancers_exclc44<-Dat$overall_pan_exclc44_cancer
	Dat$non_smoking_cancers_exclc44[Pos_case]<-NA	
	return(Dat)	
}

