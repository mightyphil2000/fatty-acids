# install.packages("devtools")
# library(devtools)
# librzary(reinstallr)
# reinstallr(path = '~/fatty-acids/colocalisation/scripts/regional_association_plots_functions.R')
library(plyr) 
library(coloc)
library(hyprcoloc)
# devtools::install_github("MRCIEU/TwoSampleMR")
library(gassocplot)
# install_github("jrs95/gassocplot")

regional_association_plot<-function(gene=NULL,study=NULL,gtex_tissues=NULL,fix_charge=FALSE,region=NULL,gene_tab=NULL,ref=NULL,fa.tab=NULL,gtex_data=NULL,eqtlgen_data=NULL,bbj_eqtl_data=NULL,Plot_height=NULL,gwis_file=gwis_file1,crc_data=NULL,lun_data=NULL,cancer_data=NULL,Top.marker=NULL,ld_eas_pop=NULL,fix_log10=TRUE){

	# replace lines 8/79 with format_data3()
	# Plot<-format_data3()
	# for(j in 1:length(gene_tab$gene)){		
	# tissues=c("colon sigmoid","whole blood","lung","adipose subcutaneous","artery coronary","liver","colon transverse","adipose visceral omentum")

	# snp<-gene_tab$SNP[gene_tab$gene == gene]
	
	Plot<-format_data3(crc_data=crc_data,lun_data=lun_data,cancer_data=cancer_data,gtex_data=gtex_data,fa.tab=fa.tab,ref=ref,eqtlgen_data=eqtlgen_data,gene_tab=gene_tab,bbj_eqtl_data=bbj_eqtl_data,ld_eas_pop=ld_eas_pop,gene=gene,region=region,gtex_tissues=gtex_tissues,fix_charge=fix_charge,studies=study,gwis_file=gwis_file,fix_log10=fix_log10)		




# 	if(is.null(ld_eas_pop)){
# 		Region_object_list<-define_region(gene=gene,region=region,ref=ref) #region refers to how many base pairs to plot . region=NULL plots 1million base pairs (the default). 
# 		SNPlist<-data.frame(Region_object_list[3],stringsAsFactors=F)
# 		SNPlist<-SNPlist$SNP
# 	}
	
# 	if(!is.null(ld_eas_pop)){
# 		Region_object_list<-define_region_eas(gene=gene,region=region,ref=ref,ld_eas_pop=ld_eas_pop)
# 		SNPlist<-unlist(Region_object_list[3])		
# 	}

# 	ref3<-data.frame(Region_object_list[1],stringsAsFactors=F)
# 	Gen.table<-data.frame(Region_object_list[2],stringsAsFactors=F)
	
# 	fa.tab2<-fa.tab[fa.tab$study %in% study,]				
# 	# fa.tab2<-format_fastudies(Dat=fa.tab2)

# 	if(all(study == "CHARGE") & fix_charge){
# 		fa.tab_fix<-fa.tab2[fa.tab2$trait == "AA:DGLA",]
	
# 		# fa.tab_fix<-fa.tab2[fa.tab2$study == "CHARGE" & fa.tab2$trait == "AA:DGLA",]
# 		fa.tab_fix$trait<-"DPAn6:AA_deleteme"
# 		fa.tab2<-rbind(fa.tab2,fa.tab_fix)
# 	}			

# 	tissues_keep<-tissues[tissues %in% gtex_tissues ]
	
# 	Traits1<-format_data1(Fa.tab=fa.tab2,ref.dat=ref3,gwis_test=TRUE,gwis_file=gwis_file,gene=gene,fix_log10=fix_log10)
# 	# study="CHARGE"
# 	# unique(Traits1$marker)
# 	# dups<-unique(fa.tab2$SNP[duplicated(fa.tab2$SNP)])
# # fa.tab2[fa.tab2$SNP == dups[1],]
		
# 	if(!is.null(gtex_data)){
# 		Traits2<-format_data1(Fa.tab=gtex_data,ref.dat=ref3,gene=gene,gene_tab=gene_tab,tissues=tissues_keep,gtex_test=TRUE,fix_log10=fix_log10)
# 	}

# 	if(!is.null(eqtlgen_data)){
# 		Traits3<-format_data1(Fa.tab=eqtlgen_data,ref.dat=ref3,gene=gene,gene_tab=gene_tab,eqtlgen_test=TRUE,fix_log10=fix_log10)
# 	}
		
# 	if(!is.null(bbj_eqtl_data)){
# 		Traits4<-format_data1(Fa.tab=bbj_eqtl_data,ref.dat=ref3,gene=gene,gene_tab=gene_tab,bbj_eqtl_test=TRUE,fix_log10=fix_log10)
# 	}

# 	if(!is.null(crc_data)){
# 		Traits5<-crc_data 		
# 		# format_crc(Dat=crc_data,ref.dat=ref3)
# 	}

# 	if(!is.null(lun_data)){
# 		Traits6<-lun_data
# 		# format_lun(Dat=lun_data,ref.dat=ref3)
# 	}

# 	# Traits7<-format_data1(Fa.tab=fa.tab2,ref.dat=ref3,gwis_test=TRUE,gwis_file=gwis_file,study="Framingham",gene=gene)
	
# 	trait_list<-ls()[grep("Traits[0-9]",ls())] 
	
# 	Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
# 		!is.null(nrow(eval(parse(text=trait_list[x]))))))
# 	rm(list=trait_list[!Test]) #remove the objects with no data
# 	trait_list<-ls()[grep("Traits[0-9]",ls())] 

# 	trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))		

# 	Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,study=study)
	ld.matrix<-data.frame(Plot[1])			
	Z.matrix<-data.frame(Plot[2])
	Markers<-data.frame(Plot[3])

	# write.table(Markers$marker,"~/fatty-acids/colocalisation/data/snplist_eas.txt",col.names=F,row.names=F,quote=F)
	Trait_names<-unlist(Plot[4])	
	Trait_names_simple<-unlist(Plot[names(Plot) == "Traits2"])
	#fix trait names for plots of crc, lc, ea and gene expression in East Asians or Europeans 
	# Trait_names<-gsub("Blood1","Blood",Trait_names)
	# Trait_names<-gsub("Blood2","Blood",Trait_names)
	# Trait_names<-gsub("Blood3","Blood",Trait_names)
	# Trait_names<-trimws(unlist(lapply(1:length(Trait_names),FUN=function(x)
	# 	unlist(strsplit(Trait_names,split="\n")[x])[1]
	# 	)))

	Markers$Markers.pos<-as.integer(Markers$Markers.pos)
	basepairs<-max(Markers$Markers.pos)-min(Markers$Markers.pos)
	# Z.matrix<-Z.matrix/2 #if the Y axis is flooring the P value can fix by dividing the Z score by 2

	extra_info<-""
	if(fix_charge) extra_info<-"fix_charge"

	Trait_title<-unlist(strsplit(Trait_names[1],split=" "))[1]
	Trait_title<-gsub(":","to",Trait_title)

	plot.title<-make_title(Dir="~/fatty-acids/colocalisation/results/plots/",Info=c(study,Trait_title,gene,basepairs,extra_info),type=".png",plot_strategy=NULL,nsnps=nrow(Markers),Trait_names=Trait_names_simple,gwis_file=gwis_file,ld_eas_pop=ld_eas_pop)

	if(is.null(Plot_height)){
		Plot_height<-plot_height(Dat=length(Trait_names))
	}

	names(Markers)<-gsub("Markers.","",names(Markers))

	if(is.null(Top.marker)){
		png(plot.title, width=500,height=Plot_height)
			stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names)
		dev.off()
	}

	
	if(!is.null(Top.marker)){
		plot.title<-gsub(".png",paste0("_",Top.marker,".png"),plot.title)
		png(plot.title, width=500,height=Plot_height)
		# for(i in 1:length(Trait_names)){
			# print(i)
			# plot.title2<-gsub(".png",paste0("_",i,".png"),plot.title)
			# print(plot.title2)
			# print(Trait_names[i])
			# png(plot.title2, width=500,height=400)
				stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names,top.marker=Top.marker)
			dev.off()	
		# }
		
	}
	


	
	# out.file<-gsub(".png","_Nsnps.txt",plot.title)
	# write.table(paste0(gene," genomic region in ",study," (N SNPs =", nrow(Markers),")"),out.file,col.names=F,row.names=F,quote=F)
			
}



load_data<-function(cancer1=NULL,Crc_dat=NULL,cancer2=NULL,Lun_dat=NULL,cancer3=NULL,Cancer_dat=NULL,gtex_file=NULL,eqtlgen_file=NULL,gwis_file=NULL,ref_dat=NULL,ref_turn_off=FALSE,bbj_eqtl_file=NULL,UseBioMart=FALSE){
	
	# cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",#genetic associations for lung cancer from TRICL consortium in 1 mb region centred around the index SNP for each region. This file was created using ~/fatty-acids/scripts/extract_SNPs_cancer.R. 
		
	# cancer1_name="Lung cancer",
	
	# cancer2="~/MR_FattyAcids/data/summary data/colorectal_cancer/061119/crc_snps.Rdata",#genetic associations for colorectal cancer from GECCO/CORECT consortium in 1 mb region centred around the index SNP for each region. This file was created using ~/fatty-acids/scripts/extract_SNPs_cancer.R. 
	
	# cancer2_name="Colorectal cancer",
	
	# gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",#genetic associations for gene expression from gtex version8 in 1 mb region centred around the index SNP for each region in estimated in all GTEx participants, including non-Europeans (84.6% white/European, 12.9% African American). The file can also be found here /projects/MRC-IEU/users/ph14916/gtex2. This file was created using ~/fatty-acids/colosalition/scripts/extract_regions.sh and bash_input_gtex2.txt . build 38
	
	# eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", #genetic associations for gene expression from eqtlgen in a 1mb region centred around the index SNP for each region. This file was created using extract_regions_eqtlgen.sh and extract_regions_eqtlgen2.sh
	
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata", #genetic associations from gwis in a 1mb regino centred around the index SNP for each region. Assumed to contain three data frames called Charge1, Shin1 and Framingham1. The file can also be found here /projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/. This file was created using the script ~/fatty-acids/colocalisation/scripts/extract_snp_data.R 
	
	# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz", #reference data. all SNPs in random 10,000 European participants in UK Biobank and their chromosomal coordinates in GRCh37. I obtained this file from Tom Richardson. The file can also be found here /projects/MRC-IEU/users/ph14916/ukb_bed 
	
	# ref_turn_off=FALSE

	if(UseBioMart){
		library(biomaRt)
		Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
		Attr<-listAttributes(Mart)
	}

	if(!ref_turn_off){
		ref <- read.table(ref_dat,stringsAsFactors=F,head=F)
		if(ref_dat=="~/fatty-acids/colocalisation/data/ukb_fads_ref_hg19.txt") {
			ref<-format_ref(ref=ref)

		}

	}

	if(!is.null(cancer1)){
		load(cancer1)
		Lun_dat<-Lun	
		Lun_dat<-format_lun(Dat=Lun_dat,ref.dat=ref)
	}

	if(!is.null(cancer2)){
		if(sum(grep(".txt",cancer2))==1){
			CRC<-read.table(cancer2,sep="\t",stringsAsFactors=F,head=T)
			Crc_dat<-format_crc_accc(Dat=CRC,ref.dat=ref)
		}else{
			load(cancer2)
			Crc_dat<-format_crc_gecco(Dat=CRC,ref.dat=ref)
			head(Crc_dat)
		}
	}
		
	if(!is.null(cancer3)){
		load(cancer3)
		Dat_name<-unlist(strsplit(cancer3,split="/"))
		Dat_name<-Dat_name[5]
		Dat_name<-unlist(strsplit(Dat_name,split=".RData"))
		cancer_dat<-eval(parse(text=Dat_name))
		Cancer_dat<-format_cancer(Dat=cancer_dat,trait="outcome",marker="rsid",eaf="eaf",effect_allele="Effect.Allele",other_allele="Other.Allele",beta="lnor",se="se",pval="p",cases="ncase",controls="ncontrol",ref.dat=ref)			

	}
	

	load(gwis_file)
	if(sum(grep("eastasians",gwis_file))==0){
		Charge1$study <- "CHARGE"
		Shin1$study <- "Shin"
		Tin1$study <- "Framingham"
		fa.tab<-do.call(rbind,list(Charge1,Shin1,Tin1))
	}
	if(sum(grep("eastasians",gwis_file))==1){
		ChargeEA1$study <- "ChargeEA"
		Dor4$study <- "Dorajoo/SCHS"
		Dor1<-format_schs(Dat=Dor4)
		fa.tab<-do.call(rbind,list(ChargeEA1,Dor1))	
	}
	
	fa.tab<-format_fastudies(Dat=fa.tab)
	# if(is.null(bbj_eqtl_file)){
	# bbj_eqtl_data<-NULL
	
	gtex_data<-NULL	
	if(!is.null(gtex_file)){
		load(gtex_file)
		gtex_data<-format_gtex(df4)
	}

	eqtlgen_data<-NULL	
	if(!is.null(eqtlgen_file)){
		eqtlgen_data<-read.table(eqtlgen_file,sep="\t",stringsAsFactors=F,head=T)
		eqtlgen_data$tissue<-"blood"
		eqtlgen_data$study<-"eQTLGen"
			
	}
	
	bbj_eqtl_data<-NULL
	if(!is.null(bbj_eqtl_file)){
		bbj_eqtl_data<-load_bbjeqtl(File=bbj_eqtl_file)
	}

	

	gene_tab<-make_gene_table(snps=c("rs174528","rs174528","rs3734398","rs603424","rs6584379"),genes=c("FADS1","FADS2","ELOVL2","SCD","HIF1AN"),Mart=NULL ) #function assumes that these are the SNPs and genes of interest
# =Mart
	return(list(gene_tab,ref,fa.tab,gtex_data,eqtlgen_data,bbj_eqtl_data,gwis_file,Lun_dat,Crc_dat,Cancer_dat))
}





format_schs<-function(Dat=NULL){
			Pos<-which(!is.na(Dat$effect_allele_freq))
			Dat$eaf[Pos]<-Dat$effect_allele_freq[Pos]
			Dat$SNP[Pos]<-Dat$snp[Pos]
			Dat$pval[Pos]<-Dat$p[Pos]
			Dat$samplesize[Pos]<-Dat$n[Pos]
			Dat<-Dat[,!names(Dat) %in% c("snp","effect_allele_freq","p","n")]
			return(Dat)
		}

load_bbjeqtl<-function(File=bbj_eqtl_file){
		load(File)
		bbj_bcells$tissue<-"B cells"
		bbj_blood$tissue<-"Blood"
		bbj_cd4tcells$tissue<-"CD4 T cells"
		bbj_cd8tcells$tissue<-"CD8 T cells"
		bbj_mono$tissue<-"Monocytes"
		bbj_nkcells$tissue<-"NK cells"
		List<-list(bbj_bcells,bbj_blood,bbj_cd4tcells,bbj_cd8tcells,bbj_mono,bbj_nkcells)
		bbj_eqtl<-do.call(rbind,List)		
		bbj_eqtl<-bbj_eqtl[,names(bbj_eqtl) != "file"]
		Ens<-bbj_eqtl$gene
		Ens<-unlist(strsplit(Ens,split="\\."))
		bbj_eqtl$gene<-Ens[seq(1,length(Ens),by=2)]
		bbj_eqtl$study<-"BBJ"
		bbj_eqtl$se<-bbj_eqtl$beta/bbj_eqtl$t.stat
		return(bbj_eqtl)
	}


format_gtex<-function(Dat=NULL){
	Dat$study<-"GTEx"
	Ens<-Dat$gene_id
	Ens<-unlist(strsplit(Ens,split="\\."))
	Dat$Ens<-Ens[seq(1,length(Ens),by=2)]
	return(Dat)
}


make_gene_table<-function(snps=NULL,genes=NULL,Mart=NULL){
	if(!is.null(Mart)){
		gene_tab<-data.frame(matrix(c(snps,genes),ncol=2,nrow=4),stringsAsFactors=F)
		names(gene_tab)<-c("SNP","gene")
		# load("/Users/ph14916/fatty-acids/colocalisation/data/Ens_gene.Rdata")
		Ens_gene<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="external_gene_name",values=gene_tab$gene,mart=Mart) #unique(Fa.tab$gene_id)
		# save(Ens_gene,file="/Users/ph14916/fatty-acids/colocalisation/data/Ens_gene.Rdata")
		gene_tab<-merge(gene_tab,Ens_gene,by.x="gene",by.y="external_gene_name")
	}
	# lapply(1:ncol(gene_tab1),FUN=function(x)
	# 	class(gene_tab[,x])
	# )
	# gene_tab1

	if(is.null(Mart)){
		genes<-c("ELOVL2","FADS1","FADS2", 
			"SCD","HIF1AN","SEC31B","NDUFB8", #genes in chromosome 10 SCD region 
			"VPS37C","BSCL2","FADS3","TMEM258", #rs2727271 and rs2524299 are eQTLs for these non FADS1/2 genes in gtex. in eqtl gene also eqtl for MYRF
			"DAGLA","FEN1","RAB3IL1","BEST1","FTH1","SYT7","MYRF") #other genes in 500kb window around FADS1/2
		# "MYRE"
		snps<-c("rs3734398","rs174528","rs174528","rs603424","rs6584379","rs6584379","rs6584379","rs174528","rs174528","rs174528","rs174528",
			"rs174528","rs174528","rs174528","rs174528","rs174528","rs174528","rs174528")
		ens<-c("ENSG00000197977" ,"ENSG00000149485","ENSG00000134824","ENSG00000099194","ENSG00000166135","ENSG00000075826","ENSG00000166136", #ENSG for NDUFB8
			"ENSG00000167987","ENSG00000168000","ENSG00000221968","ENSG00000134825",#ENSG ids for VPS37C, BSCL2, FADS3 and TMEM258
			"ENSG00000134780","ENSG00000168496","ENSG00000167994","ENSG00000167995","ENSG00000167996","ENSG00000011347","ENSG00000124920" #ENSG IDs for DAGLA, FEN1, RAB3IL1, BEST1, FTH1, SYT7 and MYRF
			)

		gene_tab<-data.frame(do.call(cbind,list(genes,snps,ens)),stringsAsFactors=F)
		names(gene_tab)<-c("gene","SNP","ensembl_gene_id")
	}
	return(gene_tab)
}


define_region<-function(gene=NULL,region=NULL,ref=NULL){
	# SNPlist<-readLines(paste("~/fatty-acids/colocalisation/data/",index_snp,".txt",sep=""))
	
	# SNPlist<-readLines(paste("~/fatty-acids/colocalisation/data/",index_snp,"_r_matrix_ukb.frq",sep=""))
	# Gen.table<-read.table(paste("~/fatty-acids/colocalisation/data/",index_snp,"_r_matrix_ukb.ld.gz",sep=""),sep="\t",head=F,stringsAsFactors=F)		

	genes2<-c("ELOVL2","FADS1","FADS2","SCD")
	genes3<-c("elovl2","fads","fads","scd")
	genes_df<-data.frame(cbind(genes2,genes3),stringsAsFactors=F)
	index_gene<-genes_df$genes3[genes_df$genes2==gene]

	if(gene %in% c("HIF1AN","SEC31B","NDUFB8")) {
		index_gene<-"scd"
	}

	if(gene %in% c("VPS37C","BSCL2","FADS3","TMEM258",
		"DAGLA","FEN1","RAB3IL1","BEST1","FTH1","SYT7","MYRF")){
		index_gene<-"fads"
	}

	# index_gene<-"fads"
	Gen.table<-read.table(paste0("~/fatty-acids/colocalisation/data/",index_gene,"_r_matrix_ukb.ld.gz"),sep="\t",head=F,stringsAsFactors=F)		
	SNPlist<-read.table(paste0("~/fatty-acids/colocalisation/data/",index_gene,"_r_matrix_ukb.frq"),head=T,stringsAsFactors=F)

	row.names(Gen.table)<-SNPlist$SNP
	colnames(Gen.table)<-SNPlist$SNP	
	# Pos<-which(row.names(Gen.table) %in% c("rs968567","rs174546"))
	# Gen.table[Pos,Pos]
	ref3<-ref[ref$V2 %in% SNPlist$SNP,]
	if(!is.null(region)){
		region2<-region/2
		Max<-max(ref3$V4)
		Max<-Max-region2
		Min<-min(ref3$V4)
		Min<-Min+region2
		ref4<-ref3[ref3$V4 > Min & ref3$V4<Max,]
		SNPlist<-SNPlist[SNPlist$SNP %in% ref4$V2,]
		Gen.table<-Gen.table[row.names(Gen.table) %in% SNPlist$SNP,colnames(Gen.table) %in% SNPlist$SNP]
		ref3<-ref4		
		# max(ref3$V4)-min(ref3$V4)
	}

	return(list(ref3,Gen.table,SNPlist))
}



define_region_eas<-function(gene=NULL,region=NULL,ref=NULL,ld_eas_pop=ld_eas_pop){
	# SNPlist<-readLines(paste("~/fatty-acids-mr/coloc/",index_snp,".txt",sep=""))
	# Gen.table<-read.table(paste("~/fatty-acids-mr/coloc/",index_snp,"_r_matrix.ld",sep=""),sep="\t",head=F,stringsAsFactors=F)	
	
	gene_eas_tab<-data.frame(cbind(c("FADS1","FADS2","ELOVL2","SCD"),c("fads","fads","elovl2","scd")),stringsAsFactors=F)
	gene_eas<-gene_eas_tab$X2[gene_eas_tab$X1==gene]
	if(gene %in% c("HIF1AN","SEC31B","NDUFB8")) {
		gene_eas<-"SCD"
	}
	
	if(is.null(ld_eas_pop)) stop("ld_eas_pop is NULL")
	Gen.table<-read.table(paste("~/fatty-acids/colocalisation/data/r2_matrix_",ld_eas_pop,"_",gene_eas,".txt",sep=""),sep="\t",head=T,stringsAsFactors=F)
	
	row.names(Gen.table)<-Gen.table$RS_number
	Gen.table<-Gen.table[,names(Gen.table) != "RS_number"   ]
	Gen.table<-sqrt(Gen.table) #regional association plot input assumes correlation is an unsquared correlation coefficient. For east asian LD files downloaded from ldlink the correlation is an R2

	SNPlist<-names(Gen.table)
	ref3<-ref[ref$V2 %in% SNPlist,]
	
	if(!is.null(region)){
		region2<-region/2
		Max<-max(ref3$V4)
		Max<-Max-region2
		Min<-min(ref3$V4)
		Min<-Min+region2
		ref4<-ref3[ref3$V4 > Min & ref3$V4<Max,]
		SNPlist<-SNPlist[SNPlist %in% ref4$V2]
		Gen.table<-Gen.table[row.names(Gen.table) %in% SNPlist,colnames(Gen.table) %in% SNPlist]
		ref3<-ref4
	}
	return(list(ref3,Gen.table,SNPlist))
}

format_fastudies<-function(Dat=NULL){
	File<-c("AA_to_DGLA_adjSNP.tab", "GLA_to_LA_adjSNP.tab","AA_to_DGLA.tab","ADA_to_AA.tab","DGLA_to_GLA.tab","DHA_to_DPA_n3.tab","DPA_n3_to_EPA.tab","GLA_to_LA.tab","OA_to_SA.tab","POA_to_PA.tab","DPA_n6_to_AA.tab","score_lnD5D_pooled_allchr_qc1.tab", "score_lnD6D_pooled_allchr_qc1.tab" ,
		"DHA_to_EPA.tab")
	Trait<-c("AA:DGLAadj","GLA:LAadj","AA:DGLA","ADA:AA","DGLA:GLA","DHA:DPAn3","DPAn3:EPA","GLA:LA","OA:SA","POA:PA","DPAn6:AA", "AA:DGLA lnD5Dpooled" ,"GLA:LA lnD6Dpooled","DHA:EPA")

	Trait_file<-data.frame(matrix(c(File,Trait),ncol=2,nrow=length(File)),stringsAsFactors=F)
	names(Trait_file)<-c("file","trait")
	
	schs_file<-c("_case.tab","_control.tab","_pooled.tab")
	schs_trait<-c(" cases"," controls", " pooled")
	
	File1<-unlist(lapply(1:length(schs_file),FUN=function(x) gsub(".tab",schs_file[x],Trait_file$file)))
	Trait1<-unlist(lapply(1:length(schs_trait),FUN=function(x) paste0(Trait_file$trait,schs_trait[x])))
	Trait_file1<-data.frame(cbind(File1,Trait1),stringsAsFactors=F)
	names(	Trait_file1)<-c("file","trait")
	Trait_file<-rbind(Trait_file,Trait_file1)

	# dim(unique(Trait_file))
	Dat$file<-gsub(".gz","",Dat$file)
	Dat2<-merge(Dat,Trait_file,by="file")
	maf<-Dat2$eaf
	maf[which(maf>0.5)]<-1-maf[which(maf>0.5)]
	Dat2$maf<-maf	
	return(Dat2)
}

format_data1<-function(Fa.tab=NULL,ref.dat=NULL,fatty_acids=FALSE,gene="",gene_tab=NULL,tissues=NULL,gtex_test=FALSE,eqtlgen_test=FALSE,gwis_test=FALSE,bbj_eqtl_test=FALSE,gwis_file=NULL,fix_log10=NULL){

	if(fatty_acids){ #this section was writen for the non_GWIS fatty acid traits  
		Files<-unique(Fa.tab$file)
		Files_split<-strsplit(Files,split="/")
		K<-length(unlist(Files_split[1]))
		Files_split<-unlist(Files_split)
		Files_split2<-Files_split[seq(K,length(Files_split),by=K)]
		L_dat<-NULL
		for(i in 1:length(Files)){
			print(i)
			Dat<-Fa.tab[Fa.tab$file == Files[i],]
			Dat$filename<-Files_split2[i]
			L_dat[[i]]<-Dat
		}
		Fa.tab<-do.call(rbind,L_dat)
		if(any(names(Fa.tab) == "info")){
			Fa.tab<-Fa.tab[Fa.tab$info > 0.8, ]
		}

		Key<-read.table("~/fatty-acids/data/filekeyv2.txt",sep="\t",stringsAsFactors=F,head=T)
		Key<-Key[,c("pmid","author","consortium","trait","chain","chain.length","SigSNPs","sample_size.analysis","population","sex","year","filename","imputed","id.analysis","omega","unsaturation","saturation","type")]
		Key$filename<-gsub("\xca","",Key$filename)
		Key$filename<-gsub("tbl.gz","tbl",Key$filename)
		Key$filename<-trimws(Key$filename)
		Key$trait<-gsub("eicosenoic acid \\(20:1)", "eicosenoic acid (20:1n9)",Key$trait)
		Key<-Key[!is.na(Key$filename),]
		Key$filename2<-gsub(".tab","_elolv5.txt",Key$filename)
		# Files[!Files %in% Key$filename2]
		
		Fa.tab.m<-merge(Fa.tab,Key,by="filename")
	}

	if(gtex_test){
		Fa.tab.m<-merge(Fa.tab,gene_tab[,c("gene","ensembl_gene_id")],by.x="Ens",by.y="ensembl_gene_id")
		Fa.tab.m<-Fa.tab.m[Fa.tab.m$gene == gene,]
		Fa.tab.m<-Fa.tab.m[Fa.tab.m$tissue %in% tissues,]				
		# }
		if(nrow(Fa.tab.m)>0){
			Fa.tab.m<-Fa.tab.m[Fa.tab.m$ma_count>4,]
			Eqtl<-paste(Fa.tab.m$gene," expression in ",Fa.tab.m$tissue,sep="")
			Eqtl2<-paste(Fa.tab.m$gene," expression in ",Fa.tab.m$tissue," (GTEx) mac(median=",median(Fa.tab.m$ma_count),", min=",min(Fa.tab.m$ma_count),", max=",max(Fa.tab.m$ma_count),")",sep="")			
			Fa.tab.m$trait<-Eqtl
			Fa.tab.m$trait2<-Eqtl2
			Chr<-unlist(strsplit(Fa.tab.m$variant_id,split="_"))
			Chr<-Chr[seq(1,length(Chr),by=5)]
			Chr<-unlist(strsplit(Chr,split="chr"))
			Fa.tab.m$chr<-Chr[seq(2,length(Chr),by=2)]
			Pos<-unlist(strsplit(Fa.tab.m$variant_id,split="_"))
			Fa.tab.m$pos<-Pos[seq(2,length(Pos),by=5)]
			Fa.tab.m$pos<-as.numeric(Fa.tab.m$pos)
			names(Fa.tab.m)[names(Fa.tab.m) == "slope"]<-"beta" 
			names(Fa.tab.m)[names(Fa.tab.m) == "slope_se"]<-"se" 
			
		}
	}

	if(eqtlgen_test | bbj_eqtl_test){		
		Fa.tab$Ens<-Fa.tab$gene
		Fa.tab<-Fa.tab[,names(Fa.tab) != "gene"]
		Fa.tab.m<-merge(Fa.tab,gene_tab[,c("gene","ensembl_gene_id")],by.x="Ens",by.y="ensembl_gene_id")	
		Fa.tab.m<-Fa.tab.m[Fa.tab.m$gene ==gene,]
		
		# Fa.tab.m<-merge(Fa.tab,Test,by.x="Ens",by.y="ensembl_gene_id")
		if(nrow(Fa.tab.m)>0){			
			Fa.tab.m$trait<-paste0(Fa.tab.m$gene," expression in ",Fa.tab.m$tissue)
			Fa.tab.m$trait2<-paste0(Fa.tab.m$gene," expression in ",Fa.tab.m$tissue," (",Fa.tab.m$study,")")		
		}
		
	}
	# else{
	if(any(names(Fa.tab) %in% c("filename","file"))){
		names(Fa.tab)[names(Fa.tab) == "file"]<-"filename"
		Files<-unique(Fa.tab$filename)
		if(sum(grep("/",Files))!=0){
			Files_split<-strsplit(Files,split="/")
			K<-length(unlist(Files_split[1]))
			Files_split<-unlist(Files_split)
			Files_split2<-Files_split[seq(K,length(Files_split),by=K)]
			L_dat<-NULL
			for(i in 1:length(Files)){
				print(i)
				Dat<-Fa.tab[Fa.tab$file == Files[i],]
				Dat$filename<-Files_split2[i]
				L_dat[[i]]<-Dat
			}
			Fa.tab<-do.call(rbind,L_dat)
		}
		Test1<-sum(grep("TRICL",Fa.tab$filename[1])==1)!=0
		if(sum(grep("TRICL",Fa.tab$filename[1])==1)!=0){
			Fa.tab$trait<-NA	
			Fa.tab$trait[Fa.tab$filename == "Onco_TRICL_032116_Adeno.csv.tab"]<-"Lung adenocarcinoma"
			Fa.tab$trait[Fa.tab$filename == "Onco_TRICL_032116_Ever.csv.tab"]<-"Lung cancer in ever smokers"
			Fa.tab$trait[Fa.tab$filename ==  "Onco_TRICL_032116_Never.csv.tab"] <-"Lung cancer in never smokers"
			Fa.tab$trait[Fa.tab$filename ==  "Onco_TRICL_032116_Overall.csv.tab"] <-"Lung cancer"
			Fa.tab$trait[Fa.tab$filename ==  "Onco_TRICL_032116_Small.csv.tab"] <-"Small cell lung cancer"
			Fa.tab$trait[Fa.tab$filename ==  "Onco_TRICL_032116_Squam.csv.tab"  ] <-"Squamous cell lung cancer"
			Fa.tab.m<-Fa.tab
		}		
	}

	# This section formats the GWIS ratios
	if(gwis_test){
	# if(!any(c(fatty_acids,Test1,Test2,Test3))){
		Fa.tab.m<-Fa.tab
		# names(Fa.tab.m)[names(Fa.tab.m) == "samplesize"]<-"n"	
		Fa.tab.m$trait[Fa.tab.m$trait=="AA:DGLA"]<-"AA:DGLA / D5D"
		Fa.tab.m$trait[Fa.tab.m$trait=="GLA:LA"]<-"GLA:LA / D6D"
		if(sum(grep("deleteme",unique(Fa.tab.m$trait)))!=0){
			Fa.tab.m$trait[Fa.tab.m$trait == "DPAn6:AA_deleteme"]<- "DPAn6:AA_deleteme / ELOVL2"
		}else{
			Fa.tab.m$trait[Fa.tab.m$trait == "DPAn6:AA"]<- "DPAn6:AA / ELOVL2"
		}
		Fa.tab.m$trait[Fa.tab.m$trait=="DHA:DPAn3"]<-"DHA:DPAn3 / ELOVL2"
		Fa.tab.m$trait[Fa.tab.m$trait=="ADA:AA"]<-"ADA:AA / ELOVL2/5"
		Fa.tab.m$trait[Fa.tab.m$trait=="DPAn3:EPA"]<-"DPAn3:EPA / ELOVL2/5"
		Fa.tab.m$trait[Fa.tab.m$trait=="DGLA:GLA"]<-"DGLA:GLA / ELOVL5"
		Fa.tab.m$trait[Fa.tab.m$trait=="POA:PA"]<-"POA:PA / SCD"
		Fa.tab.m$trait[Fa.tab.m$trait=="OA:SA"]<-"OA:SA / SCD"
		Fa.tab.m$trait2<-paste(Fa.tab.m$trait," (",Fa.tab.m$study,")",sep="")
		
		# D5D<-c("AA:DGLA")
		# D6D<-c("GLA:LA")
		# ELOVL2<-c("DPAn6:AA","DHA:DPAn3")
		# ELOVL2or5<-c("ADA:AA","DPAn3:EPA")
		# ELOVL5<-c("DGLA:GLA")
		# SCD<-c("POA:PA","OA:SA")

	}

	names(Fa.tab.m)[names(Fa.tab.m) == "snp"]<-"marker"
	names(Fa.tab.m)[names(Fa.tab.m) == "SNP"]<-"marker"

	Fa.tab.m$z<-Fa.tab.m$beta/Fa.tab.m$se
	
	if(eqtlgen_test & gene == "FADS2" & fix_log10){ #max -log10 P values larger than 1000 and therefore don\t get plotted properly in regional association plots
		Fa.tab.m$z<-Fa.tab.m$z/2
		Fa.tab.m$trait	<- paste0(Fa.tab.m$trait, " (Zdec)")
		# print(unique(gsub("eQTLGen","eQTLGen Zdec",Fa.tab.m$trait2)))
		# print("gene=")
		# print(gene)
		# print("eqtlgen_test=")
		#  print(eqtlgen_test)

	}
 	
 	imputed<-FALSE
 	Charge<-FALSE
	aa_dgla<-FALSE
 	
 	study<-unique(Fa.tab$study)
 	if(!is.null(study)){
		if(any(study == "CHARGE")){
			Charge<-TRUE	
		}
	}
 	 	
 	if(!is.null(gwis_file)){
		if(sum(grep("notimputed",gwis_file))==0) {
			imputed	<-TRUE	
			# print("gwis_file=")
			# print(gwis_file)
		}
	}
	
	if(sum(grep("AA:DGLA",Fa.tab$trait)) != 0) {
		aa_dgla <-TRUE
		print("trait=")
		print(unique(Fa.tab$trait))
	}
	# print("howzit")
	# 		print(gene)
	# 		print(Charge)
	# 		print(aa_dgla)
	# 		print(imputed)
	# 		print(fix_log10)

	if(gene == "FADS1" & Charge & aa_dgla & imputed & fix_log10){ #max -log10 P values larger than 1000 and therefore don\t get plotted properly in regional association plots
			
			Fa.tab.m$z[Fa.tab.m$study == "CHARGE"]<-Fa.tab.m$z[Fa.tab.m$study == "CHARGE"]/1.25
			# Pos<-which(rnorm(n=1000)>0)
			# Fa.tab.m$z[Pos]<-Fa.tab.m$z[Pos]/1000
			Fa.tab.m$trait2[Fa.tab.m$study == "CHARGE"]	<-gsub("CHARGE","CHARGE Zdec",Fa.tab.m$trait2[Fa.tab.m$study == "CHARGE"]	)
			print(unique(Fa.tab.m$trait2))
	}
					
	Fa.tab.m<-Fa.tab.m[!is.na(Fa.tab.m$z),]
	
	
	if(gtex_test){
		if(nrow(Fa.tab.m)>0){
			ID.temp<-paste(Fa.tab.m$trait,Fa.tab.m$marker)
			Dups<-Fa.tab.m$marker[duplicated(ID.temp)]
			Dups<-unique(Dups)
			Dup.ids<-unique(ID.temp[duplicated(ID.temp)])
			# for(i in 1:length(Dup.ids)){
			# 	print(i)
			# 	print(Dup.ids[i])
			# 	print(Fa.tab.m[which(ID.temp == Dup.ids[i]),"variant_id"])
			# }
			#drop duplicates, which seem to correspond to CNVs or multiallelic SNPS. The duplicates get created because I merged a file with rsids and position with the gtex file on position. The position is duplicated in gtex when there is more than 1 allelic form for a position. There were 3 duplicated rsids / positions on this basis. below are 3 illustrative examples

			# [1] "FADS1 expression in colon sigmoid rs369461948"
			# [1] "chr11_61499277_C_CGTT_b38" "chr11_61499277_C_T_b38"  
			# [1] "FADS1 expression in liver rs34057719"
			# [1] "chr11_61499277_C_CGTT_b38" "chr11_61499277_C_T_b38"   
			# [1] "FADS1 expression in colon sigmoid rs573394124"
			# [1] "chr11_61481705_T_TG_b38" "chr11_61481705_T_TC_b38"
	
			Fa.tab.m2<-Fa.tab.m[!Fa.tab.m$marker %in% Dups,]			
		}
	}else{
		Fa.tab.m2<-merge(Fa.tab.m[,!names(Fa.tab.m) %in% c("chr","bp","pos")],ref.dat,by.x="marker",by.y="V2")
	}

	
	if(nrow(Fa.tab.m2)>0){
		names(Fa.tab.m2)[names(Fa.tab.m2) == "V1"]<-"chr"
		names(Fa.tab.m2)[names(Fa.tab.m2) == "V4"]<-"pos"
		Fa.tab.m2<-Fa.tab.m2[order(Fa.tab.m2$chr,Fa.tab.m2$pos),]
		Fa.tab.m2$marker<-as.character(Fa.tab.m2$marker)
		# if(any(duplicated(Fa.tab.m2$marker))) stop("eish duplicated markers present!")
		# print(unique(Fa.tab.m2$trait2))
		return(Fa.tab.m2)				
	}
	if(nrow(Fa.tab.m)==0){
		return("no data")
	}
}

format_data2<-function(trait_list=NULL,ld.matrix=NULL,snps=NULL,study=NULL){
		Temp<-do.call(rbind.fill,trait_list)		
		N<-length(unique(Temp$trait2))
		Temp<-Temp[which(Temp$z!="NaN"),]
		if(any(is.na(Temp$z))) stop("warning Z scores missing")
		Table_traits<-table(Temp$marker)
		Names_table<-names(Table_traits) 
		SNPs_keep<-Names_table[Table_traits==N]
		Temp<-Temp[Temp$marker %in% SNPs_keep,]
		
		if("ma_count" %in% names(Temp)){
			Min<-unlist(lapply(unique(Temp$trait),FUN=function(x) 
				min(Temp$ma_count[Temp$trait==x])))
			Max<-unlist(lapply(unique(Temp$trait),FUN=function(x) 
				max(Temp$ma_count[Temp$trait==x])))
			Med<-unlist(lapply(unique(Temp$trait),FUN=function(x) 
				median(Temp$ma_count[Temp$trait==x])))
			ID<-unique(Temp$trait)
			Dat<-data.frame(matrix(c(ID,Min,Med,Max),nrow=length(ID),ncol=4,byrow=FALSE),stringsAsFactors=FALSE)
			names(Dat)<-c("trait","Min","Med","Max")
			Dat$trait3<-paste(Dat$trait,paste0("in GTEx \n(","median MAC=",round(as.numeric(Dat$Med),0),", min=",Dat$Min,", max=",Dat$Max,")"))
			Temp<-merge(Temp,Dat,by="trait")
			Temp$trait2[!is.na(Temp$Min)] <- Temp$trait3[!is.na(Temp$Min)]
		}
		Pos.eqtlgen<-grep("eQTLGen",Temp$trait2)
		if(sum(Pos.eqtlgen)!=0){
			Min<-min(as.numeric(Temp$n[Pos.eqtlgen]))
			Max<-max(as.numeric(Temp$n[Pos.eqtlgen]))
			Med<-median(as.numeric(Temp$n[Pos.eqtlgen]))

			Temp$trait2[Pos.eqtlgen]<-paste(unique(Temp$trait[Pos.eqtlgen]),paste0("in eQTLGen", " \n(Median N=",Med,", min=",Min,", max=",Max,")"))
		}
 		# order data so that fatty acids appear in alphabetical order above the eQTLs
		fa_eqtl_traits<-unique(Temp$trait2)		
		fa_traits<-fa_eqtl_traits[grep(":",fa_eqtl_traits)]
		
		eqtl_traits<-fa_eqtl_traits[grep(":",fa_eqtl_traits,invert=T)]

		cancer_terms<-c("carcinoma","adenocarcinoma","cancer")

		Cancer_rows<-eqtl_traits[unlist(lapply(cancer_terms,FUN=function(x)
			grep(x,eqtl_traits)))]

		eqtl_traits<-eqtl_traits[!eqtl_traits %in% Cancer_rows]
		# eqtl_traits[grep("carcinoma",eqtl_traits,invert=T)]
		# eqtl_traits<-eqtl_traits[grep("cancer",eqtl_traits,invert=T)]

		fa_traits<-fa_traits[order(fa_traits)]

		# if(plot_strategy=="1trait_allstudies_alltissues"){
		
		if(length(study) > 1){
			unique(Temp$trait)
			pos<-which(Temp$trait2 %in% fa_traits)			
			# pos<-grep(unique(Temp$trait2),fa_traits)			
			Temp$trait[pos]<-Temp$trait2[pos]
			fa_traits<-unique(Temp$trait[pos])			

			Pos2<-grep("CHARGE",fa_traits)
			Pos3<-grep("CHARGE",fa_traits,invert=T)
			fa_traits<-fa_traits[c(Pos2,Pos3)]
		}

		eqtl_traits<-eqtl_traits[order(eqtl_traits)]
		eqtl_traits1<-eqtl_traits[grep("expression in blood",eqtl_traits,invert=TRUE,ignore.case=T)]
		eqtl_traits2<-eqtl_traits[grep("expression in blood",eqtl_traits,ignore.case=T)]
		eqtl_traits<-c(eqtl_traits1,eqtl_traits2)
		# sort(fa_traits)
		fa_eqtl_traits2<-c(fa_traits,eqtl_traits)
		
		Num<-1:length(fa_eqtl_traits2)
		dat_temp<-data.frame(matrix(c(fa_eqtl_traits2,Num),nrow=length(Num),ncol=2),stringsAsFactors=F)
		dat_temp$X2<-as.numeric(dat_temp$X2)
		dat_temp$X2[grep("AA:DGLA",dat_temp$X1)]<-1
		dat_temp$X2[grep("GLA:LA",dat_temp$X1)]<-2
		Test<-sum(grep("deleteme",dat_temp$X1))!=0
		if(Test){
			dat_temp$X2[dat_temp$X1 %in% c("DPAn6:AA_deleteme / ELOVL2","DHA:DPAn3 / ELOVL2")]<-c(3,4)
		}

		
		dat_temp$X2[grep("ADA:AA",dat_temp$X1)]<-5
		dat_temp$X2[grep("DPAn3:EPA",dat_temp$X1)]<-6
		dat_temp$X2[grep("DGLA:GLA",dat_temp$X1)]<-7
		dat_temp$X2[grep("POA:PA",dat_temp$X1)]<-8
		dat_temp$X2[grep("OA:SA",dat_temp$X1)]<-9
		dat_temp$X3<-dat_temp$X2
		dat_temp$X3[grep("pooled",dat_temp$X1,invert=T)]<-1
		dat_temp$X3[grep("pooled",dat_temp$X1)]<-2
		dat_temp$X3[grep("controls",dat_temp$X1)]<-3
		dat_temp$X3[grep("cases",dat_temp$X1)]<-4

		dat_temp1<-dat_temp[dat_temp$X1 %in% fa_traits,]
		dat_temp2<-dat_temp[!dat_temp$X1 %in% fa_traits,]
		
		dat_temp1<-dat_temp1[order(dat_temp1$X2,dat_temp1$X3),]
		dat_temp<-rbind(dat_temp1,dat_temp2)
		
		# Cancer_names<-c("cancer","carcinoma")
		Pos<-unique(unlist(lapply(cancer_terms,FUN=function(i)
			grep(i,fa_eqtl_traits))))
		if(sum(Pos)!=0){
			cancer_traits<-fa_eqtl_traits[Pos]
			cancer_traits<-cancer_traits[order(cancer_traits)]
			dat_temp<-dat_temp[!dat_temp$X1 %in% fa_eqtl_traits[Pos],]
			can_tab<-data.frame(cbind(cancer_traits,1:length(cancer_traits)),stringsAsFactors=F)
			can_tab$X3<-1:nrow(can_tab)
			names(can_tab)<-c("X1","X2","X3")			
			dat_temp<-rbind(can_tab,dat_temp)			
		}

		# manually set order for plot with crc, lc and ea and most relevant tissue for gene expression in east asians
		# dat_temp0<-dat_temp[dat_temp$X1  ==  "GLA:LA lnD6Dpooled (Dorajoo/SCHS)",]
		# Pos<-which(dat_temp$X1 == "Colorectal cancer (ACCC)")
		# Pos<-c(Pos,grep("Blood1",dat_temp$X1))
		# dat_temp1<-dat_temp[Pos,]
		# Pos<-which(dat_temp$X1 == "Esophageal squamous cell carcinoma BJ/N-UGC")
		# Pos<-c(Pos,grep("Blood2",dat_temp$X1))
		# dat_temp2<-dat_temp[Pos,]
		# Pos<-which(dat_temp$X1 == "Lung cancer BJ")
		# Pos<-c(Pos,grep("Blood3",dat_temp$X1))
		# dat_temp3<-dat_temp[Pos,]
		# dat_temp10<-do.call(rbind,list(dat_temp0,dat_temp1,dat_temp2,dat_temp3))
		# dat_temp<-dat_temp10


		# manually set order for plot with crc, lc and ea and most relevant tissue for gene expression in europeans
		# dat_temp0<-dat_temp[dat_temp$X1  ==  "AA:DGLA / D5D (CHARGE)",]
		# Pos<-which(dat_temp$X1 == "Colorectal cancer GECCO/CORECT/CCFR")
		# Pos<-c(Pos,grep("colon",dat_temp$X1))
		# dat_temp1<-dat_temp[Pos,]
		# Pos<-which(dat_temp$X1 == "Esophageal adenocarcinoma EAS/UKB")
		# Pos<-c(Pos,grep("esophagus",dat_temp$X1))
		# dat_temp2<-dat_temp[Pos,]
		# Pos<-which(dat_temp$X1 == "Lung cancer ILCCO/UKB")
		# Pos<-c(Pos,grep("lung",dat_temp$X1))
		# dat_temp3<-dat_temp[Pos,]
		# dat_temp10<-do.call(rbind,list(dat_temp0,dat_temp1,dat_temp2,dat_temp3))
		# dat_temp<-dat_temp10

		# if(sum(grep("carcinoma",fa_eqtl_traits))!=0){
		# 	cancer_traits<-fa_eqtl_traits[grep("carcinoma",fa_eqtl_traits)]
		# 	can_tab<-data.frame(cbind(cancer_traits,1:length(cancer_traits)),stringsAsFactors=F)
		# 	can_tab$X3<-1:nrow(can_tab)
		# 	names(can_tab)<-c("X1","X2","X3")			
		# 	dat_temp<-rbind(can_tab,dat_temp)			
		# }
		
		dat_temp$X2<-1:nrow(dat_temp) 						
		
		Temp<-merge(Temp,dat_temp,by.x="trait2",by.y="X1")
		Temp$trait2<-gsub("_"," ",Temp$trait2)
		# Temp$trait2<-gsub("Blood1","Blood",Temp$trait2)
		# Temp$trait2<-gsub("Blood2","Blood",Temp$trait2)
		# Temp$trait2<-gsub("Blood3","Blood",Temp$trait2)
		# Temp<-Temp[order(Temp$marker),]
		Temp<-Temp[order(Temp$pos),]
		if(any(is.na(Temp$pos))) stop("missing base pair positions")
		snporder<-unique(Temp$marker)

		Temp1<-split(Temp,f=Temp$X2)
		# FA1<-data.frame(Temp1[1])
		# Can1<-data.frame(Temp1[2])
		# Can2<-data.frame(Temp1[3])
		# any(duplicated(Can2$X3.marker))
		# Can1$X2.marker[!Can1$X2.marker %in% FA1$X1.marker]

		names(Temp1)<-paste("trait",1:length(Temp1),sep="")
		list2env(Temp1,envir=.GlobalEnv)
		# snporder<-snporder[snporder %in% row.names(ld.matrix)]
		# length(snporder)
		# dim(ld.matrix)

		ld.matrix2<-ld.matrix[row.names(ld.matrix) %in% snporder,colnames(ld.matrix) %in% snporder]
		ld.matrix2<-ld.matrix2[match(snporder,row.names(ld.matrix2)),match(snporder,colnames(ld.matrix2))]
		# head(ld.matrix2)
	
		# ld.matrix2<-ld.matrix2[order(row.names(ld.matrix2)),order(colnames(ld.matrix2))]
		# ld.matrix2<-ld.matrix[row.names(ld.matrix) %in% SNPs_keep,colnames(ld.matrix) %in% SNPs_keep]
		all(snporder == row.names(ld.matrix2))
		# length(trait1$marker)
		# length(snporder)
		# row.names(ld.matrix2)[!row.names(ld.matrix2) %in% trait1$marker]
		# which(Temp$marker == "rs198461")
		# if(any(row.names(ld.matrix2) != trait1$marker)) stop("marker order mismatched")
		
	    Z.list<-lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$z",sep=""))))
	    b.list<-lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$beta",sep=""))))
	    se.list<-lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$se",sep=""))))
	    pvalue.list<-lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$pval",sep=""))))
	    cases.list<-lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$cases",sep=""))))	     
	    controls.list<-lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$controls",sep=""))))	     
	    samplesize.list<-lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$samplesize",sep=""))))	     
		Z.matrix<-matrix(unlist(Z.list),ncol=length(Z.list),nrow=length(SNPs_keep))

		Markers<- trait1[,c("marker","chr","pos")]
		# Markers<- trait1[,c("marker","chr","pos")]
		
		# Markers<-unique(Temp[,c("marker","chr","pos")]) #this doesn't work because not all studies have used same build
		
		if(any(row.names(ld.matrix2) != Markers$marker)) stop("marker order mismatched")
		Names.list<-unique(unlist(lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$trait2",sep=""))))))
		Names.list2<-unique(unlist(lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$trait",sep=""))))))

		return(list(ld.matrix=ld.matrix2,Z.matrix=Z.matrix,Markers=Markers,Traits=Names.list,beta=b.list,se=se.list,cases=cases.list,controls=controls.list,samplesize=samplesize.list,Traits2=Names.list2,pvalues=pvalue.list))
}


make_title<-function(Dir=NULL,Info=NULL,gwis_file=NULL,type=".png",plot_strategy=NULL,nsnps=NULL,Trait_names=NULL,ld_eas_pop=NULL){

	Info2<-paste(unique(unlist(lapply(1:length(strsplit(Info,split=" ")),FUN=function(x) 
		unlist(strsplit(Info,split=" ")[x])[1]))),collapse="_")
	Info2<-gsub("_NA","",Info2)
	if(!is.null(gwis_file)){
		if(sum(grep("notimputed",gwis_file,invert=T))!=0){
			Info2<-paste0(Info2,"_imputed")
		}
		if(sum(grep("notimputed",gwis_file))!=0){
			Info2<-paste0(Info2,"_notimputed")
		}
	}
	
	if(!is.null(Trait_names)){
		if(sum(grep("adj",Trait_names))!=0){
			Info2<-paste0(Info2,"_adjSNP")
		}
	}

	if(!is.null(nsnps)){		
		Info2<-paste0(Info2,"_nsnps",nsnps)
	}

	Info2<-gsub("/","",Info2)

	Trait_names2<-Trait_names[grep("expression",Trait_names,invert=T)]

	Trait_names2<-strsplit(Trait_names2,split=" ")
	# Length_split<-unique(unlist(lapply(1:length(Trait_names2),FUN=function(x) 
		# length(unlist(Trait_names2[x])))))
	# if(length(Length_split) != 1) warning("trait names have inconsistent lengths after applying strsplit, which could cause problems with the make_title function")
	Trait_names2<-unlist(Trait_names2)
	# Trait_names2<-gsub("/","_",Trait_names2)
	# Index<-round(abs(rnorm(10000)[1])*1e6,0)
	Index<-paste(unlist(lapply(1:length(Trait_names2),FUN=function(x) nchar(Trait_names2[x]))),collapse="")
	Index<-paste(unlist(Index),collapse="")

	plot.title<-paste(Dir,paste(Info2,collapse="_"),type,sep="")
	plot.title<-gsub(":","to",plot.title)
	# if(plot_strategy == "allcancer_1trait_1study"){
		plot.title<-gsub(".png","_eQTLs_cancer.png",plot.title)
	# }
	
	# if(plot_strategy == "allcancer_1trait_1study_notissues"){
	# 	plot.title<-gsub(".png","_cancer.png",plot.title)
	# }

	if(sum(grep("cancer", Trait_names))!=0){		
		Cancer_names<-gsub(" ","_",paste(Trait_names[grep("cancer", Trait_names)],collapse="_"))
		plot.title<-gsub(".png",paste0("_",Cancer_names,".png"),plot.title)
	}

	if(!is.null(ld_eas_pop)){
		plot.title<-gsub(".png",paste0("_ldpop_",ld_eas_pop,".png"),plot.title)	
	}


	Index<-substr(Index,1,6)
	plot.title<-gsub(".png",paste0("_index",Index,".png"),plot.title)
	plot.title<-gsub("Colorectal_cancer_GECCO/CORECT/CCFR","GECCO_CORECT_CCFR",plot.title)
	# plot.title<-gsub(" ","_",plot.title)
	# plot.title<-gsub("\\(","",plot.title)
	# plot.title<-gsub(")","",plot.title)

	
	return(plot.title)
	# pdf.title<-paste("~/fatty-acids-mr/coloc/AA_to_DGLA_fads1_liver.pdf",sep="")
				# pdf(pdf.title, onefile=FALSE)
}

#     gene       SNP ensembl_gene_id
# 1 ELOVL2 rs3734398 ENSG00000197977
# 2  FADS1  rs174528 ENSG00000149485
# 3  FADS2  rs174528 ENSG00000134824
# 4    SCD  rs603424 ENSG00000099194

plot_height<-function(Dat=length(Trait_names)){
		# H1<-c(14,8,4)
		# H2<-c(3000,1285,400)
		# Height_table<-data.frame(cbind(H1,H2))
		# names(Height_table)<-c("Ntraits","Height")
		# Plot_height<- Height_table$Height[Height_table$Ntraits == Dat]
		Plot_height<-3000*(Dat/14)
		if(Dat==14) Plot_height<-3000
		return(Plot_height)
	}


format_cancer<-function(Dat=NULL,trait=NULL,trait2=NULL,marker=NULL,eaf=NULL,effect_allele=NULL,other_allele=NULL,beta=NULL,se=NULL,pval=NULL,cases=NULL,controls=NULL,ref.dat=NULL,ref.dat.marker="V2",ref.pos="V4"){

	Dat$trait<-Dat[,trait]
	Dat$trait2<-paste(Dat[,trait],Dat$study)
	
	Name_cols<-c("marker","effect_allele","other_allele","beta","se","eaf","pval","cases","controls")

	for(i in 1:length(Name_cols)){
		print(i)
		if(!is.na(eval(parse(text=Name_cols[i])))){
			names(Dat)[names(Dat) == eval(parse(text=Name_cols[i]))]<-Name_cols[i]
		}
	}

	any(duplicated(paste(Dat$marker,Dat$ID)))
	Dat<-Dat[!duplicated(paste(Dat$marker,Dat$ID)),]
	Dat$beta<-as.numeric(Dat$beta)
	Dat$se<-as.numeric(Dat$se)
	Dat$z<-Dat$beta/Dat$se
	Cols.keep<-names(Dat)[!names(Dat) %in% c("V1","V2","V3","V4")]
	Dat<-Dat[,Cols.keep]
	Dat<-merge(Dat,ref.dat,by.x="marker",by.y=ref.dat.marker)
	Dat$chr<-Dat$V1
	Dat$pos<-Dat[,ref.pos]		
	Dat<-Dat[order(Dat$chr,Dat$pos),]	
	Dat$cases<-as.numeric(Dat$cases)
	Dat$controls<-as.numeric(Dat$controls)
	Dat$samplesize<-Dat$cases+Dat$controls
	Dat$s<-Dat$cases/Dat$samplesize
	Dat$eaf<-as.numeric(Dat$eaf)
	maf<-Dat$eaf
	maf[which(maf>0.5)]<-1-maf[which(maf>0.5)]
	Dat$maf<-maf
	return(Dat)
}


format_crc_gecco<-function(Dat=NULL,ref.dat=NULL,ref.dat.marker="V2",ref.pos="V4"){
# code for situation where don't know how many cancers are being included or don't want to make assumptions about numbers of cancers
# length(Cancer_objects)
# eval(Cancer_objects[1])
# Cancers<-lapply(1:length(Cancer_objects),FUN=function(x)
#     	eval(parse(text=Cancer_objects[x])))
# names(Cancers)<-paste("Cancer",1:length(Cancers),sep="")
# list2env(Cancers,envir=.GlobalEnv)

# assuming only 2 cancers present
	Dat$trait<-"Colorectal cancer"
	Dat$trait2<-"Colorectal cancer GECCO/CORECT/CCFR"
	names(Dat)[names(Dat) == "rs_number"]<-"marker"
	names(Dat)[names(Dat) == "Freq1"]<-"eaf"
	names(Dat)[names(Dat) == "Allele1"]<-"effect_allele"
	names(Dat)[names(Dat) == "Allele2"]<-"other_allele"
	names(Dat)[names(Dat) == "Effect"]<-"beta"
	names(Dat)[names(Dat) == "StdErr"]<-"se"
	names(Dat)[names(Dat) == "P.value"]<-"pval"
	# snps<-unique(Dat$marker[duplicated(Dat$marker)])
	# Dat[Dat$marker %in% snps,]
	Dat<-Dat[!duplicated(Dat$marker),] #duplicated SNPs seem to correspond to same base pair position but with different alleles. Perhaps these are non biallelic SNPs / ie. are triallelic

	Dat$z<-Dat$beta/Dat$se	
	Dat<-merge(Dat,ref.dat,by.x="marker",by.y=ref.dat.marker)
	# names(Dat)[names(Dat) == "V1"]<-"chr"
	Dat$chr<-Dat$V1
	Dat$pos<-Dat[,ref.pos]
	# names(Dat)[names(Dat) == "V4"]<-"pos"
	Dat<-Dat[order(Dat$chr,Dat$pos),]	
	Dat$cases<-82546
	Dat$controls<-211703 
	Dat$samplesize<-Dat$cases+Dat$controls
	Dat$s<-Dat$cases/Dat$samplesize
	maf<-Dat$eaf
	maf[which(maf>0.5)]<-1-maf[which(maf>0.5)]
	Dat$maf<-maf
	return(Dat)
}

format_lun<-function(Dat=NULL,ref.dat=NULL,ref.dat.marker="V2",ref.pos="V4"){
# code for situation where don't know how many cancers are being included or don't want to make assumptions about numbers of cancers
# length(Cancer_objects)
# eval(Cancer_objects[1])
# Cancers<-lapply(1:length(Cancer_objects),FUN=function(x)
#     	eval(parse(text=Cancer_objects[x])))
# names(Cancers)<-paste("Cancer",1:length(Cancers),sep="")
# list2env(Cancers,envir=.GlobalEnv)

# assuming only 2 cancers present
	Dat$trait<-"Lung cancer"
	names(Dat)[names(Dat) == "effect_allele_freq"]<-"eaf"
	names(Dat)[names(Dat) == "p"]<-"pval"
	names(Dat)[names(Dat) == "snp"]<-"marker"
	Dat<-Dat[Dat$file == "Onco_TRICL_032116_Overall.csv.tab",]
	if(any(duplicated(Dat$marker))) stop("duplicated SNPs present")
	
	Dat$z<-Dat$beta/Dat$se
	Dat$trait2<-Dat$trait
	head(Dat)
	Dat<-merge(Dat,ref.dat,by.x="marker",by.y=ref.dat.marker)
	Dat$chr<-Dat$V1
	Dat$pos<-Dat[,ref.pos]	
	Dat<-Dat[order(Dat$chr,Dat$pos),]
	return(Dat)
}


format_crc_accc<-function(Dat=NULL,ref.dat=NULL,ref.dat.marker="V2",ref.pos="V4"){
	SNP<-gregexpr(":",Dat$snp)
	Test<-unlist(lapply(1:length(SNP),FUN=function(x)
		length(unlist(SNP[x]))))
	Pos<-which(Test==3)
	Dat<-Dat[Pos,]
	SNP<-unlist(strsplit(Dat$snp,split=":"))
	Dat$chr<-SNP[seq(1,length(SNP),by=4)]
	Dat$bp<-SNP[seq(2,length(SNP),by=4)]	
	Dat<-merge(Dat,ref.dat,by.x=c("chr","bp"),by.y=c("V1",ref.pos))	
	Dat$trait<-"Colorectal cancer"
	Dat$trait2<-"Colorectal cancer (ACCC)"
	names(Dat)[names(Dat) == "Allele1"]<-"effect_allele"
	names(Dat)[names(Dat) == "Allele2"]<-"other_allele"
	names(Dat)[names(Dat) == "Freq1"]<-"eaf"
	names(Dat)[names(Dat) == "Effect"]<-"beta"
	names(Dat)[names(Dat) == "StdErr"]<-"se"
	names(Dat)[names(Dat) == "P.value"]<-"pval"
	names(Dat)[names(Dat) == "HetPVal"]<-"Phet"
	names(Dat)[names(Dat) == ref.dat.marker]<-"marker"

	Dat$study<-"ACCC" #The Asian Colorectal Cancer Consortium?
	Dat$population<-"East Asian"
	Dat$pmid<-"31826910"
	# Dat$effect_allele_confirmed<-TRUE
	# Dat$UKbiobank<-FALSE
	Dat$cases<-23572
	Dat$controls<-48700 
	Dat$samplesize<-Dat$cases+Dat$controls
	Dat$s<-Dat$cases/Dat$samplesize
	maf<-Dat$eaf
	maf[which(maf>0.5)]<-1-maf[which(maf>0.5)]
	Dat$maf<-maf
	Dat$z<-Dat$beta/Dat$se
	names(Dat)[names(Dat)=="bp"]<-"pos"
	Dat<-Dat[order(Dat$chr,Dat$pos),]
	Dups<-Dat$marker[duplicated(Dat$marker)]
	Dat<-Dat[!Dat$marker %in% Dups,] #duplicates seem to correspond to cnvs	
	return(Dat)
}



	# D5D	AA_to_DGLA
	# D5D	EPA_to_ETA
	# D6D	GLA_to_LA
	# D6D	SDA_to_ALA
	# ELOVL2	DPAn6_to_AA
	# ELOVL2	DHA_to_DPAn3
	# ELOVL5	DGLA_to_GLA
	# ELOVL5	ETA_to_SDA
	# SCD	POA_to_PA
	# SCD	OA_to_SA
	# ELOVL2/5	ADA_to_AA
	# ELOVL2/5	DPAn3_to_EPA



# region = NULL #how many base pairs should be included in the regional association plot? If set to NULL, the entire available region (1 million base pairs).  

# fix_charge=TRUE#Add DPAn6_to_AA to CHARGE. This trait is missing from CHARGE but not missing from Shin and Framingham. The purpose is to make the dimensions of the plots consistent, which makes comparing them side by side easier. The plot with DPAn6_to_AA_deleteme or CHARGE_deleteme should be deleted/ignore/covered up with a white square (sorry a bit hacky but I don't know how to change the underlying plotting code from gassocplot to do something less hacky). Note that the colocalisation results including fix_charge=TRUE are not going to be accurate when trait is DPAn6_to_AA. Therefore the function does not return colocalisation results when fix_charge=TRUE

# plot_strategy="alltraits_1study_alltissues", #Which coloc strategy to use. 
	# "everything" will include all traits, all studies and all tissues in the colocalisation analysis and in the regional association plot. 
	# "alltraits_1study_alltissues" includes 1 study (CHARGE, Framingham or Shin), all fatty acid traits and all gene expression traits (all tissues) per analysis 
	# "1trait_allstudies_alltissues" includes 1trait, all studies and all tissues per analysis.  
	# allcancer_1trait_alltissues includes all cancers, 1 fatty acid trait and gene expression across all tissues
	# allcancer_1trait_1study includes all cancers and 1 fatty acid trait from 1 study and all tissue eQTLs. This is used for example to make regional association plots of lung cancer, colorectal cancer and AA:DGLA (CHARGE consortium), together with various tissue eQTLs

# gene_table:	

   # gene       SNP ensembl_gene_id
# 1 ELOVL2 rs3734398 ENSG00000197977
# 2  FADS1  rs174528 ENSG00000197977
# 3  FADS2  rs174528 ENSG00000134824
# 4    SCD  rs603424 ENSG00000099194

format_data3<-function(crc_data=NULL,lun_data=NULL,cancer_data=NULL,gtex_data=NULL,fa.tab=NULL,ref=NULL,eqtlgen_data=NULL,gene_tab=NULL,bbj_eqtl_data=NULL,ld_eas_pop=NULL,gene=NULL,region=NULL,gtex_tissues=NULL,fix_charge=FALSE,studies=NULL,gwis_file=NULL,fix_log10=TRUE,pre_formatted_data=FALSE){	

	# for(j in 1:length(gene_tab$gene)){		
	# tissues=c("colon sigmoid","whole blood","lung","adipose subcutaneous","artery coronary","liver","colon transverse","adipose visceral omentum")

	# print(gene)
	# snp<-gene_tab$SNP[gene_tab$gene == gene]
	
	if(is.null(ld_eas_pop)){
		Region_object_list<-define_region(gene=gene,region=region,ref=ref) #region refers to how many base pairs to plot . region=NULL plots 1million base pairs (the default). 
		SNPlist<-data.frame(Region_object_list[3],stringsAsFactors=F)
		SNPlist<-SNPlist$SNP
	}
	
	if(!is.null(ld_eas_pop)){
		Region_object_list<-define_region_eas(gene=gene,region=region,ref=ref,ld_eas_pop=ld_eas_pop)
		SNPlist<-unlist(Region_object_list[3])		
	}

	ref3<-data.frame(Region_object_list[1],stringsAsFactors=F)
	Gen.table<-data.frame(Region_object_list[2],stringsAsFactors=F)
	
	if(!is.null(studies)){
		fa.tab2<-fa.tab[fa.tab$study %in% studies,]				
	}else{
		fa.tab2<-fa.tab
	}
	
	# fa.tab2<-format_fastudies(Dat=fa.tab2)

	# print(study)
	
	if(all(unique(fa.tab$study) == "CHARGE") & fix_charge){
		fa.tab_fix<-fa.tab2[fa.tab2$trait == "AA:DGLA",]
	
		# fa.tab_fix<-fa.tab2[fa.tab2$study == "CHARGE" & fa.tab2$trait == "AA:DGLA",]
		fa.tab_fix$trait<-"DPAn6:AA_deleteme"
		fa.tab2<-rbind(fa.tab2,fa.tab_fix)
	}			

	# tissues_keep<-tissues[tissues %in% gtex_tissues ]
	Traits1<-format_data1(Fa.tab=fa.tab2,ref.dat=ref3,gwis_test=TRUE,gwis_file=gwis_file,gene=gene,fix_log10=fix_log10)
	
	# head(Traits1)
	# study="CHARGE"
	# unique(Traits1$marker)
	# dups<-unique(fa.tab2$SNP[duplicated(fa.tab2$SNP)])
# fa.tab2[fa.tab2$SNP == dups[1],]
		
	if(!is.null(gtex_data)){
		Traits2<-format_data1(Fa.tab=gtex_data,ref.dat=ref3,gene=gene,gene_tab=gene_tab,tissues=gtex_tissues,gtex_test=TRUE,fix_log10=fix_log10)
	}

	if(!is.null(eqtlgen_data)){
		Traits3<-format_data1(Fa.tab=eqtlgen_data,ref.dat=ref3,gene=gene,gene_tab=gene_tab,eqtlgen_test=TRUE,fix_log10=fix_log10)
	}
		
	if(!is.null(bbj_eqtl_data)){
		Traits4<-format_data1(Fa.tab=bbj_eqtl_data,ref.dat=ref3,gene=gene,gene_tab=gene_tab,bbj_eqtl_test=TRUE,fix_log10=fix_log10)
	}

	if(!is.null(crc_data)){
		Traits5<-crc_data 
		# format_crc(Dat=crc_data,ref.dat=ref3)
	}

	if(!is.null(lun_data)){
		Traits6<-lun_data
		# format_lun(Dat=lun_data,ref.dat=ref3)
	}

	if(!is.null(cancer_data)){
		# if(any(names(cancer_data) == "fisherp")){
			# cancer_data$z<-qnorm(cancer_data$fisherp/2,lower.tail=F)		
			# print("fisherp!")
		# }
		cancer_data_list<-split(cancer_data,f=cancer_data$study)

		names(cancer_data_list)<-paste("Traits",7:c(6+length(cancer_data_list)),sep="")
		list2env(cancer_data_list,env = environment())		
	}

	# if(!is.null(crc_data.eas)){
	# 	Traits7<-crc_data.eas
	# 	# format_lun(Dat=lun_data,ref.dat=ref3)
	# }

	# if(!is.null(fa.tab.eas)){
	# 	Traits8<-format_data1(Fa.tab=fa.tab.eas,ref.dat=ref3,gwis_test=TRUE,gwis_file=gwas_file.eas,gene=gene,fix_log10=fix_log10)
	# }

	# Traits7<-format_data1(Fa.tab=fa.tab2,ref.dat=ref3,gwis_test=TRUE,gwis_file=gwis_file,study="Framingham",gene=gene)	
	trait_list<-ls()[grep("Traits[0-9]",ls())] 
	
	Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
		!is.null(nrow(eval(parse(text=trait_list[x]))))))
	rm(list=trait_list[!Test]) #remove the objects with no data
	# ,envir=.GlobalEnv
	trait_list<-ls()[grep("Traits[0-9]",ls())] 
	trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))
	Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,study=studies)
	
	if(pre_formatted_data){
		return(trait_list2)
	}else{
		return(Plot)
	}

}

# format_data3<-function(crc_data1=NULL,gtex_data1=NULL,fa.tab1=NULL,ref1=NULL,eqtlgen_data1=NULL,ld_eas_pop=NULL,gene=NULL,region=NULL,tissues=NULL){	
# 	if(is.null(ld_eas_pop)){
# 		Region_object_list<-define_region(gene=gene,region=region,ref=ref1) #region refers to how many base pairs to plot . region=NULL plots 1million base pairs (the default). 
# 		SNPlist<-data.frame(Region_object_list[3],stringsAsFactors=F)
# 		SNPlist<-SNPlist$SNP
# 	}

# 	if(!is.null(ld_eas_pop)){
# 		Region_object_list<-define_region_eas(gene=gene,region=region,ref=ref1,ld_eas_pop=ld_eas_pop)
# 		SNPlist<-unlist(Region_object_list[3])		
# 	}

# 	ref3<-data.frame(Region_object_list[1],stringsAsFactors=F)
# 	Gen.table<-data.frame(Region_object_list[2],stringsAsFactors=F)


# 	fa.tab2<-format_fastudies(Dat=fa.tab1)					
# 	Traits<-unique(fa.tab2$trait)
# 	Traits1<-format_data1(Fa.tab=fa.tab2,ref.dat=ref3,gwis_test=TRUE)					
# 	# length(unique(Traits2$marker))
# 	# unique(Traits3[,c("trait","trait2")])
# 	Traits2<-format_data1(Fa.tab=gtex_data1,ref.dat=ref3,gene=gene,gene_tab=gene_tab1,tissues=tissues,gtex_test=TRUE)
# 	Traits3<-format_data1(Fa.tab=eqtlgen_data1,ref.dat=ref3,gene=gene,gene_tab=gene_tab1,eqtlgen_test=TRUE)
# 	Traits4<-format_crc(Dat=crc_data1,ref.dat=ref3)
# 	trait_list<-ls()[grep("Traits[0-9]",ls())] 
# 	Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
# 		!is.null(nrow(eval(parse(text=trait_list[x]))))))
# 	rm(list=trait_list[!Test]) #remove the objects with no data
# 	trait_list<-ls()[grep("Traits[0-9]",ls())] 
# 	trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))				
# 	Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,Studies=unique(fa.tab2$study))		
# 	return(Plot)		
# }




get_nice_names<-function(dat=NULL,Trait_names=NULL){
	# traits colocalise
	Trait_names<-gsub("\n","",Trait_names)
	Pos<-dat$traits!="None"
	Trait_id<-gsub("T","",dat$traits[Pos])
	Traits_match_id<-lapply(1:length(Trait_id),FUN=function(x) Trait_names[as.numeric(unlist(strsplit(Trait_id[x],split=", ")))])
	Traits_match_id2<-unlist(lapply(1:length(Traits_match_id),FUN=function(i) paste(unlist(Traits_match_id[i]),collapse=", ")))	
	dat$traits_nice_names<-NA 
	dat$traits_nice_names[which(Pos)]<-Traits_match_id2
	# dropped traits
	Pos<-dat$traits=="None"
	Trait_id<-gsub("T","",dat$dropped_trait[Pos])
	Traits_match_id<-lapply(1:length(Trait_id),FUN=function(x) Trait_names[as.numeric(unlist(strsplit(Trait_id[x],split=", ")))])
	Traits_match_id2<-unlist(lapply(1:length(Traits_match_id),FUN=function(i) paste(unlist(Traits_match_id[i]),collapse=", ")))	
	dat$dropped_trait_nice_names<-NA 
	dat$dropped_trait_nice_names[which(Pos)]<-Traits_match_id2
	return(dat)
}


get_gene_name<-function(){
	X<-Trait_names[grep("expression",Trait_names)]
	Pos<-unlist(gregexpr(" expression",X))
	gene<-unique(trimws(substring(X,1,Pos)))
	return(gene)
}

coloc_abf_function<-function(Dat=NULL,Traits=NULL,type_1=NULL,type_2=NULL,sdY=1,s=NULL){
	i<-match(Traits,Dat$Traits)
	if(type_1 == "quant"){
		dataset1 <- list(beta=unlist(Dat$beta[i[1]]),varbeta=unlist(Dat$se[i[1]])^2, type=type_1,sdY=sdY)	
	}
	if(type_1 == "cc"){
		dataset1 <- list(beta=unlist(Dat$beta[i[1]]), varbeta=unlist(Dat$se[i[1]])^2,type=type_1,s=s)
	}		
	if(type_2 == "cc"){
		dataset2 <- list(beta=unlist(Dat$beta[i[2]]), varbeta=unlist(Dat$se[i[2]])^2,type=type_2,s=s)
	}
	if(type_2 == "quant"){
		dataset2 <- list(beta=unlist(Dat$beta[i[2]]),varbeta=unlist(Dat$se[i[2]])^2, type=type_2,sdY=sdY)
	}
	result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
	#List into data frame.
	Df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
	Traits<-gsub("\n","",Traits)
	Df$traits<-paste(Traits,collapse=" & ")

	#Label the columns in the data frame.
	names(Df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf","traits")
	return(Df)		
}


hyprcoloc_function<-function(Dat=NULL,gene=NULL,gwis_file=NULL,binary.outcomes=rep(0,length(Dat$Traits))){	
	# Hypercolocalisation 	
	names(Dat)<-NULL
	ld.matrix<-data.frame(Dat[1])
	Markers<-data.frame(Dat[3])
	Trait_names<-unlist(Dat[4])
	B.matrix<-data.frame(Dat[5])		
	SE.matrix<-data.frame(Dat[6])
	traits <- paste0("T", 1:dim(B.matrix)[2])
	# traits <- Trait_names
	row.names(B.matrix)<-Markers$marker
	row.names(SE.matrix)<-Markers$marker
	row.names(ld.matrix)<-Markers$marker	
	names(B.matrix)<-traits
	B.matrix<-as.matrix(B.matrix)	
	names(SE.matrix)<-traits
	SE.matrix<-as.matrix(SE.matrix)	
	names(ld.matrix)<-traits
	ld.matrix<-as.matrix(ld.matrix)	
	
	res <- hyprcoloc(B.matrix, SE.matrix, trait.names=traits, snp.id=Markers$marker,prior.1 = 1e-04,#default 
       prior.2 = 0.98,#default 
       binary.outcomes=binary.outcomes,ld.matrix = ld.matrix)
      
       # ld.matrix
	# , sensitivity = TRUE
	res_dat<-res$results
	if(is.null(gene)){
		gene<-get_gene_name()	
	}
	res_dat$gene<-gene
	Markers$pos<-as.numeric(Markers$pos)
	res_dat$min<-min(Markers$pos)	
	res_dat$max<-max(Markers$pos)
	res_dat$bp<-max(Markers$pos)-min(Markers$pos)
	res_dat$gwis_fileimputed<-FALSE
	if(!is.null(gwis_file)){
		if(sum(grep("notimputed",gwis_file,invert=T))!=0){
			res_dat$gwis_fileimputed<-TRUE
		}
	}
	res_dat<-get_nice_names(dat=res_dat,Trait_names=Trait_names)	
	res_dat$nsnps<-length(Markers$marker)
	# res_dat$trait_group<-Traits[i]
	# Index<-paste(gene_tab$gene[j],Traits[i])
	# head(res_dat)
	# res_list[[Index]]<-res_dat
	return(res_dat)
}

make_file_name<-function(Dir=NULL,Method=NULL,Info=NULL){
	Notimputed=sum(grep("notimputed",gwis_file1))==1
	Pos<-1:length(Info)
	Pos2<-unlist(lapply(c(Trait),FUN=function(x)
		grep(x,Info)))
	Pos<-Pos[which(!Pos %in% Pos2)]
	Info<-paste(Info[Pos],collapse="_")
	Info<-unlist(strsplit(Info,split="expression"))[1]

	Info<-gsub(" ","_",Info)
	Info<-gsub("/","_",Info)
	if(!Notimputed) File_name<-paste0(Dir,Method,"_",Info,"_",gsub(":","_",Trait),"_","_imputed.txt")
	if(Notimputed) File_name<-paste0(Dir,Method,"_",Info,"_",gsub(":","_",Trait),"_","_notimputed.txt")
	File_name<-gsub(" ","_",File_name)
	File_name<-gsub("__","_",File_name)
	if(sum(grep("eastasians",gwis_file1))==1) File_name<-gsub(".txt","_east_asian.txt",File_name)
	return(File_name)
}


b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}

format_ref<-function(ref=NULL){
	V1<-unlist(strsplit(ref$V1,split="chr"))
	V1<-V1[V1!=""]
	ref$V1<-as.numeric(V1)			
	V2<-ref$V4
	V4<-ref$V2
	ref$V2<-V2
	ref$V4<-V4
	return(ref)
}

# gtex ALT allele is effect allele
# variant_id:  variant ID in the format {chr}_{pos_first_ref_base}_{ref_seq}_{alt_seq}_b37

call_coloc_abf_res<-function(Dat=NULL,trait1_cancer=FALSE){
	coloc_abf_res<-NULL
	for(i in 2:length(unlist(Dat$Traits))){
	# for(i in 4:length(unlist(Dat$Traits))){
	# for( i in 1){
		print(i)	
		# Traits<-Dat$Traits[c(2,i)] #modification for lung cancer
		# Traits<-Dat$Traits[c(3,i)]
		Traits<-Dat$Traits[c(1,i)]

		print('--------')
		print(Traits)

		if(trait1_cancer){
			s<-unique(unlist(Dat$cases[1]))/unique(unlist(Dat$controls[1]))
			coloc_abf_res[[i]]<-coloc_abf_function(Dat=Dat,Traits=Traits,type_1="cc",type_2="quant",s=s)
		}

		# if(Traits[1] == "Coloretcal cancer") {
		# 	s<-unique(crc_data1$s)
		# 	coloc_abf_res[[i]]<-coloc_abf_function(Dat=Dat,Traits=Traits,type_1="cc",type_2="quant",s=s)
		# }
		# if(Traits[1] == "Lung cancer") {
		# 	s<-29000/85449 
		# 	type_1<-"cc"		
		# 	coloc_abf_res[[i]]<-coloc_abf_function(Dat=Dat,Traits=Traits,type_1=type_1,type_2="quant",s=s)
		if(!trait1_cancer){
			coloc_abf_res[[i]]<-coloc_abf_function(Dat=Dat,Traits=Traits,type_1="quant",type_2="quant")
		}
	}
	return(coloc_abf_res)
}
