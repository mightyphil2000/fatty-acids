# install.packages("devtools")
# library(devtools)
# install_github("jrs95/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
# install_github("jrs95/gassocplot")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(hyprcoloc)
library(plyr)
library(gassocplot)
library(biomaRt)

Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
Attr<-listAttributes(Mart)


coloc_results<-colocalisation(
	cancer1="~/fatty-acids/data/tricl_snps_ukb.Rdata",#genetic associations for lung cancer from TRICL consortium in 1 mb region centred around the index SNP for each region. This file was created using ~/fatty-acids/scripts/extract_SNPs_cancer.R. 
	cancer1_name="Lung cancer",
	cancer2="~/MR_FattyAcids/data/summary data/colorectal_cancer/061119/crc_snps.Rdata",#genetic associations for colorectal cancer from GECCO/CORECT consortium in 1 mb region centred around the index SNP for each region. This file was created using ~/fatty-acids/scripts/extract_SNPs_cancer.R. 
	cancer2_name="Colorectal cancer",
	gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",#genetic associations for gene expression from gtex version8 in 1 mb region centred around the index SNP for each region in estimated in all GTEx participants, including non-Europeans (84.6% white/European, 12.9% African American). The file can also be found here /projects/MRC-IEU/users/ph14916/gtex2. This file was created using ~/fatty-acids/colosalition/scripts/extract_regions.sh and bash_input_gtex2.txt . build 38
	eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", #genetic associations for gene expression from eqtlgen in a 1mb region centred around the index SNP for each region. This file was created using extract_regions_eqtlgen.sh and extract_regions_eqtlgen2.sh
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_imputed.Rdata",
	# "~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios.Rdata", #genetic associations from gwis in a 1mb regino centred around the index SNP for each region. Assumed to contain three data frames called Charge1, Shin1 and Framingham1. The file can also be found here /projects/MRC-IEU/users/ph14916/fatty_acids_summary/gwis/. This file was created using the script ~/fatty-acids/colocalisation/scripts/extract_snp_data.R 
	ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz", #reference data. all SNPs in random 10,000 European participants in UK Biobank and their chromosomal coordinates in GRCh37. I obtained this file from Tom Richardson. The file can also be found here /projects/MRC-IEU/users/ph14916/ukb_bed 
	ref_turn_off=TRUE, #don't load the ref dataset e.g. because it is already loaded in memory. this can be convenient because it takes a long time to load the reference data
	tissues=c("colon sigmoid","whole blood","lung","adipose subcutaneous","artery coronary","liver","colon transverse","adipose visceral omentum"), #the tissues to include. Must correspond to file name in gtex. possible tissues: "colon sigmoid","whole blood","lung","adipose subcutaneous","artery coronary","liver","colon transverse","adipose visceral omentum"
	coloc_strategy="alltraits_1study_alltissues", #Which coloc strategy to use. 
	# "everything" will include all traits, all studies and all tissues in the colocalisation analysis and in the regional association plot. 
	# "alltraits_1study_alltissues" includes 1 study (CHARGE, Framingham or Shin), all fatty acid traits and all gene expression traits (all tissues) per analysis 
	# "1trait_allstudies_alltissues" includes 1trait, all studies and all tissues per analysis.  
	# allcancer_1trait_alltissues includes all cancers, 1 fatty acid trait and gene expression across all tissues
	# allcancer_1trait_1study includes all cancers and 1 fatty acid trait from 1 study and all tissue eQTLs. This is used for example to make regional association plots of lung cancer, colorectal cancer and AA:DGLA (CHARGE consortium), together with various tissue eQTLs
	turn_off_plot=FALSE, #for convenience can turn the plot function off and only do colocalisation. 
	region=NULL, #how many base pairs should be included in the regional association plot? If set to NULL, the entire available region (1 million base pairs).  
	fix_charge=TRUE #Add DPAn6_to_AA to CHARGE. This trait is missing from CHARGE but not missing from Shin and Framingham. The purpose is to make the dimensions of the plots consistent, which makes comparing them side by side easier. The plot with DPAn6_to_AA_deleteme or CHARGE_deleteme should be deleted/ignore/covered up with a white square (sorry a bit hacky but I don't know how to change the underlying plotting code from gassocplot to do something less hacky). Note that the colocalisation results including fix_charge=TRUE are not going to be accurate when trait is DPAn6_to_AA. Therefore the function does not return colocalisation results when fix_charge=TRUE
	)

res_dat<-do.call(rbind,coloc_results[[1]])
Col1<-which(names(res_dat) == "region")
Col2<-which(names(res_dat) != "region")
res_dat<-res_dat[,c(Col1,Col2)]
res_dat2<-res_dat[res_dat$traits != "None",]

write.table(res_dat2,paste("~/fatty-acids/results/coloc/hyprcoloc_",coloc_results[[2]],"_",coloc_results[[3]],".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(res_dat,paste("~/fatty-acids/results/coloc/hyprcoloc_",coloc_results[[2]],"_",coloc_results[[3]],"_incl_dropped_traits.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)



colocalisation<-function(gtex=NULL,eqtlgen=NULL,gwis=NULL,ref_dat=NULL,tissues=NULL,coloc_strategy=NULL,turn_off_plot=FALSE,ref_turn_off=FALSE,region=NULL,fix_charge=FALSE){

	List1<-ls()
	load(cancer1)
	load(cancer2)
	List2<-ls()
	Cancer_objects<-List2<-List2[!List2 %in% List1]

	load(gtex_file)
	eqtlgen_data<-read.table(eqtlgen_file,sep="\t",stringsAsFactors=F,head=T)
	dim(eqtlgen_data)
	load(gwis_file)

	gtex2<-format_gtex(df4)
	# gtex3<-gtex2[gtex2$Ens == "ENSG00000149485" & gtex2$file == "Whole_Blood.allpairs.txt",]
	# # head(gtex3)
	# gtex3$Z<-gtex3$slope/gtex3$slope_se
	# gtex4<-gtex3[order(abs(gtex3$Z),decreasing=T),]
	# gtex4[grep("61595564",gtex4$variant_id),]

	Charge1$study <- "CHARGE"
	Shin1$study <- "Shin"
	Tin1$study <- "Framingham"
	fa.tab<-do.call(rbind,list(Charge1,Shin1,Tin1))
	
	if(!ref_turn_off){
		ref <- read.table(ref_dat,stringsAsFactors=F,head=F)
	}
	
	# ref[ref$V2 %in% c("rs174574","rs968567"),]
	# ref$V4[ref$V2 == "rs12210577"]-500000
	# ref$V4[ref$V2 == "rs12210577"]+500000
	# Tissues<-c("Liver.allpairs.txt","Adipose_Subcutaneous.allpairs.txt","Adipose_Visceral_Omentum.allpairs.txt","Whole_Blood.allpairs.txt")

	gene_tab<-make_gene_table(snps=c("rs174528","rs174528","rs3734398","rs603424"),genes=c("FADS1","FADS2","ELOVL2","SCD")) #function assumes that these are the SNPs and genes of interest

	res_list<-NULL
	
	reg_plot<-function(j=NULL)
	j<-3
	for(j in 1:length(gene_tab$gene)){
		
		print(gene_tab$gene[j])
		snp<-gene_tab$SNP[j]
		Region_object_list<-define_region(index_snp=snp,region=region,ref=ref) #region refers to how many base pairs to plot . region=NULL plots 1million base pairs (the default). 
		ref3<-data.frame(Region_object_list[1],stringsAsFactors=F)
		Gen.table<-data.frame(Region_object_list[2],stringsAsFactors=F)
		SNPlist<-unlist(Region_object_list[3])

		if(coloc_strategy=="alltraits_1study_alltissues"){
			Studies<-unique(fa.tab$study)
			for(i in 1:length(Studies)){				
				print(i)
				print(Studies[i])
				fa.tab2<-format_fastudies(Dat=fa.tab)
				# unique(Fa.tab$trait)
				if(Studies[i] == "CHARGE" & fix_charge){
					fa.tab_fix<-fa.tab2[fa.tab2$study == "CHARGE" & fa.tab2$trait == "AA:DGLA",]
					fa.tab_fix$trait<-"DPAn6:AA_deleteme"
					fa.tab2<-rbind(fa.tab2,fa.tab_fix)
				}			

				tissues_keep<-tissues[tissues %in% c("whole blood","adipose subcutaneous","liver","adipose visceral omentum")]
				Traits1<-format_data1(Fa.tab=fa.tab2[fa.tab2$study == Studies[i],],ref.dat=ref3,gwis_test=TRUE)
				Traits2<-format_data1(Fa.tab=gtex2,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,tissues=				tissues_keep,gtex=TRUE)
				Traits3<-format_data1(Fa.tab=eqtlgen_data,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,eqtlgen=TRUE)
				trait_list<-ls()[grep("Traits[0-9]",ls())] 
				Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
					!is.null(nrow(eval(parse(text=trait_list[x]))))))
				rm(list=trait_list[!Test]) #remove the objects with no data
				trait_list<-ls()[grep("Traits[0-9]",ls())] 

				trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))
				Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,coloc_strategy=coloc_strategy,study=Studies[i])
				ld.matrix<-data.frame(Plot[1])
				Z.matrix<-data.frame(Plot[2])
				Markers<-data.frame(Plot[3])
				Trait_names<-unlist(Plot[4])
				Markers$pos<-as.integer(Markers$pos)
				basepairs<-max(Markers$pos)-min(Markers$pos)
				# Z.matrix<-Z.matrix/2 #if the Y axis is flooring the P value can fix by dividing the Z score by 2

				extra_info<-""
				if(fix_charge) extra_info<-"fix_charge"
				if(!turn_off_plot){
					plot.title<-make_title(Dir="~/fatty-acids/colocalisation/results/plots/",Info=c(Studies[i],gene_tab$gene[j],basepairs,extra_info),gwis_file=gwis_file,type=".png")
					png(plot.title, width=500,height=3000)
						stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names)
					dev.off()
				}
				
				top.marker=""


				out.file<-gsub(".png","_Nsnps.txt",plot.title)
				write.table(paste0(print(Studies[i])," (N SNPs =", nrow(Markers),")"),out.file,col.names=F,row.names=F,quote=F)
			

				
				# Hypercolocalisation 
				# B.matrix<-data.frame(Plot[5])
				# SE.matrix<-data.frame(Plot[6])
				# traits <- paste0("T", 1:dim(B.matrix)[2])
				# # traits <- Trait_names
				# row.names(B.matrix)<-Markers$marker
				# row.names(SE.matrix)<-Markers$marker
				# names(B.matrix)<-traits
				# B.matrix<-as.matrix(B.matrix)
				# names(SE.matrix)<-traits
				# SE.matrix<-as.matrix(SE.matrix)
				# res <- hyprcoloc(B.matrix, SE.matrix, trait.names=traits, snp.id=Markers$marker)	
				# res_dat<-res$results
				# res_dat$region<-gene_tab$gene[j]
				# res_dat$min<-min(Markers$pos)
				# res_dat$max<-max(Markers$pos)
				# res_dat$bp<-basepairs	
				# res_dat<-get_nice_names(dat=res_dat,Trait_names=Trait_names)
				# res_dat$study_group<-Studies[i]
				# Index<-paste(gene_tab$gene[j],Studies[i])
				# res_list[[Index]]<-res_dat
			
			}			
		}


		# coloc_strategy<-"allcancer_1trait_1study"
		if(coloc_strategy=="allcancer_1trait_1study"){
			Traits<-unique(fa.tab$file)
			Studies<-unique(fa.tab$study)
			for(i in 1:length(Traits)){
				for(k in 1:length(Studies)){
					print(i)
					print(Traits[i])
					print(i)
					print(Studies[i])			
					fa.tab2<-format_fastudies(fa.tab)
					Traits<-unique(fa.tab2$trait)
					Traits1<-format_data1(Fa.tab=fa.tab2[fa.tab2$trait == Traits[i] & fa.tab2$study == Studies[k],],ref.dat=ref3)
					Traits2<-format_data1(Fa.tab=gtex2,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,tissues=tissues,gtex=TRUE)
					Traits3<-format_data1(Fa.tab=eqtlgen_dat,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,eqtlgen=TRUE)
					Traits4<-format_cancer()
					trait_list<-ls()[grep("Traits[0-9]",ls())] 
					Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
						!is.null(nrow(eval(parse(text=trait_list[x]))))))
					rm(list=trait_list[!Test]) #remove the objects with no data
					trait_list<-ls()[grep("Traits[0-9]",ls())] 
					trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))				
					Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,Studies=unique(fa.tab2$study),coloc_strategy=coloc_strategy)
					ld.matrix<-data.frame(Plot[1])
					Z.matrix<-data.frame(Plot[2])
					Markers<-data.frame(Plot[3])
					Trait_names<-unlist(Plot[4])
					Markers$pos<-as.integer(Markers$pos)
					basepairs<-max(Markers$pos)-min(Markers$pos)
					# Z.matrix<-Z.matrix/2 #if the Y axis is flooring the P value can fix by dividing the Z score by 2
					extra_info<-""
					if(fix_charge) extra_info<-"fix_charge"
					if(!turn_off_plot){
						
						# top_marker<-"rs2727270"				
						# top_marker<-"rs2280018"
						# top_marker<-"rs968567"
						# top_marker<-"rs2524299"
						# top_marker<-"rs174546"
						plot.title<-make_title(Dir="~/fatty-acids/results/coloc/plots/",Info=c(Traits[i],gene_tab$gene[j],basepairs,extra_info,top_marker),type=".png")						
						if(!is.null(top_marker)){
							png(plot.title, width=500,height=1500)
							stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names,top.marker=top_marker)
							dev.off()
						}else{
							png(plot.title, width=500,height=1500)
								stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names)
							dev.off()
						}
					}
					# top.marker=""

					# Hypercolocalisation 
					B.matrix<-data.frame(Plot[5])
					SE.matrix<-data.frame(Plot[6])
					traits <- paste0("T", 1:dim(B.matrix)[2])
					# traits <- Trait_names
					row.names(B.matrix)<-Markers$marker
					row.names(SE.matrix)<-Markers$marker
					names(B.matrix)<-traits
					B.matrix<-as.matrix(B.matrix)
					names(SE.matrix)<-traits
					SE.matrix<-as.matrix(SE.matrix)
					res <- hyprcoloc(B.matrix, SE.matrix, trait.names=traits, snp.id=Markers$marker)	
					res_dat<-res$results
					res_dat$region<-gene_tab$gene[j]
					res_dat$min<-min(Markers$pos)
					res_dat$max<-max(Markers$pos)
					res_dat$bp<-basepairs	
					res_dat<-get_nice_names(dat=res_dat,Trait_names=Trait_names)
					res_dat$trait_group<-Traits[i]
					Index<-paste(gene_tab$gene[j],Traits[i])
					res_list[[Index]]<-res_dat
				}
			}						
		}
		# coloc_strategy<-"allcancer_1trait_1study_notissues"
		if(coloc_strategy=="allcancer_1trait_1study_notissues"){
			Traits<-unique(fa.tab$file)
			Studies<-unique(fa.tab$study)
			for(i in 1:length(Traits)){
				for(k in 1:length(Studies)){
					print(i)
					print(Traits[i])
					print(i)
					print(Studies[i])			
					fa.tab2<-format_fastudies(fa.tab)
					Traits<-unique(fa.tab2$trait)
					trait_list<-ls()[grep("Traits[0-9]",ls())]
					rm(list=trait_list) 
					Traits1<-format_data1(Fa.tab=fa.tab2[fa.tab2$trait == Traits[i] & fa.tab2$study == Studies[k],],ref.dat=ref3)
					# Traits2<-format_data1(Fa.tab=gtex2,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,tissues=tissues,gtex=TRUE)
					# Traits3<-format_data1(Fa.tab=eqtlgen_dat,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,eqtlgen=TRUE)
					Traits2<-format_cancer()
					trait_list<-ls()[grep("Traits[0-9]",ls())] 
					Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
						!is.null(nrow(eval(parse(text=trait_list[x]))))))
					rm(list=trait_list[!Test]) #remove the objects with no data
					trait_list<-ls()[grep("Traits[0-9]",ls())] 
					trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))				
					Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,Studies=unique(fa.tab2$study),coloc_strategy=coloc_strategy)
					ld.matrix<-data.frame(Plot[1])
					Z.matrix<-data.frame(Plot[2])
					Markers<-data.frame(Plot[3])
					Trait_names<-unlist(Plot[4])
					Markers$pos<-as.integer(Markers$pos)
					basepairs<-max(Markers$pos)-min(Markers$pos)
					# Z.matrix<-Z.matrix/2 #if the Y axis is flooring the P value can fix by dividing the Z score by 2
					extra_info<-""
					if(fix_charge) extra_info<-"fix_charge"
					if(!turn_off_plot){
					# 	plot.title<-make_title(Dir="~/fatty-acids/results/coloc/plots/",Info=c(Traits[i],gene_tab$gene[j],basepairs,extra_info),type=".png")
					# 	png(plot.title, width=500,height=1500)
					# 		stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names)
					# 	dev.off()
					# }

					# top_marker<-"rs2727270"				
					# top_marker<-"rs2280018"
					# top_marker<-"rs968567"
					# top_marker<-"rs2524299"
					# top_marker<-"rs174546"
					# top_marker<-NULL
					plot.title<-make_title(Dir="~/fatty-acids/results/coloc/plots/",Info=c(Traits[i],gene_tab$gene[j],basepairs,extra_info,top_marker),type=".png")						
						if(!is.null(top_marker)){
							png(plot.title, width=500,height=600)
							stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names,top.marker=top_marker)
							dev.off()
						}else{
							png(plot.title, width=500,height=600)
								stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names)
							dev.off()
						}
					}

					# top.marker=""

					# Hypercolocalisation 
					B.matrix<-data.frame(Plot[5])
					SE.matrix<-data.frame(Plot[6])
					traits <- paste0("T", 1:dim(B.matrix)[2])
					# traits <- Trait_names
					row.names(B.matrix)<-Markers$marker
					row.names(SE.matrix)<-Markers$marker
					names(B.matrix)<-traits
					B.matrix<-as.matrix(B.matrix)
					names(SE.matrix)<-traits
					SE.matrix<-as.matrix(SE.matrix)
					res <- hyprcoloc(B.matrix, SE.matrix, trait.names=traits, snp.id=Markers$marker)	
					res_dat<-res$results
					res_dat$region<-gene_tab$gene[j]
					res_dat$min<-min(Markers$pos)
					res_dat$max<-max(Markers$pos)
					res_dat$bp<-basepairs	
					res_dat<-get_nice_names(dat=res_dat,Trait_names=Trait_names)
					res_dat$trait_group<-Traits[i]
					Index<-paste(gene_tab$gene[j],Traits[i])
					res_list[[Index]]<-res_dat
				}
			}						
		}



		}


		if(coloc_strategy=="1trait_allstudies_alltissues"){
			Traits<-unique(fa.tab$file)
			for(i in 1:length(Traits)){
				print(i)
				print(Traits[i])				
				fa.tab2<-format_fastudies(fa.tab)
				Traits<-unique(fa.tab2$trait)
				if(Traits[i] == "DPAn6:AA" & fix_charge){
					fa.tab_fix<-fa.tab2[fa.tab2$trait == "AA:DGLA" & fa.tab2$study == "CHARGE",]
					fa.tab_fix$trait<-"DPAn6:AA"
					fa.tab_fix$study<-"CHARGE_deleteme"
					fa.tab2<-rbind(fa.tab2,fa.tab_fix)
						
				}

				Traits1<-format_data1(Fa.tab=fa.tab2[fa.tab2$trait == Traits[i],],ref.dat=ref3)
				Traits2<-format_data1(Fa.tab=gtex2,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,tissues=tissues,gtex=TRUE)
				Traits3<-format_data1(Fa.tab=eqtlgen_dat,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,eqtlgen=TRUE)
				trait_list<-ls()[grep("Traits[0-9]",ls())] 
				Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
					!is.null(nrow(eval(parse(text=trait_list[x]))))))
				rm(list=trait_list[!Test]) #remove the objects with no data
				trait_list<-ls()[grep("Traits[0-9]",ls())] 

				trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))
			
				Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,Studies=unique(fa.tab2$study),coloc_strategy=coloc_strategy)
				ld.matrix<-data.frame(Plot[1])
				Z.matrix<-data.frame(Plot[2])
				Markers<-data.frame(Plot[3])
				Trait_names<-unlist(Plot[4])
				Markers$pos<-as.integer(Markers$pos)
				basepairs<-max(Markers$pos)-min(Markers$pos)
				# Z.matrix<-Z.matrix/2 #if the Y axis is flooring the P value can fix by dividing the Z score by 2
				extra_info<-""
				if(fix_charge) extra_info<-"fix_charge"
				if(!turn_off_plot){
					plot.title<-make_title(Dir="~/fatty-acids/results/coloc/plots/",Info=c(Traits[i],gene_tab$gene[j],basepairs,extra_info),type=".png")
					png(plot.title, width=500,height=1500)
						stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names)
					dev.off()
				}

				 # top.marker=""

				# Hypercolocalisation 
				B.matrix<-data.frame(Plot[5])
				SE.matrix<-data.frame(Plot[6])
				traits <- paste0("T", 1:dim(B.matrix)[2])
				# traits <- Trait_names
				row.names(B.matrix)<-Markers$marker
				row.names(SE.matrix)<-Markers$marker
				names(B.matrix)<-traits
				B.matrix<-as.matrix(B.matrix)
				names(SE.matrix)<-traits
				SE.matrix<-as.matrix(SE.matrix)
				res <- hyprcoloc(B.matrix, SE.matrix, trait.names=traits, snp.id=Markers$marker)	
				res_dat<-res$results
				res_dat$region<-gene_tab$gene[j]
				res_dat$min<-min(Markers$pos)
				res_dat$max<-max(Markers$pos)
				res_dat$bp<-basepairs	
				res_dat<-get_nice_names(dat=res_dat,Trait_names=Trait_names)
				res_dat$trait_group<-Traits[i]
				Index<-paste(gene_tab$gene[j],Traits[i])
				res_list[[Index]]<-res_dat
			}						
		}

		
		if(coloc_strategy =="everything"){
			fa.tab2<-format_fastudies(fa.tab)
			Traits1<-format_data1(Fa.tab=fa.tab2,ref.dat=ref3)
			Traits2<-format_data1(Fa.tab=gtex2,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,tissues=tissues,gtex=TRUE)
			Traits3<-format_data1(Fa.tab=eqtlgen_dat,ref.dat=ref3,gene=gene_tab$gene[j],gene_tab=gene_tab,eqtlgen=TRUE)
			trait_list<-ls()[grep("Traits[0-9]",ls())] 
			Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
				!is.null(nrow(eval(parse(text=trait_list[x]))))))
			rm(list=trait_list[!Test]) #remove the objects with no data
			trait_list<-ls()[grep("Traits[0-9]",ls())] 
			trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))
			Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,coloc_strategy=coloc_strategy)
			ld.matrix<-data.frame(Plot[1])
			Z.matrix<-data.frame(Plot[2])
			Markers<-data.frame(Plot[3])
			Trait_names<-unlist(Plot[4])
			Markers$pos<-as.integer(Markers$pos)
			basepairs<-max(Markers$pos)-min(Markers$pos)
			if(!turn_off_plot){
				plot.title<-make_title(Dir="~/fatty-acids/results/coloc/plots/",Info=c("everything",gene_tab$gene[j],basepairs),type=".png")
				png(plot.title, width=500,height=4500)
					stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=Trait_names)
				dev.off()
			}
			# Z.matrix<-Z.matrix/2 #if the Y axis is flooring the P value can fix by dividing the Z score by 2
			# Hypercolocalisation 
			B.matrix<-data.frame(Plot[5])
			SE.matrix<-data.frame(Plot[6])
			traits <- paste0("T", 1:dim(B.matrix)[2])
			# traits <- Trait_names
			row.names(B.matrix)<-Markers$marker
			row.names(SE.matrix)<-Markers$marker
			names(B.matrix)<-traits
			B.matrix<-as.matrix(B.matrix)
			names(SE.matrix)<-traits
			SE.matrix<-as.matrix(SE.matrix)
			res <- hyprcoloc(B.matrix, SE.matrix, trait.names=traits, snp.id=Markers$marker)	
			res_dat<-res$results
			res_dat$region<-gene_tab$gene[j]
			res_dat$max<-max(Markers$pos)
			res_dat$min<-min(Markers$pos)
			res_dat$bp<-basepairs	
			res_dat<-get_nice_names(dat=res_dat,Trait_names=Trait_names)
			res_list[[j]]<-res_dat
		}		
	}	
	if(is.null(region)){
		region1<-"1mb"
	}else{
		region1<-paste(region,"bp",sep="")
	}
	if(!fix_charge){
		return(list(res_list,coloc_strategy,region1))
	}
	if(fix_charge){
		return("colocalisation results contaminated by fake data for CHARGE. Set fix_charge to FALSE to get colocalisation results")
	}

}


# stack_assoc_plot(markers=Markers, corr=ld.matrix, z=Z.matrix,traits=c(unique(Trait1$trait),unique(Trait2$trait)),top.marker="rs968567")

# Row<-which(row.names(ld.matrix) == "rs61896141" )
# Col<-which(colnames(ld.matrix) == "rs174535")
# r.corr<-ld.matrix[Row,Col]
# r2<-r.corr^2
# rs174546
# rs174547 (top SNP ARA in Guan)
# rs102275 (top SNP after adjustment for rs174547)

format_data1<-function(Fa.tab=NULL,ref.dat=NULL,fatty_acids=FALSE,gene="",gene_tab=NULL,tissues=NULL,gtex=FALSE,eqtlgen=FALSE,gwis_test=FALSE){
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

	if(gtex){
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

	if(eqtlgen){

		Fa.tab$Ens<-Fa.tab$gene
		Fa.tab<-Fa.tab[,names(Fa.tab) != "gene"]
		# if(is.null(gene)){
		# 	Test<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="ensembl_gene_id",values=Ens ,mart=Mart) #unique(Fa.tab$gene_id)
		# }
		# if(!is.null(gene)){
		# 	Test<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="external_gene_name",values=gene,mart=Mart) #unique(Fa.tab$gene_id)
		# }
		Fa.tab.m<-merge(Fa.tab,gene_tab[,c("gene","ensembl_gene_id")],by.x="Ens",by.y="ensembl_gene_id")	
		Fa.tab.m<-Fa.tab.m[Fa.tab.m$gene ==gene,]
		# Fa.tab.m<-merge(Fa.tab,Test,by.x="Ens",by.y="ensembl_gene_id")
		if(nrow(Fa.tab.m)>0){
			Fa.tab.m$trait<-paste(Fa.tab.m$gene," expression in blood",sep="")
			Fa.tab.m$trait2<-paste(Fa.tab.m$gene," expression in blood (eQTLGen)"," median=",median(Fa.tab.m$n),min(Fa.tab.m$n),sep="")
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
	if(eqtlgen & gene == "FADS2"){ #max -log10 P values larger than 1000 and therefore don\t get plotted properly in regional association plots
		Fa.tab.m$z<-Fa.tab.m$z/2
	}
	
	Fa.tab.m<-Fa.tab.m[which(Fa.tab.m$z!="NaN"),]
	
	if(gtex){
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
		return(Fa.tab.m2)
	}
	if(nrow(Fa.tab.m)==0){
		return("no data")
	}
}


# dat<-list(cancer1)
format_cancer<-function(){
	# code for situation where don't know how many cancers are being included or don't want to make assumptions about numbers of cancers
	# length(Cancer_objects)
	# eval(Cancer_objects[1])
	# Cancers<-lapply(1:length(Cancer_objects),FUN=function(x)
	#     	eval(parse(text=Cancer_objects[x])))
	# names(Cancers)<-paste("Cancer",1:length(Cancers),sep="")
	# list2env(Cancers,envir=.GlobalEnv)

	# assuming only 2 cancers present
	
	cancer1<-eval(parse(text=Cancer_objects[1]))
	cancer2<-eval(parse(text=Cancer_objects[2]))
	cancer1$trait<-"Colorectal cancer"
	names(cancer1)[names(cancer1) == "rs_number"]<-"marker"
	names(cancer1)[names(cancer1) == "Freq1"]<-"eaf"
	names(cancer1)[names(cancer1) == "Allele1"]<-"effect_allele"
	names(cancer1)[names(cancer1) == "Allele2"]<-"other_allele"
	names(cancer1)[names(cancer1) == "Effect"]<-"beta"
	names(cancer1)[names(cancer1) == "StdErr"]<-"se"
	names(cancer1)[names(cancer1) == "P.value"]<-"pval"
	# snps<-unique(cancer1$marker[duplicated(cancer1$marker)])
	# cancer1[cancer1$marker %in% snps,]
	cancer1<-cancer1[!duplicated(cancer1$marker),] #duplicated SNPs seem to correspond to same base pair position but with different alleles. Perhaps these are non biallelic SNPs / ie. are triallelic
	
	cancer2$trait<-"Lung cancer"
	names(cancer2)[names(cancer2) == "effect_allele_freq"]<-"eaf"
	names(cancer2)[names(cancer2) == "p"]<-"pval"
	names(cancer2)[names(cancer2) == "snp"]<-"marker"
	cancer2<-cancer2[cancer2$file == "/projects/MRC-IEU/research/data/evidencehub/summary/gwas/dev/release_candidate/data/results/mrbase/passqc/Onco_TRICL_032116_Overall.csv.tab",]
	if(any(duplicated(cancer2$marker))) stop("duplicated SNPs present")
	cancer_dat<-rbind.fill(cancer1,cancer2)

	cancer_dat$z<-cancer_dat$beta/cancer_dat$se
	cancer_dat$trait2<-cancer_dat$trait
	cancer_dat<-merge(cancer_dat,ref3,by.x="marker",by.y="V2")
	names(cancer_dat)[names(cancer_dat) == "V1"]<-"chr"
	names(cancer_dat)[names(cancer_dat) == "V4"]<-"pos"
	cancer_dat<-cancer_dat[order(cancer_dat$chr,cancer_dat$pos),]
	return(cancer_dat)
}

# unique(Traits2$marker) == Traits1$marker
# cancer_dat[cancer_dat$marker == "rs2524299",]
# Temp[Temp$marker=="rs2524299",]
format_data2<-function(trait_list=NULL,ld.matrix=NULL,snps=NULL,Studies=NULL,coloc_strategy=NULL,study=NULL){
		Temp<-do.call(rbind.fill,trait_list)
		N<-length(unique(Temp$trait2))
		Temp<-Temp[which(Temp$z!="NaN"),]
		if(any(is.na(Temp$z))) stop("warning Z scores missing")
		Table_traits<-table(Temp$marker)
		Names_table<-names(Table_traits) 
		SNPs_keep<-Names_table[Table_traits==N]
		Temp<-Temp[Temp$marker %in% SNPs_keep,]
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

		Pos.eqtlgen<-grep("eQTLGen",Temp$trait2)
		if(sum(Pos.eqtlgen)!=0){
			Min<-min(Temp$n[Pos.eqtlgen])			
			Max<-max(Temp$n[Pos.eqtlgen])
			Med<-median(Temp$n[Pos.eqtlgen])

			Temp$trait2[Pos.eqtlgen]<-paste(unique(Temp$trait[Pos.eqtlgen]),paste0("in eQTLGen", " \n(Median N=",Med,", min=",Min,", max=",Max,")"))
		}
 		# order data so that fatty acids appear in alphabetical order above the eQTLs
		fa_eqtl_traits<-unique(Temp$trait)
		fa_traits<-fa_eqtl_traits[grep(":",fa_eqtl_traits)]
		
		eqtl_traits<-fa_eqtl_traits[grep(":",fa_eqtl_traits,invert=T)]
		eqtl_traits<-eqtl_traits[grep("cancer",eqtl_traits,invert=T)]
		fa_traits<-fa_traits[order(fa_traits)]
		cancer_traits<-fa_eqtl_traits[grep("cancer",fa_eqtl_traits)]
		
		if(coloc_strategy=="1trait_allstudies_alltissues"){
			Temp$trait[which(Temp$trait == fa_traits)]<-Temp$trait2[which(Temp$trait == fa_traits)]
			fa_traits<-paste(fa_traits," (",Studies[order(Studies)],")",sep="")
		}
		
		if(coloc_strategy=="everything"){
			# Temp$trait[which(Temp$trait %in% fa_traits)]<-Temp$trait2[which(Temp$trait %in% fa_traits)]
			fa_traits<-unique(Temp$trait2[which(Temp$trait %in% fa_traits)])
			# fa_traits<-paste(fa_traits," (",Studies[order(Studies)],")",sep="")
		}

		eqtl_traits<-eqtl_traits[order(eqtl_traits)]
		eqtl_traits1<-eqtl_traits[grep("expression in blood",eqtl_traits,invert=TRUE)]
		eqtl_traits2<-eqtl_traits[grep("expression in blood",eqtl_traits)]
		eqtl_traits<-c(eqtl_traits1,eqtl_traits2)
		fa_eqtl_traits<-c(fa_traits,eqtl_traits)
		if(coloc_strategy ==  "allcancer_1trait_1study" | coloc_strategy == "allcancer_1trait_1study_notissues" ){
			fa_eqtl_traits<-c(fa_traits,eqtl_traits,cancer_traits)
		}

		Num<-1:length(fa_eqtl_traits)
		dat_temp<-data.frame(matrix(c(fa_eqtl_traits,Num),nrow=length(Num),ncol=2),stringsAsFactors=F)
		dat_temp$X2<-as.numeric(dat_temp$X2)
		dat_temp$X2[dat_temp$X1 == "AA:DGLA / D5D"]<-1
		dat_temp$X2[dat_temp$X1 == "GLA:LA / D6D"]<-2
		Test<-sum(grep("deleteme",dat_temp$X1))!=0
		if(Test){
			dat_temp$X2[dat_temp$X1 %in% c("DPAn6:AA_deleteme / ELOVL2","DHA:DPAn3 / ELOVL2")]<-c(3,4)
		}
		if(coloc_strategy=="alltraits_1study_alltissues" & !Test){
			if(study != "CHARGE"){
			dat_temp$X2[dat_temp$X1 %in% c("DPAn6:AA / ELOVL2","DHA:DPAn3 / ELOVL2")]<-c(3,4)
			}
			if(study == "CHARGE"){
				dat_temp$X2[dat_temp$X1 == "DHA:DPAn3 / ELOVL2"]<-3

			}
		}
		dat_temp$X2[dat_temp$X1 %in% c("ADA:AA / ELOVL2/5","DPAn3:EPA / ELOVL2/5")]<-c(5,6)
		dat_temp$X2[dat_temp$X1 =="DGLA:GLA / ELOVL5"]<-7
		dat_temp$X2[dat_temp$X1 %in% c("POA:PA / SCD","OA:SA / SCD")]<-c(8,9)
		dat_temp<-dat_temp[order(dat_temp$X2),]

		if(coloc_strategy=="everything"){
			Pos<-unlist(lapply(c("CHARGE","Framingham","Shin"),FUN=function(x)
				grep(x,dat_temp$X1)))
			dat_temp1<-dat_temp[Pos,]
			Pos<-unlist(lapply(c("AA:DGLA","GLA:LA","DPAn6:AA","DHA:DPAn3","ADA:AA","DPAn3:EPA","DGLA:GLA","POA:PA","OA:SA"),FUN=function(x)
				grep(x,dat_temp1$X1)))
			dat_temp2<-dat_temp1[Pos,]
			Pos<-unlist(lapply(eqtl_traits,FUN=function(x)
				grep(x,dat_temp$X1)))
			dat_temp<-rbind(dat_temp2,dat_temp[Pos,])
		}

		dat_temp$X2<-1:nrow(dat_temp)

		# "AA:DGLA / D5D","GLA:LA / D6D","DPAn6:AA","DHA:DPAn3","DPAn6:AA"
		if(coloc_strategy=="everything"){
			Pos<-unlist(lapply(c("CHARGE","Framingham","Shin"),FUN=function(x)
				grep(x,Temp$trait2)))
			Temp$trait[Pos]<-Temp$trait2[Pos]
			Temp<-merge(Temp,dat_temp,by.x="trait",by.y="X1")
		}else{			
			Temp<-merge(Temp,dat_temp,by.x="trait",by.y="X1")
		}
		Temp$trait2<-gsub("_"," ",Temp$trait2)
		
		# any(duplicated(paste(Temp$trait2,Temp$marker)))
		# table(Temp$marker)
	
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
		
		ld.matrix2<-ld.matrix[row.names(ld.matrix) %in% snporder,colnames(ld.matrix) %in% snporder]
		ld.matrix2<-ld.matrix2[match(snporder,row.names(ld.matrix2)),match(snporder,colnames(ld.matrix2))]

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
		Z.matrix<-matrix(unlist(Z.list),ncol=length(Z.list),nrow=length(SNPs_keep))
		Markers<- trait1[,c("marker","chr","pos")]
		# Markers<-unique(Temp[,c("marker","chr","pos")]) #this doesn't work because not all studies have used same build
		
		if(any(row.names(ld.matrix2) != Markers$marker)) stop("marker order mismatched")
		Names.list<-unique(unlist(lapply(1:length(Temp1),FUN=function(x)
	    	eval(parse(text=paste("trait",x,"$trait2",sep=""))))))

		return(list(ld.matrix2,Z.matrix,Markers,Names.list,b.list,se.list))
}


format_gtex<-function(Dat=NULL){
	Dat$study<-"GTEx"
	Ens<-Dat$gene_id
	Ens<-unlist(strsplit(Ens,split="\\."))
	Dat$Ens<-Ens[seq(1,length(Ens),by=2)]
	return(Dat)
}


format_fastudies<-function(Dat=NULL){
	File<-c("AA_to_DGLA.tab","ADA_to_AA.tab","DGLA_to_GLA.tab","DHA_to_DPA_n3.tab","DPA_n3_to_EPA.tab","GLA_to_LA.tab","OA_to_SA.tab","POA_to_PA.tab","DPA_n6_to_AA.tab" )
	Trait<-c("AA:DGLA","ADA:AA","DGLA:GLA","DHA:DPAn3","DPAn3:EPA","GLA:LA","OA:SA","POA:PA","DPAn6:AA" )
	Trait_file<-data.frame(matrix(c(File,Trait),ncol=2,nrow=length(File)),stringsAsFactors=F)
	names(Trait_file)<-c("file","trait")
	Dat$file<-gsub(".gz","",Dat$file)
	Dat2<-merge(Dat,Trait_file,by="file")
	return(Dat2)
}

define_region<-function(index_snp=NULL,region=NULL,ref=NULL){
	SNPlist<-readLines(paste("~/fatty-acids-mr/coloc/",index_snp,".txt",sep=""))
	Gen.table<-read.table(paste("~/fatty-acids-mr/coloc/",index_snp,"_r_matrix.ld",sep=""),sep="\t",head=F,stringsAsFactors=F)	
	row.names(Gen.table)<-SNPlist
	colnames(Gen.table)<-SNPlist
	ref3<-ref[ref$V2 %in% SNPlist,]
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

make_title<-function(Dir=NULL,Info=NULL,gwis_file=NULL,type=".png"){
	Info<-paste(unique(unlist(lapply(1:length(strsplit(Info,split=" ")),FUN=function(x) 
		unlist(strsplit(Info,split=" ")[x])[1]))),collapse="_")
	Info<-gsub("_NA","",Info)
	if(!is.null(gwis_file)){
		grep("imputed",gwis_file)
		Info<-paste0(Info,"_imputed")
	}

	plot.title<-paste(Dir,paste(Info,collapse="_"),type,sep="")
	plot.title<-gsub(":","to",plot.title)
	if(coloc_strategy == "allcancer_1trait_1study"){
		plot.title<-gsub(".png","_eQTLs_cancer.png",plot.title)
	}
	if(coloc_strategy == "allcancer_1trait_1study_notissues"){
		plot.title<-gsub(".png","_cancer.png",plot.title)
	}

	# plot.title<-gsub(" ","_",plot.title)
	# plot.title<-gsub("\\(","",plot.title)
	# plot.title<-gsub(")","",plot.title)
	
	return(plot.title)
	# pdf.title<-paste("~/fatty-acids-mr/coloc/AA_to_DGLA_fads1_liver.pdf",sep="")
				# pdf(pdf.title, onefile=FALSE)
}

make_title2<-function(Dir=NULL,Info=NULL,type=".png"){
	if(length(Info) == 2){
		file.name<-paste(Dir,paste(Info,collapse=""),type,sep="")
	}
	# if(length(Info) == 3){
	return(file.name)
}

make_gene_table<-function(snps=NULL,genes=NULL){
	gene_tab<-data.frame(matrix(c(snps,genes),ncol=2,nrow=4),stringsAsFactors=F)
	names(gene_tab)<-c("SNP","gene")
	Ens_gene<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="external_gene_name",values=gene_tab$gene,mart=Mart) #unique(Fa.tab$gene_id)
	gene_tab<-merge(gene_tab,Ens_gene,by.x="gene",by.y="external_gene_name")
	return(gene_tab)
}

get_nice_names<-function(dat=NULL,Trait_names=NULL){
	# traits colocalise
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


harmonise_dat<-function(Dat=NULL,ref_dat=NULL){
	names(ref_dat)[names(ref_dat) == "effect_allele"]<-"ref_ea"
	names(ref_dat)[names(ref_dat) == "other_allele"]<-"ref_oa"
	names(ref_dat)[names(ref_dat) == "eaf"]<-"ref_eaf"
	Dat<-merge(Dat,ref_dat,by="SNP")
	Alleles<-paste(Dat$effect_allele,Dat$other_allele,sep="")
	Pos<-Alleles %in% c("GC","CG","TA","AT")
	Dat2<-Dat[Pos,] #palindromic SNPs
	Dat3<-Dat[!Pos,] #non palindromic SNPs
	
	# harmonise non palindromic SNPs
	ea<-Dat3$effect_allele
	oa<-Dat3$other_allele
	eaf<-Dat3$eaf
	beta<-Dat3$beta
	if(any(ea != Dat3$ref_ea & ea != Dat3$ref_oa)) stop("some SNPs are on different strands to the reference")
	Pos<-which(ea != Dat3$ref_ea) #positions where effect allele is different from effect allele in reference set
	Dat3$effect_allele[Pos]<-oa[Pos]
	Dat3$other_allele[Pos]<-ea[Pos]
	Dat3$eaf[Pos]<-1-eaf[Pos]
	Dat3$beta[Pos]<-beta[Pos]*-1
	all(Dat3$effect_allele == Dat3$ref_ea)
	
	# deal with palindromic SNPs
	# check that palindromic SNPs "appear to be on same strand". 
	ea<-Dat2$effect_allele
	if(any(ea != Dat2$ref_ea & ea != Dat2$ref_oa)) warning("some palindromic SNPs appear to be on different strands to the reference, indicating allele errors. these SNPs will be dropped")
	Pos<-ea != Dat2$ref_ea & ea != Dat2$ref_oa
	Dat2_1<-Dat2[!Pos,]# drop SNPs that are not palindromic in both the test and reference datasets. E.g. one rs11231053 is C/G in CHARGE GLA:LA but is A/G UK biobank bim file and dbSNP.  
	# ukb[which(ukb$V2 == "rs11231053"),]

	# code all palindromic SNPs so that effect allele is the minor allele. The effect allele is always ≤0.5 in the reference dataset (i.e. CHS LD file). Drop palindromic SNPs with MAF close to 0.5 defined as >0.42
	if(any(ref_dat$eaf>0.5)) stop("some effect alleles are not minor allele in reference dataset")
	if(any(Dat2_1$ref_eaf>0.5)) stop("some effect alleles of palindromic SNPs are not the minor allele")
	ea<-Dat2_1$effect_allele
	oa<-Dat2_1$other_allele
	eaf<-Dat2_1$eaf
	beta<-Dat2_1$beta

	Pos<-eaf>0.5
	Dat2_1$eaf[Pos]<-1-eaf[Pos]
	Dat2_1$beta[Pos]<-beta[Pos]*-1
	Dat2_1$effect_allele[Pos]<-oa[Pos]
	Dat2_1$other_allele[Pos]<-ea[Pos]
	Dat2_1[,c("eaf","ref_eaf")]
	# exclude palindromic SNPs with MAF >0.42
	Pos<-Dat2_1$ref_eaf>=0.43 | Dat2_1$eaf>=0.43
	# Dat2_1[which(Pos),]
	# Lun3[Lun3$SNP %in% c("rs7112985","rs7946441"),]
	# Crc2[Crc2$SNP %in% c("rs7112985","rs7946441"),]
	# chs2[chs2$SNP %in% c("rs7112985","rs7946441"),]
	# d5d2[d5d2$SNP %in% c("rs7112985","rs7946441"),]
	# d6d2[d6d2$SNP %in% c("rs7112985","rs7946441"),]
	Dat2_2<-Dat2_1[!Pos,]

	# join non palindromic and palindromic SNPs together
	Har<-rbind(Dat2_2,Dat3)
	# make sure order of SNPs in Har is same as in ref_dat
	Pos<-match(ref_dat$SNP,Har$SNP)
	Pos<-Pos[!is.na(Pos)]
	Har<-Har[Pos,]
	all(ref_dat$SNP[ref_dat$SNP %in% Har$SNP] == Har$SNP)
	return(Har)
}






# res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, bb.selection = "align"); might make more sense because AA:DGLA might not colocalaise with FADS1 because it has a different causal variant but all the traits definitely have a causal variant in region somewhere. On other hand expression in different tissues might no have causal variant in regin. 
			# The regional selection criterion is computed from a collection of hypotheses which assume that all traits do not colocalize because one of the traits does not have a causal variant in the region. The alignment selection criterion, however, is computed from hypotheses which assume that all traits do not colocalize because one of the traits has a causal variant elsewhere in the region.


