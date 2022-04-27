# source("~/fatty-acids/colocalisation/scripts/regional_association_plots_functions.R")
data_list<-load_data(
	bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata", 
	ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  
	ref_turn_off=FALSE)

gene_tab1<-data.frame(data_list[1],stringsAsFactors=F)
ref1<-data.frame(data_list[2],stringsAsFactors=F)
fa.tab1<-data.frame(data_list[3],stringsAsFactors=F)
bbj_eqtl_data1<-data.frame(data_list[4],stringsAsFactors=F)


# tissues=c("colon sigmoid","whole blood","lung","adipose subcutaneous","artery coronary","liver","colon transverse","adipose visceral omentum")

gene<-"FADS1"
study<-"Dorajoo/SCHS"
snp<-gene_tab$SNP[gene_tab$gene == gene]
Region_object_list<-define_region(index_snp=snp,region=region,ref=ref) #region refers to how many base pairs to plot . region=NULL plots 1million base pairs (the default). 
ref3<-data.frame(Region_object_list[1],stringsAsFactors=F)
Gen.table<-data.frame(Region_object_list[2],stringsAsFactors=F)
SNPlist<-unlist(Region_object_list[3])

fa.tab<-fa.tab[fa.tab$study == study,]				
fa.tab2<-format_fastudies(Dat=fa.tab)

Traits1<-format_data1(Fa.tab=fa.tab2,ref.dat=ref3,gwis_test=TRUE,gwis_file=gwis_file,study=study,gene=gene)
Traits2<-format_data1(Fa.tab=bbj_eqtl_data,ref.dat=ref3,gene=gene,gene_tab=gene_tab,bbj_eqtl_test=TRUE)

trait_list<-ls()[grep("Traits[0-9]",ls())] 

Test<-unlist(lapply(1:length(trait_list),FUN=function(x) 
	!is.null(nrow(eval(parse(text=trait_list[x]))))))
rm(list=trait_list[!Test]) #remove the objects with no data
trait_list<-ls()[grep("Traits[0-9]",ls())] 

trait_list2<-lapply(1:length(trait_list),FUN=function(x) eval(parse(text=trait_list[x])))
Plot<-format_data2(trait_list=trait_list2,ld.matrix=Gen.table,snps=SNPlist,plot_strategy=plot_strategy,study=study)
ld.matrix<-data.frame(Plot[1])
Z.matrix<-data.frame(Plot[2])
Markers<-data.frame(Plot[3])

snplist<-paste(Markers$marker,collapse="\n")

paste0("curl -k -H \\"Content-Type: application/json\\" -X POST -d \\'{\\"snps\\": \\"rs3\nrs4\\", \\"pop\\": \\"JPT\\",\\"r2_d\\": \\"d\\"}\\' \\'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?token=a086cd05a12d\\'")

paste("curl -k -H "Content-Type:") 
	application/json\" -X POST -d \'{\"snps\": \"rs3\nrs4\", \"pop\": \"JPT\",\"r2_d\": \"d\"}\' \'https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix?token=a086cd05a12d\'")


load_data<-function(gwis_file=NULL,ref_dat=NULL,ref_turn_off=FALSE,bbj_eqtl_file=NULL,UseBioMart=FALSE){
	if(UseBioMart){
		library(biomaRt)
		Mart <- useMart(host="www.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
		Attr<-listAttributes(Mart)
	}

	ChargeEA1$study <- "ChargeEA"
	Dor4$study <- "Dorajoo/SCHS"
	Dor1<-format_schs(Dat=Dor4)
	fa.tab<-do.call(rbind,list(ChargeEA1,Dor1))	
	
	bbj_eqtl_data<-load_bbjeqtl(File=bbj_eqtl_file)

	if(!ref_turn_off){
		ref <- read.table(ref_dat,stringsAsFactors=F,head=F)
	}
	
	gene_tab<-make_gene_table(snps=c("rs174528","rs174528","rs3734398","rs603424"),genes=c("FADS1","FADS2","ELOVL2","SCD"),Mart=NULL ) #function assumes that these are the SNPs and genes of interest
# =Mart
	return(list(gene_tab,ref,fa.tab,bbj_eqtl_data))
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
		genes<-c("ELOVL2","FADS1","FADS2","SCD")
		snps<-c("rs3734398","rs174528","rs174528","rs603424")
		ens<-c("ENSG00000197977" ,"ENSG00000149485","ENSG00000134824","ENSG00000099194")
		gene_tab<-data.frame(do.call(cbind,list(genes,snps,ens)),stringsAsFactors=F)
		names(gene_tab)<-c("gene","SNP","ensembl_gene_id")
	}
	return(gene_tab)
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

	dim(unique(Trait_file))
	Dat$file<-gsub(".gz","",Dat$file)
	Dat2<-merge(Dat,Trait_file,by="file")
	return(Dat2)
}
