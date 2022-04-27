# ALT is effect allele in gtex
# maf not necessarily eaf
#  "plus genomic strand" is reference strand

rm(list=ls())
# source("~/fatty-acids-mr/mr/plot_mr_results_functions.R")
library(ggforestplot)
library(ggplot2)
# install.packages("meta")
library(meta)
library(plyr)

snp_tab<-make_snp_tab()
gene_tab<-make_gene_tab()


# ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz"
# ref1 <- read.table(ref_dat,stringsAsFactors=F,head=F)

gtex_dat<-format_gtex()
eqtlgen_dat<-format_eqtlgen()
eqtl_bbj<-format_bbj()
fa_tab_eur<-format_fa_eur()
fa_tab_eas<-format_fa_eas()
 # fa_tab_eas
plot_dat<-do.call(rbind.fill,list(eqtlgen_dat,gtex_dat,fa_tab_eur))
plot_eas<-do.call(rbind.fill,list(eqtl_bbj,fa_tab_eas))

plot_eas2<-format_data_eas(Dat=plot_eas)
unique(plot_eas2$SNP)
forestplot(df = plot_eas2,
			 logodds = FALSE,
			 name=tissue,
						  estimate=b_sd,
						  se=se_sd,
						  shape=SNP,
						  colour = gene,
						   xlab = "")

plot_dat2<-format_data(Dat=plot_dat,type=FALSE,adjust=TRUE)
plot_dat2[plot_dat2$SNP =="rs174546" ,]
plot_dat[plot_dat$SNP =="rs174546" & plot_dat$trait == "GLA:LA" ,c("SNP","effect_allele","trait","beta","se",tissue","file")]

C allele rs17 higher GLA:L"
unique(plot_dat2$SNP)
plot_dat2$tissue<-paste0("gene expression\n",plot_dat2$tissue)

plot_dat2$tissue[grep("CHARGE",plot_dat2$tissue)]<-"enzyme activity\nblood (CHARGE)"

plot<-forestplot(df = plot_dat2,
			 logodds = FALSE,
			 name=tissue,
						  estimate=b_sd,
						  se=se_sd,
						  # shape="",
						  colour = gene,
						   xlab = "")
pdf("~/fatty-acids/mr/results/plots/fads_instruments_pub_v1.pdf")
	plot
dev.off()

plot_dat3<-format_data(Dat=plot_dat,snpFADS=FALSE, gene_expression=FALSE,Charge=FALSE,adjust=FALSE)
plot_dat3$study[grep("adj",plot_dat3$file)]<-"CHARGE\n(adjusted)"
plot_dat3$study[plot_dat3$study == "CHARGE"]<-"CHARGE\n(unadjusted)"
plot_dat3$study[plot_dat3$study == "Framingham"]<-"Framingham\n(unadjusted)"
plot_dat3$study[plot_dat3$study == "Shin"]<-"Shin\n(unadjusted)"
plot_dat3a <- plot_dat3[plot_dat3$study == "CHARGE\n(adjusted)",]
plot_dat3b <- plot_dat3[plot_dat3$study != "CHARGE\n(adjusted)",]
plot_dat3<-rbind(plot_dat3a,plot_dat3b)

plot<-forestplot(df = plot_dat3,
			 logodds = FALSE,
			 name=study,
						  estimate=b_sd,
						  se=se_sd,
						  shape=SNP,
						  colour = gene,
						   xlab = "")
pdf("~/fatty-acids/mr/results/plots/fads_instruments_cfs_v1.pdf")
	plot
dev.off()

plot_dat4<-format_data(Dat=plot_dat,type=FALSE,adjust=FALSE,snpFADS=FALSE)

head(plot_dat4$trait2)

Pos<-plot_dat4$adj == "adjusted"
plot_dat4$trait2[Pos]<-paste0(plot_dat4$trait2[Pos]," (adjusted)")
unique(plot_dat4$trait2)
plot_dat4[plot_dat4$adj=="adjusted",]
plot_dat4$snp_trait
plot<-forestplot(df = plot_dat3,
			 logodds = FALSE,
			 name=study,
						  estimate=b_sd,
						  se=se_sd,
						  shape=SNP,
						  colour = gene,
						   xlab = "")

pdf("~/fatty-acids/mr/results/plots/fads_instruments_v11.pdf")
	plot
dev.off()




				 # labs(title=Title.plot,size=1)+
				 # theme(plot.title = element_text(size = text.title))+
				 theme(text = element_text(size=text.names))


make_gene_tab<-function(){
	genes<-c("FADS1","FADS2")
	trait<-c("AA:DGLA","GLA:LA")
	ens<-c("ENSG00000149485","ENSG00000134824")
	gene_tab<-data.frame(do.call(cbind,list(genes,ens,trait)),stringsAsFactors=F)
	names(gene_tab)<-c("gene","ensembl_gene_id","trait")
	return(gene_tab)
}   
   # gene       SNP ensembl_gene_id

# 1 ELOVL2 rs3734398 ENSG00000197977
# 2  FADS1  rs174528 ENSG00000197977
# 3  FADS2  rs174528 ENSG00000134824
# 4    SCD  rs603424 ENSG00000099194


make_snp_tab<-function(){
	# build 19 / GCRh 37
	snps<-c("rs61896141","rs174546","rs968567","rs2727271","rs2524299")
	pos<-c(61556039,61569830,61595564,61603358,61604782)
	chr<-rep(11,length(snps))
	snp_tab<-data.frame(do.call(cbind,list(snps,pos,chr)),stringsAsFactors=F)
	names(snp_tab)<-c("SNP","pos","chr")
	return(snp_tab)
}


make_snp_eas<-function(){
	# build 19 / GCRh 37
	# snps<-c("rs61896141","rs174546","rs968567","rs2727271","rs2524299")
	snps<-c("rs174584","rs174546","rs174572")
	# pos<-c(61556039,61569830,61595564,61603358,61604782)
	chr<-rep(11,length(snps))
	snp_tab<-data.frame(do.call(cbind,list(snps,chr)),stringsAsFactors=F)
	# "pos",
	names(snp_tab)<-c("SNP","chr")
	return(snp_tab)
}


format_eqtlgen<-function(){
	eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt"
	Dat<-read.table(eqtlgen_file,sep="\t",stringsAsFactors=F,head=T)
	Dat$tissue<-"blood"
	Dat$study<-"eQTLGen"
	gene_tab<-make_gene_tab()
	names(Dat)[names(Dat)=="gene"]<-"ensembl_gene_id"
	Dat<-merge(Dat,gene_tab,by="ensembl_gene_id")
	Dat<-Dat[,c("snp","gene","effect_allele","other_allele","effect_allele_freq","beta","se","p","chr","pos","tissue","study","n","ensembl_gene_id")]
	snp_tab<-make_snp_tab()
	names(Dat)[names(Dat) == "snp"]<-"SNP"
	names(Dat)[names(Dat) == "effect_allele_freq"]<-"eaf"
	Dat<-Dat[Dat$SNP %in% snp_tab$SNP, ]
	names(Dat)[names(Dat) == 'p']<-"pval"
	names(Dat)[names(Dat) == "n"]<-"samplesize"
	Dat$population<-"European"
	Dat$trait <- Dat$gene
	Dat$type<-"gene expression"
	return(Dat)
}


format_gtex<-function(){
	gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata"
	load(gtex_file)
	Dat<-df4
	tissues=c("colon sigmoid","whole blood","lung","adipose subcutaneous","liver","colon transverse","adipose visceral omentum")
	gene_tab<-make_gene_tab()
	snp_tab<-make_snp_tab()
	Dat<-merge(Dat,snp_tab,by="SNP")
	Ens<-Dat$gene_id
	Ens<-unlist(strsplit(Ens,split="\\."))
	Dat$ensembl_gene_id<-Ens[seq(1,length(Ens),by=2)]
	Dat<-merge(Dat,gene_tab,by="ensembl_gene_id")
	Dat<-Dat[Dat$tissue %in% tissues,c("SNP","slope","slope_se","tissue","ensembl_gene_id","maf","variant_id","gene")]
	Dat$tissue[Dat$tissue == "whole blood"]<-"blood"
	alleles<-unlist(strsplit(Dat$variant_id,split="_"))
	Dat$ref<-alleles[seq(3,length(alleles),by=5)]
	Dat$alt<-alleles[seq(4,length(alleles),by=5)]
	names(Dat)[names(Dat) == "alt"]<-"effect_allele"
	names(Dat)[names(Dat) == "ref"]<-"other_allele"
	names(Dat)[names(Dat) == "slope"]<-"beta"
	names(Dat)[names(Dat) == "slope_se"]<-"se"
	Dat$eaf<-NA # eaf not reported by gtex but alleles reported for + strand
	Dat$study<-"GTEx"
	Dat$trait <- Dat$gene
	Dat$population<-"European"
	Dat$samplesize<-948 
	Dat$type<-"gene expression"
	Dat$tissue[Dat$tissue == "adipose subcutaneous"]<-"adipose \nsubcutaneous"
	Dat$tissue[Dat$tissue == "adipose visceral omentum"]<-"adipose \nvisceral omentum"
	return(Dat)
}

format_fa_eur<-function(
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata"){	
	load(gwis_file)
	Charge1$study <- "CHARGE"
	Shin1$study <- "Shin"
	Tin1$study <- "Framingham"
	fa.tab<-do.call(rbind,list(Charge1,Shin1,Tin1))
	snp_tab<-make_snp_tab()
	fa.tab<-merge(fa.tab,snp_tab,by="SNP")
	fa.tab$trait<-NA
	fa.tab$trait[grep("AA_to_DGLA",fa.tab$file)]<-"AA:DGLA"
	fa.tab$trait[grep("GLA_to_LA",fa.tab$file)]<-"GLA:LA"
	fa.tab<-fa.tab[!is.na(fa.tab$trait),]
	gene_tab<-make_gene_tab()
	fa.tab<-merge(fa.tab,gene_tab,by="trait")
	fa.tab$population<-"European"
	fa.tab$tissue<-"blood"
	fa.tab$type<-"enzyme activity" 
	# fa.tab$type[fa.tab$study == "CHARGE"] <-"enzyme activity (CHARGE)"
	# fa.tab$type[fa.tab$file == "AA_to_DGLA_adjSNP.tab" ] <-"enzyme activity (CHARGE \nadj rs174547)"
	# fa.tab$type[fa.tab$file == "GLA_to_LA_adjSNP.tab" ] <-"enzyme activity (CHARGE \nadj rs174547)"
	# fa.tab$type[fa.tab$study == "Framingham"] <-"enzyme activity (Fram)"
	# fa.tab$type[fa.tab$study == "Shin"] <-"enzyme activity (Shin)"
	# fa.tab$tissue[fa.tab$study %in% c("CHARGE","Framingham")]<-"whole blood"
	# fa.tab<-fa.tab[fa.tab$study == "CHARGE",]
	# fa.tab<-fa.tab[fa.tab$file %in% c("GLA_to_LA_adjSNP.tab", "AA_to_DGLA.tab"),]
	return(fa.tab)
}

format_fa_eas<-function(){
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata"
	load(gwis_file)
	ChargeEA1$study <- "ChargeEA"
	Dor4$study <- "Dorajoo/SCHS"
	Dor1<-format_schs(Dat=Dor4)
	fa.tab<-do.call(rbind,list(ChargeEA1,Dor1))	
	snp_tab<-make_snp_eas()	
	fa.tab<-merge(fa.tab,snp_tab,by="SNP")
	fa.tab$trait<-NA	
	fa.tab$trait[grep("AA_to_DGLA",fa.tab$file)]<-"AA:DGLA"
	fa.tab$trait[fa.tab$file=="score_lnD5D_pooled_allchr_qc1.tab"]<-"AA:DGLA"
	fa.tab$trait[fa.tab$file=="score_lnD6D_pooled_allchr_qc1.tab"]<-"GLA:LA"
	fa.tab$trait[grep("GLA_to_LA",fa.tab$file)]<-"GLA:LA"
	fa.tab[!is.na(fa.tab$trait),c("trait","file")]
	fa.tab<-fa.tab[!is.na(fa.tab$trait),]
	gene_tab<-make_gene_tab()
	fa.tab<-merge(fa.tab,gene_tab,by="trait")
	fa.tab$population<-"East Asian"
	fa.tab$tissue<-"blood"
	fa.tab<-fa.tab[fa.tab$file %in% c("score_lnD5D_pooled_allchr_qc1.tab","score_lnD6D_pooled_allchr_qc1.tab"),]
	return(fa.tab)
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

b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}

# Dat=plot_dat
format_data<-function(Dat=NULL,Charge=TRUE,gene_expression=TRUE,adjust=TRUE,type=TRUE,SNP=c("rs174546","rs968567"),snpFADS=TRUE){
	Dat<-Dat[Dat$SNP %in% SNP ,]
	Dat$maf[is.na(Dat$maf)]<-Dat$eaf[is.na(Dat$maf)]
	Pos<-Dat$maf>=0.5	
	Dat$maf[Pos]<-1-Dat$maf[Pos]
	Dat$z<-Dat$beta/Dat$se
	Dat$b_sd<-b_sd(z=Dat$z,maf=Dat$maf,n=Dat$samplesize)
	Dat$se_sd<-Dat$b_sd/Dat$z

	Dat.h<-harmonise_dat(Dat)
	Dat<-Dat.h

	if(snpFADS){
		Dat1<-Dat[Dat$SNP == "rs174546" & Dat$gene == "FADS1",]
		Dat2<-Dat[Dat$SNP == "rs968567" & Dat$gene == "FADS2",]
		Dat<-rbind(Dat1,Dat2)
		Dat$snpgene<-NA
		Dat$snpgene[Dat$gene == "FADS1"]<-"rs174546 FADS1"
		Dat$snpgene[Dat$gene == "FADS2"]<-"rs968567 FADS2"
	}

	if(Charge){
		Dat<-Dat[Dat$study!="Framingham",]
		Dat<-Dat[Dat$study!="Shin",]
	}

	if(!gene_expression){
		Dat<-Dat[Dat$study %in% c("CHARGE","Shin","Framingham"),]		
	}
	if(adjust){
		Dat1<-Dat[Dat$file!="AA_to_DGLA_adjSNP.tab" | is.na(Dat$file),]
		# if(type){
			# Dat$type[Dat$type=="enzyme activity (CHARGE \nadj rs174547)"]<-"enzyme activity (CHARGE)"
			Dat<-Dat[Dat$file!= "GLA_to_LA.tab" | is.na(Dat$file),]			
		
		# }
		# Dat<-Dat[Dat$file!="GLA_to_LA.tab" | is.na(Dat$file),]
	}
	
	Dat$snp_trait<-paste0(Dat$SNP," ",Dat$trait,"\n(",Dat$gene,")")
	Dat$adj<-"unadjusted"
	Dat$adj[grep("adj",Dat$file)]<-"adjusted"

	Dat$trait2<-Dat$trait
	Dat$trait2[Dat$trait2 == "AA:DGLA"]<-"FADS1 enzyme activity"
	Dat$trait2[Dat$trait2 == "GLA:LA"]<-"FADS2 enzyme activity"
	Dat$trait2[Dat$trait2 == "FADS1"]<-"FADS1 gene expression"
	Dat$trait2[Dat$trait2 == "FADS2"]<-"FADS2 gene expression"

	Dat$tissue[Dat$tissue == "blood" & Dat$study == "CHARGE"]<-"blood (CHARGE)"
	Dat$tissue[Dat$tissue == "blood" & Dat$study == "GTEx"]<-"blood (GTEx)"
	Dat$tissue[Dat$tissue == "blood" & Dat$study == "eQTLGen"]<-"blood (eQTLGen)"

	Dat1<-Dat[grep("blood",Dat$tissue),]
	Dat2<-Dat1[Dat1$study == "CHARGE",]
	Dat1<-Dat1[Dat1$study != "CHARGE",]
	Dat3<-Dat[grep("blood",Dat$tissue,invert=T),]
	Dat<-rbind(Dat2,Dat1)
	Dat<-rbind(Dat,Dat3)
	Dat$trait2<-Dat$trait
	Dat$trait<-paste0(Dat$gene," (",Dat$type,")")
	return(Dat)
}


harmonise_dat<-function(Dat=NULL){
	ref_dat<-Dat[!is.na(Dat$eaf),c("effect_allele","other_allele","eaf","SNP")]
	Dat[Dat$study == "BBJ",]
	ref_dat<-ref_dat[!duplicated(ref_dat$SNP),]
	ea<-ref_dat$effect_allele
	oa<-ref_dat$other_allele
	eaf<-ref_dat$eaf
	Pos<-which(ref_dat$eaf>0.5)
	ref_dat$effect_allele[Pos]<-oa[Pos]
	ref_dat$other_allele[Pos]<-ea[Pos]
	ref_dat$eaf[Pos]<-1-eaf[Pos]
	
	if(!all(ref_dat$eaf<0.5)) stop("eaf not minor allele")

	names(ref_dat)[names(ref_dat) == "effect_allele"]<-"ref_ea"
	names(ref_dat)[names(ref_dat) == "other_allele"]<-"ref_oa"
	names(ref_dat)[names(ref_dat) == "eaf"]<-"ref_eaf"
	Dat<-merge(Dat,ref_dat,by="SNP")
	Alleles<-paste(Dat$effect_allele,Dat$other_allele,sep="")
	Pos<-Alleles %in% c("GC","CG","TA","AT")
	if(any(Pos)) stop("palindromic SNPs present")
	Dat2<-Dat[Pos,] #palindromic SNPs
	Dat3<-Dat[!Pos,] #non palindromic SNPs
	
	# harmonise non palindromic SNPs
	ea<-Dat3$effect_allele
	oa<-Dat3$other_allele
	eaf<-Dat3$eaf
	beta<-Dat3$beta
	b_sd<-Dat3$b_sd
	if(any(ea != Dat3$ref_ea & ea != Dat3$ref_oa)) stop("some SNPs are on different strands to the reference")

	Pos<-which(ea != Dat3$ref_ea) #positions where effect allele is different from effect allele in reference set
	Dat3$effect_allele[Pos]<-oa[Pos]
	Dat3$other_allele[Pos]<-ea[Pos]
	Dat3$eaf[Pos]<-1-eaf[Pos]
	Dat3$beta[Pos]<-beta[Pos]*-1
	Dat3$b_sd[Pos]<-b_sd[Pos]*-1
	all(Dat3$effect_allele == Dat3$ref_ea)
	
	return(Dat3)
	# deal with palindromic SNPs
	# check that palindromic SNPs "appear to be on same strand". 
	# ea<-Dat2$effect_allele
	# if(any(ea != Dat2$ref_ea & ea != Dat2$ref_oa)) warning("some palindromic SNPs appear to be on different strands to the reference, indicating allele errors. these SNPs will be dropped")
	# Pos<-ea != Dat2$ref_ea & ea != Dat2$ref_oa
	# Dat2_1<-Dat2[!Pos,]# drop SNPs that are not palindromic in both the test and reference datasets. E.g. one rs11231053 is C/G in CHARGE GLA:LA but is A/G UK biobank bim file and dbSNP.  
	# # ukb[which(ukb$V2 == "rs11231053"),]

	# # code all palindromic SNPs so that effect allele is the minor allele. The effect allele is always ≤0.5 in the reference dataset (i.e. CHS LD file). Drop palindromic SNPs with MAF close to 0.5 defined as >0.42
	# if(any(ref_dat$eaf>0.5)) stop("some effect alleles are not minor allele in reference dataset")
	# if(any(Dat2_1$ref_eaf>0.5)) stop("some effect alleles of palindromic SNPs are not the minor allele")
	# ea<-Dat2_1$effect_allele
	# oa<-Dat2_1$other_allele
	# eaf<-Dat2_1$eaf
	# beta<-Dat2_1$beta

	# Pos<-eaf>0.5
	# Dat2_1$eaf[Pos]<-1-eaf[Pos]
	# Dat2_1$beta[Pos]<-beta[Pos]*-1
	# Dat2_1$effect_allele[Pos]<-oa[Pos]
	# Dat2_1$other_allele[Pos]<-ea[Pos]
	# Dat2_1[,c("eaf","ref_eaf")]
	# # exclude palindromic SNPs with MAF >0.42
	# Pos<-Dat2_1$ref_eaf>=0.43 | Dat2_1$eaf>=0.43
	# # Dat2_1[which(Pos),]
	# # Lun3[Lun3$SNP %in% c("rs7112985","rs7946441"),]
	# # Crc2[Crc2$SNP %in% c("rs7112985","rs7946441"),]
	# # chs2[chs2$SNP %in% c("rs7112985","rs7946441"),]
	# # d5d2[d5d2$SNP %in% c("rs7112985","rs7946441"),]
	# # d6d2[d6d2$SNP %in% c("rs7112985","rs7946441"),]
	# Dat2_2<-Dat2_1[!Pos,]

	# # join non palindromic and palindromic SNPs together
	# Har<-rbind(Dat2_2,Dat3)
	# # make sure order of SNPs in Har is same as in ref_dat
	# Pos<-match(ref_dat$SNP,Har$SNP)
	# Pos<-Pos[!is.na(Pos)]
	# Har<-Har[Pos,]
	# all(ref_dat$SNP[ref_dat$SNP %in% Har$SNP] == Har$SNP)
	# return(Har)
}

format_bbj<-function(){
	bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata"
	load(bbj_eqtl_file)
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
	
	names(bbj_eqtl)[names(bbj_eqtl) == "ALT"]<-"effect_allele"
	names(bbj_eqtl)[names(bbj_eqtl) == "REF"]<-"other_allele"
	bbj_eqtl$population<-"East Asia"
	snp_tab<-make_snp_eas()
	gene_tab<-make_gene_tab()
	names(bbj_eqtl)[names(bbj_eqtl)=="gene"]<-"ensembl_gene_id"
	names(bbj_eqtl)[names(bbj_eqtl)=="p.value"]<-"pval"
	bbj_eqtl<-merge(bbj_eqtl,gene_tab,by="ensembl_gene_id")
	bbj_eqtl<-bbj_eqtl[bbj_eqtl$SNP %in% snp_tab$SNP, ]
	bbj_eqtl$trait<-bbj_eqtl$gene
	# bbj_eqtl$tissue
	bbj_eqtl$samplesize<-NA
	bbj_eqtl$samplesize[bbj_eqtl$tissue=="CD4 T cells"]<-103 
	bbj_eqtl$samplesize[bbj_eqtl$tissue=="Blood"]<-98
	bbj_eqtl$samplesize[bbj_eqtl$tissue=="Monocytes"]<-105
	bbj_eqtl$samplesize[bbj_eqtl$tissue=="NK cells"]<-104
	bbj_eqtl$samplesize[bbj_eqtl$tissue=="CD8 T cells"]<-103
	bbj_eqtl$samplesize[bbj_eqtl$tissue=="B cells"]<-104
	# bbj_eqtl$samplesize[bbj_eqtl$tissue=="B cells"]<-104
	bbj_eqtl[bbj_eqtl$SNP=="rs968567",]
	if(sum(which(bbj_eqtl$SNP=="rs968567")) == 0) message("rs968567 not present in bbj")
	bbj_eqtl$eaf<-NA
	bbj_eqtl$eaf[which(bbj_eqtl$SNP=="rs174546" & bbj_eqtl$effect_allele == "T")]<-0.327
	bbj_eqtl$eaf[which(bbj_eqtl$SNP=="rs174584" & bbj_eqtl$effect_allele == "A")]<-0.341 

	# snp<-c("rs174546")
	# allele1<-c("C")
	# freq1<-c(0.673)
	# allele2<-c("T")
	# freq2<-c(0.327)
	

	return(bbj_eqtl)
}

# Dat<-plot_eas
format_data_eas<-function(Dat=NULL,SNP=c("rs174546","rs968567","rs174584")){
	Dat<-Dat[Dat$SNP %in% SNP ,]
	Dat$maf<-Dat$eaf
	Pos<-which(Dat$maf>=0.5	)
	Dat$maf[Pos]<-1-Dat$maf[Pos]
	Dat$z<-Dat$beta/Dat$se
	Dat$b_sd<-b_sd(z=Dat$z,maf=Dat$maf,n=Dat$samplesize)
	Dat$se_sd<-Dat$b_sd/Dat$z

	Dat.h<-harmonise_dat(Dat)
	# Dat.h[,c("effect_allele","b_sd","eaf")]
	# Dat[,c("effect_allele","b_sd","eaf")]
	Dat<-Dat.h

	# if(snpFADS){
	# 	Dat1<-Dat[Dat$SNP == "rs174546" & Dat$gene == "FADS1",]
	# 	Dat2<-Dat[Dat$SNP == "rs968567" & Dat$gene == "FADS2",]
	# 	Dat<-rbind(Dat1,Dat2)
	# 	Dat$snpgene<-NA
	# 	Dat$snpgene[Dat$gene == "FADS1"]<-"rs174546 FADS1"
	# 	Dat$snpgene[Dat$gene == "FADS2"]<-"rs968567 FADS2"
	# }

	# if(Charge){
	# 	Dat<-Dat[Dat$study!="Framingham",]
	# 	Dat<-Dat[Dat$study!="Shin",]
	# }

	# if(!gene_expression){
	# 	Dat<-Dat[Dat$study %in% c("CHARGE","Shin","Framingham"),]		
	# }
		
	Dat$snp_trait<-paste0(Dat$SNP," ",Dat$trait,"\n(",Dat$gene,")")
	Dat$trait2<-Dat$trait
	Dat$trait2[Dat$trait2 == "AA:DGLA"]<-"FADS1 enzyme activity"
	Dat$trait2[Dat$trait2 == "GLA:LA"]<-"FADS2 enzyme activity"
	Dat$trait2[Dat$trait2 == "FADS1"]<-"FADS1 gene expression"
	Dat$trait2[Dat$trait2 == "FADS2"]<-"FADS2 gene expression"

	Dat$tissue[Dat$tissue == "Blood"]<-"blood"
	Dat$tissue[Dat$tissue == "blood" & Dat$study == "Dorajoo/SCHS"]<-"blood (Dorajoo/SCHS)"
	Dat$tissue[Dat$tissue == "blood" & Dat$study == "BBJ"]<-"blood (BBJ)"
	
	Dat1<-Dat[grep("blood",Dat$tissue),]
	Dat2<-Dat1[Dat1$study == "Dorajoo/SCHS",]
	Dat1<-Dat1[Dat1$study != "Dorajoo/SCHS",]
	Dat3<-Dat[grep("blood",Dat$tissue,invert=T),]
	Dat<-rbind(Dat2,Dat1)
	Dat<-rbind(Dat,Dat3)
	Dat$trait2<-Dat$trait
	Dat$trait<-paste0(Dat$gene," (gene expression)")
	Dat$trait[Dat$trait2 %in% c("AA:DGLA","GLA:LA" )]<-paste0(Dat$gene[Dat$trait2 %in% c("AA:DGLA","GLA:LA" )]," (enzyme activity)")
	return(Dat)
}

# make_ref_dat<-function(){
# 	snp<-c("rs174546")
# 	allele1<-c("C")
# 	freq1<-c(0.673)
# 	allele2<-c("T")
# 	freq2<-c(0.327)
	

	# command line way to obtain sub population allele frequencyes
	# https://www.internationalgenome.org/category/variants/
# grep CEU integrated_call_samples.20101123.ALL.panel | cut -f1 > CEU.samples.list

# vcf-subset -c CEU.samples.list ALL.chr13.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz | fill-an-ac |
#     bgzip -c > CEU.chr13.phase1.vcf.gz
#     </pre>
# }	

# unique(plot_dat[plot_dat$study == "CHARGE",c("SNP","effect_allele","other_allele","eaf")])















