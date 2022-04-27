source("~/fatty-acids/colocalisation/scripts/regional_association_plots_functions.R")
source("~/fatty-acids/mr/scripts/mr_functions.R")

devtools::install_github("NightingaleHealth/ggforestplot")
library(plyr)
library(ggforestplot)
# library(ggplot2)

data_list<-load_data(
	cancer1="~/fatty-acids/colocalisation/data/tricl_snps_ukb.Rdata",
	cancer2="~/MR_FattyAcids/data/summary_data/colorectal_cancer/061119/crc_snps.Rdata",
	# bbj_eqtl_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_eqtlbbj.Rdata",
	gtex_file="~/fatty-acids/colocalisation/data/GTEx_Analysis_v8_eQTLs_GTEx_Analysis_v8_eQTL_all_associations_alltissues_allregions_ukbsnps.Rdata",
	eqtlgen_file="~/fatty-acids/colocalisation/data/eqtlgen_data.txt", 
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_imputed.Rdata",
	gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed.Rdata",
	# gwis_file="~/fatty-acids/colocalisation/data/snplist_coloc_ukb_gwis_ratios_notimputed_eastasians.Rdata",
	ref_dat="~/fatty-acids/colocalisation/data/UKBB_10K.bim.gz",  
	ref_turn_off=FALSE) #don't load ref_dat to save time, e.g. because already loaded

gene_tab1<-data.frame(data_list[1],stringsAsFactors=F)
ref1<-data.frame(data_list[2],stringsAsFactors=F)
fa.tab1<-data.frame(data_list[3],stringsAsFactors=F)
gtex_data1<-data.frame(data_list[4],stringsAsFactors=F)
eqtlgen_data1<-data.frame(data_list[5],stringsAsFactors=F)
bbj_eqtl_data1<-data.frame(data_list[6],stringsAsFactors=F)
gwis_file1<-unlist(data_list[7])
lun_data1<-data.frame(data_list[8],stringsAsFactors=F)
crc_data1<-data.frame(data_list[9],stringsAsFactors=F)

gtex_data2<-effect_allele_gtex_data(dat=gtex_data1)
gtex_data3<-format_ens(dat=gtex_data2)
gtex_data4<-format_exposure2(dat=gtex_data3,beta="slope",se="slope_se",pval="pval_nominal",samplesize = NULL,effect_allele = "effect_allele",other_allele="other_allele",eaf ="maf",rsid = "SNP",ID = "tissue",exposure = "Ens",snps<-c("rs174546","rs968567"),standardise_beta=FALSE)

crc1<-read.table("~/MR_FattyAcids/data/summary_data/Final/colorectal.txt",sep="\t",head=T,stringsAsFactors=F)
crc1$ID<-60 
crc1<-crc1[!duplicated(crc1$rsid),]
crc2<-read.table("~/MR_FattyAcids/data/summary_data/Final/colorectal_accc.txt",sep="\t",head=T,stringsAsFactors=F)
crc2$ID<-3     
crc<-rbind.fill(crc1,crc2)
outcome_dat <-format_outcomes2(dat=crc)


mr_res<-tsmr2(outcome_dat=outcome_dat,exposure_dat=gtex_data4,action=3)
mr_res1<-mr_res[mr_res$SNP == "rs968567" & mr_res$population == "European" & !mr_res$tissue %in% c("artery coronary","lung","liver") & mr_res$gene == "FADS2", ]


mr_res1[,c("tissue","gene","OR","LCI","UCI","pval","beta.exposure","effect_allele.exposure")]
mr_res1$weight<-1/mr_res1$se/5




exp<-read.table("~/fatty-acids/mr/data/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",head=T,stringsAsFactors=F)
exp[exp$population=="East Asian",]
exposure_dat<-format_exposure2(dat=exp,standardise_beta=TRUE)
# exposure_dat[exposure_dat == "rs174546",]
# dat[,c("beta.exposure","other_allele.exposure","eaf.exposure","population.x")]

mr_res<-tsmr2(outcome_dat=outcome_dat,exposure_dat=exposure_dat)
mr_res2<-mr_res[mr_res$SNP == "rs968567" & mr_res$exposure == "GLA:LA" & mr_res$population.x == "European",]

mr_res3<-rbind.fill(mr_res1,mr_res2)
mr_res3$tissue[mr_res3$exposure == "GLA:LA"]<-"blood"
mr_res3$tissue[mr_res3$tissue=="whole blood"]<-"blood"
mr_res3$exposure[mr_res3$exposure == "ENSG00000134824"]<-"FADS2"

mr_res3[,c("exposure","tissue","beta.exposure","beta.outcome","b","se")]

Plot1<-mr_res3[mr_res3$tissue == "blood",]
Plot2<-Plot1[Plot1$exposure=="GLA:LA",]
Plot3<-Plot1[Plot1$exposure!="GLA:LA",]
Plot4<-mr_res3[mr_res3$tissue != "blood",]
plot.dat<-do.call(rbind,list(Plot2,Plot3,Plot4))

plot.dat$exposure[plot.dat$exposure=="FADS2"]<-"FADS2\ngene expression"
plot.dat$weight<-1/plot.dat$se.exposure/5
P1<-forest_plot_1_to_many(mr_res = plot.dat,b = "beta.exposure",se = "se.exposure",TraitM = "exposure",col1_width = 1,col1_title = "N",exponentiate = FALSE,trans = "identity",ao_slc = FALSE,lo =NULL,up =NULL, by="tissue",xlab="",xlab_top = "",addcols = NULL,addcol_widths = NULL,addcol_titles = "", subheading_size = 12,shape_points = 15,colour_scheme = "black",col_text_size = 5,weight ="weight",row_line_colour="black",x_axis_label_size=10)

png("~/fatty-acids/mr/results/plots/rs968567_fads2_D6D.png", width = 960, height = 960)
	print(P1) 
dev.off()



plot.dat$weight<-1/plot.dat$se/7
P2<-forest_plot_1_to_many(mr_res = plot.dat,b = "b",se = "se",TraitM = "exposure",col1_width = 1,col1_title = "N",exponentiate = FALSE,trans = "identity",ao_slc = FALSE,lo =NULL,up =NULL, by="tissue",xlab="",xlab_top = "",addcols = NULL,addcol_widths = NULL,addcol_titles = "", subheading_size = 12,shape_points = 15,colour_scheme = "black",col_text_size = 5,weight ="weight",row_line_colour="black",x_axis_label_size=10)
png("~/fatty-acids/mr/results/plots/rs968567_fads2_D6D_colorectalcancer.png", width = 960, height = 960)
	print(P2) 
dev.off()
# forest_plot_1_to_many(mr_res = plot.dat,b = "beta.exposure",se = "se.exposure",TraitM = "tissue",col1_width = 2,col1_title = "N",exponentiate = FALSE,trans = "identity",ao_slc = FALSE,lo =NULL,up =NULL, by="exposure",xlab="",xlab_top = "",addcols = NULL,addcol_widths = NULL,addcol_titles = "", subheading_size = 10,shape_points = 15,colour_scheme = "black",col_text_size = 5,weight ="weight",row_line_colour="black",x_axis_label_size=8)


plot.dat$weight<-1/plot.dat$se.exposure/5
P1<-forestplot(df = plot.dat,
			logodds = FALSE,
			name=tissue,
				  estimate=beta.exposure,
				  se=se.exposure,
				  shape=NULL,
				  colour = exposure,
				   xlab = "")+
				  ggplot2::theme(legend.position="none")
png("~/fatty-acids/mr/results/plots/rs968567_fads2_D6D_ggforest.png", width = 400, height = 400)
	print(P1) 
dev.off()

P2<-forestplot(df = plot.dat,
			logodds = TRUE,
			name=tissue,
				  estimate=b,
				  se=se,
				  shape=NULL,
				  colour = exposure,
				  xlab = "")+
				  ggplot2::theme(legend.position="none")

			# 	 theme(plot.title = element_text(size = ""))+
			# theme(text = element_text(size=14)))

png("~/fatty-acids/mr/results/plots/rs968567_fads2_D6D_colorectalcancer_ggforest.png", width = 400, height = 400)
	print(P2) 
dev.off()

