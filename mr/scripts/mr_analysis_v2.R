source("~/fatty-acids/mr/scripts/mr_functions.R")
library(plyr)
library(TwoSampleMR)

crc1<-read.table("~/MR_FattyAcids/data/summary_data/Final/colorectal.txt",sep="\t",head=T,stringsAsFactors=F)
crc1$ID<-60 

crc1<-crc1[!duplicated(crc1$rsid),]

dups<-crc1$rsid[duplicated(crc1$rsid)]
temp<-crc1[crc1$rsid %in% dups,]
temp<-temp[order(temp$rsid),]


crc2<-read.table("~/MR_FattyAcids/data/summary_data/Final/colorectal_accc.txt",sep="\t",head=T,stringsAsFactors=F)
crc2$ID<-3     

crc<-rbind.fill(crc1,crc2)

outcome_dat <-format_outcomes2(dat=crc)


exp<-read.table("~/fatty-acids/mr/data/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",head=T,stringsAsFactors=F)
exp[exp$population=="East Asian",]
exposure_dat<-format_exposure2(dat=exp,standardise_beta=TRUE)
# exposure_dat[exposure_dat == "rs174546",]
# dat[,c("beta.exposure","other_allele.exposure","eaf.exposure","population.x")]

mr_res<-tsmr2(outcome_dat=outcome_dat,exposure_dat=exposure_dat)
mr_res1<-format_mr_res(dat=mr_res,snp_keep="rs968567",keep_adjrs174546=FALSE)
mr_res1[,c("exposure","outcome","population.y")]
mr_res1[,c("SNP","exposure","OR","LCI","UCI","pval")]

dat<-mr_res1

res_out

res_exp<-format_plot_data(mr_res=mr_res1,beta="beta.exposure",se="se.exposure",modify_weight=5,beta_reverse=FALSE)
res_out<-format_plot_data(mr_res=mr_res1,beta="beta.outcome",se="se.outcome",modify_weight=20,exp=TRUE,beta_reverse=FALSE)
res_mr<-format_plot_data(mr_res=mr_res1,beta="b",se="se",modify_weight=10,beta_reverse=FALSE,exp="TRUE")

res_out[,c("exposure","OR","lci","uci","pval.outcome","population","effect_allele.outcome")]
unique(res_exp[,c("exposure","beta.exposure","se","lci","uci","pval.exposure","population","eaf.outcome","effect_allele.outcome")])
unique(res_mr[,c("exposure","OR","lci","uci","pval","population")])

res_exp$beta.exposure
res_exp$se.exposure

Plot1<-forest_plot_1_to_many(mr_res = res_exp,b = "beta.exposure",se = "se.exposure",TraitM = "population",col1_width = 2,col1_title = "N",
		exponentiate = FALSE,trans = "identity",ao_slc = FALSE,lo =0,up =1, by="exposure",xlab="",xlab_top = "",addcols = NULL,addcol_widths = NULL,addcol_titles = "", subheading_size = 1,shape_points = 15,colour_scheme = "black",col_text_size = 15,weight ="weight",row_line_colour="black",x_axis_label_size=8)

Plot2<-forest_plot_1_to_many(mr_res = res_out,b = "beta.outcome",se = "se.outcome",TraitM = "ncase.outcome",col1_width = 2,col1_title = "",
		exponentiate = TRUE,trans = "identity",ao_slc = FALSE,lo =0.95,up =1.20, by="exposure",xlab="",xlab_top = "",addcols = NULL,addcol_widths = NULL,addcol_titles = "", subheading_size = 1,shape_points = 15,colour_scheme = "black",col_text_size = 15,weight ="weight",row_line_colour="black",x_axis_label_size=8)

Plot3<-forest_plot_1_to_many(mr_res = res_mr,b = "b",se = "se",TraitM = "ncase.outcome",col1_width = 2,col1_title = "",
		exponentiate = TRUE,trans = "identity",ao_slc = FALSE,lo =0.95,up =1.20, by="exposure",xlab="",xlab_top = "",addcols = NULL,addcol_widths = NULL,addcol_titles = "", subheading_size = 1,shape_points = 15,colour_scheme = "black",col_text_size = 15,weight ="weight",row_line_colour="black",x_axis_label_size=8)

png("~/fatty-acids/mr/results/plots/d5d_d6d_schs_charge.png", width = 250, height = 500)
	print(Plot1) 
dev.off()

png("~/fatty-acids/mr/results/plots/variant_crc_accc_gecco.png", width = 250, height = 500)
	print(Plot2) 
dev.off()


png("~/fatty-acids/mr/results/plots/mr_crc_accc_gecco.png", width = 250, height = 500)
	print(Plot3) 
dev.off()


