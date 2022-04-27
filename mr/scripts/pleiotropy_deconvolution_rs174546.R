library(ggforestplot)
library(ggplot2)

source("~/fatty-acids/mr/scripts/mr_functions.R")
a<-read.table("~/fatty-acids/mr/results/pleiotropy_deconvolution_rs174546.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)

Pos<-which(a$unit=="mg/dL decrease")
uci<-1/a$or_can_exp_lci[Pos]
lci<-1/a$or_can_exp_uci[Pos]
a$or_can_exp_lci[Pos]<-lci 
a$or_can_exp_uci[Pos]<-uci
a$unit<-gsub(" decrease","",a$unit)

a$unit_value[a$unit_value=="100x10^9"]<-100e9
a$unit_value<-as.numeric(a$unit_value)
Pos<-a$unit!="sd" & !is.na(a$unit)
a[Pos,c("unit","unit_value","sd_value_opengwas")]
a$scale_factor<-1
a$scale_factor[Pos]<-a$unit_value[Pos]/a$sd_value_opengwas[Pos]
a$scale_factor[is.na(a$unit)]<-NA

a$ln_or_can_exp<-log(a$or_can_exp)
a$se_can_exp<-(log(a$or_can_exp_uci)-log(a$or_can_exp_lci))/(1.96*2)
a$ln_or_can_exp_sd<-a$ln_or_can_exp
a$ln_or_can_exp_sd<-a$ln_or_can_exp/a$scale_factor
a$se_can_exp_sd<-a$se_can_exp/a$scale_factor

a$or_can_exp_sd<-exp(a$ln_or_can_exp_sd)
a$lci_can_exp_sd<-exp(a$ln_or_can_exp_sd-1.96*a$se_can_exp_sd)
a$uci_can_exp_sd<-exp(a$ln_or_can_exp_sd+1.96*a$se_can_exp_sd)

a$lnor_decon<-a$ln_or_can_exp_sd*a$b_exp_sd 
a$se_decon<-a$se_can_exp_sd*a$b_exp_sd 
a$or_decon<-exp(a$lnor_decon)  
a$lci_decon<-exp(a$lnor_decon-1.96*a$se_decon)  
a$uci_decon<-exp(a$lnor_decon+1.96*a$se_decon)  

a$lnor_can_ea<-log(a$or_can_ea)
a$se_can_ea<-(log(a$or_can_ea_uci)-log(a$or_canc_ea_lci))/(1.96*2)

a$p_diff<-z_test_diff2(lnor_total=a$lnor_can_ea,se_total=a$se_can_ea,lnor_decon=a$lnor_decon,se_decon=a$se_decon)

a$lnor_decon_expected<-a$lnor_can_ea/a$b_exp_sd
a$se_decon_expected<-a$se_can_ea/a$b_exp_sd
a$or_decon_expected<-exp(a$lnor_decon_expected)
a$lci_decon_expected<-exp(a$lnor_decon_expected-1.96*a$se_decon_expected)
a$uci_decon_expected<-exp(a$lnor_decon_expected+1.96*a$se_decon_expected)

write.table(a,"~/fatty-acids/mr/results/pleiotropy_deconvolution_rs174546_results.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

head(a)

a1<-a[a$cancer %in% c("colorectal cancer","lung cancer"),]
b1<-a1[,c("cancer","exposure","lnor_can_ea","se_can_ea")]
names(b1)<-c("cancer","exposure","lnor","se")
b1$effect<-"total effect"
b2<-a1[,c("cancer","exposure","lnor_decon","se_decon")]
names(b2)<-c("cancer","exposure","lnor","se")
b2$effect<-"indirect effect"

b<-rbind(b1,b2)
b1<-b[b$cancer == "colorectal cancer",]
b2<-b[b$cancer == "lung cancer",]

Plot<-forestplot(df = b1,
	name=exposure,
	  logodds = TRUE,
	  estimate=lnor,
	  se=se,
	  shape=NULL,
	  colour = effect,
	   xlab = "OR (95% CI) for cancer per copy of the C allele")+
theme(plot.title = element_text(size = 10))