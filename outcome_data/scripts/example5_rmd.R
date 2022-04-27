library(CheckSumStats)
source("fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")

crc<-preformat_accc_3()
crc2<-preformat_accc_3_2(dat=crc)

write.table(crc2,"~/MR_FattyAcids/data/crc_test_dat.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

scp ~/MR_FattyAcids/data/crc_test_dat.txt ~/MR_FattyAcids/data/glioma_test_dat.txt /Users/ph14916/mrQC/inst/extdata


File<-system.file("extdata", "crc_test_dat.txt", package = "CheckSumStats")
crc<-read.table(File,sep="\t",head=TRUE,stringsAsFactors=FALSE)
dat<-format_data(dat=crc,outcome="Colorectal cancer",population="East Asian",ncase=23572,ncontrol=48700,study="ACCC",rsid="rsid",effect_allele="Allele1",other_allele="Allele2",lnor="Effect",lnor_se="StdErr",eaf="Freq1",p="P.value")
Pred<-predict_lnor_sh(dat=dat)
lm(Pred$lnor_pred ~ Pred$lnor)$coefficients
#> (Intercept)   Pred$lnor 
#> 0.002193062 0.231316576 
Plot4<-make_plot_pred_effect(dat=data.frame(Pred))
Plot4

png("~/mrQC/man/figures/README-example5.png",width=800)
	Plot4
dev.off()

summary(dat$Nstudies)
#> Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 1.00   11.75   14.00   11.69   14.00   15.00 
Pos<-Pred$Nstudies>=14
lm(Pred$lnor_pred[Pos] ~ Pred$lnor[Pos])$coefficients
#> (Intercept) Pred$lnor[Pos] 
#> 6.634408e-05   8.379027e-01 

Pos<-Pred$Nstudies<14
lm(Pred$lnor_pred[Pos] ~ Pred$lnor[Pos])$coefficients
#> (Intercept) Pred$lnor[Pos] 
#> 4.587284e-05   2.140147e-01 



preformat_accc_3_2<-function(dat=NULL){
	dat<-dat[crc$region == "FADS",]
	test<-sample(x=1:nrow(dat),size=100)
	dat<-dat[test,]
	Dir<-dat$Direction
	Dir1<-Dir[1]
	Dir2<-gsub("\\?","",Dir)
	dat$Nstudies<-nchar(Dir2)
	return(dat)
}
	

	# dat<-dat[order(dat$Nstudies),]
	# Nstudies<-dat$Nstudies
	# Nstudies2<-as.character(dat$Nstudies)
	# Nstudies2[Nstudies>=14]<-">14 studies"
	# Nstudies2[Nstudies<14 & Nstudies >= 11]<-"10-14 studies"
	# Nstudies2[Nstudies<=10 & Nstudies >5]<-"5-10 studies"
	# Nstudies2[Nstudies<=5]<-"<5 studies"