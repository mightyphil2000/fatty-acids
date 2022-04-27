library(TwoSampleMR)
ao<-available_outcomes()
source("~/fatty-acids/mr/scripts/mr_functions.R")

exp<-read.table("~/fatty-acids/mr/data/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",head=T,stringsAsFactors=F)
exp[exp$exposure == "AA:DGLA" & exp$population=="East Asian",]
exposure_dat<-format_exposure2(dat=exp,standardise_beta=TRUE)
exposure_dat<-exposure_dat[exposure_dat$exposure=="AA:DGLA" & exposure_dat$population=="European" & exposure_dat$SNP == "rs174546",]

id<-ao$id[which(ao$consortium == "GSCAN")]
out_dat <- extract_outcome_data(
    snps = "rs174546",
    outcomes = id
)

head(out_dat)

dat <- harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = out_dat
)
res<-mr(dat)
ao1<-ao[ao$id %in% id,]
ao1[,c("trait","sample_size")]

lnor<-res$b[res$outcome == "smoking initiation || id:ieu-b-4877"]
se<-res$se[res$outcome == "smoking initiation || id:ieu-b-4877"]

or<-round(exp(lnor),3)
lci<- round(exp(lnor - 1.96*se),3)
uci<- round(exp(lnor + 1.96*se),3)

c(or,lci,uci)

or2<-round(exp(lnor2),2)
lci2<- round(exp(lnor2 - 1.96*se2),2)
uci2<- round(exp(lnor2 + 1.96*se2),2)


matrix(c(or,lci,uci),nrow=2,ncol=3)
matrix(c(or2,lci2,uci2),nrow=2,ncol=3)
# , method_list=c("mr_egger_regression", "mr_ivw"))


# dat<-dat[!dat$remove,]
# dups<-unique(dat$SNP[duplicated(dat$SNP)])
# dat[dat$SNP %in% dups,c("SNP","beta.outcome","se.outcome")]

# median(exp_dat$samplesize.exposure)
dat=out_dat[out_dat$outcome == "Cigarettes smoked per day || id:ieu-b-142", ]

b_sd<-function(dat=NULL){
	# exp_dat<-exp_dat[exp_dat$exposure == "Cigarettes smoked per day || id:ieu-a-961",]
	z<-abs(dat$beta.outcome/dat$se.outcome)
	p<-dat$eaf.outcome
	# Pos<-p>0.5
	# p[Pos]<-1-p[Pos]
	n<-dat$samplesize.outcome
	b_sd<-z/sqrt(2*p*(1-p)*(n+z^2)) 
	se_sd<-b_sd/z
	return(list(b_sd,se_sd))
}



b<-ieugwasr::associations(id="ieu-b-142", variants="rs11264100",proxies=0) 
data.frame(b)

A	G	0.876	35.8	2.22E-09	
beta=-0.022165931	
SE=0.003704627	
N=335394