library(TwoSampleMR)
ao<-available_outcomes()
id<-ao$id[which(ao$consortium == "GSCAN")]
exp_dat <- extract_instruments(outcomes=id)
# exp_dat$beta.exposure<-unlist(b_sd()[1])
# exp_dat$se.exposure<-unlist(b_sd()[2])

out_dat <- extract_outcome_data(
    snps = exp_dat$SNP,
    outcomes = c("ieu-a-985","ieu-a-986","ieu-a-987")
)

dups<-unique(out_dat$SNP[duplicated(out_dat$SNP)])
# Temp<-out_dat[out_dat$SNP %in% dups,]
# Temp[order(Temp$SNP),c("SNP","beta.outcome","se.outcome","effect_allele.outcome","other_allele.outcome")]
# exp_dat[exp_dat$SNP=="rs58379124",]
# Temp[Temp$SNP == "rs58379124",]

dat <- harmonise_data(
    exposure_dat = exp_dat, 
    outcome_dat = out_dat
)

dat<-dat[!dat$remove,]
dat2<-dat[!duplicated(paste(dat$SNP,dat$id.exposure,dat$id.outcome)),]
Pos<-grep("exposure",names(dat2))
Pos1<-which(names(dat2)=="SNP")
Pos<-c(Pos1,Pos)
exposure_dataset<-dat2[,Pos]
exposure_dataset<-exposure_dataset[exposure_dataset$exposure== "Cigarettes smoked per day || id:ieu-b-142",]
write.table(exposure_dataset,"~/fatty-acids/mr/data/instruments_cig.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
res<-mr(dat2)
write.table(res,"~/fatty-acids/mr/results/mr_res_gscan_cancer.txt",sep="\t",col.names=TRUE,row.names=TRUE)

scale_factor<-1/0.017
res$or<-round(exp(res$b),2)
res$lci<-round(exp(res$b-1.96*res$se),2)
res$uci<-round(exp(res$b+1.96*res$se),2)

res$or2<-round(exp(res$b/scale_factor),2)
res$lci2<-round(exp(res$b/scale_factor-1.96*res$se/scale_factor),2)
res$uci2<-round(exp(res$b/scale_factor+1.96*res$se/scale_factor),2)


res1<-res[res$method == "Inverse variance weighted",]
res1[,c("id.outcome","exposure","or","lci","uci")]
res1[,c("id.outcome","exposure","or2","lci2","uci2")]
res1

res2<-res1[res1$id.outcome == "ieu-a-987",]
res3<-res1[res1$id.outcome == "ieu-a-985",]
res4<-res1[res1$id.outcome == "ieu-a-986",]


lnor<-res3$b[res3$exposure == "Cigarettes smoked per day || id:ieu-b-142"]
se<-res3$se[res3$exposure == "Cigarettes smoked per day || id:ieu-b-142"]
# ieu-a-987 lung cancer 
# "ieu-a-985 lung cancer in ever smokers
#ieu-a-986 lung cancer in never smokers 


# lnor<-res$b[res$method == "Inverse variance weighted"]
lnor2<-lnor/(1/0.017) #scaled to reflect effect of PUFA desaturase activity on smoking 0.003 (SE=0.0058)
# se<-res$se[res$method == "Inverse variance weighted"]
se2<-se/(1/0.017)

# lnor<-res$b[res$method == "MR Egger"]
# lnor2<-lnor/250 #scaled to reflect effect of PUFA desaturase activity on smoking 0.003 (SE=0.0058)
# se<-res$se[res$method == "MR Egger"]
# se2<-se/250

or<-round(exp(lnor),2)
lci<- round(exp(lnor - 1.96*se),2)
uci<- round(exp(lnor + 1.96*se),2)

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
b_sd<-function(){
	# exp_dat<-exp_dat[exp_dat$exposure == "Cigarettes smoked per day || id:ieu-a-961",]
	z<-abs(exp_dat$beta.exposure/exp_dat$se.exposure)
	p<-exp_dat$eaf
	# Pos<-p>0.5
	# p[Pos]<-1-p[Pos]
	n<-257120
	b_sd<-z/sqrt(2*p*(1-p)*(n+z^2)) 
	se_sd<-b_sd/z
	return(list(b_sd,se_sd))
}

