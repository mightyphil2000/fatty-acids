# rerun Sean Harrison's function for SNPs with imputation info/r2 scores <0.8
source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

Dat<-collate_dat(postqc=FALSE)

Dat$info<-as.numeric(Dat$info)
Dat$info1<-as.numeric(Dat$info1)
Dat$info2<-as.numeric(Dat$info2)
Dat$info3<-as.numeric(Dat$info3)

L<-NULL
for(i in 1:nrow(Dat)){
	# print(i)
	L[[i]]<-mean(c(Dat$info[i],Dat$info1[i],Dat$info2[i],Dat$info3[i]), na.rm=T)
}

Dat$mean_info<-unlist(L)

Dat1<-Dat[!is.na(Dat$mean_info),]
Dat1<-Dat1[!is.na(Dat1$eaf),]
Dat1$z<-Dat1$lnor/Dat1$se
maf<-Dat1$eaf
Pos<-which(maf>0.5)
maf[Pos]<-1-maf[Pos]
Dat1$maf<-maf

IDS<-unique(Dat1$ID)
IDS1<-IDS[1:16]
IDS2<-IDS[17:26]
IDS3<-IDS[26:40]
IDS4<-IDS[41:length(IDS)]


Data<-Dat1[Dat1$ID %in% IDS1,]
Datb<-Dat1[Dat1$ID %in% IDS2,]
Datc<-Dat1[Dat1$ID %in% IDS3,]
Datd<-Dat1[Dat1$ID %in% IDS4,]

List<-ls()
List<-List[!List %in% c("Data","Datb","Datc","Datd")]
rm(list=List)

save.image(file="~/fatty-acids/outcome_data/data/input_predlnor_sh_impinfo.Rdata")
load("~/fatty-acids/outcome_data/data/input_predlnor_sh_impinfo.Rdata")

# IDS<-unique(Data$ID)
# dat2<-Data[Data$ID == IDS[1],]
# Res<-pred_lnor_sh(dat2=dat2[1:10,])
# Res1<-pred_lnor_sh(dat2=dat2[1:10,])

# X<-unlist(Res[1])
# Y<-unlist(Res1[1])
# X
# Y

pred_lnor_sh(dat2=Datb)
pred_lnor_sh(dat2=Datc)
pred_lnor_sh(dat2=Datd)

any(is.na(Data$eaf))
Datc<-Datc[is.na(Datc$eaf]
any(is.na(Datc$eaf))

Data$lnor_sh<-unlist(log_or)
Data$se_sh<-unlist(log_or_se)
save(Data,file="~/Data_info_lnorsh.Rdata")

Datb$lnor_sh<-unlist(log_or)
Datb$se_sh<-unlist(log_or_se)
save(Datb,file="~/Datb_info_lnorsh.Rdata")

Datc$lnor_sh<-unlist(log_or)
Datc$se_sh<-unlist(log_or_se)
save(Datc,file="~/Datc_info_lnorsh.Rdata")

Datd$lnor_sh<-unlist(log_or)
Datd$se_sh<-unlist(log_or_se)
save(Datd,file="~/Datd_info_lnorsh.Rdata")

