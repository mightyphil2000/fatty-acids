load("~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")

load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
source("~/fatty-acids/mr/scripts/mr_functions.R")

sort(disc.tab9$cancer)
crc1<-mr_res1[mr_res1$cancer %in% c("Distal colorectal cancer","Proximal colorectal cancer"),]
crc2<-mr_res1[mr_res1$cancer %in% c("Rectal cancer","Colon cancer"),]
crc1<-mr_res1[mr_res1$cancer %in% c("Distal colorectal cancer","Proximal colorectal cancer"),]
crc2<-mr_res1[mr_res1$cancer %in% c("Rectal cancer","Colon cancer"),]
luc<-mr_res1[mr_res1$cancer %in% c("Lung adenocarcinoma","Small cell lung carcinoma","Squamous cell lung cancer"),]


# skin cancer
skc<-mr_res1[mr_res1$cancer %in% c("Basal cell carcinoma","Melanoma","Squamous cell carcinoma" ),]

skc<-skc[order(as.numeric(skc$cases),decreasing=TRUE),]
skc<-skc[!duplicated(skc$cancer),]

# bcc versus scc
skc1<-skc[skc$cancer %in% c("Basal cell carcinoma","Squamous cell carcinoma"),]

skc2<-skc[skc$cancer %in% c("Basal cell carcinoma","Melanoma"),] 

corr<-corr_function(test=skc1)
M<-matrix(c(1,corr,corr,1),nrow=2,ncol=2)
skc1$se_decoupled<-decoupling(s=as.numeric(skc1$se),C=M)

P1<-z_test_diff2(betas=as.numeric(skc1$b),ses=as.numeric(skc1$se_decoupled))
P2<-z_test_diff2(betas=as.numeric(skc2$b),ses=as.numeric(skc2$se)) #melanoma entirely independent of bcc in cases and controls

skc1[,c("outcome","cancer","cases","study.abbreviation","b","se_decoupled","se")]

skc2[,c("outcome","cancer","cases","study.abbreviation","b","se")]

# distal crc versus proximal crc
corr<-corr_function(test=crc1)
M<-matrix(c(1,corr,corr,1),nrow=2,ncol=2)
crc1$se_decoupled<-decoupling(s=as.numeric(crc1$se),C=M)
P1<-z_test_diff2(betas=as.numeric(crc1$b),ses=as.numeric(crc1$se_decoupled))

# rectal versus colon cancer
corr<-corr_function(test=crc2)
M<-matrix(c(1,corr,corr,1),nrow=2,ncol=2)
crc2$se_decoupled<-decoupling(s=as.numeric(crc2$se),C=M)
P2<-z_test_diff2(betas=as.numeric(crc2$b),ses=as.numeric(crc2$se_decoupled))

crc1[,c("se","se_decoupled")]

# squamous versus small cell versus adeno
luc1<-luc[c(1,2),]
luc2<-luc[c(1,3),]
luc3<-luc[c(2,3),]
corr1<-corr_function(test=luc1)
corr2<-corr_function(test=luc2)
corr3<-corr_function(test=luc3)

  1       2    3 
1 1     corr1 corr2
2 corr1   1   corr3
3 corr2  corr3  1

M<-matrix(c(1,corr1,corr2,corr1,1,corr3,corr2,corr3,1),nrow=3,ncol=3)
luc$se_decoupled<-decoupling(s=as.numeric(luc$se),C=M)

Q_res<-Q_test(se=as.numeric(luc$se_decoupled),beta=as.numeric(luc$b))

Q<-unlist(Q_res[1])
Q.p<-unlist(Q_res[2])



P1<-z_test_diff2(betas=as.numeric(crc1$b),ses=as.numeric(crc1$se_decoupled))



#  dcrc pcrc 
# dcrc 1    corr
# pcrc corr    1

corr_function<-function(test=NULL)
{
	nk1<-as.numeric(test$cases[1])
	nk0<-as.numeric(test$controls[1])
	nk=nk1+nk0
	nl1<-as.numeric(test$cases[2])
	nl0<-as.numeric(test$controls[2])
	# test[2,]
	nl<-nl1+nl0
	# nkl0<-min(nk0,nl0)
	nkl0<-min(c(nl0,nk0))
	nkl1<-0  #assume cases are independent

	corr<-(nkl0*sqrt((nk1*nl1)/(nk0*nl0))+nkl1*sqrt((nk0*nl0)/(nk1*nl1)))/sqrt(nk*nl)
	return(corr)
}


z_test_diff2<-function(betas=NULL,ses=NULL)
{		
	b1<-betas[1]
	b2<-betas[2]
	se1<-ses[1]
	se2<-ses[2]
	diff<-b1-b2
	diff_se = sqrt(se1^2 + se2^2)
	Z <- diff/diff_se
	# pnorm(abs(1.96),lower.tail=FALSE)*2
	P<-pnorm(abs(Z),lower.tail=FALSE)*2
	return(P)
}

decoupling(s=mr_res$se,C=M)


Q_test<-function(se=NULL,beta=NULL)
{
	w<-1/se^2
	b.fixed<-sum(beta*w)/(sum(w))
	se.fixed<-sqrt(sum(w)^-1)
	z<-abs(b.fixed/se.fixed)
	p.fixed<-pnorm(z,lower.tail=F)*2
	nstudies.fixed<-length(beta)
	Q<-sum((b.fixed-beta)^2*w)
	df.Q<-length(beta)-1		
	Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)
	return(list(Q=Q,Q.p=Q.p))
}