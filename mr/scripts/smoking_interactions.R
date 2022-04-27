library(foreign)
library(multcomp)


# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/
load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/UKBB_pan_inclc44_cancer_smoking_rs174546.Rdata")

Dat1<-Dat[Dat$smoking =="Ever" ,]
Dat2<-Dat[Dat$smoking =="Never" ,]


head(Dat)

# rs174546_2 coded to reflect D5D raising allele
# rs174546 reflects D5D lowering allele
# rs174546_3 scaled to reflect SD change in D5D

mylogit <- glm(cancer ~ rs174546_3*smoking +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")

mylogit1 <- glm(cancer ~ rs174546_3 +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
mylogit2 <- glm(cancer ~ rs174546_3 +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat2, family = "binomial")

Int<-summary(mylogit)
Ever<-summary(mylogit1)
Never<-summary(mylogit2)

 3.616e-02-3.341e-02

lnor_ever <- Ever$coefficients[2,"Estimate"]
se_ever<-Ever$coefficients[2,"Std. Error"]
lci_ever<-lnor_ever-1.96*se_ever
uci_ever<-lnor_ever+1.96*se_ever

OR_ever<-round(exp(lnor_ever),2)
LCI_ever<-round(exp(lci_ever),2)
UCI_ever<-round(exp(uci_ever),2)
P_ever<-Ever$coefficients[2,"Pr(>|z|)"]

lnor_never<-Never$coefficients[2,"Estimate"]
se_never<-Never$coefficients[2,"Std. Error"]
lci_never<-lnor_never-1.96*se_never
uci_never<-lnor_never+1.96*se_never

OR_never<-round(exp(lnor_never),2)
LCI_never<-round(exp(lci_never),2)
UCI_never<-round(exp(uci_never),2)
P_never<-Never$coefficients[2,"Pr(>|z|)"]

c(OR_ever,LCI_ever,UCI_ever,P_ever)
c(OR_never,LCI_never,UCI_never,P_never)

exp(3.616e-02 )

summary(mylogit)
summary(glht(mylogit, linfct = c("rs174546_2 + rs174546_2:smokingNever = 0")))
?glht

Never<- 0.002368
Never.se<- 0.009982 

nlci<-Never-1.96*Never.se
nuci<-Never+1.96*Never.se

Ever<- 3.110e-02
Ever.se<- 1.031e-02 

elci<-Ever-1.96*Ever.se
euci<-Ever+1.96*Ever.se


c(round(exp(Never),2),round(exp(nlci),2),round(exp(nuci),2))
c(round(exp(Never/0.86),2),round(exp(nlci/0.86),2),round(exp(nuci/0.86),2))

c(round(exp(Ever),2),round(exp(elci),2),round(exp(euci),2))
c(round(exp(Ever/0.86),2),round(exp(elci/0.86),2),round(exp(euci/0.86),2))



   0.009982 
2.873e-02 + -3.110e-02

1/exp(-3.110e-02/0.86)

-3.110e-02+2.873e-02
1/exp(-0.002368)

1/exp(-0.0315424/0.86)
1/exp(-0.0021293/0.86)
1/exp(2.873e-02/0.86)

summary(glht(mod, linfct = c("genderMale + genderMale:agehiTRUE = 0")))


cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/smoking/