cd ~/fatty-acids/colocalisation/data/
 
load("snplist_coloc_ukb_gwis_ratios_imputed.Rdata")

Charge1[Charge1$SNP %in% c("rs603424","rs6584379") & Charge1$file =="POA_to_PA.tab.gz",]

b<-c(0.001859806,0.001090576)
se<-c(0.0001956959, 0.0001890183)
z = b/se
p<-c(0.200787,0.220472)
n<-8631
b_sd = z/sqrt(2*p*(1-p)*(n+z^2)) 
var<-1
r2<-2*b_sd^2*p*(1-p)/var
r2_sum<-sum(r2)


k<-2
F<-r2_sum*(n-1-k)/((1-r2_sum)*k )

# median detectable odds ratio assuming 80% power
#median n cases= 4523
rsq<-c(0.01,0.014) 
b1<-log(1.8)
b1<-log(1.65)
sig<-0.05 #alpha
n<-4523*2 
ratio<-1
pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))

# Present results only for those cancer where power to detect odds ratio 1.5 at least 50% for all enzymes (i.e for SCD which has lowest r2)
# Hodgkin’s lymphoma (smallest cancer presented in April meeting)
Ncases<-3077
Ncontrols<-13680
n<-Ncontrols+Ncases
ratio<-Ncases/Ncontrols
rsq<-c(0.01,0.014) 
b1<-log(1.5)
sig<-0.05 #alpha
pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))

# B- cell childhood acute lymphoblastic leukemia
Ncases<-2442
Ncontrols<-14609
n<-Ncontrols+Ncases
ratio<-Ncases/Ncontrols
rsq<-c(0.01,0.014) 
b1<-log(1.5)
sig<-0.05 #alpha
pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))

# Neuroblastoma
Ncases<-2101
Ncontrols<-4202
n<-Ncontrols+Ncases
ratio<-Ncases/Ncontrols
rsq<-c(0.01,0.014) 
b1<-log(1.5)
sig<-0.05 #alpha
pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))


# bladder cancer
Ncases<-1799
Ncontrols<-4745
n<-Ncontrols+Ncases
ratio<-Ncases/Ncontrols
rsq<-c(0.01,0.014) 
b1<-log(1.5)
sig<-0.05 #alpha
pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))

# Multiple myeloma
Ncases<-1714
Ncontrols<-10391
n<-Ncontrols+Ncases
ratio<-Ncases/Ncontrols
rsq<-c(0.01,0.014) 
b1<-log(1.5)
sig<-0.05 #alpha
pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))


# meningioma
Ncases<-1606
Ncontrols<-9823
n<-Ncontrols+Ncases
ratio<-Ncases/Ncontrols
rsq<-c(0.01,0.014) 
b1<-log(1.5)
sig<-0.05 #alpha
pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))
