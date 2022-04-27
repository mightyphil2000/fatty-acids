setwd("C:/Users/stephen burgess/Dropbox (Cambridge University)/People/Fatima Batool")

fads_ld  = read.table("fattyacids_mvmr/fads_r_matrix_ukb.ld", header=FALSE, stringsAsFactors=FALSE)
dim(fads_ld)
fads_frq  = read.table("fattyacids_mvmr/fads_r_matrix_ukb.frq", header=TRUE, stringsAsFactors=FALSE)
dim(fads_frq)
load("fattyacids_mvmr/data_for_cismvmr.RData")
str(Dat)
table(Dat)
 
beta1_prune_0.2 = NULL
beta2_prune_0.2 = NULL
beta1_prune_0.3 = NULL
beta2_prune_0.3 = NULL
beta1_prune_0.4 = NULL
beta2_prune_0.4 = NULL
beta1_prune_PCA = NULL
beta2_prune_PCA = NULL
beta1_prune_FIVW = NULL
beta2_prune_FIVW = NULL
se1_prune_0.2 = NULL
se2_prune_0.2 = NULL
se1_prune_0.3 = NULL
se2_prune_0.3 = NULL
se1_prune_0.4 = NULL
se2_prune_0.4 = NULL
se1_prune_PCA = NULL
se2_prune_PCA = NULL
se1_prune_FIVW = NULL
se2_prune_FIVW = NULL


for (j in 1:7) {

if (j == 1) { D5D <- Dat[Dat$trait=="AA:DGLA / D5D",]; D6D <- Dat[Dat$trait=="GLA:LA / D6D",] }
if (j == 2) { D5D <- Dat[Dat$trait=="FADS1 expression in adipose subcutaneous",];     D6D <- Dat[Dat$trait=="FADS2 expression in adipose subcutaneous",] }
if (j == 3) { D5D <- Dat[Dat$trait=="FADS1 expression in adipose visceral omentum",]; D6D <- Dat[Dat$trait=="FADS2 expression in adipose visceral omentum",] }
if (j == 4) { D5D <- Dat[Dat$trait=="FADS1 expression in blood",];                    D6D <- Dat[Dat$trait=="FADS2 expression in blood",] }
if (j == 5) { D5D <- Dat[Dat$trait=="FADS1 expression in whole blood",];              D6D <- Dat[Dat$trait=="FADS2 expression in whole blood",] }
if (j == 6) { D5D <- Dat[Dat$trait=="FADS1 expression in colon sigmoid",];            D6D <- Dat[Dat$trait=="FADS2 expression in colon sigmoid",] }
if (j == 7) { D5D <- Dat[Dat$trait=="FADS1 expression in liver",];                    D6D <- Dat[Dat$trait=="FADS2 expression in liver",] }

outcome <- Dat[Dat$trait=="Colorectal cancer",]

## Get common variants in three datasets
commonvariants <- Reduce(intersect, list(D5D$marker,
                             D6D$marker, 
                             outcome$marker))
length(commonvariants)

## sort common variants
cv <- commonvariants [order(commonvariants)]
## Extract data for these common variants only. Note this will also give sorted data 
D5D <- D5D[D5D$marker%in%cv, ] #D5D sorting is different than following two.
D6D <- D6D[D6D$marker%in%cv, ]
outcome <- outcome[outcome$marker%in%cv, ] 

## order data
D5D <-   D5D[ order(D5D$marker), ] 
D6D <-   D6D[ order(D6D$marker), ] 
outcome <-  outcome[ order(outcome$marker), ] 
## writing cv on a file

rownames(fads_ld) <- colnames(fads_ld) <-  fads_frq$SNP
dim(fads_ld)
length(cv%in%fads_frq$SNP)
cv.cor.index <- which(fads_frq$SNP%in%cv)
fads_ld_966 <- fads_ld[cv.cor.index,cv.cor.index]

fads_ld_966.order = fads_ld_966[order(rownames(fads_ld_966)), order(rownames(fads_ld_966))]
which(colnames(fads_ld_966.order)!=D5D$marker)

effect.assoc = D5D$effect_allele
effect.cor   = fads_frq$A1[cv.cor.index][order(rownames(fads_ld_966))]

other.assoc = D5D$other_allele
other.cor   = fads_frq$A2[cv.cor.index][order(rownames(fads_ld_966))]
marker.cor = fads_frq$SNP[cv.cor.index][order(rownames(fads_ld_966))]

# which(effect.assoc!=effect.cor&effect.assoc!=other.cor) ## CHECK
# which(other.assoc!=effect.cor&other.assoc!=other.cor)   ## CHECK
# which(marker.cor!=D5D$marker)

flip.cor = ifelse(effect.assoc==effect.cor, 1, -1)

fads_ld_966.order.align = fads_ld_966.order*flip.cor%o%flip.cor

## Preparing input 
bx1 <- D5D$beta
bx2 <- D6D$beta
by <- outcome$beta
sey <- outcome$se
bx <- cbind(bx1, bx2)
sex <- cbind(D5D$se, D6D$se)
cormat <- as.matrix(fads_ld_966.order.align)
p <- length(bx1)

## variant pruning for ivw
D5D_edit = D5D
D6D_edit = D6D
cormat <- as.matrix(fads_ld_966.order.align)

varset = NULL; thres=0.2

while (sum(!is.na(D5D_edit$beta))>1) {
 D5Dz = abs(D5D_edit$beta)/D5D_edit$se
 varset = c(varset, D5D$marker[which.max(D5Dz)])
 D5D_edit$beta[which(abs(cormat[which.max(D5Dz),])>thres)] <- NA
 D6D_edit$beta[which(abs(cormat[which.max(D5Dz),])>thres)] <- NA
 D6Dz = abs(D6D_edit$beta)/D6D_edit$se
 varset = c(varset, D5D$marker[which.max(D6Dz)])
 D5D_edit$beta[which(abs(cormat[which.max(D6Dz),])>thres)] <- NA
 D6D_edit$beta[which(abs(cormat[which.max(D6Dz),])>thres)] <- NA
}
 
## Preparing input 
bx1 <- D5D$beta[which(D5D$marker%in%varset)]
bx2 <- D6D$beta[which(D5D$marker%in%varset)]
by <- outcome$beta[which(D5D$marker%in%varset)]
sey <- outcome$se[which(D5D$marker%in%varset)]
bx <- cbind(bx1, bx2)
sex <- cbind(D5D$se[which(D5D$marker%in%varset)], D6D$se[which(D5D$marker%in%varset)])
cormat <- as.matrix(fads_ld_966.order.align[which(D5D$marker%in%varset), which(D5D$marker%in%varset)])
p <- length(which(D5D$marker%in%varset))
#calling
library(MendelianRandomization)
mvivw = mr_mvivw(mr_mvinput(bx,sex,by,sey, corr=cormat))
beta1_prune_0.2[j] = mvivw$Estimate[1]
beta2_prune_0.2[j] = mvivw$Estimate[2]
se1_prune_0.2[j] = mvivw$StdError[1]
se2_prune_0.2[j] = mvivw$StdError[2]

## variant pruning for ivw
D5D_edit = D5D
D6D_edit = D6D
cormat <- as.matrix(fads_ld_966.order.align)

varset = NULL; thres=0.3

while (sum(!is.na(D5D_edit$beta))>1) {
 D5Dz = abs(D5D_edit$beta)/D5D_edit$se
 varset = c(varset, D5D$marker[which.max(D5Dz)])
 D5D_edit$beta[which(abs(cormat[which.max(D5Dz),])>thres)] <- NA
 D6D_edit$beta[which(abs(cormat[which.max(D5Dz),])>thres)] <- NA
 D6Dz = abs(D6D_edit$beta)/D6D_edit$se
 varset = c(varset, D5D$marker[which.max(D6Dz)])
 D5D_edit$beta[which(abs(cormat[which.max(D6Dz),])>thres)] <- NA
 D6D_edit$beta[which(abs(cormat[which.max(D6Dz),])>thres)] <- NA
}
 
## Preparing input 
bx1 <- D5D$beta[which(D5D$marker%in%varset)]
bx2 <- D6D$beta[which(D5D$marker%in%varset)]
by <- outcome$beta[which(D5D$marker%in%varset)]
sey <- outcome$se[which(D5D$marker%in%varset)]
bx <- cbind(bx1, bx2)
sex <- cbind(D5D$se[which(D5D$marker%in%varset)], D6D$se[which(D5D$marker%in%varset)])
cormat <- as.matrix(fads_ld_966.order.align[which(D5D$marker%in%varset), which(D5D$marker%in%varset)])
p <- length(which(D5D$marker%in%varset))
#calling
library(MendelianRandomization)
mvivw = mr_mvivw(mr_mvinput(bx,sex,by,sey, corr=cormat))
beta1_prune_0.3[j] = mvivw$Estimate[1]
beta2_prune_0.3[j] = mvivw$Estimate[2]
se1_prune_0.3[j] = mvivw$StdError[1]
se2_prune_0.3[j] = mvivw$StdError[2]

## variant pruning for ivw
D5D_edit = D5D
D6D_edit = D6D
cormat <- as.matrix(fads_ld_966.order.align)

varset = NULL; thres=0.4

while (sum(!is.na(D5D_edit$beta))>1) {
 D5Dz = abs(D5D_edit$beta)/D5D_edit$se
 varset = c(varset, D5D$marker[which.max(D5Dz)])
 D5D_edit$beta[which(abs(cormat[which.max(D5Dz),])>thres)] <- NA
 D6D_edit$beta[which(abs(cormat[which.max(D5Dz),])>thres)] <- NA
 D6Dz = abs(D6D_edit$beta)/D6D_edit$se
 varset = c(varset, D5D$marker[which.max(D6Dz)])
 D5D_edit$beta[which(abs(cormat[which.max(D6Dz),])>thres)] <- NA
 D6D_edit$beta[which(abs(cormat[which.max(D6Dz),])>thres)] <- NA
}
 
## Preparing input 
bx1 <- D5D$beta[which(D5D$marker%in%varset)]
bx2 <- D6D$beta[which(D5D$marker%in%varset)]
by <- outcome$beta[which(D5D$marker%in%varset)]
sey <- outcome$se[which(D5D$marker%in%varset)]
bx <- cbind(bx1, bx2)
sex <- cbind(D5D$se[which(D5D$marker%in%varset)], D6D$se[which(D5D$marker%in%varset)])
cormat <- as.matrix(fads_ld_966.order.align[which(D5D$marker%in%varset), which(D5D$marker%in%varset)])
p <- length(which(D5D$marker%in%varset))
#calling
library(MendelianRandomization)
mvivw = mr_mvivw(mr_mvinput(bx,sex,by,sey, corr=cormat))
beta1_prune_0.4[j] = mvivw$Estimate[1]
beta2_prune_0.4[j] = mvivw$Estimate[2]
se1_prune_0.4[j] = mvivw$StdError[1]
se2_prune_0.4[j] = mvivw$StdError[2]



###############

bx1 <- D5D$beta
bx2 <- D6D$beta
by <- outcome$beta
sey <- outcome$se
bx <- cbind(bx1, bx2)
sex <- cbind(D5D$se, D6D$se)
rho <- as.matrix(fads_ld_966.order.align)

Phi = as.matrix(fads_ld_966.order.align)
K       = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2/sum((prcomp(Phi, scale=FALSE)$sdev^2)))>0.99)[1]
       # K is number of principal components to include in analysis
       # this code includes principal components to explain 99% of variance in the risk factor
betaXG1 = as.numeric(bx1%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
betaXG2 = as.numeric(bx2%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
betaYG0 = as.numeric(by%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
sebetaYG0 = as.numeric(sey%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
Omegapca = sey%o%sey*rho
pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omegapca%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]

mvivw = mr_mvivw(mr_mvinput(cbind(betaXG1, betaXG2),cbind(rep(1, length(betaXG1)), rep(1, length(betaXG1))),
                          betaYG0, rep(1, length(betaXG1)), corr=pcOmega))

beta1_prune_PCA[j] = mvivw$Estimate[1]
beta2_prune_PCA[j] = mvivw$Estimate[2]
se1_prune_PCA[j] = mvivw$StdError[1]
se2_prune_PCA[j] = mvivw$StdError[2]

bx1 <- D5D$beta
bx2 <- D6D$beta
by <- outcome$beta
sey <- outcome$se
bx <- cbind(bx1, bx2)
sex <- cbind(D5D$se, D6D$se)
rho <- as.matrix(fads_ld_966.order.align)

fivw.out = FIVW.cor(rho, bx, by, sey, p, K)

beta1_prune_FIVW[j] = fivw.out$fmvivw[1,1]
beta2_prune_FIVW[j] = fivw.out$fmvivw[2,1]
se1_prune_FIVW[j] = fivw.out$fmvivw[1,2]
se2_prune_FIVW[j] = fivw.out$fmvivw[2,2]

}

  
#################




FIVW.cor <- function(rho, bx, by, sey, p, r){
eveccovg <- sqrt(p)*(eigen(rho)$vectors[,1:r]) #1:p, 1:r
#dim(eveccovg) <- c(p,r)
factors <- (rho%*%eveccovg)/p #num [1:p, 1:r]
evec <- eigen(t(factors)%*%factors)$vectors #num [1:r, 1:r]
eval <- eigen(t(factors)%*%factors)$values #num [1:r, 1:r]
factors <- factors%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec))) # num [1:p, 1:r]
lambda <- t(rho)%*%factors # num [1:p, 1:r]
evec <- eigen((t(lambda)%*%lambda)/p)$vectors # num [1:r, 1:r]
eval <- eigen((t(lambda)%*%lambda)/p)$values # num [1:r]
lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec))) # num [1:p, 1:r]
 
#MVIVW using approximate factor
Omega <-  (sey%o%sey)*rho   
Omegatr <- t(lambda)%*%Omega%*%lambda
betaxtr <- t(lambda)%*%bx
betaytr <- t(lambda)%*%by
 
L1 <- t(betaxtr)%*%solve(Omegatr)%*%betaxtr
L2 <-  t(betaxtr)%*%solve(Omegatr)%*%betaytr
fmvivw <- solve(L1)%*%L2
sefmvivw <-  sqrt(diag(solve(t(betaxtr)%*%solve(Omegatr)%*%betaxtr)))   
res <- cbind(coeff=fmvivw[,1], se=sefmvivw )
res <- list("est.no.factor"=r,  "loadings"=lambda, "fmvivw"=res)
return(res)
}

cbind(beta1_prune_0.2, se1_prune_0.2, beta1_prune_0.3, se1_prune_0.3, beta1_prune_0.4, se1_prune_0.4); 
cbind(beta1_prune_0.2, se1_prune_0.2, beta1_prune_PCA, se1_prune_PCA, beta1_prune_FIVW, se1_prune_FIVW);

cbind(beta2_prune_0.2, se2_prune_0.2, beta2_prune_0.3, se2_prune_0.3, beta2_prune_0.4, se2_prune_0.4); 
cbind(beta2_prune_0.2, se2_prune_0.2, beta2_prune_PCA, se2_prune_PCA, beta2_prune_FIVW, se2_prune_FIVW);

cbind(round(beta1_prune_FIVW, 3), round(se1_prune_FIVW, 3), round(beta2_prune_FIVW, 3), round(se2_prune_FIVW, 3))