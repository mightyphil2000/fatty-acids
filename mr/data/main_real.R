## For real data 
##rm(list=ls())
library(MASS)
library(MendelianRandomization)

## following two function takes cor of gv (gr) as an in put 
est.factor <-  function(rho){
pca <- prcomp(rho,scale=FALSE)
pca.sd <- pca$sdev^2
pca.sum <- sum(pca.sd)
ratio <- pca.sd/pca.sum
ratio.cum <- cumsum(ratio)
r <- which(ratio.cum > 0.99)[1] 
return(r)
}

## Approximate factor model with MV-IVW estimates
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

## Extension of Burgess (2017) IVW method to multivariate risk factors
# input beta matrix for all risk factors
# input direct correlation matrix gr
IVWPCA <- function(rho, bx, sey, r){
  matsum <- matrix(0, nrow=nrow(bx), ncol=nrow(bx))
  for( i in 1:ncol(bx)){
    mat1 <- (bx[,i])%o%(bx[,i])  
    matsum <-   matsum + mat1  
  }
  
Phimulti <- matsum/((sey)%o%(sey))*rho
Omega <-  (sey%o%sey)*rho
#summary(prcomp(Phimulti, scale=FALSE))
#K = est.factor(gr) 
K = r
bx0 = t(prcomp(Phimulti, scale=FALSE)$rotation[,1:K])%*%bx
by0 = t(prcomp(Phimulti, scale=FALSE)$rotation[,1:K])%*%by
Omegapca = sey%o%sey*rho
pcOmega = t(prcomp(Phimulti, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phimulti, scale=FALSE)$rotation[,1:K]
L1 <-  solve(t(bx0)%*%solve(pcOmega)%*%bx0)
L2 <- t(bx0)%*%solve(pcOmega)%*%by0
beta_IVWcorrel.pc =L1%*%L2
se_IVWcorrel.fixed.pc =  sqrt(diag(solve(t(bx0)%*%solve(pcOmega)%*%bx0)))    
res <- cbind("coeff"=beta_IVWcorrel.pc, "se"=se_IVWcorrel.fixed.pc)
return(res)
}

### Classical MVIVW 
ivw <- function(rho, bx, sex, by, sey){
mrout <- mr_mvivw(mr_mvinput(bx, sex, by, sey, correl=rho), model="fixed")
 res <- cbind("coeff"= mrout$Estimate, "ses"=mrout$StdError)
return(res)
}

## Ready to Call
