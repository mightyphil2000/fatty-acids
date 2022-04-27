
standardise_beta_gla_la_charge<-function(dat=NULL){
	z<-dat$beta/dat$se
	p<-dat$eaf
	n<-7596 #In Guan CHARGE GWAS, sample size in GLA analysis; n=8631 in analysis of LA
 	dat$beta<-z/sqrt(2*p*(1-p)*(n+z^2)) 
 	dat$se<-dat$beta/z
 	return(dat)
}

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