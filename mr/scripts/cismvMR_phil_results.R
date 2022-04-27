setwd("~/fatty-acids/mr")
source("scripts/cismvMR_functions.R")

fads_ld  = read.table("data/fads_r_matrix_ukb.ld", header=FALSE, stringsAsFactors=FALSE)
dim(fads_ld)
fads_frq  = read.table("data/fads_r_matrix_ukb.frq", header=TRUE, stringsAsFactors=FALSE)
dim(fads_frq)
load("data/data_for_cismvmr.RData")
sort(unique(Dat$trait))
unique(Dat$trait2)
Dat2<-Dat[Dat$trait=="GLA:LA / D6D",]
Dat<-Dat[Dat$trait!="GLA:LA / D6D",]
Dat2<-standardise_beta_gla_la_charge(dat=Dat2)
Dat<-rbind(Dat2,Dat)
 
beta1_prune_FIVW = NULL
beta2_prune_FIVW = NULL
se1_prune_FIVW = NULL
se2_prune_FIVW = NULL
nsnps<-NULL

z1<-beta1_prune_FIVW/se1_prune_FIVW
z2<-beta2_prune_FIVW/se2_prune_FIVW

p1<-2*pnorm(abs(z1),lower.tail=FALSE)
p2<-2*pnorm(abs(z2),lower.tail=FALSE)

exposure1<-c("AA:DGLA / D5D","FADS1 expression in adipose subcutaneous","FADS1 expression in adipose visceral omentum","FADS1 expression in blood","FADS1 expression in whole blood","FADS1 expression in colon sigmoid","FADS1 expression in liver")
exposure2<-c("GLA:LA / D6D","FADS2 expression in adipose subcutaneous","FADS2 expression in adipose visceral omentum","FADS2 expression in blood","FADS2 expression in whole blood","FADS2 expression in colon sigmoid","FADS2 expression in liver")

Res<-data.frame(matrix(c(exposure1,exposure2,beta1_prune_FIVW,se1_prune_FIVW,z1,p1,beta2_prune_FIVW,se2_prune_FIVW,z2,p2,nsnps),byrow=FALSE,nrow=length(beta1_prune_FIVW),ncol=11))
names(Res)<-c("exposure1","exposure2","beta1_prune_FIVW","se1_prune_FIVW","z1","p1","beta2_prune_FIVW","se2_prune_FIVW","z2","p2","nsnps")
Res$outcome<-"Colorectal cancer"
Res$panel<-"UKB"

write.table(Res,"~/fatty-acids/mr/results/cismvMR_phil_results.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)



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
	rho <- as.matrix(fads_ld_966.order.align)
	p <- length(bx1)

	Phi = as.matrix(fads_ld_966.order.align)
	K       = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2/sum((prcomp(Phi, scale=FALSE)$sdev^2)))>0.99)[1]

	fivw.out = FIVW.cor(rho, bx, by, sey, p, K)

	beta1_prune_FIVW[j] = fivw.out$fmvivw[1,1]
	beta2_prune_FIVW[j] = fivw.out$fmvivw[2,1]
	se1_prune_FIVW[j] = fivw.out$fmvivw[1,2]
	se2_prune_FIVW[j] = fivw.out$fmvivw[2,2]
	nsnps[j] = p
}

 