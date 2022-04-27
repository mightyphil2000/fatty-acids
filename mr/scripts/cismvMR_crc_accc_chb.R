setwd("~/fatty-acids/mr")
source("scripts/cismvMR_functions.R")

# chb (chinese from Beijing reference panel) from 1000 genomes n=103
# ACCC samples come from south (23%) and central & northern China. 
fads_ld  = read.table("data/fads_r_matrix_chb.ld", header=FALSE, stringsAsFactors=FALSE)
head(fads_ld)
fads_frq  = read.table("data/fads_r_matrix_chb.frq", header=TRUE, stringsAsFactors=FALSE)
dim(fads_frq)
head(fads_frq)

# colorectal cancer in ACCC study
# data_for_cismvmr_colorectal_accc.RData
# Esophageal squamous cell carcinoma in meta analysis of BJ and N-UGC studies
# data_for_cismvmr_esophageal_scc.RData
# Lung cancer in Biobank Japan
# data_for_cismvmr_lung_cancer_bj.RData

load("data/data_for_cismvmr_colorectal_accc.RData")
sort(unique(Dat$trait))

beta1_prune_FIVW = NULL
beta2_prune_FIVW = NULL
se1_prune_FIVW = NULL
se2_prune_FIVW = NULL
nsnps<-NULL

z1<-beta1_prune_FIVW/se1_prune_FIVW
z2<-beta2_prune_FIVW/se2_prune_FIVW

p1<-2*pnorm(abs(z1),lower.tail=FALSE)
p2<-2*pnorm(abs(z2),lower.tail=FALSE)

exposure1<-c("AA:DGLA lnD5Dpooled","FADS1 expression in B cells","FADS1 expression in CD4 T cells","FADS1 expression in Monocytes","FADS1 expression in Blood","FADS1 expression in CD8 T cells","FADS1 expression in NK cells")

exposure2<-c("GLA:LA lnD6Dpooled","FADS2 expression in B cells","FADS2 expression in CD4 T cells","FADS2 expression in Monocytes","FADS2 expression in Blood","FADS2 expression in CD8 T cells","FADS2 expression in NK cells")

Res<-data.frame(matrix(c(exposure1,exposure2,beta1_prune_FIVW,se1_prune_FIVW,z1,p1,beta2_prune_FIVW,se2_prune_FIVW,z2,p2,nsnps),byrow=FALSE,nrow=length(beta1_prune_FIVW),ncol=11))
names(Res)<-c("exposure1","exposure2","beta1_prune_FIVW","se1_prune_FIVW","z1","p1","beta2_prune_FIVW","se2_prune_FIVW","z2","p2","nsnps")
Res$outcome<-"Colorectal cancer"
Res$panel<-"CHB"

write.table(Res,"~/fatty-acids/mr/results/cismvMR_crc_accc_chb.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# j<-1
for (j in 1:7) {

	if (j == 1) { D5D <- Dat[Dat$trait=="AA:DGLA lnD5Dpooled",]; D6D <- Dat[Dat$trait=="GLA:LA lnD6Dpooled",] }
	if (j == 2) { D5D <- Dat[Dat$trait=="FADS1 expression in B cells",];     D6D <- Dat[Dat$trait=="FADS2 expression in B cells",] }
	if (j == 3) { D5D <- Dat[Dat$trait=="FADS1 expression in CD4 T cells",]; D6D <- Dat[Dat$trait=="FADS2 expression in CD4 T cells",] }
	if (j == 4) { D5D <- Dat[Dat$trait=="FADS1 expression in Monocytes",];                    D6D <- Dat[Dat$trait=="FADS2 expression in Monocytes",] }
	if (j == 5) { D5D <- Dat[Dat$trait=="FADS1 expression in Blood",];              D6D <- Dat[Dat$trait=="FADS2 expression in Blood",] }
	if (j == 6) { D5D <- Dat[Dat$trait=="FADS1 expression in CD8 T cells",];            D6D <- Dat[Dat$trait=="FADS2 expression in CD8 T cells",] }
	if (j == 7) { D5D <- Dat[Dat$trait=="FADS1 expression in NK cells",];                    D6D <- Dat[Dat$trait=="FADS2 expression in NK cells",] }

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

 