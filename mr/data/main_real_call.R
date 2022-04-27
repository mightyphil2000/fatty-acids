## This file is to call the functions in main_real
## After data preparation/manipulate, prepare the input, then source main_real.R and main_real_call.R


## Estimate factors
r <- est.factor(cormat )
r
## Factor IVW using cor
fcor_out<- FIVW.cor(cormat, bx, by, sey , p, r)
cat( "betas afm: \n",  fcor_out$fmvivw[,1], "\n")
cat( "se afm: \n", fcor_out$fmvivw[,2], "\n")
# fcor_out$fmvivw[,1]*(sd_x/sd_y)


#PCA IVW
#r <- FIVW_out$est.no.factor
#r <- est.factor(gr) 
pca_out <- IVWPCA(cormat, bx,  sey, r)
cat( "betas pca: \n", pca_out[,1], "\n")
cat( "se pca: \n", pca_out[,2], "\n")
# pca_out[,1]*(sd_x/sd_y)

#IVW
ivw_out <- ivw(cormat, bx, sex, by, sey)
cat( "betas ivw: \n", ivw_out[,1], "\n")
cat( "se ivw: \n", ivw_out[,2], "\n")
#stand_ivw_beta[i,] <- ivw_beta[i,]*(sd_x/sd_y)


round(fcor_out$fmvivw[,1], 4)
round(fcor_out$fmvivw[,2], 4)
round(pca_out[,1], 4)
round(pca_out[,2], 4)
round(ivw_out[,1], 4)
round(ivw_out[,2], 4)