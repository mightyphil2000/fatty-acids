source("~/fatty-acids/mr/scripts/mr_functions.R")
library(ieugwasr)
ao<-ieugwasr::gwasinfo()
# snps_gscan the effect of rs174546 and rs2524299 on smoking in GSCAN 
# gscan_cancer the effect of smoking on cancer
# snp1_cancer the effect of rs174546 and rs2524299 on cancer (the total effect )
# we want to work out the indirect effect of rs174546 and rs2524299 on cancer mediated by smoking using the product of coefficients method . This is the smoking-cancer effect * the snp-smoking effect
# given the causal effect of smoking on cancer, what is the effect of the SNPs on cancer mediated by their effect on smoking? We essentially scale the SNP-cancer effect to reflect the scale of the SNP-smoking change and the magnitude of the smoking-cancer effect. In theory, this effect is the indirect effect and should be the same as the observed SNP-cancer effect if the relationship to cancer is entirely mediated by smoking

snps_gscan<-read.table("~/fatty-acids/mr/results/snplookups_rs174546_rs2524299_gscan.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)
gscan_cancer<-read.table("~/fatty-acids/mr/results/mr_res_gscan_cancer.txt",sep="\t",stringsAsFactors=FALSE,head=TRUE)
dec_res<-gscan_decon(method="Inverse variance weighted")
write.table(dec_res,"~/fatty-acids/mr/results/deconv_gscan_cancer.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# ao<-ieugwasr::gwasinfo()
# out_dat duplicates lung cancer data in snp1_cancer and snp2_cancer
# out_dat <- extract_outcome_data(snps = c("rs174546","rs2524299"),outcomes = c("ieu-a-985","ieu-a-986","ieu-a-987")
# )

# out_dat$or<-exp(out_dat$beta)
# out_dat$lci<-exp(out_dat$beta-1.96*out_dat$se)
# out_dat$uci<-exp(out_dat$beta+1.96*out_dat$se)



# out_dat2<-merge(ao[,c("note","id")],out_dat,by.x="id",by.y="id.outcome")
# out_dat2[out_dat2$SNP == "rs174546",c("SNP","outcome","or","lci","uci","note")]




