source("~/fatty-acids/mr/scripts/mr_functions.R")

res_csi<-format_res_csi()
res_csi$id.outcome
res_csi[,c()]
snp_csi<-format_snp_csi()
decon_res<-csi_decon()

write.table(decon_res,"~/fatty-acids/mr/results/decon_csi_cancer.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

# snp_csi_effect<-snp_csi$beta*-1
	# snp_csi_effect_lci<-snp_csi$uci*-1
	# snp_csi_effect_uci<-snp_csi$lci*-1 

	# res_csi$b is effect of csi on cancer
	# snp_csi_effect is effect of SNP on csi 




# # effect of aadgla on csi
# ac<-read.table("~/fatty-acids/mr/results/mr_results_aadgla_csi.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)
# # change in csi per Sd increase in AA:DGLA
# ac$b
# ac$lci <- ac$b-1.96*ac$se
# ac$uci <- ac$b+1.96*ac$se

# scale_factor1<-1/ac$lci
# scale_factor2<-1/ac$b
# scale_factor3<-1/ac$uci
# change in log odds per 1-unit increase in csi 

# rs174546 effect on CSI 



# use product of coefficients method to scale result, which is equivalent to estimating the indirectione effect ie effect of SNP on cancer mediated by CSI
#express per copy of the minor allele, which is also the allele associated with higher csi