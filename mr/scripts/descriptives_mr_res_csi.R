load("~/fatty-acids/mr/results/res_csi.Rdata")

res_csi[res_csi$outcome2 %in% c("Lung cancer","Colorectal cancer","Basal cell carcinoma","Overall cancer"), ]
res_csi$b<-as.numeric(res_csi$b)
res_csi$se<-as.numeric(res_csi$se)
res_csi$or<-round(exp(res_csi$b),2)
res_csi$lci<-round(exp(res_csi$b-1.96*res_csi$se),2)
res_csi$uci<-round(exp(res_csi$b+1.96*res_csi$se),2)


# effect of aadgla on csi
ac<-read.table("~/fatty-acids/mr/results/mr_results_aadgla_csi.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)
# change in csi per Sd increase in AA:DGLA
ac$b
ac$lci <- ac$b-1.96*ac$se
ac$uci <- ac$b+1.96*ac$se

scale_factor<-1/ac$uci
# change in log odds per 1-unit increase in csi 
res_csi$or2<-round(exp(res_csi$b/scale_factor),2)
res_csi$lci2<-round(exp(res_csi$b/scale_factor-1.96*res_csi$se/scale_factor),2)
res_csi$uci2<-round(exp(res_csi$b/scale_factor+1.96*res_csi$se/scale_factor),2)

res_csi[res_csi$outcome2 %in% c("Lung cancer","Colorectal cancer","Basal cell carcinoma","Overall cancer"), ]
