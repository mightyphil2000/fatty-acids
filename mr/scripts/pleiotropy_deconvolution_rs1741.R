source("~/fatty-acids/mr/scripts/mr_functions.R")
a<-read.table("~/fatty-acids/mr/results/pleiotropy_deconvolution_rs1741.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE)

a$ln_or_can_exp<-log(a$or_can_exp)
a$se_can_exp<-(log(a$or_can_exp_uci)-log(a$or_can_exp_lci))/(1.96*2)
a$ln_or_can_exp_sd<-a$ln_or_can_exp
a$se_can_exp_sd<-a$se_can_exp

a$lnor_decon<-a$ln_or_can_exp_sd*a$b_exp_sd 
a$se_decon<-a$se_can_exp_sd*a$b_exp_sd 
a$or_decon<-exp(a$lnor_decon)  
a$lci_decon<-exp(a$lnor_decon-1.96*a$se_decon)  
a$uci_decon<-exp(a$lnor_decon+1.96*a$se_decon)  

a$lnor_can_ea<-log(a$or_can_ea)
a$se_can_ea<-(log(a$or_can_ea_uci)-log(a$or_canc_ea_lci))/(1.96*2)

a$p_diff<-z_test_diff2(lnor_total=a$lnor_can_ea,se_total=a$se_can_ea,lnor_decon=a$lnor_decon,se_decon=a$se_decon)
a[,c("or_can_ea","or_decon","p_diff")]

a$lnor_decon_expected<-a$lnor_can_ea/a$b_exp_sd
a$se_decon_expected<-a$se_can_ea/a$b_exp_sd
a$or_decon_expected<-exp(a$lnor_decon_expected)
a$lci_decon_expected<-exp(a$lnor_decon_expected-1.96*a$se_decon_expected)
a$uci_decon_expected<-exp(a$lnor_decon_expected+1.96*a$se_decon_expected)

write.table(a,"~/fatty-acids/mr/results/pleiotropy_deconvolution_rs1741_results.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
