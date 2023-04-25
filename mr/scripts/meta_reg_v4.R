library(plyr)
library(ggforestplot)
library(ggplot2)
library(meta)
library(metafor)
library(gsheet)

source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")

mr_res<-format_metareg_v2()
sort(mr_res$cancer)
# mr_res<-format_metareg_v2(restrict_to_european_studies=TRUE)

# mr_res<-mr_res[mr_res$population != "European",]
M<-create_correlation_matrix(study=mr_res)
M<-M$matrix

rownames(M)<-mr_res$study.abbreviation
colnames(M)<-mr_res$study.abbreviation
# corr_results_list2<-res$corr_results_list2
mr_res$se_decoupled<-decoupling(s=mr_res$se,C=M)
Weights<-1/mr_res$se^2
Weights_decoupled<-1/mr_res$se_decoupled^2

#########
#smoking#
######### 

table(mr_res$smoking1)
dim(mr_res)
# recode endometrial as smoking related cancer

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
m<-extract_meta_reg_results(dat=Model)
m$model<-"smoking1"

lc<-including_laryngeal_cancer()
b1<-lc$lnor
se1<-lc$se
w1<-1/se1^2
sm1<-1
can1<-"laryngeal cancer"
Model1<-summary(rma.uni(yi=as.numeric(c(mr_res$b,b1)),sei=as.numeric(c(mr_res$se,se1)),weights=c(Weights,w1),mods=c(mr_res$smoking1,sm1),intercept=TRUE,slab=c(mr_res$cancer,can1),method="REML",weighted=TRUE))
m1<-extract_meta_reg_results(dat=Model1)
m1$model<-"incl_laryngeal_cancer"

#lnor=0.0492, se=0.0167 0.0033 with corrected effect sizes
#lnor=0.0610, se=0.0232; 0.0085 with uncorrected effect sizes

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
md<-extract_meta_reg_results(dat=Model_decoupled)
md$model<-"smk_decoupled_se"

mr_res2<-mr_res
mr_res2$smoking1[mr_res2$cancer == "Endometrial cancer"]<-1
Model2<-summary(rma.uni(yi=as.numeric(mr_res2$b),sei=as.numeric(mr_res2$se),weights=Weights,mods=mr_res2$smoking1,intercept=TRUE,slab=mr_res2$cancer,method="REML",weighted=TRUE))
m2<-extract_meta_reg_results(dat=Model2)
m2$model<-"recode_endometrial"

m_smk<-do.call(rbind,list(m,md,m1,m2))


caption<-paste0("R2=",round(Model$R2,2),"% p=",round(Model$pval[2],4))

plot_dat<-plot_meta_reg_smoking()

# lnor<-plot_dat$b[is.na(plot_dat$weight)]
# se<-plot_dat$se[is.na(plot_dat$weight)]
# or<-round(exp(lnor),2)
# lci<-round(exp(lnor-1.96*se),2)
# uci<-round(exp(lnor+1.96*se),2)
# matrix(c(or,lci,uci),nrow=2,ncol=3)

table(plot_dat$Colour1)
# names(plot_dat)
# plot_dat<-plot_dat[,names(plot_dat)!="Colour1"]
P1<-forestplot(df = plot_dat,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = NULL,xlab = "") +
		geom_point(shape=plot_dat$Shape,size=1/plot_dat$se/10,fill=c(plot_dat$Colour),colour = plot_dat$Colour)+
		theme(legend.position = "none",text = element_text(size=20))

		# +labs(title=caption)+theme(size=1,plot.title = element_text(hjust = 1))


png("~/fatty-acids/mr/results/plots/ggforest_metareg_smoking_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()



#######################
# inflammatory cancers#
#######################
table(mr_res$infl)
Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$infl,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
m<-extract_meta_reg_results(dat=Model)
m$model<-"infl"

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$infl,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
md<-extract_meta_reg_results(dat=Model_decoupled)
md$model<-"infl_decoupled_se"

plot_dat<-plot_meta_reg_infl()
# lnor<-plot_dat$b[is.na(plot_dat$weight)]
# se<-plot_dat$se[is.na(plot_dat$weight)]
# or<-round(exp(lnor),2)
# lci<-round(exp(lnor-1.96*se),2)
# uci<-round(exp(lnor+1.96*se),2)
# matrix(c(or,lci,uci),nrow=2,ncol=3)

P1<-forestplot(df = plot_dat,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = NULL,xlab = "") +
		geom_point(shape=plot_dat$Shape,size=1/plot_dat$se/10,fill=plot_dat$Colour,colour = plot_dat$Colour)+
		theme(legend.position = "none")+
		theme(text = element_text(size=20))

length(plot_dat$Shape)
png("~/fatty-acids/mr/results/plots/ggforest_metareg_infl_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()



###########################################
# inflammatory plus infectious age cancers#
###########################################

table(mr_res$infl_agent)
table(mr_res$infl)

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$infl_agent,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

m1<-extract_meta_reg_results(dat=Model)
m1$model<-"infl_plus_agent"


Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$infl_agent,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

m1md<-extract_meta_reg_results(dat=Model_decoupled)
m1md$model<-"infl_plus_agent_decoupled_se"
m_infl<-do.call(rbind,list(m,md,m1,m1md))


plot_dat<-plot_meta_reg_infl_agent()
# lnor<-plot_dat$b[is.na(plot_dat$weight)]
# se<-plot_dat$se[is.na(plot_dat$weight)]
# or<-round(exp(lnor),2)
# lci<-round(exp(lnor-1.96*se),2)
# uci<-round(exp(lnor+1.96*se),2)
# matrix(c(or,lci,uci),nrow=2,ncol=3)

P1<-forestplot(df = plot_dat,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = NULL,xlab = "") +
		geom_point(shape=plot_dat$Shape,size=1/plot_dat$se/10,fill=plot_dat$Colour,colour = plot_dat$Colour)+
		theme(legend.position = "none")+
		theme(text = element_text(size=20))

length(plot_dat$Shape)
png("~/fatty-acids/mr/results/plots/ggforest_metareg_infl_agents_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()


########################################
# digestive system versus other systems#
########################################
# mouth and throat cancer classified as digestive system cancer

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$system2_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
m<-extract_meta_reg_results(dat=Model)
m$model<-"digestive_system"


table(mr_res$system,mr_res$system2_num)


Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$system2_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
md<-extract_meta_reg_results(dat=Model_decoupled)
md$model<-"digestive_system_decoupled_se"
m_dig<-do.call(rbind,list(m,md))

plot_dat<-plot_meta_reg_digestive()
# lnor<-plot_dat$b[is.na(plot_dat$weight)]
# se<-plot_dat$se[is.na(plot_dat$weight)]
# or<-round(exp(lnor),2)
# lci<-round(exp(lnor-1.96*se),2)
# uci<-round(exp(lnor+1.96*se),2)
# matrix(c(or,lci,uci),nrow=2,ncol=3)

P1<-forestplot(df = plot_dat,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = NULL,xlab = "") +
		geom_point(shape=plot_dat$Shape,size=1/plot_dat$se/10,fill=plot_dat$Colour,colour = plot_dat$Colour)+
		theme(legend.position = "none")+
		theme(text = element_text(size=20))

length(plot_dat$Shape)
png("~/fatty-acids/mr/results/plots/ggforest_metareg_digestive_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()


# mouth and throat cancer classified as respiratory cancer
Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$system_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$system_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))



#######################################
# external versus internal environment#
#######################################

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$external_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$external_num,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))


########################################################
# Cancer level characteristics from SEER and Vogelstein#
########################################################

mr_res<-format_metareg_v2(vogelstein=TRUE)
dim(mr_res)
# mr_res[,c("cancer","stem_cell_divisions")]
# mr_res<-mr_res[mr_res$population != "European",]
M<-create_correlation_matrix(study=mr_res)
M<-M$matrix
rownames(M)<-mr_res$study.abbreviation
colnames(M)<-mr_res$study.abbreviation
# corr_results_list2<-res$corr_results_list2
mr_res$se_decoupled<-decoupling(s=mr_res$se,C=M)
Weights<-1/mr_res$se^2
Weights_decoupled<-1/mr_res$se_decoupled^2

#########################################
#lifetime number of stem cell divisions #
#########################################

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$stem_cell_divisions,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
# p_lsc<-Model$pval[2]
m<-extract_meta_reg_results(dat=Model)
m$model<-"stem_cell"

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$stem_cell_divisions,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

md<-extract_meta_reg_results(dat=Model_decoupled)
md$model<-"stem_cell_decoupled_se"
m_sc<-do.call(rbind,list(m,md))

caption<-paste0("R2=",Model$R2,"% p=",round(Model$pval[2],2))
Lab<-toupper(substring(mr_res$cancer,1,3))
P1<-ggplot(mr_res, aes(x=stem_cell_divisions, y=b)) + 
  geom_point(size=1/mr_res$se/5) + geom_text(label=Lab,hjust=0.1,vjust=2.2,size=10)+scale_x_continuous(trans='log10')+ theme(text = element_text(size=25))+labs(x = "Cumulative number of stem cell divisions per lifetime (log10 scale)",y="log odds ratio for cancer per 1-SD change in genetically elevated PUFA desaturase activity",title=caption)+theme(plot.title = element_text(hjust = 1))

png("~/fatty-acids/mr/results/plots/ggforest_metareg_stem_cell_division_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()

write.table(data.frame(matrix(c(Lab,mr_res$cancer),nrow=length(Lab),ncol=2)),"~/fatty-acids/mr/results/plots/Label_key_ggforest_metareg_stem_cell_division_mr_results_rep_v3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)



#####################################
# cancer level characteristics SEERS#
#####################################

mr_res<-format_metareg_v2(seer=TRUE)
M<-create_correlation_matrix(study=mr_res)
M<-M$matrix
rownames(M)<-mr_res$study.abbreviation
colnames(M)<-mr_res$study.abbreviation
# corr_results_list2<-res$corr_results_list2
mr_res$se_decoupled<-decoupling(s=mr_res$se,C=M)
Weights<-1/mr_res$se^2
Weights_decoupled<-1/mr_res$se_decoupled^2
###################
# cancer incidence#
###################
# rma.uni
Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$incidence,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
m<-extract_meta_reg_results(dat=Model)
m$model<-"incidence"

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$incidence,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

md<-extract_meta_reg_results(dat=Model_decoupled)
md$model<-"incidence_decoupled_se"

caption<-paste0("R2=",Model$R2,"% p=",round(Model$pval[2],2))
Lab<-toupper(substring(mr_res$cancer,1,3))

# scale_x_continuous(trans='log10')+
P1<-ggplot(mr_res, aes(x=incidence, y=b)) + 
  geom_point(size=1/mr_res$se/5) + geom_text(label=Lab,size=10,vjust=-0.9,hjust=0.1,)+ theme(text = element_text(size=25))+labs(x = "Incidence per 100,000 per year in SEER",y="log odds ratio for cancer per 1-SD change in genetically elevated PUFA desaturase activity",title=caption)+theme(plot.title = element_text(hjust = 1))

png("~/fatty-acids/mr/results/plots/ggforest_metareg_incidence_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()

write.table(data.frame(matrix(c(Lab,mr_res$cancer),nrow=length(Lab),ncol=2)),"~/fatty-acids/mr/results/plots/Label_key_ggforest_metareg_incidence_mr_results_rep_v3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

# ####################
# # cancer incidence#
# ###################
# Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$incidence,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

# Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$incidence,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

################
# survival time#
################

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$survival_time,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

m1<-extract_meta_reg_results(dat=Model)
m1$model<-"survival_time"

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$survival_time,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
m1d<-extract_meta_reg_results(dat=Model_decoupled)
m1d$model<-"survival_time_decoupled_se"

caption<-paste0("R2=",Model$R2,"% p=",round(Model$pval[2],2))
Lab<-toupper(substring(mr_res$cancer,1,3))


P1<-ggplot(mr_res, aes(x=survival_time, y=b)) + 
  geom_point(size=1/mr_res$se/5) + geom_text(label=Lab,size=10,vjust=-0.9,hjust=0.1)+ theme(text = element_text(size=25))+labs(x = "5 year survival rate (%) in SEER",y="log odds ratio for cancer per 1-SD change in genetically elevated PUFA desaturase activity",title=caption)+theme(plot.title = element_text(hjust = 1))

png("~/fatty-acids/mr/results/plots/ggforest_metareg_survival_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()

write.table(data.frame(matrix(c(Lab,mr_res$cancer),nrow=length(Lab),ncol=2)),"~/fatty-acids/mr/results/plots/Label_key_ggforest_metareg_survival_mr_results_rep_v3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

######################
#Median age diagnosis#
######################

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$median.age.diagnosis,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
m2<-extract_meta_reg_results(dat=Model)
m2$model<-"median_age_diagnosis"

Model_decoupled<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se_decoupled),weights=Weights_decoupled,mods=mr_res$median.age.diagnosis,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))
m2d<-extract_meta_reg_results(dat=Model_decoupled)
m2d$model<-"median_age_diagnosis_decoupled_se"
m_seer<-do.call(rbind,list(m,md,m1,m1d,m2,m2d))

caption<-paste0("R2=",Model$R2,"% p=",round(Model$pval[2],2))
Lab<-toupper(substring(mr_res$cancer,1,3))

P1<-ggplot(mr_res, aes(x=median.age.diagnosis, y=b)) + 
  geom_point(size=1/mr_res$se/5) + geom_text(label=Lab,size=10,vjust=-0.9,hjust=0.1)+ theme(text = element_text(size=25))+labs(x = "Median age at diagnosis in SEER",y="log odds ratio for cancer per 1-SD change in genetically elevated PUFA desaturase activity",title=caption)+theme(plot.title = element_text(hjust = 1))

png("~/fatty-acids/mr/results/plots/ggforest_metareg_age_mr_results_rep_v3.png", width = 1200, height = 1600)
	print(P1) 
dev.off()

write.table(data.frame(matrix(c(Lab,mr_res$cancer),nrow=length(Lab),ncol=2)),"~/fatty-acids/mr/results/plots/Label_key_ggforest_metareg_age_mr_results_rep_v3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


#####
# collate all results
##
a<-do.call(rbind,list(m_smk,m_infl,m_dig,m_sc,m_seer))
a<-a[,c(which(names(a) == "model"),1:(ncol(a)-1))]
write.table(a,"~/fatty-acids/mr/results/meta_reg_results.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

