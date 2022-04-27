# library(devtools)
# library(TwoSampleMR)
library(plyr)
library(ggforestplot)
library(ggplot2)
library(meta)
library(metafor)

source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")


mr_res<-format_metareg(prune_blood_cancers=FALSE,drop_skin_cancers=FALSE)

mr_res<-format_metareg2()
mr_res$OR<-round(exp(as.numeric(mr_res$b)),2)
mr_res$LCI<-round(exp(as.numeric(mr_res$b)-1.96*as.numeric(mr_res$se)),2)
mr_res$UCI<-round(exp(as.numeric(mr_res$b)+1.96*as.numeric(mr_res$se)),2)

mr_res<-mr_res[order(as.numeric(mr_res$cases),decreasing=T),c("outcome","OR","LCI","UCI","pval","cases")]
mr_res<-mr_res[order(as.numeric(mr_res$power05),decreasing=T),c("outcome","OR","LCI","UCI","pval","cases")]

mr_res[as.numeric(mr_res$cases)>10000,]

table(mr_res$smoking1)


# P1<-
forestplot(df = mr_res,logodds = TRUE,name=plot_name,
				  estimate=b, se=se,
				  shape=NULL, colour = smoking,xlab = "")+
			theme(plot.title = element_text(size = ""))+
			theme(text = element_text(size=14))

dev.off()
# size=mr_res$weight
mr_res$Shape="15"
forestplot(df = mr_res,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = smoking,xlab = "")+
ggplot2::scale_shape_manual(
    values = c(22L))+ 
 theme(legend.position = "none")+
 theme(text = element_text(size=5))

 # geom_point(shape=15,size=mr_res$weight)
 ggplot2::scale_shape_manual(
    values = c(22L))+ 
 theme(legend.position = "none")+
 theme(text = element_text(size=5))+

# values = c(23L, 21L, 21L, 21L, 21L),
#     labels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS")

+
		geom_point(shape=15,size=5)+
		theme(plot.title = element_text(size = ""),text = element_text(size=5))



# metareg()
# ?metareg
# ?rma.uni

# Res1<-mr_res[mr_res$infl==1,]
# Res2<-mr_res[mr_res$infl==0,]

Res1<-mr_res[mr_res$smoking1==1,]
Res2<-mr_res[mr_res$smoking1==0,]
metagen(TE=mr_res$b,seTE=mr_res$se,studlab=mr_res$cancer,sm="OR",comb.random=TRUE,method.tau="REML")
# metagen(TE=mr_res$b,seTE=mr_res$se,studlab=mr_res$cancer,sm="OR",comb.random=TRUE,method.tau="DL")

Q_all<-metagen(TE=mr_res$b,seTE=mr_res$se,studlab=mr_res$cancer,sm="OR",comb.fixed=TRUE)
Q_all<-Q_all$Q

Q_smoke<-metagen(TE=Res1$b,seTE=Res1$se,studlab=Res1$cancer,sm="OR",comb.fixed=TRUE)
Q_smoke<-Q_smoke$Q
Q_nsmoke<-metagen(TE=Res2$b,seTE=Res2$se,studlab=Res2$cancer,sm="OR",comb.fixed=TRUE)
Q_nsmoke<-Q_nsmoke$Q

Q_all
Q_smoke
Q_nsmoke
Q_int<-Q_all-Q_smoke-Q_nsmoke

pchisq(q=Q_int,df=1,lower.tail=F)


# metagen(TE=Res1$b,seTE=Res1$se,studlab=Res1$cancer,sm="OR",comb.random=TRUE,method.tau="REML")
# metagen(TE=Res2$b,seTE=Res2$se,studlab=Res2$cancer,sm="OR",comb.random=TRUE,method.tau="REML")


# mr_res$weight<-1/mr_res$se
# mr_res$weight<-mr_res$cases
table(mr_res$smoking1)
table(mr_res$smoking2)
mr_res<-mr_res[mr_res$smoking2 %in% c(1,2),]
mr_res$cancer[mr_res$smoking2 ==2]
mr_res$cancer[mr_res$smoking2 ==0]
mr_res$b
laryngeal_cancer_lnor <-(log(0.71)*-1)/0.72 #lnor per 1 SD increase FADS1/2 activity 25194280, rs174549 r2=0.9786 with rs174546 in CHB+JPT (north east asians)
laryngeal_cancer_se <-((log(0.80)-log(0.63))/(1.96*2))/0.72 
pnorm(laryngeal_cancer_lnor/laryngeal_cancer_se,lower.tail=F)*2

exp(laryngeal_cancer_lnor)
exp(laryngeal_cancer_lnor-1.96*laryngeal_cancer_se)
exp(laryngeal_cancer_lnor+1.96*laryngeal_cancer_se)
 

b1<-laryngeal_cancer_lnor
se1<-laryngeal_cancer_se
Weights<-1/mr_res$se^2
Weights1<-1/se1^2
summary(mr_res$b)
table(mr_res$smoking1)

mr_res1<-mr_res[mr_res$cancer!="Colorectal cancer" & mr_res$cancer!="Lung cancer",]
Combinep_fisher(Datp=as.numeric(mr_res$pval[mr_res$smoking1==1]))
Combinep_fisher(Datp=as.numeric(mr_res$pval[mr_res$smoking1==0]))
Combinep_fisher(Datp=as.numeric(mr_res$pval))

Combinep_fisher<-function(Datp=NULL){
	Chi<--2*sum(log(Datp))
	df<-length(Datp)*2
	p_fisher<-pchisq(Chi, df,lower.tail=F)
	return(p_fisher)
}


Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="DL",weighted=TRUE))

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="FE",weighted=TRUE))

Model<-summary(rma.uni(yi=as.numeric(c(mr_res$b,b1)),sei=as.numeric(c(mr_res$se,se1)),weights=c(Weights,Weights1),mods=c(mr_res$smoking1,1),intercept=TRUE,slab=c(mr_res$cancer,"Laryngeal cancer"),method="REML",weighted=TRUE))

Model1<-summary(rma.uni(yi=as.numeric(c(mr_res$b,b1)),sei=as.numeric(c(mr_res$se,se1)),weights=c(Weights,Weights1),mods=c(mr_res$smoking1,1),intercept=TRUE,slab=c(mr_res$cancer,"Laryngeal cancer"),method="DL",weighted=TRUE))

Z<-0.0595/0.0091
Z<-0.0629/0.0090 
pnorm(Z,lower.tail=F)*2

round(c(exp(Model$b[2]),exp(Model$ci.lb[2]),exp(Model$ci.ub[2])),2)
round(c(exp(Model$b[1]),exp(Model$ci.lb[1]),exp(Model$ci.ub[1])),2)

Dat1<-c("Other cancers",Model$b[1],Model$se[1])
Dat2<-c("Smoking related cancers",Model$b[2],Model$se[2])
Dat<-data.frame(matrix(c(Dat2,Dat1),nrow=2,ncol=3,byrow=TRUE),stringsAsFactors=F)
names(Dat)<-c("plot_name","b","se")
Dat$smoking[1] <- "Smoking related cancers (overall effect)"
Dat$smoking[2] <- "Other cancers (overall effect)"
Dat$b<-as.numeric(Dat$b)
Dat$se<-as.numeric(Dat$se)
Dat$Shape<-"diamond"
Dat2<-rbind.fill(Dat[1,],Res1)
Dat3<-rbind.fill(Dat[2,],Res2)
mr_res2<-rbind(Dat2,Dat3)
mr_res2$Colour[is.na(mr_res2$Colour)]<-0
mr_res2$Colour<-mr_res2$Colour+1
mr_res2$Shape
mr_res2$cases[mr_res2$smoking=="Smoking related cancers (overall effect)"]<-sum(mr_res2$cases[mr_res2$smoking=="Increases risk"])
mr_res2$cases[mr_res2$smoking=="Other cancers (overall effect)"]<-sum(mr_res2$cases[mr_res2$smoking %in% c("No relationship","Decreases risk","Unclear relationship")])


mr_res2[,c("plot_name","b","se","smoking","cases","se")]
P1<-forestplot(df = mr_res2,logodds = TRUE,name=plot_name,
	  estimate=b, se=se,shape=NULL,
		colour = NULL,xlab = "")+
geom_point(shape=mr_res2$Shape,size=1/mr_res2$se/10,fill=c(mr_res2$Colour),colour = mr_res2$Colour)+
theme(legend.position = "none")+
theme(text = element_text(size=20))
png("~/fatty-acids/mr/results/plots/ggforest_metareg_smoking.png", width = 1200, height = 1600)
	print(P1) 
dev.off()


ggplot2::scale_shape_manual(
    values = c(22L),colour="black")+ 
 theme(legend.position = "none")+
 theme(text = element_text(size=5))


geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+
					theme(plot.title = element_text(size = ""),text = element_text(size=20))