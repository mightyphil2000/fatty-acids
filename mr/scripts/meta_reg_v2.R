library(plyr)
library(ggforestplot)
library(ggplot2)
library(meta)
library(metafor)

source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")


mr_res<-format_metareg(prune_blood_cancers=FALSE,drop_skin_cancers=FALSE)

Weights<-1/mr_res$se^2

mr_res[order(mr_res$site),c("cancer","Cancer.Group","site","study.abbreviation.y","cases","controls")]
names(mr_res)
names(test)[names(test)=="study.abbreviation.y"]<-"study.abbreviation"
test<-mr_res[c(1,9),]

Model<-summary(rma.uni(yi=as.numeric(mr_res$b),sei=as.numeric(mr_res$se),weights=Weights,mods=mr_res$smoking1,intercept=TRUE,slab=mr_res$cancer,method="REML",weighted=TRUE))



# Plot

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


1/1.37


######################################
# include laryngeal squamous carcinoma#
########################################
laryngeal_cancer_lnor <-(log(0.71)*-1)/0.72 #lnor per 1 SD increase FADS1/2 activity 25194280, rs174549 r2=0.9786 with rs174546 in CHB+JPT (north east asians)
 
lnor<-log(0.71) #originally reported odds ratio
lci<-log(0.63)
uci<-log(0.80) 
se<-(uci-lci)/(1.96*2)
lnor_d6d<-lnor/0.72*-1 #originally reported odds ratio refers to allele that lowers d6d
se_d6d<-se/0.72

laryngeal_cancer_se <-((log(0.80)-log(0.63))/(1.96*2))/0.72 #0.72 is estimated effect of allele on D6D in East Asians 
pnorm(laryngeal_cancer_lnor/laryngeal_cancer_se,lower.tail=F)*2


exp(laryngeal_cancer_lnor)
exp(laryngeal_cancer_lnor-1.96*laryngeal_cancer_se)
exp(laryngeal_cancer_lnor+1.96*laryngeal_cancer_se)
 

# b1<-laryngeal_cancer_lnor
b1<-lnor_d6d
# se1<-laryngeal_cancer_se
se1<-se_d6d
Weights<-1/mr_res$se^2
Weights1<-1/se1^2

summary(mr_res$b)
table(mr_res$smoking1)

B<-as.numeric(c(mr_res$b,b1))
SE<-as.numeric(c(mr_res$se,se1))
W<-as.numeric(c(Weights,Weights1))
Smoking1<-c(mr_res$smoking1,1)
table(Smoking1)
Cancer<-c(mr_res$cancer,"laryngeal squamous carcinoma")

Model<-summary(rma.uni(yi=B,sei=SE,weights=W,mods=Smoking1,intercept=TRUE,slab=Cancer,method="REML",weighted=TRUE))

############################################
# exclude cancers with fewer than 500 cases#
############################################

mr_res1<-mr_res[mr_res$cases>500,]
B<-as.numeric(mr_res1$b)
SE<-as.numeric(mr_res1$se)
W<-1/mr_res1$se^2
Smoking1<-mr_res1$smoking1
Cancer<-mr_res1$cancer

table(Smoking1)

Model<-summary(rma.uni(yi=B,sei=SE,weights=W,mods=Smoking1,intercept=TRUE,slab=Cancer,method="REML",weighted=TRUE))


# inflammatory conditions
# mr_res<-format_metareg(prune_blood_cancers=FALSE,drop_skin_cancers=FALSE)
mr_res<-format_metareg2()
B<-as.numeric(mr_res$b)
SE<-as.numeric(mr_res$se)
W<-1/mr_res$se^2
Infl<-mr_res$infl
mr_res$cancer[Infl==1]
mr_res$cancer[Infl==0]
table(Infl)
Cancer<-mr_res$cancer
Weights<-1/mr_res$se^2
Model<-summary(rma.uni(yi=B,sei=SE,weights=Weights,mods=Infl,intercept=TRUE,slab=Cancer,method="REML",weighted=TRUE))
# Model<-summary(rma.uni(yi=B[Infl==1],sei=SE[Infl==1],weights=Weights[Infl==1],intercept=TRUE,slab=Cancer[Infl==1],method="REML",weighted=TRUE))
# Model<-summary(rma.uni(yi=B[Infl==0],sei=SE[Infl==0],weights=Weights[Infl==0],intercept=TRUE,slab=Cancer[Infl==0],method="REML",weighted=TRUE))
round(c(exp(Model$b[2]),exp(Model$ci.lb[2]),exp(Model$ci.ub[2])),2)
round(c(exp(Model$b[1]),exp(Model$ci.lb[1]),exp(Model$ci.ub[1])),2)

exp(Model$b)
exp(0.0638)




metagen(TE=B[Infl==1],seTE=SE[Infl==1],studlab=Cancer[Infl==1],sm="OR",comb.random=TRUE,method.tau="REML")
metagen(TE=B[Infl==0],seTE=SE[Infl==0],studlab=Cancer[Infl==0],sm="OR",comb.random=TRUE,method.tau="REML")



# chronic inflammatory conditions + infectous agents
mr_res<-format_metareg3()
mr_res<-mr_res[mr_res$cases>500,]
B<-as.numeric(mr_res$b)
SE<-as.numeric(mr_res$se)
W<-1/mr_res$se^2
Infl<-mr_res$infl
mr_res$cancer[Infl==1]
mr_res$cancer[Infl==0]
table(Infl)
Cancer<-mr_res$cancer
Weights<-1/mr_res$se^2

Model<-summary(rma.uni(yi=B,sei=SE,weights=Weights,mods=Infl,intercept=TRUE,slab=Cancer,method="REML",weighted=TRUE))
round(c(exp(Model$b[2]),exp(Model$ci.lb[2]),exp(Model$ci.ub[2])),2)
round(c(exp(Model$b[1]),exp(Model$ci.lb[1]),exp(Model$ci.ub[1])),2)
