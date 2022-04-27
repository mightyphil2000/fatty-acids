library(plyr)
library(ggforestplot)
library(ggplot2)
library(meta)
library(metafor)

source("~/fatty-acids/mr/scripts/mr_functions.R")
source("~/fatty-acids/mr/scripts/mr_plot_functions.R")


mr_res<-format_metareg(prune_blood_cancers=FALSE,drop_skin_cancers=FALSE)

load("~/fatty-acids/mr/results/res_csi.Rdata")

head(mr_res)

Dat<-merge(res_csi,mr_res,by.x="outcome2",by.y="cancer",)
Dat<-Dat[Dat$population.x == "European",]
# Dat<-merge(mr_res,res_csi,by.x="cancer",by.y="outcome2")
# Dat<-Dat[Dat$population.y == "European",]
Dat$diff<-as.numeric(Dat$ncase.outcome)-as.numeric(Dat$cases)
Dat<-Dat[abs(Dat$diff)<200,]

# Dat<-Dat[as.numeric(Dat$ncase.outcome)>500,]

Weights<-1/as.numeric(Dat$se.y)^2
# Weights1<-1/as.numeric(Dat$se.x)^2

Y<-as.numeric(Dat$b.y)
X<-as.numeric(Dat$b.x)
se.y<-as.numeric(Dat$se.y)
Cancer<-Dat$outcome2

Dat[,c("ncase.outcome","cases")]

Model<-summary(rma.uni(yi=Y,sei=se.y,weights=Weights,mods=X,intercept=TRUE,slab=Cancer,method="REML",weighted=TRUE))





rma(b.y~b.x,sei=Dat$se.y,weights=Weights,data=Dat)

Weights<-1/as.numeric(Dat$se.y)^2

Weights1<-1/as.numeric(Dat$se.x)^2
Weights2<-as.numeric(Dat$ncase.outcome)
Weights3<-as.numeric(Dat$case)
summary(lm(b.y ~ b.x,data=Dat,weights=Weights1))
summary(lm(b.y ~ b.x,data=Dat,weights=Weights2))
summary(lm(b.y ~ b.x,data=Dat,weights=Weights3))
summary(lm(b.y ~ b.x,data=Dat,weights=Weights4))




Plot
head(Dat)

Dat[Dat$effect.x>0.10,c("id.outcome.x","id.outcome.y","population.y","ncase.outcome")]




Plot
plot_weight<-1/as.numeric(Dat$se.x)
# plot_weight<-1/as.numeric(Dat$se.y)/5
# plot_weight<-as.numeric(Dat$ncase.outcome)/5000
# plot_weight<-as.numeric(Dat$case)/5000

Colour<-"black"
Title<-"title"
Title_size1<-0
Ylab<-"b.y"
Xlab<-"b.x"
Subtitle<-"subtitle"
Title_xaxis_size<-2
Subtitle_size1<-1
Dat$effect.x<-Dat$b.x
Dat$effect.y<-Dat$b.y
# Dat$effect.x<-exp(as.numeric(Dat$b.x))
# Dat$effect.y<-exp(as.numeric(Dat$b.y))

Plot<-ggplot2::ggplot(Dat, ggplot2::aes(x=effect.x, y=effect.y)) + ggplot2::geom_point(colour=Colour,size=plot_weight) +ggplot2::ggtitle(Title) +ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size1,hjust = 0))+ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = Subtitle_size1)) 