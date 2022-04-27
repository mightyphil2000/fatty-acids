source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
library(ggforestplot)
library(ggplot2)

Cancers<-c("Colorectal cancer","Lung cancer","skin cancer","all cause")
mr_res2<-mr_res1[mr_res1$study.abbreviation == "UKB",]
Pos<-unlist(lapply(Cancers,FUN=function(x) grep(x,mr_res2$outcome)))
mr_res3<-mr_res2[Pos,]

mr_res4<-mr_res3[mr_res3$exposure=="AA:DGLA",c("b","se","pval","cases","controls","exposure","outcome")]
mr_res4<-mr_res4[mr_res4$outcome != "Malignant non-melanoma skin cancer 156" ,] 
mr_res4$outcome<-c("Colorectal cancer","Lung cancer","Pan-cancer\nexcl non-melanoma SC","Malignant skin cancer","Pan-cancer")
mr_res4$outcome[mr_res4$outcome=="Pan-cancer\nexcl non-melanoma SC"]<-"Overall cancer\nexcl non-melanoma SC"
mr_res4$outcome[mr_res4$outcome=="Pan-cancer"]<-"Overall cancer"

load("~/fatty-acids/mr/results/sibling_results.Rdata")
load("~/fatty-acids/mr/results/mr_crc_lc_msc_d5d.Rdata")
Model
b<-Model$coefficients[2,1]
se<-Model$coefficients[2,2]
pval<-Model$coefficients[2,4]
mr_selected_cancers<-data.frame(matrix(c("\"Sig\" cancers",b,se,pval,  25424,270330,"AA:DGLA" ),nrow=1,ncol=7))
names(mr_selected_cancers)<-c("outcome","b","se","pval","cases","controls","exposure")
mr_res5<-rbind(mr_res4,mr_selected_cancers)
mr_res5$design<-"between unrelated\nindividuals"

Out2<-Out[,c("PHEN", "BETA_WF","SE_BETA_WF","P_BETA_WF","ncase","ncont")]
Out2$PHEN<-c("Overall cancer\nexcl non-melanoma SC","Overall cancer","Lung cancer","Colorectal cancer","Malignant skin cancer","\"Sig\" cancers")
names(Out2)<-c("outcome","b","se","pval","cases","controls")
Out2$exposure <- "AA:DGLA"
Out2$design <- "within siblings"

# Out$PHEN
Out3<-Out[,c("PHEN", "BETA_BF","SE_BETA_BF","P_BETA_BF","ncase","ncont")]
Out3$PHEN<-c("Overall cancer\nexcl non-melanoma SC","Overall cancer","Lung cancer","Colorectal cancer","Malignant skin cancer","\"Sig\" cancers")
names(Out3)<-c("outcome","b","se","pval","cases","controls")
Out3$exposure <- "AA:DGLA"
Out3$design <- "between families"

mr_res5$order[order(as.numeric(mr_res5$cases),decreasing=TRUE)]<-1:nrow(mr_res5)
Out2$order[order(as.numeric(Out2$cases),decreasing=TRUE)]<-1:nrow(Out2)
Out3$order[order(as.numeric(Out3$cases),decreasing=TRUE)]<-1:nrow(Out3)
res<-rbind(mr_res5,Out2)
# res<-rbind(res,Out3)
res<-res[order(res$order),]
res$b<-round(as.numeric(res$b),3)
res$se<-round(as.numeric(res$se),3)
res$pval<-round(as.numeric(res$pval),3)
res$or<-round(exp(res$b),3)
res$lci<-round(exp(res$b-1.96*res$se),3)
res$uci<-round(exp(res$b+1.96*res$se),3)

res[res$design == "within siblings",c("outcome","cases","or","lci","uci","pval")]
# res$order<-c(5,6,9,10,11,12,7,8,1,2,3,4)
# res<-res[order(res$order,res$cases),]
# res[,c("outcome","design")]


Plot<-plot_dat(Dat=res,colour="design")

png("~/fatty-acids/mr/results/plots/within_sibling.png",width=1000,height=1000)
	Plot
dev.off()

plot_dat<-function(Dat=NULL,text.names=20,text.title=25,Shape=NULL,colour=NULL){ 

	Title.plot<-"OR (95% confidence interval) per SD increase in FADS1/2 activity"
    Dat$Design<-Dat[,colour]
    Dat$outcome[Dat$outcome %in% "\"Sig\" cancers"]<-"Cancer (3 sites)"
    Dat$outcome2<-paste0(Dat$outcome,"\nNcases=",Dat$cases)  

    Dat$b<-as.numeric(Dat$b)
    Dat$se<-as.numeric(Dat$se)
    p<-forestplot(df = Dat,
            logodds = TRUE,
            name=outcome2,
                  estimate=b,
                  se=se,
                  shape=Shape,
                  colour = Design,
                   xlab = "")+
		    theme(legend.position = "none")+
            # labs(title=Title.plot,size=1)+
            # theme(plot.title = element_text(size = text.title))+
            theme(text = element_text(size=text.names))+
            geom_point(size=6,shape=15,colour=c("black","red","black","red","black","red","black","red","black","red","black","red"))
    return(p)
}
# 
# geom_point(shape=mr_res$shape,size=mr_res$weight,fill=c("black"))+