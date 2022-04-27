source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")

# load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
# load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")

Dat<-load_data_info()

Dat<-Dat[!Dat$study %in% c("UCSF_MAYO","UKB","UCSF_AGS + SFAGS"),]

Dat<-Dat[abs(Dat$bias)<400,]
dat<-Dat
Study<-dat$study
dat$Y<-abs(dat$bias)
# dat$Y<-dat$bias
dat$X<-dat$mean_info

Model<-summary(lm(Y~X+study,dat))

Xlab<-"Imputation info or r2 metric"
Ylab<-"% absolute bias in predicted log odds ratio"
Title<-""
Title_size<-0

png("~/fatty-acids/outcome_data/results/plots/predicted_lnor/info_scores.png",width=1000,height=500)

ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Study))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=10))
dev.off()
