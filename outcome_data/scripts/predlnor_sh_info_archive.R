Dat<-load_data_info()

Dat<-Dat[!Dat$study %in% c("UCSF_MAYO","UKB","UCSF_AGS + SFAGS"),]

Dat<-Dat[abs(Dat$bias)<400,]
dat<-Dat
Study<-dat$study
dat$Y<-abs(dat$bias)
dat$Y<-dat$bias
dat$X<-dat$mean_info

Model<-summary(lm(Y~X+study,dat))

summary(dat$mean_info)

dat[which(dat$mean_info>1.10),c("info","info1","info2","info3","mean_info")]

Xlab<-"Imputation info or r2 metric"
Ylab<-"Absolute bias in predicted log odds ratio"

ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=Study))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))

unique(Dat$study[duplicated(paste(Dat$rsid,Dat$study))])
plot(abs(Dat1$bias),Dat1$mean_info)

sort(unique(Dat$study))

Dat1<-Dat[Dat$ID == 23,]
Dat
plot(Dat1$bias,Dat1$mean_info)
Dat$study

summary(lm(abs(bias) ~ mean_info+study,Dat))