Dat<-format_digestive()
Dat[Dat$Cancer.Group=="Esophageal cancer",c("cancer","id.outcome","study.abbreviation","Cancer.Group")]
mr_res[mr_res$Cancer.Group == "Esophageal cancer",c("cancer","id.outcome","study.abbreviation","Cancer.Group")]
mr_res2[mr_res2$Cancer.Group == "Esophageal cancer",c("cancer","id.outcome","study.abbreviation","Cancer.Group")]

Dat1<-mr_res[mr_res$cancer == "Esophageal squamous cell carcinoma",c("study.abbreviation","cases","controls","b","se","pval")]
Dat2<-mr_res2[mr_res2$Cancer.Group == "Esophageal cancer",c("study.abbreviation","cases","controls","b","se","pval")]
Dat2<-Dat2[Dat2$study.abbreviation != "UKB",]
Dat<-rbind(Dat1,Dat2)

Dat<-meta_analysis(Dat)
Dat1<-mr_res[mr_res$Cancer.Group == "Lung cancer",c("cancer","OR","LCI","UCI","pval","cases","study.abbreviation")]
Dat2<-mr_res2[mr_res2$Cancer.Group == "Lung cancer",c("cancer","OR","LCI","UCI","pval","cases","study.abbreviation")]
Dat<-rbind(Dat1,Dat2)
Dat[order(Dat$cancer,Dat$cases,decreasing=T),]

meta_analysis<-function(Dat=NULL){	
		b<-Dat$b
		se<-Dat$se
		# p<-temp$p
		w<-1/se^2
		b.fixed<-sum(b*w)/(sum(w))
		se.fixed<-sqrt(sum(w)^-1)
		z<-abs(b.fixed/se.fixed)
		p.fixed<-pnorm(z,lower.tail=F)*2
		nstudies.fixed<-length(b)
		cases<-sum(Dat$cases)
		controls<-sum(Dat$controls)
		study<-"Overall fixed effect"
		Q<-sum((b.fixed-b)^2*w)
		df.Q<-length(b)-1		
		Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)		
		return(matrix(c(b.fixed,se.fixed,p.fixed,nstudies.fixed,cases,controls,study,Q,df.Q,Q.p)))
}