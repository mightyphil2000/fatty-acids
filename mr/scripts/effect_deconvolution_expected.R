load("~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")

mr_res1$b<-as.numeric(mr_res1$b)
mr_res1$se<-as.numeric(mr_res1$se)
Pos<-is.na(mr_res1$OR)
mr_res1$OR[Pos]<-exp(mr_res1$b[Pos])
mr_res1$LCI[Pos] <- exp(mr_res1$b[Pos]-1.96*mr_res1$se[
	Pos])
mr_res1$UCI[Pos] <- exp(mr_res1$b[Pos]+1.96*mr_res1$se[Pos])


mr_res1[mr_res1$outcome == "Lung cancer 75" ,c("OR","LCI","UCI")] 
mr_res1[grep("basal",mr_res1$outcome,ignore.case=TRUE),c("outcome","population","b","OR","LCI","UCI","cases","study.abbreviation","cases")]

names(mr_res1)

mr_res1[mr_res1$outcome== "Cancer (all cause)",]

