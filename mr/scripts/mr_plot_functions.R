


format_bystudy<-function(cancer=NULL,overall_effect=TRUE){
	# load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
	# mr_res2<-mr_res1[mr_res1$population == "European",]
	# mr_res2<-mr_res2[mr_res2$exposure == "AA:DGLA",]
	# mr_res3<-mr_res1[mr_res1$population != "European",]
	# mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	# mr_res1<-rbind(mr_res2,mr_res3)
	load("~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
	# mr_res1$cancer[grep(";",mr_res1$id.outcome)]
	mr_res1$cancer[is.na(mr_res1$cancer)]<-mr_res1$outcome[is.na(mr_res1$cancer)]
	mr_res1<-mr_res1[which(mr_res1$cancer == cancer),]	
	mr_res1$b<-as.numeric(mr_res1$b)
	mr_res1$se<-as.numeric(mr_res1$se)
	mr_res1$pval<-as.numeric(mr_res1$pval)
	mr_res1$OR<-round(exp(mr_res1$b),2)
	mr_res1$LCI<-round(exp(mr_res1$b-mr_res1$se*1.96),2)
	mr_res1$UCI<-round(exp(mr_res1$b+mr_res1$se*1.96),2)
	
	mr_res1<-mr_res1[order(as.numeric(mr_res1$cases),decreasing=T),]		
	mr_res1$study.abbreviation[is.na(mr_res1$study.abbreviation)]<-mr_res1$study[is.na(mr_res1$study.abbreviation)]
	mr_res1$study.abbreviation2<-mr_res1$study.abbreviation
	mr_res1$shape<-15
	mr_res1$shape[mr_res1$study.abbreviation == "Overall fixed effect"]<-23
	if(!overall_effect){
		mr_res1<-mr_res1[mr_res1$study.abbreviation !="Overall fixed effect",]
	}
	mr_res1$study.abbreviation<-paste0(mr_res1$study.abbreviation,"\ncases=",mr_res1$cases)	
	mr_res1<-mr_res1[order(as.numeric(mr_res1$cases),decreasing=T),]	
	mr_res1$weight<-1/mr_res1$se/5		

	return(mr_res1)
}



format_nulls<-function(power05=TRUE){
	load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
	mr_res1$b<-as.numeric(mr_res1$b)
	mr_res1$se<-as.numeric(mr_res1$se)
	mr_res1$pval<-as.numeric(mr_res1$pval)
	mr_res1$OR<-round(exp(mr_res1$b),2)
	mr_res1$LCI<-round(exp(mr_res1$b-mr_res1$se*1.96),2)
	mr_res1$UCI<-round(exp(mr_res1$b+mr_res1$se*1.96),2)
	# sort(unique(mr_res1$cancer))

	# mr_res1[mr_res1$id.outcome == 140,]
	Cancers<-unique(disc.tab9$cancer)
	if(power05){
		Cancers<-unique(disc.tab9$cancer[disc.tab9$power05>0.80])
	}

	# unique(mr_res1$cancer[mr_res1$cancer %in% Cancers])

	# disc.tab9$power10[disc.tab9$cancer == "Esophageal squamous cell carcinoma"]
	# length(which(disc.tab9$power10>0.8))
	Cancers2<-unique(mr_res1$outcome[mr_res1$pval<0.05/67])
	Pos<-regexpr("[:0-9:]",Cancers2)
	Cancers3<-Cancers2[Pos==-1]
	Cancers4<-Cancers2[Pos!=-1]
	Pos<-regexpr("[:0-9:]",Cancers4)
	Cancers4<-substring(Cancers4,1,Pos-2)
	Cancers2<-unique(c(Cancers3,Cancers4))
	Pos<-regexpr("[:0-9:]",mr_res1$outcome)
	mr_res2<-mr_res1[Pos!=-1,]
	mr_res3<-mr_res1[Pos==-1,]
	Pos<-regexpr("[:0-9:]",mr_res2$outcome)
	mr_res2$outcome<-substring(mr_res2$outcome,1,Pos-2)
	mr_res1<-rbind(mr_res2,mr_res3)
	mr_res1<-mr_res1[mr_res1$outcome %in% Cancers,]
	mr_res1<-mr_res1[!mr_res1$outcome %in% Cancers2,]
	mr_res1<-mr_res1[order(mr_res1$Cancer.Group),]
	mr_res2<-mr_res1[mr_res1$population == "European",]
	mr_res2<-mr_res2[mr_res2$exposure == "AA:DGLA",]
	mr_res3<-mr_res1[mr_res1$population != "European",]
	mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	mr_res1<-rbind(mr_res2,mr_res3)
	mr_res1<-mr_res1[order(as.numeric(mr_res1$cases),decreasing=T),]
	mr_res1<-mr_res1[!duplicated(mr_res1$outcome),]
	mr_res1$cancer<-mr_res1$outcome
	mr_res1$outcome<-paste0(mr_res1$outcome,"\ncases=",mr_res1$cases)
	
	mr_res1$Cancer.Group[mr_res1$Cancer.Group %in% c("Female reproductive cancers","Ovarian cancer")]<-"Female reproductive cancers"
	mr_res1$Cancer.Group[mr_res1$Cancer.Group %in% c("Colorectal cancer","Kidney cancer","Skin cancer")]<-"Other cancers"
	mr_res1<-mr_res1[order(mr_res1$Cancer.Group),]
	mr_res1<-mr_res1[order(as.numeric(mr_res1$cases),decreasing=T),]
	# mr_res1<-mr_res1[order(mr_res1$b),]
	# mr_res1<-mr_res1[mr_res1$cancer != "Cancer (excluding non-melanoma skin cancer)",]
	mr_res2<-mr_res1[mr_res1$Cancer.Group == "Other cancers",]
	mr_res3<-mr_res1[mr_res1$Cancer.Group != "Other cancers",]
	mr_res1<-rbind(mr_res3,mr_res2)
	mr_res1$weight<-1/mr_res1$se/12
	return(mr_res1)
}


format_all_discovery<-function(Power=NULL){
	load("~/fatty-acids/mr/results/mr_results_discovery_v2.Rdata")
	load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")

	mr_res_disc$OR<-round(exp(mr_res_disc$b),2)
	mr_res_disc$LCI<-round(exp(mr_res_disc$b-mr_res_disc$se*1.96),2)
	mr_res_disc$UCI<-round(exp(mr_res_disc$b+mr_res_disc$se*1.96),2)
	# sort(unique(mr_res_disc$cancer))

	# mr_res_disc[mr_res_disc$id.outcome == 140,]
	Cancers<-unique(disc.tab9$cancer)
	if(!is.null(Power)){

		Cancers<-unique(disc.tab9$cancer[disc.tab9[,Power]>0.80])
	}

	mr_res_disc<-mr_res_disc[mr_res_disc$cancer %in% Cancers,]
	
	# disc.tab9$power10[disc.tab9$cancer == "Esophageal squamous cell carcinoma"]
	# length(which(disc.tab9$power10>0.8))
	Cancers2<-unique(mr_res_disc$cancer[mr_res_disc$pval<0.05/67])
	# Pos<-regexpr("[:0-9:]",Cancers2)
	# Cancers3<-Cancers2[Pos==-1]
	# Cancers4<-Cancers2[Pos!=-1]
	# Pos<-regexpr("[:0-9:]",Cancers4)
	# Cancers4<-substring(Cancers4,1,Pos-2)
	# Cancers2<-unique(c(Cancers3,Cancers4))
	# Pos<-regexpr("[:0-9:]",mr_res_disc$outcome)
	# mr_res2<-mr_res_disc[Pos!=-1,]
	# mr_res3<-mr_res_disc[Pos==-1,]
	# Pos<-regexpr("[:0-9:]",mr_res2$outcome)
	# mr_res2$outcome<-substring(mr_res2$outcome,1,Pos-2)
	# mr_res_disc<-rbind(mr_res2,mr_res3)
	# mr_res_disc<-mr_res_disc[mr_res_disc$cancer %in% Cancers,]
	# mr_res_disc<-mr_res_disc[!mr_res_disc$cancer %in% Cancers2,]
	
	mr_res_disc<-mr_res_disc[order(mr_res_disc$Cancer.Group),]
	mr_res_disc$b_fix<-abs(mr_res_disc$b)
	mr_res_disc<-mr_res_disc[order(mr_res_disc$b_fix),]	
	mr_res2<-mr_res_disc[mr_res_disc$population == "European",]
	mr_res2<-mr_res2[mr_res2$exposure == "AA:DGLA",]
	
	# mr_fix<-mr_res2[mr_res2$cancer == "Rectal cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Rectal cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)
	
	# mr_fix<-mr_res2[mr_res2$cancer == "Lung adenocarcinoma",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Lung adenocarcinoma",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)
	
	# mr_fix<-mr_res2[mr_res2$cancer == "Kidney cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Kidney cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Basal cell carcinoma",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Basal cell carcinoma",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Leukaemia",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Leukaemia",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)
	
	# mr_fix<-mr_res2[mr_res2$cancer == "Oral cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Oral cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Invasive mucinous ovarian cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Invasive mucinous ovarian cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Neuroblastoma",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Neuroblastoma",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Low malignant potential serous ovarian cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Low malignant potential serous ovarian cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Diffuse large b cell lymphoma",]
	# mr_res2<-mr_res2[mr_res2$cancer !="Diffuse large b cell lymphoma",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_res2<-mr_res2[order(mr_res2$b,decreasing=TRUE),]	

	# mr_fix<-mr_res2[mr_res2$cancer == "Endometrioid ovarian cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Endometrioid ovarian cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)
	
	# mr_fix<-mr_res2[mr_res2$cancer == "Clear cell ovarian cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Clear cell ovarian cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Low grade serous ovarian cancer",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Low grade serous ovarian cancer",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	
	# mr_fix<-mr_res2[mr_res2$cancer == "Central nervous system and eye cancer"  ,]
	# mr_res2<-mr_res2[mr_res2$cancer != "Central nervous system and eye cancer"  ,]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Hodgkin’s lymphoma"  ,]
	# mr_res2<-mr_res2[mr_res2$cancer != "Hodgkin’s lymphoma"   ,]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Multiple myeloma"  ,]
	# mr_res2<-mr_res2[mr_res2$cancer != "Multiple myeloma"   ,]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Lymphoma",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Lymphoma",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Glioma",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Glioma",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)
	
	# mr_fix<-mr_res2[mr_res2$cancer == "Esophageal adenocarcinoma",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Esophageal adenocarcinoma",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)

	# mr_fix<-mr_res2[mr_res2$cancer == "Cancer of digestive organs",]
	# mr_res2<-mr_res2[mr_res2$cancer != "Cancer of digestive organs",]
	# mr_fix<-mr_fix[1,]
	# mr_res2<-rbind(mr_res2,mr_fix)
	# mr_res2<-mr_res2[order(mr_res2$b),]	


	mr_res3<-mr_res_disc[mr_res_disc$population != "European",]
	mr_res3<-mr_res3[mr_res3$exposure == "GLA:LA",]
	mr_res_disc<-rbind(mr_res2,mr_res3)
	mr_res_disc$cases<-as.numeric(mr_res_disc$cases)
	mr_res_disc<-mr_res_disc[order(as.numeric(mr_res_disc$cases),decreasing=T),]

	mr_res_disc<-mr_res_disc[!duplicated(mr_res_disc$cancer),]
	mr_res_disc$outcome<-paste0(mr_res_disc$cancer,"\ncases=",mr_res_disc$cases)
	
	mr_res_disc$Cancer.Group[mr_res_disc$Cancer.Group %in% c("Female reproductive cancers","Ovarian cancer")]<-"Female reproductive cancers"
	mr_res_disc$Cancer.Group[mr_res_disc$Cancer.Group %in% c("Colorectal cancer","Kidney cancer","Skin cancer")]<-"Other cancers"
	mr_res_disc<-mr_res_disc[order(as.numeric(mr_res_disc$cases),decreasing=T),]
	
	return(mr_res_disc)
}



format_dat2<-function(dat=NULL,p_threshold=NULL,sort=NULL){
	id.disc<-disc.tab9$ID
	# id.disc[grep("86",id.disc)]

	# mr_res1$id.outcome[grep("49",mr_res1$id.outcome)]
	
	id.disc<-id.disc[grep(";",id.disc,invert=T)]
	dat1<-dat[dat$exposure !="GLA:LAadj_rs174546",]
	dat2<-dat1[grep(";",dat1$id.outcome),]
	dat3<-dat1[dat1$id.outcome %in% id.disc,]
	dat<-rbind(dat2,dat3)

	dat2<-dat[dat$exposure == "AA:DGLA",]
	dat2<-dat2[grep("East Asian",dat2$population,invert=T),]
	dat3<-dat[dat$exposure == "GLA:LA",]
	dat3<-dat3[dat3$population !="European",]
	dat<-rbind(dat2,dat3)	
	dat$pval<-as.numeric(dat$pval)
	if(!is.null(p_threshold)){
		dat<-dat[dat$pval<p_threshold,]
	}
	dat<-dat[order(as.numeric(dat$cases),decreasing=T),]
	
	if(!is.null(sort)){
		L<-NULL
		Sort<-unique(dat[,sort])
		for(i in 1:length(Sort)){
			L[[i]]<-dat[dat[,sort] == Sort[i],]
		}
		dat<-do.call(rbind,L)	
		if("Other" %in% Sort ){
			dat1<-dat[dat[,sort] == "Other",]
			dat2<-dat[dat[,sort] != "Other",]
			dat<-rbind(dat2,dat1)
		}

	}
	dat$b<-as.numeric(dat$b)
	dat$se<-as.numeric(dat$se)	
	dat$cancer[is.na(dat$cancer)] <- dat$outcome[is.na(dat$cancer)]
	dat$cancer<-paste(dat$cancer,"\ncases=",dat$cases)
	return(dat)
}

replication_plots<-function(cancer=NULL){
	Dat1<-mr_res_disc[mr_res_disc$cancer == cancer,c("cancer","study","cases","b","se","OR","LCI","UCI","pval","study.abbreviation","id.outcome")]
	IDS<-unlist(strsplit(Dat1$id.outcome,split=";"))
	IDS<-trimws(IDS)
	Dat2<-mr_res_rep[which(mr_res_rep$id.outcome %in% IDS),c("cancer","study","cases","b","se","OR","LCI","UCI","pval","study.abbreviation","id.outcome")]
	Dat<-rbind(Dat1,Dat2)
	Dat$study.abbreviation [grep(";",Dat$study.abbreviation )]<-"Overall\neffect"
	# Dat$weight = Dat$cases/100
	forest_plot_1_to_many(mr_res_disc = Dat,b = "b",se = "se",TraitM = "study.abbreviation",
	   col1_width = 1,col1_title = "",exponentiate = TRUE,
	   trans = "log2",ao_slc = FALSE,
	   lo = NULL,up = NULL, by = NULL, xlab = "", addcols = "cases",
	   addcol_widths = 1, addcol_titles = "", subheading_size = 1,   shape_points = 15,
	   colour_scheme = "black",  col_text_size = 5, weight = NULL )	
}

# plot_dat(Dat=Dat,colour=study.abbreviation)

format_dat<-function(Dat=NULL,p_threshold=NULL){

	Dat<-Dat[Dat$cancer != "Cancer (excluding non-melanoma skin cancer)",]
	Dat$system[Dat$system %in% c("Endocrine","Urinary")]<-"Other"
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat<-Dat[order(Dat$system),]
	Dat1<-Dat[Dat$system == "Multiple",]
	Dat2<-Dat[Dat$system != "Multiple",]
	Dat<-rbind(Dat1,Dat2)
	Dat$cancer<-gsub("low malignant potential","LMP",Dat$cancer)
	Dat$cancer<-gsub("Low malignant potential","LMP",Dat$cancer)
	Dat$cancer2<-paste0(Dat$cancer,"\ncases=",Dat$cases)
	if(!is.null(p_threshold)){
		Dat<-Dat[Dat$pval<p_threshold,]
	}
	return(Dat)
}

plot_dat<-function(Dat=NULL,text.names=4,text.title=1,Shape=NULL,colour=NULL){
	p<-forestplot(df = Dat,
			logodds = TRUE,
			name=cancer2,
				  estimate=b,
				  se=se,
				  shape=Shape,
				  colour = colour,
				   xlab = "")+
			# labs(title=Title.plot,size=1)+
			theme(plot.title = element_text(size = text.title))+
			theme(text = element_text(size=text.names))
	return(p)
}

format_site<-function(Site=NULL){
	dat1<-mr_res_disc[mr_res_disc$site == Site,]
	dat2<-mr_res_rep[mr_res_rep$system == Site,]
	Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	return(Dat)
}

format_ovarian<-function(){
	Dat$cancer<-gsub("low malignant potential","LMP",Dat$cancer)
	Dat$cancer<-gsub("Low malignant potential","LMP",Dat$cancer)
	Dat$ovarian_subtype<-"Other"
	Dat$ovarian_subtype[Dat$cancer=="Ovarian cancer"]<-"Ovarian cancer"
	Dat$ovarian_subtype[grep("mucinous",Dat$cancer)]<-"Mucinous"
	Dat$ovarian_subtype[grep("Mucinous",Dat$cancer)]<-"Mucinous"
	Dat$ovarian_subtype[grep("serous",Dat$cancer)]<-"Serous"
	Dat$ovarian_subtype[grep("Serous",Dat$cancer)]<-"Serous"
	Dat$cancer<-gsub(" ovarian cancer","\novarian cancer",Dat$cancer)
	Dat1<-Dat[Dat$ovarian_subtype == "Ovarian cancer",]
	Dat2<-Dat[Dat$ovarian_subtype != "Ovarian cancer",]
	Dat3<-Dat2[Dat2$ovarian_subtype == "Other",]
	Dat4<-Dat2[Dat2$ovarian_subtype != "Other",]
	Dat<-rbind(Dat1,Dat4)
	Dat<-rbind(Dat,Dat3)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat<-Dat[order(Dat$ovarian_subtype),]
	Dat1<-Dat[Dat$ovarian_subtype=="Ovarian cancer",]
	Dat2<-Dat[Dat$ovarian_subtype!="Ovarian cancer",]
	Dat3<-Dat2[Dat2$ovarian_subtype=="Serous",]
	Dat4<-Dat2[Dat2$ovarian_subtype!="Serous",]
	Dat<-rbind(Dat1,Dat3)
	Dat<-rbind(Dat,Dat4)
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	# Dat[,c("cancer","cases","ovarian_subtype")]
	# Dat[order(Dat$ovarian_subtype),c("cancer","ovarian_subtype")]
	# Dat$cancer<-gsub("serous","",Dat$cancer)
	# Dat$cancer<-gsub("Serous","",Dat$cancer)
	# Dat$cancer<-gsub("mucinous","",Dat$cancer)
	# Dat$cancer<-gsub("Mucinous","",Dat$cancer)
	# Dat$cancer<-gsub("  "," ",Dat$cancer)
	# gsub("ovarian cancer","",Dat$cancer)
	# Dat$cancer<-trimws(Dat$cancer)
	# Dat$cancer<-gsub("(^[[:alpha:]])", "\\U\\1", Dat$cancer, perl=TRUE)
	return(Dat)
}

# Dat<-format_site(Site="Ovary")
format_ovarian2<-function(){
	Dat$cancer<-gsub("low malignant potential","LMP",Dat$cancer)
	Dat$cancer<-gsub("Low malignant potential","LMP",Dat$cancer)
	Dat$ovarian_subtype<-"Other"
	Dat$ovarian_subtype[Dat$cancer=="Ovarian cancer"]<-"Ovarian cancer"
	Dat$ovarian_subtype[grep("mucinous",Dat$cancer)]<-"Mucinous"
	Dat$ovarian_subtype[grep("Mucinous",Dat$cancer)]<-"Mucinous"
	Dat$ovarian_subtype[grep("serous",Dat$cancer)]<-"Serous"
	Dat$ovarian_subtype[grep("Serous",Dat$cancer)]<-"Serous"
	Dat$cancer<-gsub(" ovarian cancer","\novarian cancer",Dat$cancer)
	Dat1<-Dat[Dat$ovarian_subtype == "Ovarian cancer",]
	Dat2<-Dat[Dat$ovarian_subtype != "Ovarian cancer",]
	Dat3<-Dat2[Dat2$ovarian_subtype == "Other",]
	Dat4<-Dat2[Dat2$ovarian_subtype != "Other",]
	Dat<-rbind(Dat1,Dat4)
	Dat<-rbind(Dat,Dat3)	
	return(Dat)
}



format_digestive<-function(){
	Dat<-mr_res_disc[mr_res_disc$system == "Digestive",]
	Dat<-Dat[order(Dat$cases,decreasing=T),]	
	Dat$Cancer.Group[Dat$cancer %in% c("Pancreatic cancer","Cancer of digestive organs","Liver cancer") ]<-"Other"
	Dat<-Dat[order(Dat$Cancer.Group),]
	Dat1<-Dat[Dat$Cancer.Group == "Colorectal cancer",]
	Dat2<-Dat[Dat$Cancer.Group != "Colorectal cancer",]
	Dat3<-Dat2[Dat2$Cancer.Group != "Other",]
	Dat4<-Dat2[Dat2$Cancer.Group == "Other",]
	Dat<-rbind(Dat1,Dat3)
	Dat<-rbind(Dat,Dat4)
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	return(Dat)


}

format_skin<-function(){
	dat1<-mr_res_disc[mr_res_disc$system == "Integumentary",]
	dat2<-mr_res_rep[mr_res_rep$system == "Integumentary",]
	Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	# Dat<-Dat[order(Dat$study.abbreviation),]
	return(Dat)
}

format_skin2<-function(){
	dat1<-mr_res_disc[mr_res_disc$system == "Integumentary",]
	dat2<-mr_res_rep[mr_res_rep$system == "Integumentary",]
	Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat<-Dat[order(Dat$cancer,decreasing=T),]
	# Dat<-Dat[order(Dat$study.abbreviation),]
	Dat$cancer2<-paste0(Dat$cancer,"\n(",Dat$cases,")")
	return(Dat)
}

format_skin3<-function(){
	Dat<-mr_res_disc[mr_res_disc$system == "Integumentary",]
	# dat2<-mr_res_rep[mr_res_rep$system == "Integumentary",]
	# Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	# Dat<-Dat[order(Dat$study.abbreviation),]
	Dat$cancer<-gsub("cancer","carcinoma",Dat$cancer)
	Dat$cancer<-gsub(" carcinoma","\ncarcinoma",Dat$cancer)
	Dat$cancer<-gsub("Malignant non-melanoma skin\ncarcinoma","Malignant non-melanoma",Dat$cancer)
	Dat$cancer2<-paste0(Dat$cancer,"\n(",Dat$cases,")")

	return(Dat)
}

format_prostate<-function(){
	Dat<-mr_res_disc[mr_res_disc$site == "Prostate",]
	# dat2<-mr_res_rep[mr_res_rep$system == "Prostate",]
	# Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat$cancer<-gsub("Advanced prostate cancer","Advanced\nprostate cancer",Dat$cancer)
	Dat$cancer<-gsub("Early-onset prostate cancer","Early-onset\nprostate cancer",Dat$cancer)
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	return(Dat)
}


format_allcause<-function(){
	dat1<-mr_res_disc[mr_res_disc$system == "Multiple",]
	dat2<-mr_res_rep[mr_res_rep$system == "Multiple",]
	Dat<-rbind.fill(dat1,dat2)	
	Dat$allcause1<-Dat$study.abbreviation
 	Dat$allcause1[Dat$cancer == "Cancer (excluding non-melanoma skin cancer)"]<-"UKB2 excl\nnon-melanoma\nskin cancer"
	Dat$allcause1[Dat$allcause1 == "UKB"]<-"UKB1 incl\nnon-melanoma\nskin cancer"
	Dat$allcause1[Dat$study.abbreviation == "FinnGen"]<-c("FinnGen2","FinnGen1")
	Dat$allcause1[Dat$study.abbreviation == "UKB; FinnGen"]<-"Overall fixed effect\n(UKB1+FinnGen1)"
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat$cancer2<-paste0(Dat$allcause1,"\n(",Dat$cases,")")
	return(Dat)
}


format_respiratory<-function(){
	Dat<-mr_res_disc[mr_res_disc$system == "Respiratory",]
	# dat2<-mr_res_rep[mr_res_rep$system == "Respiratory",]
	# Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	# Dat<-Dat[!duplicated(Dat$cancer),]		
	Dat<-Dat[order(Dat$Cancer.Group),]	
	Dat$cancer[Dat$cancer == "Respiratory and intrathoracic cancer"]<-"Respiratory &\nintrathoracic cancer"
	Dat$cancer[Dat$cancer == "Oral cavity and pharyngeal cancer"]<-"Oral cavity &\npharyngeal cancer"
	Dat$cancer[Dat$cancer == "Lung cancer in never smokers"]<-"Lung cancer\n(never smokers)"
	Dat$cancer[Dat$cancer == "Lung cancer in ever smokers"]<-"Lung cancer\n(ever smokers)"
	Dat$cancer[Dat$cancer == "Squamous cell lung cancer"]<-"Squamous cell carcinoma"
	Dat$cancer[Dat$cancer == "Small cell lung carcinoma"]<-"Small cell carcinoma"
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	return(Dat)
}


format_respiratory2<-function(){
	Dat<-mr_res_disc[mr_res_disc$system == "Respiratory",]
	Dat2<-mr_res_disc[mr_res_disc$Cancer.Group == "Esophageal cancer",]
	Dat<-rbind(Dat,Dat2)
	# Dat3<-mr_res_rep[which(mr_res_rep$Cancer.Group == "Esophageal cancer"),c("cancer","study.abbreviation")]

	# dat2<-mr_res_rep[mr_res_rep$system == "Respiratory",]
	# Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	# Dat<-Dat[!duplicated(Dat$cancer),]		
	Dat<-Dat[order(Dat$Cancer.Group),]	
	Dat$cancer[Dat$cancer == "Respiratory and intrathoracic cancer"]<-"Respiratory &\nintrathoracic cancer"
	Dat$cancer[Dat$cancer == "Oral cavity and pharyngeal cancer"]<-"Oral cavity &\npharyngeal cancer"
	Dat$cancer[Dat$cancer == "Lung cancer in never smokers"]<-"Lung cancer\n(never smokers)"
	Dat$cancer[Dat$cancer == "Lung cancer in ever smokers"]<-"Lung cancer\n(ever smokers)"
	Dat$cancer[Dat$cancer == "Squamous cell lung cancer"]<-"Squamous cell carcinoma"
	Dat$cancer[Dat$cancer == "Small cell lung carcinoma"]<-"Small cell carcinoma"
	Dat1<-Dat[Dat$Cancer.Group=="Esophageal cancer",]
	Dat2<-Dat[Dat$Cancer.Group!="Esophageal cancer",]	
	Dat<-rbind(Dat2,Dat1)
	Dat$cancer[Dat$cancer == "Esophageal adenocarcinoma"]<-"Esophageal\nadenocarcinoma"
	Dat$cancer[Dat$cancer == "Esophageal squamous cell carcinoma"]<-"Esophageal\nsquamous cell carcinoma"
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	return(Dat)
}

format_blood<-function(){
	dat1<-mr_res_disc[mr_res_disc$system == "Blood",]
	dat2<-mr_res_rep[mr_res_rep$system == "Blood",]
	Dat<-rbind.fill(dat1,dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat<-Dat[!duplicated(Dat$cancer),]		
	Dat$Cancer.Group2<-Dat$Cancer.Group
	Dat$Cancer.Group2[Dat$Cancer.Group2 %in% c("Lymphoma","Non-follicular lymphoma","Hodgkin’s lymphoma", "Multiple myeloma")]<-"Other lymphoma or myeloma"
	Dat$Cancer.Group2[Dat$Cancer.Group2 %in% c("Lymphoid leukaemia","Myeloid leukemia")]<-"Leukemia"
	Dat<-Dat[order(Dat$Cancer.Group2),]	
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	return(Dat)

}

format_blood2<-function(){
	Dat<-mr_res_disc[mr_res_disc$system == "Blood",]
	Dat<-Dat[order(Dat$cases,decreasing=T),]	
	Dat$Cancer.Group2<-Dat$Cancer.Group
	Dat$Cancer.Group2[Dat$Cancer.Group2 %in% c("Lymphoma","Non-follicular lymphoma","Hodgkin’s lymphoma", "Multiple myeloma")]<-"Other lymphoma or myeloma"
	Dat$Cancer.Group2[Dat$Cancer.Group2 %in% c("Lymphoid leukaemia","Myeloid leukemia")]<-"Leukemia"
	Dat<-Dat[order(Dat$Cancer.Group2),]	
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	return(Dat)

}

format_breast<-function(){
	Dat<-mr_res_disc[mr_res_disc$Cancer.Group == "Breast cancer",]	
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat$cancer2<-paste0(Dat$cancer,"\n(cases=",Dat$cases,")")
	return(Dat)
}



format_colorectal<-function(){
	Dat1<-mr_res_disc[which(mr_res_disc$Cancer.Group =="Colorectal cancer"),]
	Dat2<-mr_res_rep[which(mr_res_rep$Cancer.Group =="Colorectal cancer"),]
	Dat<-rbind.fill(Dat1,Dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	Dat<-Dat[!duplicated(Dat$cancer),]
	Dat$weight<-1/Dat$se/5
	Dat$crc_group<-"By site"
	Dat$crc_group[Dat$cancer =="Colorectal cancer"]<-""
	Dat$crc_group[grep("male",Dat$cancer)]<-"By sex"

	# Dat$cancer[Dat$cancer=="Colorectal cancer in males"  ]<-"Colorectal cancer\nin males" 
	# Dat$cancer[Dat$cancer=="Colorectal cancer in females"  ]<-"Colorectal cancer\nin females" 
	return(Dat)
}

format_colorectal2<-function(){
	Dat1<-mr_res_disc[which(mr_res_disc$Cancer.Group =="Colorectal cancer"),]
	Dat2<-mr_res_rep[which(mr_res_rep$Cancer.Group =="Colorectal cancer"),]
	Dat<-rbind.fill(Dat1,Dat2)
	Dat<-Dat[order(Dat$cases,decreasing=T),]
	# Dat<-Dat[!duplicated(Dat$cancer),]
	Dat$weight<-1/Dat$se/5
	Dat$crc_group<-"By location"
	Dat$crc_group[Dat$cancer =="Colorectal cancer"]<-""
	Dat$cancer[Dat$cancer =="Colorectal cancer"]<-Dat$study.abbreviation[Dat$cancer =="Colorectal cancer"]
	Dat$cancer[Dat$cancer == "GECCO/CORECT/CCFR; ACCC; FinnGen" ]<-"Overall fixed effect (all studies)"	
	Dat$crc_group[grep("male",Dat$cancer)]<-"By sex"

	# Dat$cancer[Dat$cancer=="Colorectal cancer in males"  ]<-"Colorectal cancer\nin males" 
	# Dat$cancer[Dat$cancer=="Colorectal cancer in females"  ]<-"Colorectal cancer\nin females" 
	return(Dat)
}






# meta_analysis<-function(dat1=NULL){	
# 		b<-dat1$lnor
# 		se<-dat1$se
# 		# p<-temp$p
# 		w<-1/se^2
# 		b.fixed<-sum(b*w)/(sum(w))
# 		se.fixed<-sqrt(sum(w)^-1)
# 		z<-abs(b.fixed/se.fixed)
# 		p.fixed<-pnorm(z,lower.tail=F)*2
# 		nstudies.fixed<-length(b)
# 		cancer<-unique(dat1$cancer)
# 		if(length(cancer)!=1) cancer<-unique(dat1$Cancer.Group)
# 		if(length(cancer)!=1) 	cancer<-unique(paste(dat1$system,"system cancers"))
# 		if(length(cancer)!=1) stop("length of cancer not 1")
# 		ids.fixed<-IDS[i]
# 		cases<-sum(dat1$cases)
# 		controls<-sum(dat1$controls)
# 		study<-"Overall fixed effect"
# 		Q<-sum((b.fixed-b)^2*w)
# 		df.Q<-length(b)-1		
# 		Q.p<-pchisq(Q, df.Q, lower.tail = FALSE)}
# }


format_all_discovery_v2<-function(){
	load("~/fatty-acids/mr/results/mr_results_rep_v3.Rdata")
	# mr_res1$id.outcome[grep("86",mr_res1$id.outcome)]
	mr_res1$pval<-as.numeric(mr_res1$pval)
	mr_res1$b<-as.numeric(mr_res1$b)
	mr_res1$se<-as.numeric(mr_res1$se)
	mr_res1$OR<-round(exp(mr_res1$b),2)
	mr_res1$LCI<-round(exp(mr_res1$b-mr_res1$se*1.96),2)
	mr_res1$UCI<-round(exp(mr_res1$b+mr_res1$se*1.96),2)
	mr_res1<-mr_res1[!is.na(mr_res1$cases),]
	mr_res1$cases<-as.numeric(mr_res1$cases)
	mr_res1<-mr_res1[mr_res1$cases>1000,]
	mr_res1<-mr_res1[order(mr_res1$cases,decreasing=TRUE),]
	mr_res1$cancer[is.na(mr_res1$cancer)]<-mr_res1$outcome[is.na(mr_res1$cancer)]
	mr_res1<-mr_res1[!duplicated(mr_res1$cancer),]
	mr_res1<-mr_res1[grep("males",mr_res1$cancer,invert=TRUE),]
	mr_res1<-mr_res1[grep("smokers",mr_res1$cancer,invert=TRUE),]
	return(mr_res1)
}

format_plot_dat<-function(){
	plot_dat$cancer[plot_dat$cancer=="Cancer (all cause)"]<-"Cancer (all sites)"
	plot_dat$cancer[plot_dat$cancer=="Cancer (excluding non-melanoma skin cancer)"]<-"Cancer (excl nmsc)"
	plot_dat$cancer[plot_dat$cancer=="Low grade & low malignant potential serous ovarian cancer"]<-"LG & LMP serous ovarian cancer"
	plot_dat$cancer[plot_dat$cancer=="Low malignant potential serous ovarian cancer"]<-"LMP serous ovarian cancer"
	plot_dat$cancer[plot_dat$cancer=="Low malignant potential mucinous ovarian cancer"]<-"LMP mucinous ovarian cancer"
	plot_dat$cancer[plot_dat$cancer=="Low malignant 
	potential ovarian cancer"]<-"LMP ovarian cancer"
	plot_dat$cancer[plot_dat$cancer=="Central nervous system and eye cancer"]<-"CNS & eye cancer"
	plot_dat$cancer[plot_dat$cancer=="Malignant non-melanoma skin cancer"  ]<-"Non-melanoma"
	plot_dat$cancer[plot_dat$cancer=="Malignant skin cancer"  ]<-"Overall skin cancer"
	plot_dat$cancer<-gsub("Low malignant potential","LMP",plot_dat$cancer)
	plot_dat$cancer<-gsub("Low grade","LG",plot_dat$cancer)
	plot_dat$cancer<-gsub("Noncardia","NC",plot_dat$cancer)


	

	plot_dat$cancer[plot_dat$cancer=="Respiratory and intrathoracic cancer"  ]<-"Respiratory & intrathoracic cancer"  
	plot_dat$cancer[plot_dat$cancer=="Esophageal squamous cell carcinoma"  ]<-"Esophageal SCC" 
	plot_dat$cancer<-trimws(gsub("cancer","",plot_dat$cancer))
	
	plot_dat$cancer[plot_dat$cancer=="Hodgkin’s lymphoma"]<-"Hodgkin lymphoma"
	plot_dat$cancer[plot_dat$cancer=="Non-hodgkin's lymphoma"]<-"Non-Hodgkin lymphoma"
	plot_dat$cancer[plot_dat$cancer=="Diffuse large b cell lymphoma"]<-"Diffuse large B cell lymphoma"

	plot_dat$outcome<-paste0(plot_dat$cancer,"\ncases=",plot_dat$cases)

	# plot_dat$cancer[nchar(plot_dat$cancer)>30]
	# plot_dat$cases[plot_dat$cancer=="Lung cancer"]
	plot_dat<-plot_dat[order(plot_dat$cases,decreasing=TRUE),]
	plot_dat$system[plot_dat$system=="Reproductive"]<-"AReproductive"
	plot_dat$system[plot_dat$system=="Digestive"]<-"BDigestive"
	plot_dat$system[plot_dat$system=="Blood"]<-"CBlood"
	plot_dat$system[plot_dat$system=="Respiratory"]<-"DRespiratory"
	plot_dat$system[plot_dat$system=="Integumentary"]<-"EIntegumentary"
	# plot_dat$system[plot_dat$system=="Integumentary"]<-"EIntegumentary"
	# plot_dat$system[plot_dat$system=="Nervous"]<-"FNervous"
	plot_dat$system[plot_dat$system %in% c("Endocrine","Multiple","Urinary","Nervous")]<-"GOther"	
	# plot_dat<-plot_dat[!plot_dat$cancer %in% c("Leukaemia","Oral","Invasive mucinous ovarian","Neuroblastoma","LMP serous ovarian","Diffuse large b cell lymphoma"),]

	return(plot_dat)
}

plot_meta_reg_smoking<-function(){
	mod1<-c("Other cancers",Model$b[1],Model$se[1])
	mod2<-c("Smoking related cancers",Model$b[2],Model$se[2])
	mod<-data.frame(matrix(c(mod2,mod1),nrow=2,ncol=3,byrow=TRUE),stringsAsFactors=F)
	names(mod)<-c("plot_name","b","se")
	mod$smoking[1] <- "Smoking related cancers (overall effect)"
	mod$smoking[2] <- "Other cancers (overall effect)"
	mod$b<-as.numeric(mod$b)
	mod$se<-as.numeric(mod$se)
	mod$Shape<-"diamond"
	mr_res2<-mr_res[mr_res$smoking1==1,]
	mr_res3<-mr_res[mr_res$smoking1==0,]
	mr_res2<-rbind.fill(mod[1,],mr_res2)
	mr_res3<-rbind.fill(mod[2,],mr_res3)
	mr_res<-rbind(mr_res2,mr_res3)
	mr_res$Colour[is.na(mr_res$Colour)]<-0
	mr_res$Colour<-mr_res$Colour+1
	return(mr_res)
}


plot_meta_reg_infl<-function(){
	mod1<-c("Other cancers",Model$b[1],Model$se[1],0)
	mod2<-c("Chronic inflammation cancers",Model$b[2],Model$se[2],1)
	mod<-data.frame(matrix(c(mod2,mod1),nrow=2,ncol=4,byrow=TRUE),stringsAsFactors=F)
	names(mod)<-c("plot_name","b","se","infl")
	# mod$smoking[1] <- "Smoking related cancers (overall effect)"
	# mod$smoking[2] <- "Other cancers (overall effect)"
	mod$b<-as.numeric(mod$b)
	mod$se<-as.numeric(mod$se)
	mod$Shape<-"diamond"
	mr_res<-mr_res[order(mr_res$cases,decreasing=TRUE),]
	mr_res2<-mr_res[mr_res$infl==1,]
	mr_res3<-mr_res[mr_res$infl==0,]
	mr_res2<-rbind.fill(mod[1,],mr_res2)
	mr_res3<-rbind.fill(mod[2,],mr_res3)
	mr_res<-rbind(mr_res2,mr_res3)
	mr_res$infl<-as.numeric(mr_res$infl)
	mr_res$Colour<-mr_res$infl+1

	# mr_res$Colour[is.na(mr_res$Colour)]<-0
	# mr_res$Colour<-mr_res$Colour+1
	return(mr_res)
}


plot_meta_reg_infl_agent<-function(){
	mod1<-c("Other cancers",Model$b[1],Model$se[1],0)
	mod2<-c("Chronic inflammation plus infectious agents",Model$b[2],Model$se[2],1)
	mod<-data.frame(matrix(c(mod2,mod1),nrow=2,ncol=4,byrow=TRUE),stringsAsFactors=F)
	names(mod)<-c("plot_name","b","se","infl_agent")
	# mod$smoking[1] <- "Smoking related cancers (overall effect)"
	# mod$smoking[2] <- "Other cancers (overall effect)"
	mod$b<-as.numeric(mod$b)
	mod$se<-as.numeric(mod$se)
	mod$Shape<-"diamond"
	mr_res<-mr_res[order(mr_res$cases,decreasing=TRUE),]
	mr_res2<-mr_res[mr_res$infl_agent==1,]
	mr_res3<-mr_res[mr_res$infl_agent==0,]
	mr_res2<-rbind.fill(mod[1,],mr_res2)
	mr_res3<-rbind.fill(mod[2,],mr_res3)
	mr_res<-rbind(mr_res2,mr_res3)
	mr_res$infl_agent<-as.numeric(mr_res$infl_agent)
	mr_res$Colour<-mr_res$infl_agent+1

	# mr_res$Colour[is.na(mr_res$Colour)]<-0
	# mr_res$Colour<-mr_res$Colour+1
	return(mr_res)
}


plot_meta_reg_digestive<-function(){
	mod1<-c("Other cancers",Model$b[1],Model$se[1],0)
	mod2<-c("Digestive system cancers",Model$b[2],Model$se[2],1)
	mod<-data.frame(matrix(c(mod2,mod1),nrow=2,ncol=4,byrow=TRUE),stringsAsFactors=F)
	names(mod)<-c("plot_name","b","se","system2_num")
	# mod$smoking[1] <- "Smoking related cancers (overall effect)"
	# mod$smoking[2] <- "Other cancers (overall effect)"
	mod$b<-as.numeric(mod$b)
	mod$se<-as.numeric(mod$se)
	mod$Shape<-"diamond"
	mr_res<-mr_res[order(mr_res$cases,decreasing=TRUE),]
	mr_res2<-mr_res[mr_res$system2_num==1,]
	mr_res3<-mr_res[mr_res$system2_num==0,]
	mr_res2<-rbind.fill(mod[1,],mr_res2)
	mr_res3<-rbind.fill(mod[2,],mr_res3)
	mr_res<-rbind(mr_res2,mr_res3)
	mr_res$system2_num<-as.numeric(mr_res$system2_num)
	mr_res$Colour<-mr_res$system2_num+1

	# mr_res$Colour[is.na(mr_res$Colour)]<-0
	# mr_res$Colour<-mr_res$Colour+1
	return(mr_res)
}


format_plot<-function(dat=NULL,cancer=NULL){
	if(!is.null(cancer))
	{
		dat<-dat[dat$outcome2 %in%cancer,]
	}
	dat$PUFA<-NA
	dat$PUFA[grep("n3",dat$exposure)]<-"Omega 3"
	dat$PUFA[dat$exposure=="Omega-3 fatty acids"]<-"Omega 3"
	dat$PUFA[grep("n6",dat$exposure)]<-"Omega 6"
	dat$PUFA[dat$exposure=="Omega-6 fatty acids"]<-"Omega 6"
	dat$PUFA[grep("n3 or n6",dat$exposure)]<-"other"
	dat$PUFA[is.na(dat$PUFA)] <- "other"

	dat1<-dat[grep(":",dat$exposure),]
	dat1<-dat1[dat1$exposure!= "Other polyunsaturated fatty acids than 18:2",]
	dat2<-dat[grep(":",dat$exposure,invert=TRUE),]
	dat2<-rbind(dat2,dat[dat$exposure== "Other polyunsaturated fatty acids than 18:2",])
	Pos<-gregexpr(":",dat1$exposure)
	dat1$chain.length<-as.numeric(substr(dat1$exposure,start=unlist(Pos)-2,stop=unlist(Pos)-1))
	dat2$chain.length<-NA
	dat<-rbind(dat1,dat2)

	# dat<-dat[dat$population =="European",]
	# dat<-dat[dat$exposure=="Arachidonic acid",]
	dat$b<-as.numeric(dat$b)
	dat$se<-as.numeric(dat$se)
	dat$ncase.outcome<-as.numeric(dat$ncase.outcome)
	dat<-dat[order(dat$ncase.outcome,decreasing=TRUE),]
	# dat$ncase.outcome<-as.numeric(dat$ncase.outcome)
	dat$shape<-15
	dat$weight<-1/dat$se/5
	# dat$name <- paste(dat$outcome2,"\n",dat$ncase.outcome)

	dat$outcome_name<-paste(dat$outcome2,"\nN. cases=",dat$ncase.outcome)
	# dat<-dat[dat$outcome2 !="Overall cancer (excluding non-melanoma skin cancer)",]
	dat$outcome_name<-gsub("\\(excluding non-melanoma skin cancer)","\n(excl. nm-skin cancer)",dat$outcome_name)
	dat$outcome_name<-gsub("Respiratory and","Respiratory &\n",dat$outcome_name)
	dat<-dat[!dat$exposure %in% c("Ratio of omega-6 fatty acids to omega-3 fatty acids","Other polyunsaturated fatty acids than 18:2" ,"Dihomo-linolenic acid (20:3n3 or n6)", "Tetradecadienoic acid (14:2n9)"),]
	exp_dat<-find_exposure_info()
	exp_dat<-exp_dat[,c("id.exposure","consortium","exposure")]
	dat<-merge(dat,exp_dat,by="exposure")

	load("~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity_v3.Rdata")
	het1<-het[,c("outcome","id.outcome","id.exposure","method","Q","Q_pval","Q_df")]
	het1$fads<-FALSE
	load("~/fatty-acids/mr/results/secondary_pufas_inclfads_heterogeneity_v3.Rdata")
	
	het2<-het[,c("outcome","id.outcome","id.exposure","method","Q","Q_pval","Q_df")]
	het2$fads<-TRUE
	het<-rbind(het1,het2)
	het<-het[het$method !="MR Egger",]
	het<-merge(het,exp_dat,by="id.exposure")
	het[het$exposure=="Eicosapentaenoic acid (20:5n3)" , ]

	if(cancer == "Basal cell carcinoma")
	{
		het_bcc_ala<-het[het$exposure == "Alpha-linolenic acid (18:3n3)" & het$id.outcome == 1,]
	}
	# unique(het$id.outcome)
	# unique(dat$id.outcome)
	# het[het$id.outcome %in% c("75","149"),]
	if(cancer %in% c("Lung cancer","Basal cell carcinoma"))
	{	
		load("~/fatty-acids/mr/results/secondary_pufas_inclfads_heterogeneity_fisherp_v3.Rdata")
		het1<-het_fisherp
		het1$fads<-TRUE
		load("~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity_fisher_p_v3.Rdata")
		het2<-het_fisherp
		het2$fads<-FALSE
		head(het1)
		het<-rbind(het1,het2)
		
		het<-het[het$outcome == cancer & het$method == "Inverse variance weighted",]
		# dat2<-dat
		# dat<-dat2
		# dat2[dat2$exposure  == "Alpha-linolenic acid (18:3n3)",]
		
		het[unique(het$exposure)  == "Alpha-linolenic acid (18:3n3)",]

		dat<-merge(dat,het,by.x=c("outcome2","exposure","fads"),by.y=c("outcome","exposure","fads"),all.x=TRUE)

		if(cancer == "Basal cell carcinoma")
		{
			dat$Q_pval_fisher[dat$exposure == "Alpha-linolenic acid (18:3n3)" & dat$id.outcome == "1" & dat$nsnp == 2]<-het_bcc_ala$Q_pval 
		}
	
		names(dat)[names(dat) == "Q_pval"]<-"Q_pval_studies"
		names(dat)[names(dat) == "Q_pval_fisher"]<-"Q_pval"
		names(dat)[names(dat) =="nsnp"]<-"nsnp_old"

		names(dat)[names(dat) =="nsnp_avg"]<-"nsnp"
		
		dat$nsnp[which(dat$nsnp==20.5)]<-21
		dat$nsnp<-round(as.numeric(dat$nsnp),0)
		dat$nsnp[is.na(dat$nsnp)]<-dat$nsnp_old[is.na(dat$nsnp)]
		dat[,c("exposure","nsnp","fads","nsnp_old")]

	}else{

		dat<-merge(dat,het,by=c("exposure","id.outcome","fads"),all.x=TRUE)
	}

	Dups<-unique(dat$exposure[duplicated(dat$exposure)])
	dat<-dat[dat$exposure %in% Dups,]



	# dat<-dat[!dat$exposure %in% c("Omega-3 fatty acids","Omega-6 fatty acids"),]
	dat<-dat[dat$PUFA !="other",]
	dat$exposure2<-paste0(dat$exposure,"\n(",dat$consortium,")")
	dat$exposure3<-paste0(dat$exposure,"\n(",dat$consortium," nsnps=",dat$nsnp,")")
	
	
	dat$pval<-formatC(as.numeric(dat$pval), format = "e", digits = 2)
	# dat$pval<-paste(dat$pval,c("A","B","C"))
	dat$Q_pval<-formatC(as.numeric(dat$Q_pval), format = "e", digits = 2)
	dat$Q_pval<-gsub(" ","",dat$Q_pval)
	dat$pval_name<-paste0(dat$pval,"\n")
	dat$Q_pval_name<-paste0(dat$Q_pval)
	dat$Q_pval_name[dat$Q_pval_name=="NA"]<-c("NA1","NA2","NA3")
	dat$nsnps_name<-paste0(substring(dat$exposure,1,1),substring(dat$exposure,10,10),"   ",dat$nsnp)
	
	dat$FADS<-dat$fads
	# dat$FADS<-""
	dat$FADS[dat$fads]<-"included"
	dat$FADS[!dat$fads]<-"excluded"
	dat$FADS2[dat$FADS=="excluded"]<-"B"
	dat$FADS2[dat$FADS=="included"]<-"A"

	dat1<-dat[dat$exposure %in% c("Omega-6 fatty acids","Omega-3 fatty acids"),]
	dat2<-dat[!dat$exposure %in% c("Omega-6 fatty acids","Omega-3 fatty acids"),]
	exp<-unlist(strsplit(dat2$exposure,split="\\("))
	exp<-trimws(exp[seq(1,length(exp),by=2)])
	dat2$exposure_clean<-exp
	dat1$exposure_clean<-dat1$exposure
	dat<-rbind(dat1,dat2)
	names(dat)[names(dat) == "consortium.x"]<-"consortium"
	dat$exposure_study<-paste0(dat$exposure_clean," (",dat$consortium,")")
	
	# dat<-dat[order(dat$exposure),]
	dat1<-dat[dat$PUFA == "Omega 6",]
	dat<-dat[dat$PUFA != "Omega 6",]
	dat2<-dat1[dat1$exposure=="Linoleic acid (18:2n6)",]
	dat1<-dat1[dat1$exposure!="Linoleic acid (18:2n6)",]		
	dat<-dat[order(dat$PUFA,dat$chain.length,dat$exposure,dat$FADS2),]
	dat1<-dat1[order(dat1$PUFA,dat1$chain.length,dat1$exposure,dat1$FADS2),]
	dat2<-dat2[order(dat2$PUFA,dat2$chain.length,dat2$exposure,dat2$FADS2),]
	dat1<-rbind(dat2,dat1)
	dat<-rbind(dat,dat1)
	# dat$exposure[dat$exposure=="Alpha-linolenic acid (18:3n3)"   ][1]<-paste0("Omega 3 PUFAs\n\n","Alpha-linolenic acid (18:3n3)"  )
	# dat$exposure[dat$exposure=="Linoleic acid (18:2n6)"   ][1]<-paste0("Omega 6 PUFAs\n\n","Linoleic acid (18:2n6)"   )
	dat$or<-round(exp(dat$b),2)
	dat$lci<-round(exp(dat$b-1.96*dat$se),2)
	dat$uci<-round(exp(dat$b+1.96*dat$se),2)

	return(dat)
}



find_exposure_info<-function()
{
	load("~/fatty-acids/mr/data/instruments.Rdata")
	exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")
	exp_dat<-unique(exp_dat[,c("exposure","consortium","id.exposure")])
	return(exp_dat)
}

format_plot2<-function(dat=NULL,cancer=NULL){
	if(!is.null(cancer))
	{
		dat<-dat[dat$outcome %in%cancer,]
	}
	
	dat$PUFA<-NA
	dat$PUFA[grep("n3",dat$exposure)]<-"Omega 3"
	dat$PUFA[dat$exposure=="Omega-3 fatty acids"]<-"Omega 3"
	dat$PUFA[grep("n6",dat$exposure)]<-"Omega 6"
	dat$PUFA[dat$exposure=="Omega-6 fatty acids"]<-"Omega 6"
	dat$PUFA[grep("n3 or n6",dat$exposure)]<-"other"
	dat$PUFA[is.na(dat$PUFA)] <- "other"

	dat1<-dat[grep(":",dat$exposure),]
	dat1<-dat1[dat1$exposure!= "Other polyunsaturated fatty acids than 18:2",]
	dat2<-dat[grep(":",dat$exposure,invert=TRUE),]
	dat2<-rbind(dat2,dat[dat$exposure== "Other polyunsaturated fatty acids than 18:2",])
	Pos<-gregexpr(":",dat1$exposure)
	dat1$chain.length<-as.numeric(substr(dat1$exposure,start=unlist(Pos)-2,stop=unlist(Pos)-1))
	dat2$chain.length<-NA
	dat<-rbind(dat1,dat2)

	# dat<-dat[dat$population =="European",]
	# dat<-dat[dat$exposure=="Arachidonic acid",]
	dat$b<-as.numeric(dat$b)
	dat$se<-as.numeric(dat$se)
	#dat$ncase.outcome<-as.numeric(dat$ncase.outcome)
	#dat<-dat[order(dat$ncase.outcome,decreasing=TRUE),]
	# dat$ncase.outcome<-as.numeric(dat$ncase.outcome)
	dat$shape<-15
	dat$weight<-1/dat$se/5
	# dat$name <- paste(dat$outcome2,"\n",dat$ncase.outcome)

	#dat$outcome_name<-paste(dat$outcome,"\nN. cases=",dat$ncase.outcome)
	# dat<-dat[dat$outcome2 !="Overall cancer (excluding non-melanoma skin cancer)",]
	#dat$outcome_name<-gsub("\\(excluding non-melanoma skin cancer)","\n(excl. nm-skin cancer)",dat$outcome_name)
	#dat$outcome_name<-gsub("Respiratory and","Respiratory &\n",dat$outcome_name)
	dat<-dat[!dat$exposure %in% c("Ratio of omega-6 fatty acids to omega-3 fatty acids","Other polyunsaturated fatty acids than 18:2" ,"Dihomo-linolenic acid (20:3n3 or n6)", "Tetradecadienoic acid (14:2n9)"),]
	exp_dat<-find_exposure_info()
	exp_dat<-exp_dat[,c("id.exposure","consortium","exposure")]
	dat<-merge(dat,exp_dat,by="exposure")

	load("~/fatty-acids/mr/results/secondary_pufas_exclfads_heterogeneity_altmethod.Rdata")
	het1<-het[,c("outcome","id.outcome","id.exposure","method","Q","Q_pval","Q_df")]
	het1$fads<-FALSE
	load("~/fatty-acids/mr/results/secondary_pufas_inclfads_heterogeneity_altmethod.Rdata")
	
	het2<-het[,c("outcome","id.outcome","id.exposure","method","Q","Q_pval","Q_df")]
	het2$fads<-TRUE
	het<-rbind(het1,het2)
	het<-het[het$method !="MR Egger",]
	het<-merge(het,exp_dat,by="id.exposure")
	het[het$exposure=="Eicosapentaenoic acid (20:5n3)" , ]
	dat<-merge(dat,het,by=c("exposure","id.outcome","fads"),all.x=TRUE)
	Dups<-unique(dat$exposure[duplicated(dat$exposure)])
	dat<-dat[dat$exposure %in% Dups,]
	# dat<-dat[!dat$exposure %in% c("Omega-3 fatty acids","Omega-6 fatty acids"),]
	dat<-dat[dat$PUFA !="other",]
	dat$exposure2<-paste0(dat$exposure,"\n(",dat$consortium,")")
	dat$exposure3<-paste0(dat$exposure,"\n(",dat$consortium," nsnps=",dat$nsnp,")")
	
	
	dat$pval<-formatC(as.numeric(dat$pval), format = "e", digits = 2)
	# dat$pval<-paste(dat$pval,c("A","B","C"))
	dat$Q_pval<-formatC(as.numeric(dat$Q_pval), format = "e", digits = 2)
	dat$Q_pval<-gsub(" ","",dat$Q_pval)
	dat$pval_name<-paste0(dat$pval,"\n")
	dat$Q_pval_name<-paste0(dat$Q_pval)
	dat$Q_pval_name[dat$Q_pval_name=="NA"]<-c("NA1","NA2","NA3")
	dat$nsnps_name<-paste0(substring(dat$exposure,1,1),substring(dat$exposure,10,10),"   ",dat$nsnp)
	
	dat$FADS<-dat$fads
	# dat$FADS<-""
	dat$FADS[dat$fads]<-"included"
	dat$FADS[!dat$fads]<-"excluded"
	dat$FADS2[dat$FADS=="excluded"]<-"B"
	dat$FADS2[dat$FADS=="included"]<-"A"

	dat1<-dat[dat$exposure %in% c("Omega-6 fatty acids","Omega-3 fatty acids"),]
	dat2<-dat[!dat$exposure %in% c("Omega-6 fatty acids","Omega-3 fatty acids"),]
	exp<-unlist(strsplit(dat2$exposure,split="\\("))
	exp<-trimws(exp[seq(1,length(exp),by=2)])
	dat2$exposure_clean<-exp
	dat1$exposure_clean<-dat1$exposure
	dat<-rbind(dat1,dat2)
	names(dat)[names(dat) == "consortium.x"]<-"consortium"
	dat$exposure_study<-paste0(dat$exposure_clean," (",dat$consortium,")")
	
	# dat<-dat[order(dat$exposure),]
	dat1<-dat[dat$PUFA == "Omega 6",]
	dat<-dat[dat$PUFA != "Omega 6",]
	dat2<-dat1[dat1$exposure=="Linoleic acid (18:2n6)",]
	dat1<-dat1[dat1$exposure!="Linoleic acid (18:2n6)",]		
	dat<-dat[order(dat$PUFA,dat$chain.length,dat$exposure,dat$FADS2),]
	dat1<-dat1[order(dat1$PUFA,dat1$chain.length,dat1$exposure,dat1$FADS2),]
	dat2<-dat2[order(dat2$PUFA,dat2$chain.length,dat2$exposure,dat2$FADS2),]
	dat1<-rbind(dat2,dat1)
	dat<-rbind(dat,dat1)
	# dat$exposure[dat$exposure=="Alpha-linolenic acid (18:3n3)"   ][1]<-paste0("Omega 3 PUFAs\n\n","Alpha-linolenic acid (18:3n3)"  )
	# dat$exposure[dat$exposure=="Linoleic acid (18:2n6)"   ][1]<-paste0("Omega 6 PUFAs\n\n","Linoleic acid (18:2n6)"   )
	dat$or<-round(exp(dat$b),2)
	dat$lci<-round(exp(dat$b-1.96*dat$se),2)
	dat$uci<-round(exp(dat$b+1.96*dat$se),2)

	return(dat)
}