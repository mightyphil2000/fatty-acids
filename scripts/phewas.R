
setwd("~/Google Drive/MR_base/5am")

SNP<-"rs174546"
ao<-available_outcomes()
ao1<-ao[!is.na(ao$ncase),]
ao2<-ao[is.na(ao$ncase),]
ao1<-ao[ao1$ncase>500,]
ao2<-ao2[!is.na(ao2$sample_size), ]
ao2<-ao2[ao2$sample_size>500,]
ao3<-rbind(ao1,ao2)

IDS<-ao3$id

out_dat <- extract_outcome_data(
    snps = SNP,
    outcomes = IDS,
    proxies = TRUE, 
    rsq = 0.8,
)

save(out_dat,file="~/fatty-acids/data/phewas_rs174546.Rdata")

load("~/fatty-acids/data/phewas_rs174546.Rdata")
res<-read.table("IL23R-findings.txt",sep="\t",head=T,stringsAsFactors=F,fill=T,quote="")

data<-out_dat

phewas<-function(data=NULL){
	
	library(TwoSampleMR)
	ao<-available_outcomes()

	data$trait<-data$originalname.outcome
	for(i in 1:nrow(data)){
		print(i)
		trait<-trimws(unlist(strsplit(data$trait[i],split=":")))
		# dat$trait[grep(":",dat$trait)]
		trait<-trait[length(trait)]
		Test<-unlist(strsplit(trait,split=" "))
		if(nchar(Test[1])==3 & !Test[1] %in% c("ER+","ER-","Low")) {
			trait<-paste(Test[2:length(Test)],collapse=" ")
		}
		data$trait[i]<-trait
	}

	data<-data[data$trait!="Internalizing problems",]

	data$trait[data$trait %in% c("Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)",
		"ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)",
		"ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)")]<-c("Breast cancer","ER+ breast cancer","ER- breast cancer")

	data$trait<-gsub("(^[[:alpha:]])", "\\U\\1", data$trait, perl=TRUE)
	data$trait<-gsub("  ", " ", data$trait)
	data$trait<-gsub("and", "&", data$trait)
	data$trait<-gsub("\\(dvt)", "", data$trait)
	data$trait<-trimws(data$trait)
	names(data)[names(data) == "subcategory.outcome"]<-"Category"
	# data$Category[data$Category=="Autoimmune"]<-"Autoimmune/inflammatory"

	data$Category[data$trait=="Coxarthrosis [arthrosis of hip]"]<-"Bone"
	data$Category[data$trait %in% c("Muscle or soft tissue injuries","Uterine fibroids","Leiomyoma of uterus","Fibroblastic disorders")]<-"Bone/connective tissue"
	data$Category[data$Category=="Bone"]<-"Bone/connective tissue"
	data$Category[data$Category %in% c("Cardiovascular","Diabetes")]<-"Cardiometabolic"
	data$Category[data$Category %in% c("Neurological","Psychiatric")]<-"Psychiatric/neurological"
	data$Category[data$trait %in% c("Sleep disorders")]<-"Psychiatric/neurological"
	data$trait[data$trait=="Benign neoplasm of colon rectum anus & anal canal"]<-"Colorectal benign neoplasm"
	data$Category[data$trait %in% c("Diverticular disease/diverticulitis","Gastritis & duodenitis","Acute appendicitis")]<-"Autoimmune/inflammatory"   
	data$trait[data$trait=="Gastro-oesophageal reflux (gord) / gastric reflux" ]<-"Gastro-oesophageal reflux"
	data$trait[data$trait=="Kidney stone/ureter stone/bladder stone"] <-"Kidney stone"
	data$trait[data$trait=="Chronic obstructive airways disease/copd"] <-"Chronic obstructive airways disease"
	data$trait[data$trait=="Haemorrhage from respiratory passages"] <-"Respiratory haemorrhage"
	data<-data[!data$trait %in% c("Follow-up examination after treatment for conditions other than malignant neoplasms","Malignant neoplasm of breast"),]
	 
	data$trait[data$trait=="Injury or trauma resulting in loss of vision"]<-"Loss of vision"
	data$trait[data$trait=="Fissure & fistula of anal & rectal regions"]<-"Anal/rectal fissure"
	data$trait[data$trait=="Excessive frequent & irregular menstruation"]<-"Excessive/irregular menstruation"
	data$Category[data$Category=="Cancer"]<-"Neoplastic"

	# ids<-data$id[!data$id %in% res$id.outcome]
	# data$Trait_name[!data$id %in% ids]
	data$Category[data$Category %in% c("Eye","Kidney","Kindey","Lung")]<-"Other"
	data<-data[order(data$beta.outcome,decreasing=T),]
	data<-data[order(data$Category),]
	# data$beta<-data$beta*-1 

	data<-data[!data$trait %in% c("Retinal detachment","High cholesterol"),]
	
	res.other<-data[data$Category=="Other",]
	res.bone<-data[data$Category=="Bone/connective tissue",]
	data2<-data[!data$Category %in% c("Other","Bone/connective tissue"),]
	data2<-rbind(data2,res.bone)
	data2<-rbind(data2,res.other)
	data2[,c("trait","Category")]

	# Controls<-data[order(as.numeric(data$ncontrol.outcome),decreasing=T),c("consortium.outcome","ncase.outcome","ncontrol.outcome")]
	# Controls<-Controls[!duplicated(Controls$consortium.outcome),]
	# sum(as.numeric(Controls$ncontrol.outcome))
	# sum(as.numeric(data$ncase.outcome))
	# # Cases<-data[order(as.numeric(data$N_case),decreasing=T),c("Consortium","N_case","N_control")]
	# Cases$Consortium[is.na(Cases$Consortium)]<-"unknown"
	# Cases1<-Cases[Cases$Consortium!="Neale's lab",]
	# Cases1<-Cases1[!duplicated(Cases1$Consortium),]
	# Cases2<-Cases[Cases$Consortium=="Neale's lab",]
	# sum(as.numeric(Cases1$N_case))+sum(as.numeric(Cases2$N_case))

	# sum(as.numeric(Controls$N_control))+
	# sum(as.numeric(Cases1$N_case))+sum(as.numeric(Cases2$N_case))

	# phewas
	# phewas<-read.csv("rs11581607.csv",head=T,stringsAsFactors=F)
	# phewas<-phewas[!duplicated(paste(phewas$Trait,phewas$Beta,phewas$SE)),]
	data$trait[data$pval.outcome<0.05/nrow(data)]

	# il23<-data.frame(matrix(c("Interleukin-23 receptor","Sun et al","M & F",3301,29875488,"A","G",-0.42	,1.2E-17,0.05,"",""),nrow=1,ncol=12),stringsAsFactors=F)
	# names(il23)<-names(phewas)
	# Dat<-rbind(il23,phewas)
	data<-data[order(data$trait),]
	ID<-substr(data$trait,1,1)
	# ao<-available_outcomes()
	ao$trait<-tolower(ao$trait)
	data$trait<-tolower(data$trait)
	# data<-merge(data,ao[!duplicated(ao$trait),c("trait","category","subcategory")],by.x="trait",by.y="trait",all.x=T)
	# data$subcategory[which(data$subcategory=="NA") ]<-NA
	data$Category[which(data$Category=="") ]<-NA
	ao$subcategory[ao$subcategory=="NA"]<-NA
	ao$subcategory[ao$subcategory==""]<-NA
	ao1<-ao[!is.na(ao$subcategory),]
	miss<-data[is.na(data$Category),]
	data<-data[!is.na(data$Category),]
	miss<-merge(miss,ao1[,c("trait","subcategory")],by.x="trait",by.y="trait",all.x=T)
	# names(miss)[names(miss)=="subcategory.y"]<-"subcategory"
	# miss<-miss[,names(miss)!="subcategory.x"]
	library(plyr)
	data<-rbind.fill(miss,data)
	data1<-data[!is.na(data$subcategory),]
	data2<-data[is.na(data$subcategory),]
	for(i in 1:length(data1$trait)){
		print(data1$trait[i]) 
		Pos<-grep(data1$trait[i],data2$trait)
		data2$subcategory[Pos]<-data1$subcategory[i]
	}

	data1<-rbind(data1,data2[!is.na(data2$subcategory),])
	data2<-data2[is.na(data2$subcategory),]

	data2$subcategory2<-NA
	for(i in 1:length(data1$trait)){
		print(i)
		print(data1$trait[i])
		trait<-unlist(strsplit(data1$trait[i],split=" "))
		trait<-trait[1]
		# print(unique(data2$trait[Pos]))
		print(data1$subcategory[i])
		Pos<-grep(trait,data2$trait)
		data2$subcategory2[Pos]<-data1$subcategory[i] 
	}

	data2<-data2[,names(data2)!="subcategory"]
	names(data2)[names(data2)=="subcategory2"]<-"subcategory"
	data3$trait[which(data3$pval.outcome < 0.05/nrow(data3))]
	data3[data3$trait == "crohn's disease",]
	Temp<-data3[which(data3$Category == "Autoimmune / inflammatory" | data3$Category ==  "Inflammatory marker" | data3$Category == "Lung disease" | data3$Category == "Immune system" | data3$Category == "Immune cell-surface protein expression levels" | data3$Category == "protein"),c("trait","pval.outcome","Category","pmid.outcome","beta.outcome","effect_allele.outcome","eaf.outcome")]
	Temp<-Temp[order(as.numeric(Temp$pval.outcome)),]
	# Temp[Temp$Category == "Lung disease",]

	data3<-rbind(data1,data2)
	data3$subcategory[is.na(data3$subcategory)]<-"Other"
	length(which(is.na(data3$Category)))
	length(which(is.na(data3$subcategory)))
	# data3[data3$trait=="Interleukin-23 receptor",]
	library(ggplot2)
	library(ggrepel)
	data3$Pval<-as.numeric(data3$Pval)
	data3$trait[data3$Pval<0.05/22000]
	data3<-data3[!data3$trait %in% c("diagnoses - main icd10: k50 crohn's disease [regional enteritis]","diagnoses - main icd10: k51 ulcerative colitis","non-cancer illness code  self-reported: crohns disease","non-cancer illness code  self-reported: psoriasis"),]
	data3$trait[data3$trait %in% c( "diagnoses - main icd10: k51.9 ulcerative colitis, unspecified" ,"diagnoses - secondary icd10: k50.9 crohn's disease, unspecified","non-cancer illness code, self-reported: ulcerative colitis","non-cancer illness code, self-reported: crohns disease","non-cancer illness code, self-reported: psoriasis")]<-c("ICD10: K51.9 ulcerative colitis" ,"ICD10: K50.9 Crohn's disease","Self-reported  ulcerative colitis","Self-reported crohns disease","Self-reported psoriasis")
	         

	data3<-data3[order(data3$subcategory),]
	data3$ID<-1:nrow(data3)
	# data3[1:100,c("subcategory","ID")]
	Label<-data3$trait
	Label[as.numeric(data3$Pval)>0.05/22000]<-""
	Label<-gsub("(^[[:alpha:]])", "\\U\\1", Label, perl=TRUE)
	Label[Label %in% c("Crohn's disease","Inflammatory bowel disease","Ulcerative colitis") ]<-c("Crohn's disease (IIBDGC)","Inflammatory bowel disease (IIBDGC)","Ulcerative colitis (IIBDGC)") 
	Label[Label %in% c("ICD10: K51.9 ulcerative colitis","ICD10: K50.9 Crohn's disease","Self-reported  ulcerative colitis","Self-reported crohns disease","Self-reported psoriasis") ]<-c("ICD10: K51.9 ulcerative colitis (UK Biobank)","ICD10: K50.9 Crohn's disease (UK Biobank)","Self-reported ulcerative colitis (UK Biobank)","Self-reported Crohns disease (UK Biobank)","Self-reported psoriasis (UK Biobank)") 
	Label[Label=="Interleukin-23 receptor"]<-"Interleukin-23 receptor (Sun et al)"
	data3$subcategory2<-as.factor(data3$subcategory)
	# data3[data3$Pval<0.05/20000,c("trait","Consortium")]
	data3$Pval<--log10(data3$Pval)
	# table(data3$Consortium)

	pdf("phewas.pdf",width=15,height=8)
	ggplot(data3,aes(x=ID,y=Pval,color=subcategory))+
	       geom_point(size=2,shape=15) +
	       geom_label_repel(aes(label = Label),size = 2.5,show.legend=F,
	       			box.padding   = 0.35, 
	                  point.padding = 0.5,
	                  segment.color = 'grey50') +
	        # geom_text(aes(label=Label,hjust=1,vjust=0), position = position_dodge(width=0.9),  size=2,show.legend = FALSE) +
	         # scale_color_discrete(guide=F)+
	       # scale_x_continuous(breaks=data3$ID,labels=data3$Label) +
	       # scale_y_continuous(trans="log10")+ 
	       theme(axis.text.y = element_text(size=10),axis.text.x = element_text(angle = 45, hjust = 1,size=5),legend.title=element_text(size=10),legend.text=element_text(size=9),panel.background = element_rect(fill = "white",colour = "black",
	                                size = 0.5, linetype = "solid")) +
	       ylab("-log10 P value") +
	       xlab("trait")
	dev.off()