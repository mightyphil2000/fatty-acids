setwd("~/Google Drive/MR_base/5am")

ids<-find_id()
library(ieugwasr)

ao<-gwasinfo()
ao1<-ao[ao$id %in% ids,]
ao1<-ao1[which(ao1$subcategory != "Cancer"),]

ao1[is.na(ao1$ncase),c("trait","id","ncase","ncontrol")]
ao1<-ao1[!is.na(ao1$ncase),]

ao2<-ao[!is.na(ao$ncase),]
ao2$trait
ao2<-ao2[ao2$ncase>1000,]
ao2<-ao2[!is.na(ao2$trait),]
ao1$trait2<-ao1$trait
for(i in 1:nrow(ao1)){
	print(i)
	trait<-trimws(unlist(strsplit(ao1$trait[i],split=":")))
	trait<-trait[length(trait)]
	Test<-unlist(strsplit(trait,split=" "))
	if(nchar(Test[1])==3 & !Test[1] %in% c("ER+","ER-","Low")) {
		trait<-paste(Test[2:length(Test)],collapse=" ")
	}
	ao1$trait2[i]<-trait
}

ao1<-ao1[ao1$trait2!="Internalizing problems",]

ao1<-ao1[!ao1$trait2 %in% c("Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)",
	"ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)","ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"),]

ao2<-ao1[grep("not cancer",ao1$trait2,ignore.case=TRUE),]
ao1<-ao1[grep("cancer",ao1$trait2,ignore.case=TRUE,invert=TRUE),]
ao1<-ao1[grep("neoplasm",ao1$trait2,ignore.case=TRUE,invert=TRUE),]
ao1<-rbind(ao2,ao1)

ao1$trait2<-gsub("(^[[:alpha:]])", "\\U\\1", ao1$trait2, perl=TRUE)
ao1$trait2<-gsub("  ", " ", ao1$trait2)
ao1$trait2<-gsub("and", "&", ao1$trait2)
ao1$trait2<-gsub("\\(dvt)", "", ao1$trait2)
ao1$trait2<-trimws(ao1$trait2)
ao1<-ao1[which(as.numeric(ao1$ncase)>=1000),]

ao1$Category<-ao1$subcategory 
ao1$Category[ao1$trait2=="Coxarthrosis [arthrosis of hip]"]<-"Bone"
ao1$Category[ao1$trait2 %in% c("Muscle or soft tissue injuries","Uterine fibroids","Leiomyoma of uterus","Fibroblastic disorders")]<-"Bone/connective tissue"
ao1$Category[ao1$Category=="Bone"]<-"Bone/connective tissue"
ao1$Category[ao1$Category %in% c("Cardiovascular","Diabetes")]<-"Cardiometabolic"
ao1$Category[ao1$trait2 %in% c("Sleep disorders")]<-"Psychiatric/neurological"
ao1$Category[ao1$trait2 %in% c("Diverticular disease/diverticulitis","Gastritis & duodenitis","Acute appendicitis")]<-"Autoimmune/inflammatory"   
ao1$trait2[ao1$trait2=="Gastro-oesophageal reflux (gord) / gastric reflux" ]<-"Gastro-oesophageal reflux"
ao1$trait2[ao1$trait2=="Kidney stone/ureter stone/bladder stone"] <-"Kidney stone"
ao1$trait2[ao1$trait2=="Chronic obstructive airways disease/copd"] <-"Chronic obstructive airways disease"
ao1$trait2[ao1$trait2=="Haemorrhage from respiratory passages"] <-"Respiratory haemorrhage"

ao1$trait2[ao1$trait2=="Injury or trauma resulting in loss of vision"]<-"Loss of vision"
ao1$trait2[ao1$trait2=="Fissure & fistula of anal & rectal regions"]<-"Anal/rectal fissure"
ao1$trait2[ao1$trait2=="Excessive frequent & irregular menstruation"]<-"Excessive/irregular menstruation"


# ids<-ao1$id[!ao1$id %in% res$id.outcome]
# ao1$trait2_name[!ao1$id %in% ids]

find_id<-function(){
	res<-read.table("~/Google Drive/MR_base/5am/IL23R-findings.txt",sep="\t",head=T,stringsAsFactors=F,fill=T,quote="")
	diseases<-unlist(strsplit(res$outcome,split="\\|\\|"))
	IDS<-trimws(diseases[seq(from=2,to=length(diseases),by=2)])
	IDS<-unlist(strsplit(IDS,split="id:"))
	IDS<-IDS[seq(2,length(IDS),2)]

	IDS_numeric<-as.numeric(IDS)[!is.na(as.numeric(IDS))]

	ids1<-paste0("ieu-a-",IDS_numeric)
	ids2<-paste0("ieu-b-",IDS_numeric)

	ids3<-IDS[is.na(as.numeric(IDS))]
	ids3<-gsub(":","-",ids3)
	ids3<-gsub("UKB","ukb",ids3)
	ids<-c(ids1,ids2,ids3)
	return(ids)
}

save(list="ao1",file="~/fatty-acids/mr/data/ids_for_phewas.Rdata")

res.info<-merge(res,ao1[,c("id","ncase","ncontrol","trait2","Category","consortium")],by.x="id.outcome",by.y="id")




res$id.outcome
ao1$id

res.info$Category[res.info$Category %in% c("Eye","Kidney","Kindey","Lung")]<-"Other"
res.info<-res.info[order(res.info$b_pQTL,decreasing=T),]
res.info<-res.info[order(res.info$Category),]
res.info$b_pQTL<-res.info$b_pQTL*-1 

res.info<-res.info[!res.info$trait2 %in% c("Retinal detachment","High cholesterol"),]
lci<-exp(res.info$b_pQTL-1.96*res.info$se_pQTL)
uci<-exp(res.info$b_pQTL+1.96*res.info$se_pQTL)
min(lci)
max(uci)

res.info$pval<-formatC(res.info$p_pQTL, format = "e", digits = 2)
res.info1<-res.info[1:50,]
res.info2<-res.info[51:100,]


pdf("disease1.pdf")
forest_plot_1_to_many(res.info1,b="b_pQTL",se="se_pQTL",
    exponentiate=T,ao_slc=F,lo=0.35,up=8,
    trait2M="trait2",by=NULL,
    trans="log2",xlab="",
    col1_width=1.50, addcol_titles=c("",""),
    addcols=c("N_case","pval"),addcol_widths=c(0.35,0.50))
dev.off()

# res.info<-res.info[order(res.info$Category),]

res.other<-res.info[res.info$Category=="Other",]
res.bone<-res.info[res.info$Category=="Bone/connective tissue",]
res.info2<-res.info[!res.info$Category %in% c("Other","Bone/connective tissue"),]
res.info2<-rbind(res.info2,res.bone)
res.info2<-rbind(res.info2,res.other)
res.info2[,c("trait2","Category")]

table(res.info2$Category)
res.infoA<-res.info2[1:49,]
res.infoB<-res.info2[50:99,]
pdf("diseaseA_cat.pdf",height=8,width=7.5)
forest_plot_1_to_many(res.infoA,b="b_pQTL",se="se_pQTL",
    exponentiate=T,ao_slc=F,lo=0.10,up=2.6,
    trait2M="trait2",by="Category",
    trans="log2",xlab="",
    col1_width=1.50, addcol_titles=c("",""),
    subheading_size=7,
    addcols=c("N_case","pval"),
    addcol_widths=c(0.35,0.50),
    shape_points=15,
    colour_scheme="black")
dev.off()

Controls<-res.info[order(as.numeric(res.info$N_control),decreasing=T),c("Consortium","N_case","N_control")]
Controls<-Controls[!duplicated(Controls$Consortium),]
sum(as.numeric(Controls$N_control))
sum(as.numeric(res.info$N_case))
Cases<-res.info[order(as.numeric(res.info$N_case),decreasing=T),c("Consortium","N_case","N_control")]
Cases$Consortium[is.na(Cases$Consortium)]<-"unknown"
Cases1<-Cases[Cases$Consortium!="Neale's lab",]
Cases1<-Cases1[!duplicated(Cases1$Consortium),]
Cases2<-Cases[Cases$Consortium=="Neale's lab",]
sum(as.numeric(Cases1$N_case))+sum(as.numeric(Cases2$N_case))

sum(as.numeric(Controls$N_control))+
sum(as.numeric(Cases1$N_case))+sum(as.numeric(Cases2$N_case))

info<-read.table("outcome-info.txt",sep="\t",head=T,stringsAsFactors=F,fill=T,quote="",comment.char="",as.is=TRUE)
info$disease<-"notdisease"
info$disease[info$trait2_name %in% c("Diagnoses - main ICD10: G56 Mononeuropathies of upper limb",
	"Diagnoses - main ICD10: I48 Atrial fibrillation and flutter",
	"Diagnoses - main ICD10: I80 Phlebitis and thrombophlebitis",
	"Diagnoses - main ICD10: I83 Varicose veins of lower extremities",
	"Diagnoses - main ICD10: K57 Diverticular disease of intestine",
	"Diagnoses - main ICD10: K80 Cholelithiasis",
	"Non-cancer illness code  self-reported: psoriasis",
	"Non-cancer illness code  self-reported: osteoarthritis",
	"Non-cancer illness code  self-reported: gout",
	"Non-cancer illness code  self-reported: pneumothorax",
	"Non-cancer illness code  self-reported: polio / poliomyelitis",
	"Non-cancer illness code  self-reported: arthritis (nos)",
	"Non-cancer illness code  self-reported: hypertrophic cardiomyopathy (hcm / hocm)",
	"Non-cancer illness code  self-reported: vitiligo",
	"Underlying (primary) cause of death: ICD10: E85.4 Organ-limited amyloidosis",
	"Eye problems/disorders: Diabetes related eye disease",
	"Eye problems/disorders: Glaucoma",
	"Vascular/heart problems diagnosed by doctor: Angina",
	"Diagnoses - main ICD10: C61 Malignant neoplasm of prostate",
	"Diagnoses - main ICD10: D12 Benign neoplasm of colon  rectum  anus and anal canal",
	"Diagnoses - main ICD10: I30 Acute pericarditis",
	"Diagnoses - main ICD10: K20 Oesophagitis",
	"Diagnoses - main ICD10: K29 Gastritis and duodenitis",
	"Diagnoses - main ICD10: K35 Acute appendicitis",
	"Cancer code  self-reported: small intestine/small bowel cancer",
	"Diagnoses - main ICD10: M16 Coxarthrosis [arthrosis of hip]",
	"Diagnoses - main ICD10: M17 Gonarthrosis [arthrosis of knee]",
	"Diagnoses - main ICD10: M54 Dorsalgia",
	"Prostate cancer",
	"Diagnoses - main ICD10: D25 Leiomyoma of uterus",
	"Diagnoses - main ICD10: G47 Sleep disorders",
	"Diagnoses - main ICD10: M72 Fibroblastic disorders",
	"Diagnoses - main ICD10: N20 Calculus of kidney and ureter",
	"Diagnoses - main ICD10: N81 Female genital prolapse",
	"Diagnoses - main ICD10: R04 Haemorrhage from respiratory passages",
	"Diagnoses - main ICD10: R07 Pain in throat and chest",
	"Diagnoses - main ICD10: R10 Abdominal and pelvic pain",
	"Diagnoses - main ICD10: R11 Nausea and vomiting",
	"Diagnoses - main ICD10: R14 Flatulence and related conditions",
	"Diagnoses - main ICD10: R35 Polyuria",
	"Diagnoses - main ICD10: R55 Syncope and collapse",
	"Cancer code  self-reported: basal cell carcinoma",
	"Non-cancer illness code  self-reported: hypertension",
	"Non-cancer illness code  self-reported: pulmonary embolism +/- dvt",
	"Non-cancer illness code  self-reported: deep venous thrombosis (dvt)",
	"Non-cancer illness code  self-reported: asthma",
	"Non-cancer illness code  self-reported: chronic obstructive airways disease/copd",
	"Non-cancer illness code  self-reported: emphysema/chronic bronchitis",
	"Non-cancer illness code  self-reported: sleep apnoea",
	"Non-cancer illness code  self-reported: gastro-oesophageal reflux (gord) / gastric reflux",
	"Non-cancer illness code  self-reported: kidney stone/ureter stone/bladder stone",
	"Non-cancer illness code  self-reported: bladder problem (not cancer)",
	"Non-cancer illness code  self-reported: hyperthyroidism/thyrotoxicosis",
	"Non-cancer illness code  self-reported: hypothyroidism/myxoedema",
	"Non-cancer illness code  self-reported: migraine",
	"Non-cancer illness code  self-reported: depression",
	"Non-cancer illness code  self-reported: anxiety/panic attacks",
	"Non-cancer illness code  self-reported: mania/bipolar disorder/manic depression",
	"Non-cancer illness code  self-reported: bone disorder",
	"Non-cancer illness code  self-reported: joint disorder",
	"Non-cancer illness code  self-reported: osteoporosis",
	"Non-cancer illness code  self-reported: ankylosing spondylitis",
	"Non-cancer illness code  self-reported: iron deficiency anaemia",
	"Non-cancer illness code  self-reported: pernicious anaemia")]<-"disease"

info$disease[grep("illness",info$trait2_name)]<-"disease"
info$disease[grep("Diagnoses",info$trait2_name)]<-"disease"
info$disease[grep("disorders",info$trait2_name)]<-"disease"

# phewas
phewas<-read.csv("rs11581607.csv",head=T,stringsAsFactors=F)
# phewas<-phewas[!duplicated(paste(phewas$trait2,phewas$Beta,phewas$SE)),]
phewas$trait2[phewas$Pval<0.05/10000]

il23<-data.frame(matrix(c("Interleukin-23 receptor","Sun et al","M & F",3301,29875488,"A","G",-0.42	,1.2E-17,0.05,"",""),nrow=1,ncol=12),stringsAsFactors=F)
names(il23)<-names(phewas)
Dat<-rbind(il23,phewas)
Dat<-Dat[order(Dat$trait2),]
ID<-substr(Dat$trait2,1,1)
# ao<-available_outcomes()
ao$trait2<-tolower(ao$trait2)
Dat$trait2<-tolower(Dat$trait2)
Dat<-merge(Dat,ao[!duplicated(ao$trait2),c("trait2","category","subcategory")],by.x="trait2",by.y="trait2",all.x=T)
Dat$subcategory[which(Dat$subcategory=="NA") ]<-NA
Dat$subcategory[which(Dat$subcategory=="") ]<-NA
ao$subcategory[ao$subcategory=="NA"]<-NA
ao$subcategory[ao$subcategory==""]<-NA
ao1<-ao[!is.na(ao$subcategory),]
miss<-Dat[is.na(Dat$subcategory),]
Dat<-Dat[!is.na(Dat$subcategory),]
miss<-merge(miss,ao1[,c("trait2","subcategory")],by.x="trait2",by.y="trait2",all.x=T)
names(miss)[names(miss)=="subcategory.y"]<-"subcategory"
miss<-miss[,names(miss)!="subcategory.x"]
Dat<-rbind(miss,Dat)
Dat1<-Dat[!is.na(Dat$subcategory),]
Dat2<-Dat[is.na(Dat$subcategory),]
for(i in 1:length(Dat1$trait2)){
	print(Dat1$trait2[i]) 
	Pos<-grep(Dat1$trait2[i],Dat2$trait2)
	Dat2$subcategory[Pos]<-Dat1$subcategory[i]
}

Dat1<-rbind(Dat1,Dat2[!is.na(Dat2$subcategory),])
Dat2<-Dat2[is.na(Dat2$subcategory),]

Dat2$subcategory2<-NA
for(i in 1:length(Dat1$trait2)){
	print(i)
	print(Dat1$trait2[i])
	trait2<-unlist(strsplit(Dat1$trait2[i],split=" "))
	trait2<-trait2[1]
	# print(unique(Dat2$trait2[Pos]))
	print(Dat1$subcategory[i])
	Pos<-grep(trait2,Dat2$trait2)
	Dat2$subcategory2[Pos]<-Dat1$subcategory[i] 
}

Dat2<-Dat2[,names(Dat2)!="subcategory"]
names(Dat2)[names(Dat2)=="subcategory2"]<-"subcategory"

Dat3<-rbind(Dat1,Dat2)
Dat3$subcategory[is.na(Dat3$subcategory)]<-"Other"

# Dat3[Dat3$trait2=="Interleukin-23 receptor",]
library(ggplot2)
library(ggrepel)
Dat3$Pval<-as.numeric(Dat3$Pval)
Dat3$trait2[Dat3$Pval<0.05/22000]
Dat3<-Dat3[!Dat3$trait2 %in% c("diagnoses - main icd10: k50 crohn's disease [regional enteritis]","diagnoses - main icd10: k51 ulcerative colitis","non-cancer illness code  self-reported: crohns disease","non-cancer illness code  self-reported: psoriasis"),]
Dat3$trait2[Dat3$trait2 %in% c( "diagnoses - main icd10: k51.9 ulcerative colitis, unspecified" ,"diagnoses - secondary icd10: k50.9 crohn's disease, unspecified","non-cancer illness code, self-reported: ulcerative colitis","non-cancer illness code, self-reported: crohns disease","non-cancer illness code, self-reported: psoriasis")]<-c("ICD10: K51.9 ulcerative colitis" ,"ICD10: K50.9 Crohn's disease","Self-reported  ulcerative colitis","Self-reported crohns disease","Self-reported psoriasis")
         

Dat3<-Dat3[order(Dat3$subcategory),]
Dat3$ID<-1:nrow(Dat3)
# Dat3[1:100,c("subcategory","ID")]
Label<-Dat3$trait2
Label[as.numeric(Dat3$Pval)>0.05/22000]<-""
Label<-gsub("(^[[:alpha:]])", "\\U\\1", Label, perl=TRUE)
Label[Label %in% c("Crohn's disease","Inflammatory bowel disease","Ulcerative colitis") ]<-c("Crohn's disease (IIBDGC)","Inflammatory bowel disease (IIBDGC)","Ulcerative colitis (IIBDGC)") 
Label[Label %in% c("ICD10: K51.9 ulcerative colitis","ICD10: K50.9 Crohn's disease","Self-reported  ulcerative colitis","Self-reported crohns disease","Self-reported psoriasis") ]<-c("ICD10: K51.9 ulcerative colitis (UK Biobank)","ICD10: K50.9 Crohn's disease (UK Biobank)","Self-reported ulcerative colitis (UK Biobank)","Self-reported Crohns disease (UK Biobank)","Self-reported psoriasis (UK Biobank)") 
Label[Label=="Interleukin-23 receptor"]<-"Interleukin-23 receptor (Sun et al)"
Dat3$subcategory2<-as.factor(Dat3$subcategory)
# Dat3[Dat3$Pval<0.05/20000,c("trait2","Consortium")]
Dat3$Pval<--log10(Dat3$Pval)
# table(Dat3$Consortium)
