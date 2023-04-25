library(ggforestplot)
library(ggplot2)
library(TwoSampleMR)
source("~/fatty-acids/mr/scripts/mr_functions.R")
exp<-read.table("~/fatty-acids/mr/data/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",head=T,stringsAsFactors=F)
exp[exp$population=="East Asian",]
exposure_dat<-format_exposure2(dat=exp,standardise_beta=TRUE)

exposure_dat<-exposure_dat[exposure_dat$SNP == "rs174546",]
eur<-exposure_dat[exposure_dat$population == "European",]
eur<-eur[eur$exposure == "AA:DGLA",]
eas<-exposure_dat[exposure_dat$population == "East Asian",]
eas<-eas[eas$exposure == "GLA:LA",]

ids<-find_id()
# load("~/fatty-acids/mr/data/ids_for_phewas.Rdata")
out_dat <- extract_outcome_data(
    snps = unique(exposure_dat$SNP),
    outcomes = ids

)

# out_dat$outcome[grep("Crohn",out_dat$outcome)]
dim(out_dat)
ao<-ieugwasr::gwasinfo()

id_stroke<-ao$id[grep("stroke",ao$trait,ignore.case=TRUE)]

out_dat_stroke<-extract_outcome_data(
    snps = unique(exposure_dat$SNP),
    outcomes = id_stroke)

id_hemorrhage<-ao$id[grep("hemorrhage",ao$trait,ignore.case=TRUE)]
id_hemorrhage2<-ao$id[grep("bleeding",ao$trait,ignore.case=TRUE)]
id_hemorrhage3<-ao$id[grep("haemorrhage",ao$trait,ignore.case=TRUE)]



# Dat<-ieugwasr::associations(variants = unique(exposure_dat$SNP),id = id_inflammatory)

out_dat_hemorrhage<-extract_outcome_data(
    snps = unique(exposure_dat$SNP),
    outcomes = unique(c(id_hemorrhage,id_hemorrhage2,id_hemorrhage3)))
out_dat_inflammatory<-read.table("~/fatty-acids/mr/results/inflammation_lookups_opengwas.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE,quote="")

id_inflammatory<-out_dat_inflammatory$id

# out_dat_inflammatory<-format_outcome2(dat=out_dat_inflammatory)

out_dat_inflammatory<-extract_outcome_data(
    snps = unique(exposure_dat$SNP),
    outcomes = id_inflammatory)

# id_miss3<-ao1$id[91:length(ao1$id)][!ao1$id[91:length(ao1$id)] %in% out_dat10$id]

# id_miss<-c(id_miss1,id_miss2,id_miss3)

# ao<-available_outcomes()
# out_dat11<-extract_outcome_data(
#     snps = unique(exposure_dat$SNP),
#     outcomes = id_miss)


out_dat_all<-do.call(plyr::rbind.fill,list(out_dat,out_dat_hemorrhage,out_dat_stroke,out_dat_inflammatory))
# out_dat_all<-do.call(plyr::rbind.fill,list(out_dat,out_dat_hemorrhage,out_dat_stroke))
out_dat_all<-out_dat_all[!duplicated(out_dat_all$id.outcome),]
dim(out_dat_all)
# out_dat_all<-do.call(plyr::rbind.fill,list(out_dat1,out_dat2,out_dat3,out_dat4,out_dat5,out_dat6,out_dat7,out_dat8,out_dat9,out_dat10,out_dat_hemorrhage,out_dat_stroke,out_dat_inflammatory))

out_dat_all[,c("effect_allele.outcome","other_allele.outcome","eaf.outcome")]


dat <- harmonise_data(
    exposure_dat = eur, 
    outcome_dat = out_dat_all
)


mr_results<-mr(dat, method_list="mr_wald_ratio")

save(list=c("mr_results","dat","out_dat_all","eur","ao"),file="~/fatty-acids/mr/results/mr_results_phewas_v2.Rdata")
load("~/fatty-acids/mr/results/mr_results_phewas_v2.Rdata")
dim(mr_results)
head(mr_results)

# blood clotting disorder;
#  hemorrhagic stroke 
# gastrointestinal hemorrhage

ao1<-cleanup_ao(ao=ao)


# ao1<-cleanup_ao(ao=NULL)
mr.ao<-merge(mr_results,ao1,by.x="id.outcome",by.y="id")
mr.ao1<-mr.ao[mr.ao$batch %in% c("ukb-a","ukb-b","ukb-d","ukb-e"),  ]
# mr.ao1<-mr.ao1[!mr.ao1$id.outcome %in% c("ukb-b-19354","ukb-b-18113","ukb-b-9127" ,"ukb-b-20141" ),]
mr.ao2<-mr.ao[!mr.ao$batch %in% c("ukb-a","ukb-b","ukb-d","ukb-e"),  ]
# mr.ao3<-mr.ao[mr.ao$id.outcome %in% c("ukb-b-19354","ukb-b-18113","ukb-b-9127" ,"ukb-b-20141" ),]
# mr.ao2<-rbind(mr.ao2,mr.ao3)

mr.ao1<-CheckSumStats::transform_betas(dat=mr.ao1,effect="b",effect.se="se")
mr.ao<-rbind(mr.ao1,mr.ao2)

# mr.ao[grep("Eczema",mr.ao$trait2,ignore.case=TRUE),]
# mr.ao1<-mr.ao[1:41,]
# mr.ao2<-mr.ao[43:83,]

# mr.can<-load_cancer_results()
# as.numeric(mr.can[,c("pval")])

# mr.can[mr.can$trait2==unique(mr.can$trait2)[7],]
# mr.ao<-plyr::rbind.fill(mr.can,mr.ao)
# mr.can[mr.can$outcome == "Colorectal cancer","b"]


# sort(unique(mr_res2$outcome))

plot_dat<-format_plot_data()
save.image(file="~/fatty-acids/mr/results/mr_results_phewas_v3.Rdata")
load("~/fatty-acids/mr/results/mr_results_phewas_v3.Rdata")
plot_dat[plot_dat$outcome=="Pernicious anaemia\ncases=1020",]

plot_dat[plot_dat$outcome=="Pernicious anaemia\ncases=1020",]
plot_dat[plot_dat$Category=="FBleeding disorders",c("trait2","or","lci","uci")]

sum(plot_dat$ncase)
sum(plot_dat$ncontrol,na.rm=TRUE)
plot_dat$outcome[1:10]
plot_dat$outcome[plot_dat$Category=="BAutoimmune/inflammatory" ]
plot_dat$outcome[plot_dat$pval<0.05/(nrow(plot_dat)-6)]
plot_dat[plot_dat$pval<0.05/nrow(plot_dat[plot_dat$Category!="Cancer",]),c("outcome","or","lci","uci","pval")]
plot_dat[plot_dat$pval<0.05/nrow(plot_dat[plot_dat$Category!="Cancer",]),c("outcome","or","lci","uci","pval","Category")]
table(plot_dat$Category)

plot_dat<-plot_dat[order(plot_dat$Category),]
# plot_dat$b<-plot_dat$b*-1
Min<-min(plot_dat$lci)
Max<-max(plot_dat$uci)
# Max<-max(1/plot_dat$lci)
# Min<-min(1/plot_dat$uci)
Max<-Max+0.01
Min<-Min-0.01

# plot_dat<-plot_dat[plot_dat$Category !="ACancer",]
dim(plot_dat)
plot_dat1<-plot_dat[1:42,]
plot_dat2<-plot_dat[43:84,]
# plot_dat1[plot_dat1$outcome=="Esophageal squamous cell carcinoma\ncases=2013",c("or","lci","uci" )]
table(plot_dat1$Category)
table(plot_dat2$Category)
dim(plot_dat1)
dim(plot_dat2)
Psig<-0.05/nrow(plot_dat)
# plot_dat$outcome[which(plot_dat$pval<Psig)]
Color<-as.numeric(plot_dat$pval)
Pos1<-Color<=Psig
Pos2<-Color>Psig
Color[Pos1]<-"red"
Color[Pos2]<-"black"

# head(plot_dat)
# a<-plot_dat[plot_dat$Category == "BAutoimmune/inflammatory",]
# head(a)
# a<-a[,c("trait2","id.outcome","batch","year","group_name","author","consortium","sex","population","unit","ncase","ncontrol","sample_size","pmid","trait","note")]
# head(a)
# names(a)[names(a)=="trait"]<-"originalname"
# names(a)[names(a)=="trait2"]<-"outcome_nice_name"
# write.table("a","~/protein_prs_ukb/data/autoimmune_infl_diseases.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

P1<-forestplot(df = plot_dat1,logodds = TRUE,name=outcome,estimate=b,se=se,xlab = "", psignif = 0.05/nrow(plot_dat))+theme(legend.position = "none")+geom_point(shape=15,fill="black",colour=Color[1:42])  +xlim(Min, Max)+theme(plot.title = element_text(size = ""),text = element_text(size=12))
P1<-P1+ggforce::facet_col(
    facets = ~Category,
    scales = "free_y",
    space = "free" )


P2<-forestplot(df = plot_dat2,logodds = TRUE,name=outcome,estimate=b,se=se,shape=NULL,xlab = "", psignif = 0.05/nrow(plot_dat))+theme(legend.position = "none")+geom_point(shape=15,fill="black",colour=Color[43:84]) +xlim(Min, Max)+theme(plot.title = element_text(size = ""),text = element_text(size=12))
P2<-P2+ggforce::facet_col(
    facets = ~Category,
    scales = "free_y",
    space = "free" )

png("~/fatty-acids/mr/results/plots/disease_phewas1_inhibition_v3.png",width = 500, height = 1000)
	P1
dev.off()

png("~/fatty-acids/mr/results/plots/disease_phewas2_inhibition_v3.png",width = 500, height = 1000)
	P2
dev.off()



find_id<-function(){
	res<-read.table("~/Google Drive/MR_base/5am/IL23R-findings.txt",sep="\t",head=T,stringsAsFactors=F,fill=T,quote="")

	# head(res)
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

# save(list="ao1",file="~/fatty-acids/mr/data/ids_for_phewas.Rdata")


format_outcome2<-function(dat=NULL){
	names(dat)<-paste0(names(dat),".outcome")
	names(dat)<-gsub("trait2","originalname",names(dat))
	names(dat)[names(dat) == "ea.outcome"]<-"effect_allele.outcome"
	names(dat)[names(dat) == "nea.outcome"]<-"other_allele.outcome"
	names(dat)<-gsub("trait.outcome","outcome",names(dat))
	names(dat)<-gsub("p.outcome","pval.outcome",names(dat))
	dat<-dat[,names(dat) %in% names(out_dat_stroke)]
	# names(out_dat_stroke)[!names(out_dat_stroke) %in% names(dat)]
	dat$SNP<-"rs174546"
	dat$other_allele.outcome<-"T"
	return(dat)
}

 
cleanup_ao<-function(ao=NULL){
	if(is.null(ao)){
		ao<-ieugwasr::gwasinfo()
	}

	ao1<-ao[ao$id %in% c(ids,id_stroke,id_hemorrhage,id_hemorrhage2,id_inflammatory),]
	ao1<-ao1[which(ao1$subcategory != "Cancer"),]

	# ao1[is.na(ao1$ncase),c("trait","id","ncase","ncontrol")]
	ao1<-ao1[!is.na(ao1$ncase),]
	# ao1.1<-ao1[ao1$trait=="Gastrointestinal hemorrhage",] 
	ao1<-ao1[ao1$ncase>=1000,]
	# ao1<-rbind(ao1,ao1.1)

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

	ao2<-ao1[grep("not cancer",ao1$trait2,ignore.case=TRUE),]
	ao1<-ao1[grep("cancer",ao1$trait2,ignore.case=TRUE,invert=TRUE),]
	ao1<-ao1[grep("neoplasm",ao1$trait2,ignore.case=TRUE,invert=TRUE),]
	ao1<-ao1[grep("carcinoma",ao1$trait2,ignore.case=TRUE,invert=TRUE),]
	ao1<-ao1[ao1$trait2!=	"Malignant melanoma",]
	ao1<-rbind(ao2,ao1)
	ao1$trait2<-gsub("(^[[:alpha:]])", "\\U\\1", ao1$trait2, perl=TRUE)
	ao1$trait2<-gsub("  ", " ", ao1$trait2)
	ao1$trait2<-gsub("and", "&", ao1$trait2)
	ao1$trait2<-gsub("\\(dvt)", "", ao1$trait2)
	ao1$trait2<-trimws(ao1$trait2)
	# ao1<-ao1[which(as.numeric(ao1$ncase)>=1000),]
	ao1$Category<-ao1$subcategory 
	
	ao1$Category[ao1$trait2=="Coxarthrosis [arthrosis of hip]"]<-"Bone"
	ao1$Category[ao1$trait2 %in% c("Muscle or soft tissue injuries","Uterine fibroids","Leiomyoma of uterus","Fibroblastic disorders")]<-"Bone/connective tissue"
	ao1$Category[ao1$Category=="Bone"]<-"Bone/connective tissue"
	ao1$Category[ao1$Category %in% c("Cardiovascular","Diabetes")]<-"Cardiometabolic"
	ao1$Category[ao1$trait2 %in% c("Sleep disorders")]<-"Psychiatric/neurological"
	ao1$Category[ao1$trait2 %in% c("Diverticular disease/diverticulitis","Gastritis & duodenitis","Acute appendicitis","Eczema/dermatitis")]<-"Autoimmune/inflammatory"   
	ao1$trait2[ao1$trait2=="Gastro-oesophageal reflux (gord) / gastric reflux" ]<-"Gastro-oesophageal reflux"
	ao1$trait2[ao1$trait2=="Kidney stone/ureter stone/bladder stone"] <-"Kidney stone"
	ao1$trait2[ao1$trait2=="Chronic obstructive airways disease/copd"] <-"Chronic obstructive airways disease"
	ao1$trait2[ao1$trait2=="Haemorrhage from respiratory passages"] <-"Respiratory haemorrhage"

	ao1$trait2[ao1$trait2=="Injury or trauma resulting in loss of vision"]<-"Loss of vision"
	ao1$trait2[ao1$trait2=="Fissure & fistula of anal & rectal regions"]<-"Anal/rectal fissure"
	ao1$trait2[ao1$trait2=="Excessive frequent & irregular menstruation"]<-"Excessive/irregular menstruation"

	ao1<-ao1[!ao1$trait2 %in% c("Retinal detachment","High cholesterol"),]

	ao1$Category[unique(unlist(lapply(c("bleeding","haemorrhage","hemorrhage","anaemia", "Excessive/irregular menstruation"   ),FUN=function(x)
		grep(x,ao1$trait2,ignore.case=TRUE))))]<-"Bleeding disorders"	
	ao1$Category[grep("stroke",ao1$trait2,ignore.case=TRUE)]<-"Cardiometabolic"

	ao1$Category[unique(unlist(lapply(c("Asthma","Hayfever/allergic rhinitis","Nasal polyp","Oesophagitis","bronchitis","Chronic obstructive airways disease","Psoriasis" ,"Phlebitis & thrombophlebitis"  ),FUN=function(x)
		grep(x,ao1$trait2,ignore.case=TRUE))))]<-"Autoimmune/inflammatory"
	ao1$Category[unique(unlist(lapply(c("Arthritis" ,"Gout","Osteoporosis", "Joint disorder" ,"Fractured/broken bones in last 5 years","Internal derangement of knee","Gonarthrosis" ),FUN=function(x)
		grep(x,ao1$trait2,ignore.case=TRUE))))]<-"Bone/connective tissue"  

	ao1$Category[unique(unlist(lapply(c("Pulmonary embolism" ,"Angina","hypertension","diabetes","thrombosis","Varicose veins","Atrial fibrillation & flutter" ),FUN=function(x)
		grep(x,ao1$trait2,ignore.case=TRUE))))]<-"Cardiometabolic" 
ao1$Category[unique(unlist(lapply(c("Loss of vision","Sleep apnoea" ,"Hearing difficulty/problems" ,"Migraine","Depression","Anxiety/panic attacks"  ),FUN=function(x)
		grep(x,ao1$trait2,ignore.case=TRUE))))]<-"Psychiatric/neurological"
	ao1<-ao1[!ao1$trait2 %in% c("Ankle","Other bones","NA Arm" ,"Wrist",
		"Pain in throat & chest",
		"Hiatus hernia", "Inguinal hernia" ,"Diaphragmatic hernia"     , 
		"Syncope & collapse","Ventral hernia" ,"Cough on most days","Vaginal prolapse/uterine prolapse","Enlarged prostate", "Anal/rectal fissure","Cellulitis"      ,                    "Hyperplasia of prostate"  ,"Fracture resulting from simple fall", "Mononeuropathies of upper limb" ,   "Haemorrhoids" , "Polyuria","Female genital prolapse"),]


	ao1$Category[which(ao1$Category=="NA")]<-"Other"
	Pos<-unlist(regexpr("-",ao1$id))
	ao1$batch<-substring(ao1$id,1,Pos+1)
	ao1$Category<-gsub(" / ","/",ao1$Category)
	ao1<-ao1[order(ao1$ncase,decreasing=TRUE),]
	ao1<-ao1[!duplicated(ao1$trait2),]	

	# ao1$trait2[ao1$Category=="Autoimmune/inflammatory" ]

	# "\"Stroke, including SAH\" (no controls excluded)
	# "\"Stroke, excluding SAH\" (no controls excluded)\

	ao1<-ao1[ao1$trait2 != "Malignant melanoma",]
	ao1$trait2[ao1$trait2 == "STROKE"]<-"Stroke"
	ao1$trait2[ao1$trait2 == "Ischemic stroke (cardioembolic)"]<-"Cardioembolic stroke"
	ao1$trait2[ao1$trait2 == "Ischemic stroke (small-vessel)"]<-"Small-vessel stroke"
	ao1$trait2[grep("Stroke, including SAH",ao1$trait2)]<-"Stroke (incl SAH)"
	ao1$trait2[grep("Ischaemic Stroke, excluding all haemorrhages",ao1$trait2,ignore.case=TRUE)]<-"Ischaemic Stroke (excl haemorrhages)" 

	ao1$trait2[grep("Stroke, excluding SAH",ao1$trait2)]<-"Stroke (excl SAH)"
	ao1$trait2[ao1$trait2 ==  "Ischemic stroke (large artery atherosclerosis)"]<-"Large artery stroke)"
	ao1$trait2[grep("Vitreous bleeding",ao1$trait2,ignore.case=TRUE)]<-"Vitreous bleeding"
	
	ao1$trait2[grep("Eczema",ao1$trait2,ignore.case=TRUE)]<-"Eczema/dermatitis"
	# ao1$trait2[grep("cholelithiasis",ao1$trait2,ignore.case=TRUE)]
# 
	Pos<-ao1$trait2=="NA Yes" 
	trait<-ao1$trait[Pos]
	trait<-gsub(": Yes","",trait)
	ao1$trait2[Pos]<-trait
	ao1$trait2[grep("N93",ao1$trait2)]<-"Abnormal uterine & vaginal bleeding"
	ao1$trait2[grep("N95",ao1$trait2)]<-"Postmenopausal bleeding"
	ao1$trait2[ao1$trait2=="Other abnormal uterine & vaginal bleeding"]<-"Abnormal uterine & vaginal bleeding"

	return(ao1)
}



load_cancer_results<-function(){
	# source("~/fatty-acids/mr/scripts/mr_plot_functions.R")
	load("~/fatty-acids/mr/results/mr_results_rep_v2.Rdata")
	a1<-mr_res1[mr_res1$population == "European",]
	a2<-mr_res1[mr_res1$population == "East Asian",]
	a1<-a1[which(a1$exposure == "AA:DGLA"),]
	a2<-a2[which(a2$exposure == "GLA:LA"),]
	mr_res1<-rbind(a1,a2)
	Pos<-unique(unlist(lapply(c("colorectal","lung","esophageal squamous","malignant skin cancer","non-melanoma skin cancer","all cause"),FUN=function(x) grep(x,mr_res1$outcome,ignore.case=TRUE)) ))
	mr_res2<-mr_res1[Pos,]
	mr_res2<-mr_res2[!is.na(mr_res2$cases),]
	mr_res2$cases<-as.numeric(mr_res2$cases)
	mr_res2<-mr_res2[order(mr_res2$cases,decreasing=TRUE),]
	Pos<-regexpr("[0-9]",mr_res2$outcome)
	mr_res3<-mr_res2[which(Pos==-1),]
	mr_res4<-mr_res2[which(Pos!=-1),]
	Pos<-regexpr("[0-9]",mr_res4$outcome)
	mr_res4$outcome<-substring(mr_res4$outcome,1,Pos-2)
	# mr_res4[mr_res4$outcome == "Colorectal cancer",]
	mr_res2<-rbind(mr_res4,mr_res3)
	mr_res2<-mr_res2[!duplicated(mr_res2$outcome),]
	mr_res2<-mr_res2[grep("males",mr_res2$outcome,invert=TRUE),]
	mr_res2<-mr_res2[grep("females",mr_res2$outcome,invert=TRUE),]
	mr_res2<-mr_res2[grep("smokers",mr_res2$outcome,invert=TRUE),]
	mr_res2<-mr_res2[grep("Distal",mr_res2$outcome,invert=TRUE),]
	mr_res2<-mr_res2[grep("adenocarcinoma",mr_res2$outcome,invert=TRUE),]
	mr_res2<-mr_res2[grep("Proximal",mr_res2$outcome,invert=TRUE),]
	mr_res2<-mr_res2[grep("Small cell",mr_res2$outcome,invert=TRUE),]
	mr_res2<-mr_res2[grep("Squamous cell lung cancer",mr_res2$outcome,invert=TRUE),]
	mr_res2$Category<-"Cancer"
	names(mr_res2)[names(mr_res2) == "cases"]<-"ncase"	
	mr_res2$b<-as.numeric(mr_res2$b)
	mr_res2$se<-as.numeric(mr_res2$se)
	mr_res2$ncase<-as.numeric(mr_res2$ncase)
	mr_res2$trait2<-mr_res2$outcome
	mr_res2<-mr_res2[mr_res2$trait2 !="Cancer (excluding non-melanoma skin cancer)",]
		mr_res2<-mr_res2[mr_res2$trait2 !="Cancer (excluding non-melanoma skin cancer)",]
	mr_res2$trait2[mr_res2$trait2== "Cancer (all cause)"]<-"Overall cancer"   
	# mr_res2[mr_res2$trait2==mr_res2$trait2[7],]
	return(mr_res2)
}

format_plot_data<-function(){
	
	mr.ao<-mr.ao[order(mr.ao$ncase,decreasing=TRUE),]
	mr.ao<-mr.ao[!duplicated(mr.ao$trait2),]
	mr.ao$outcome<-paste0(mr.ao$trait2,"\ncases=",mr.ao$ncase)


	# plot_dat$outcome[grep("stroke",plot_dat$outcome)]
	mr.ao$pval<-as.numeric(mr.ao$pval)
	groups<-sort(unique(mr.ao$Category))
	mr.ao<-mr.ao[order(mr.ao$Category,mr.ao$ncase,decreasing=TRUE),]
	plot1<-mr.ao[mr.ao$Category == "Cancer",]
	plot4<-mr.ao[mr.ao$Category ==  "Autoimmune/inflammatory",]
	plot5<-mr.ao[mr.ao$Category ==  "Bleeding disorders" ,]
	# mr.ao$Category[mr.ao$Category ==  "bleeding disorder"]<-"Bleeding disorders"
	plot3<-mr.ao[mr.ao$Category ==  "Bone/connective tissue" ,]
	plot2<-mr.ao[mr.ao$Category ==  "Cardiometabolic"  ,]
	plot6<-mr.ao[mr.ao$Category ==  "Psychiatric/neurological" ,]
	plot7<-mr.ao[mr.ao$Category ==  "Other" ,]
	plot_dat<-do.call(rbind,list(plot1,plot2,plot3,plot4,plot5,plot6,plot7))
	

	plot_dat$Category[plot_dat$Category=="Cancer"]<-paste0("ACancer")
	plot_dat$Category[plot_dat$Category=="Autoimmune/inflammatory"]<-paste0("BAutoimmune/inflammatory")

	plot_dat$Category[plot_dat$Category=="Cardiometabolic"]<-paste0("ECardiometabolic")
	plot_dat$Category[plot_dat$Category=="Psychiatric/neurological"]<-paste0("DPsychiatric/neurological")
	plot_dat$Category[plot_dat$Category=="Bone/connective tissue"]<-paste0("CBone/connective tissue")
	plot_dat$Category[plot_dat$Category=="Bleeding disorders"]<-paste0("FBleeding disorders")
	plot_dat$Category[plot_dat$Category=="Other"]<-paste0("GOther")


	plot_dat$lci<-exp(plot_dat$b-1.96*plot_dat$se)
	plot_dat$uci<-exp(plot_dat$b+1.96*plot_dat$se)
	plot_dat$or<-exp(plot_dat$b)
	List<-lapply(c("or","lci","uci"),FUN=function(x)
	round(plot_dat[,x],2))
	
	plot_dat$or<-unlist(List[1])
	plot_dat$lci<-unlist(List[2])
	plot_dat$uci<-unlist(List[3])
	return(plot_dat)
}




# # ,weight=mr.ao$weight
# forestplot(
# df = df_compare_traits_groups,
# estimate = beta,
# pvalue = pvalue,
# psignif = 0.002,
# xlab = "1-SD increment in cardiometabolic trait\nper 1-SD increment in biomarker concentration",
# colour = trait
# ) +
# ggforce::facet_col(
# facets = ~group,
# scales = "free_y",
# space = "free"
# )
# forest_plot_1_to_many(mr.ao1,b="b",se="se",
#     exponentiate=T,ao_slc=F,lo=0.87,up=1.3,
#     traitM="outcome",by="Category",
#     trans="log2",xlab="",
#     col1_width=1.50, addcol_titles=c("",""),
#     subheading_size=7,
#     addcols=c("N_case","pval"),
#     addcol_widths=c(0.35,0.50),
#     shape_points=15,
#     colour_scheme="black")
# # dev.off()

# # Gastrointestinal
# # "Gastro-oesophageal reflux"
# # "Chronic obstructive airways disease"
# # grep("Pulmonary embolism" ,ao1$trait2)
# # Pos<-which(ao1$trait2 == "Wrist")
# # trait2<-ao1$trait[Pos]
# # trait2<-unlist(strsplit(trait2,split=": Yes"))
# # ao1$trait2[Pos]<-trait2
# # Temp<-ao1[which(ao1$Category=="NA"),]
# # Temp$trait2



# # ao1$Category[ao1$Category %in% c("Eye","Kidney","Kindey","Lung")]<-"Other"
# # res.info<-res.info[order(res.info$b_pQTL,decreasing=T),]
# # res.info<-res.info[order(res.info$Category),]
# # res.info$b_pQTL<-res.info$b_pQTL*-1 



# # res.info$pval<-formatC(res.info$p_pQTL, format = "e", digits = 2)
# # res.info1<-res.info[1:50,]
# # res.info2<-res.info[51:100,]


# pdf("disease1.pdf")
# forest_plot_1_to_many(res.info1,b="b_pQTL",se="se_pQTL",
#     exponentiate=T,ao_slc=F,lo=0.35,up=8,
#     trait2M="trait2",by=NULL,
#     trans="log2",xlab="",
#     col1_width=1.50, addcol_titles=c("",""),
#     addcols=c("N_case","pval"),addcol_widths=c(0.35,0.50))
# dev.off()

# # res.info<-res.info[order(res.info$Category),]

# res.other<-res.info[res.info$Category=="Other",]
# res.bone<-res.info[res.info$Category=="Bone/connective tissue",]
# res.info2<-res.info[!res.info$Category %in% c("Other","Bone/connective tissue"),]
# res.info2<-rbind(res.info2,res.bone)
# res.info2<-rbind(res.info2,res.other)
# res.info2[,c("trait2","Category")]

# table(res.info2$Category)
# res.infoA<-res.info2[1:49,]
# res.infoB<-res.info2[50:99,]
# pdf("diseaseA_cat.pdf",height=8,width=7.5)
# forest_plot_1_to_many(res.infoA,b="b_pQTL",se="se_pQTL",
#     exponentiate=T,ao_slc=F,lo=0.10,up=2.6,
#     trait2M="trait2",by="Category",
#     trans="log2",xlab="",
#     col1_width=1.50, addcol_titles=c("",""),
#     subheading_size=7,
#     addcols=c("N_case","pval"),
#     addcol_widths=c(0.35,0.50),
#     shape_points=15,
#     colour_scheme="black")
# dev.off()

# 	Controls<-res.info[order(as.numeric(res.info$N_control),decreasing=T),c("Consortium","N_case","N_control")]
# 	Controls<-Controls[!duplicated(Controls$Consortium),]
# 	sum(as.numeric(Controls$N_control))
# 	sum(as.numeric(res.info$N_case))
# 	Cases<-res.info[order(as.numeric(res.info$N_case),decreasing=T),c("Consortium","N_case","N_control")]
# 	Cases$Consortium[is.na(Cases$Consortium)]<-"unknown"
# 	Cases1<-Cases[Cases$Consortium!="Neale's lab",]
# 	Cases1<-Cases1[!duplicated(Cases1$Consortium),]
# 	Cases2<-Cases[Cases$Consortium=="Neale's lab",]
# 	sum(as.numeric(Cases1$N_case))+sum(as.numeric(Cases2$N_case))

# 	sum(as.numeric(Controls$N_control))+
# 	sum(as.numeric(Cases1$N_case))+sum(as.numeric(Cases2$N_case))

# 	info<-read.table("outcome-info.txt",sep="\t",head=T,stringsAsFactors=F,fill=T,quote="",comment.char="",as.is=TRUE)
# 	info$disease<-"notdisease"
# 	info$disease[info$trait2_name %in% c("Diagnoses - main ICD10: G56 Mononeuropathies of upper limb",
# 		"Diagnoses - main ICD10: I48 Atrial fibrillation and flutter",
# 		"Diagnoses - main ICD10: I80 Phlebitis and thrombophlebitis",
# 		"Diagnoses - main ICD10: I83 Varicose veins of lower extremities",
# 		"Diagnoses - main ICD10: K57 Diverticular disease of intestine",
# 		"Diagnoses - main ICD10: K80 Cholelithiasis",
# 		"Non-cancer illness code  self-reported: psoriasis",
# 		"Non-cancer illness code  self-reported: osteoarthritis",
# 		"Non-cancer illness code  self-reported: gout",
# 		"Non-cancer illness code  self-reported: pneumothorax",
# 		"Non-cancer illness code  self-reported: polio / poliomyelitis",
# 		"Non-cancer illness code  self-reported: arthritis (nos)",
# 		"Non-cancer illness code  self-reported: hypertrophic cardiomyopathy (hcm / hocm)",
# 		"Non-cancer illness code  self-reported: vitiligo",
# 		"Underlying (primary) cause of death: ICD10: E85.4 Organ-limited amyloidosis",
# 		"Eye problems/disorders: Diabetes related eye disease",
# 		"Eye problems/disorders: Glaucoma",
# 		"Vascular/heart problems diagnosed by doctor: Angina",
# 		"Diagnoses - main ICD10: C61 Malignant neoplasm of prostate",
# 		"Diagnoses - main ICD10: D12 Benign neoplasm of colon  rectum  anus and anal canal",
# 		"Diagnoses - main ICD10: I30 Acute pericarditis",
# 		"Diagnoses - main ICD10: K20 Oesophagitis",
# 		"Diagnoses - main ICD10: K29 Gastritis and duodenitis",
# 		"Diagnoses - main ICD10: K35 Acute appendicitis",
# 		"Cancer code  self-reported: small intestine/small bowel cancer",
# 		"Diagnoses - main ICD10: M16 Coxarthrosis [arthrosis of hip]",
# 		"Diagnoses - main ICD10: M17 Gonarthrosis [arthrosis of knee]",
# 		"Diagnoses - main ICD10: M54 Dorsalgia",
# 		"Prostate cancer",
# 		"Diagnoses - main ICD10: D25 Leiomyoma of uterus",
# 		"Diagnoses - main ICD10: G47 Sleep disorders",
# 		"Diagnoses - main ICD10: M72 Fibroblastic disorders",
# 		"Diagnoses - main ICD10: N20 Calculus of kidney and ureter",
# 		"Diagnoses - main ICD10: N81 Female genital prolapse",
# 		"Diagnoses - main ICD10: R04 Haemorrhage from respiratory passages",
# 		"Diagnoses - main ICD10: R07 Pain in throat and chest",
# 		"Diagnoses - main ICD10: R10 Abdominal and pelvic pain",
# 		"Diagnoses - main ICD10: R11 Nausea and vomiting",
# 		"Diagnoses - main ICD10: R14 Flatulence and related conditions",
# 		"Diagnoses - main ICD10: R35 Polyuria",
# 		"Diagnoses - main ICD10: R55 Syncope and collapse",
# 		"Cancer code  self-reported: basal cell carcinoma",
# 		"Non-cancer illness code  self-reported: hypertension",
# 		"Non-cancer illness code  self-reported: pulmonary embolism +/- dvt",
# 		"Non-cancer illness code  self-reported: deep venous thrombosis (dvt)",
# 		"Non-cancer illness code  self-reported: asthma",
# 		"Non-cancer illness code  self-reported: chronic obstructive airways disease/copd",
# 		"Non-cancer illness code  self-reported: emphysema/chronic bronchitis",
# 		"Non-cancer illness code  self-reported: sleep apnoea",
# 		"Non-cancer illness code  self-reported: gastro-oesophageal reflux (gord) / gastric reflux",
# 		"Non-cancer illness code  self-reported: kidney stone/ureter stone/bladder stone",
# 		"Non-cancer illness code  self-reported: bladder problem (not cancer)",
# 		"Non-cancer illness code  self-reported: hyperthyroidism/thyrotoxicosis",
# 		"Non-cancer illness code  self-reported: hypothyroidism/myxoedema",
# 		"Non-cancer illness code  self-reported: migraine",
# 		"Non-cancer illness code  self-reported: depression",
# 		"Non-cancer illness code  self-reported: anxiety/panic attacks",
# 		"Non-cancer illness code  self-reported: mania/bipolar disorder/manic depression",
# 		"Non-cancer illness code  self-reported: bone disorder",
# 		"Non-cancer illness code  self-reported: joint disorder",
# 		"Non-cancer illness code  self-reported: osteoporosis",
# 		"Non-cancer illness code  self-reported: ankylosing spondylitis",
# 		"Non-cancer illness code  self-reported: iron deficiency anaemia",
# 		"Non-cancer illness code  self-reported: pernicious anaemia")]<-"disease"

# 	info$disease[grep("illness",info$trait2_name)]<-"disease"
# 	info$disease[grep("Diagnoses",info$trait2_name)]<-"disease"
# 	info$disease[grep("disorders",info$trait2_name)]<-"disease"

# 	# phewas
# 	phewas<-read.csv("rs11581607.csv",head=T,stringsAsFactors=F)
# 	# phewas<-phewas[!duplicated(paste(phewas$trait2,phewas$Beta,phewas$SE)),]
# 	phewas$trait2[phewas$Pval<0.05/10000]

# 	il23<-data.frame(matrix(c("Interleukin-23 receptor","Sun et al","M & F",3301,29875488,"A","G",-0.42	,1.2E-17,0.05,"",""),nrow=1,ncol=12),stringsAsFactors=F)
# 	names(il23)<-names(phewas)
# 	Dat<-rbind(il23,phewas)
# 	Dat<-Dat[order(Dat$trait2),]
# 	ID<-substr(Dat$trait2,1,1)
# 	# ao<-available_outcomes()
# 	ao$trait2<-tolower(ao$trait2)
# 	Dat$trait2<-tolower(Dat$trait2)
# 	Dat<-merge(Dat,ao[!duplicated(ao$trait2),c("trait2","category","subcategory")],by.x="trait2",by.y="trait2",all.x=T)
# 	Dat$subcategory[which(Dat$subcategory=="NA") ]<-NA
# 	Dat$subcategory[which(Dat$subcategory=="") ]<-NA
# 	ao$subcategory[ao$subcategory=="NA"]<-NA
# 	ao$subcategory[ao$subcategory==""]<-NA
# 	ao1<-ao[!is.na(ao$subcategory),]
# 	miss<-Dat[is.na(Dat$subcategory),]
# 	Dat<-Dat[!is.na(Dat$subcategory),]
# 	miss<-merge(miss,ao1[,c("trait2","subcategory")],by.x="trait2",by.y="trait2",all.x=T)
# 	names(miss)[names(miss)=="subcategory.y"]<-"subcategory"
# 	miss<-miss[,names(miss)!="subcategory.x"]
# 	Dat<-rbind(miss,Dat)
# 	Dat1<-Dat[!is.na(Dat$subcategory),]
# 	Dat2<-Dat[is.na(Dat$subcategory),]
# 	for(i in 1:length(Dat1$trait2)){
# 		print(Dat1$trait2[i]) 
# 		Pos<-grep(Dat1$trait2[i],Dat2$trait2)
# 		Dat2$subcategory[Pos]<-Dat1$subcategory[i]
# 	}

# 	Dat1<-rbind(Dat1,Dat2[!is.na(Dat2$subcategory),])
# 	Dat2<-Dat2[is.na(Dat2$subcategory),]

# 	Dat2$subcategory2<-NA
# 	for(i in 1:length(Dat1$trait2)){
# 		print(i)
# 		print(Dat1$trait2[i])
# 		trait2<-unlist(strsplit(Dat1$trait2[i],split=" "))
# 		trait2<-trait2[1]
# 		# print(unique(Dat2$trait2[Pos]))
# 		print(Dat1$subcategory[i])
# 		Pos<-grep(trait2,Dat2$trait2)
# 		Dat2$subcategory2[Pos]<-Dat1$subcategory[i] 
# 	}

# 	Dat2<-Dat2[,names(Dat2)!="subcategory"]
# 	names(Dat2)[names(Dat2)=="subcategory2"]<-"subcategory"

# 	Dat3<-rbind(Dat1,Dat2)
# 	Dat3$subcategory[is.na(Dat3$subcategory)]<-"Other"

# 	# Dat3[Dat3$trait2=="Interleukin-23 receptor",]
# 	library(ggplot2)
# 	library(ggrepel)
# 	Dat3$Pval<-as.numeric(Dat3$Pval)
# 	Dat3$trait2[Dat3$Pval<0.05/22000]
# 	Dat3<-Dat3[!Dat3$trait2 %in% c("diagnoses - main icd10: k50 crohn's disease [regional enteritis]","diagnoses - main icd10: k51 ulcerative colitis","non-cancer illness code  self-reported: crohns disease","non-cancer illness code  self-reported: psoriasis"),]
# 	Dat3$trait2[Dat3$trait2 %in% c( "diagnoses - main icd10: k51.9 ulcerative colitis, unspecified" ,"diagnoses - secondary icd10: k50.9 crohn's disease, unspecified","non-cancer illness code, self-reported: ulcerative colitis","non-cancer illness code, self-reported: crohns disease","non-cancer illness code, self-reported: psoriasis")]<-c("ICD10: K51.9 ulcerative colitis" ,"ICD10: K50.9 Crohn's disease","Self-reported  ulcerative colitis","Self-reported crohns disease","Self-reported psoriasis")
	         

# 	Dat3<-Dat3[order(Dat3$subcategory),]
# 	Dat3$ID<-1:nrow(Dat3)
# 	# Dat3[1:100,c("subcategory","ID")]
# 	Label<-Dat3$trait2
# 	Label[as.numeric(Dat3$Pval)>0.05/22000]<-""
# 	Label<-gsub("(^[[:alpha:]])", "\\U\\1", Label, perl=TRUE)
# 	Label[Label %in% c("Crohn's disease","Inflammatory bowel disease","Ulcerative colitis") ]<-c("Crohn's disease (IIBDGC)","Inflammatory bowel disease (IIBDGC)","Ulcerative colitis (IIBDGC)") 
# 	Label[Label %in% c("ICD10: K51.9 ulcerative colitis","ICD10: K50.9 Crohn's disease","Self-reported  ulcerative colitis","Self-reported crohns disease","Self-reported psoriasis") ]<-c("ICD10: K51.9 ulcerative colitis (UK Biobank)","ICD10: K50.9 Crohn's disease (UK Biobank)","Self-reported ulcerative colitis (UK Biobank)","Self-reported Crohns disease (UK Biobank)","Self-reported psoriasis (UK Biobank)") 
# 	Label[Label=="Interleukin-23 receptor"]<-"Interleukin-23 receptor (Sun et al)"
# 	Dat3$subcategory2<-as.factor(Dat3$subcategory)
# 	# Dat3[Dat3$Pval<0.05/20000,c("trait2","Consortium")]
# 	Dat3$Pval<--log10(Dat3$Pval)
# 	# table(Dat3$Consortium)


# 	}
# }
