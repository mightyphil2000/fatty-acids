load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
a<-read.table("~/fatty-acids/mr/data/Cancer_datasets - studies.tsv",sep="\t",head=TRUE,stringsAsFactors=FALSE,fill=TRUE,quote="")
# load("~/MR_FattyAcids/data/summary_data/meta_data.Rdata")
which(a$ID == "39")
head(a)

a<-a[!is.na(a$ID),]
dim(a)

head(a)

a<-a[,c("cancer","study.abbreviation","study","site","system", "cell","Cancer.Group","cases","controls","total","original_name","excluded..reason.", "ID","PMID.OpenGWAS","OpenGWAS","GWAS.catalog","correspondence","population","subpopulation")]

a<-a[a$excluded..reason.=="no",]
b<-c("39","138; 30","37","49","87","5","48","155","86","150; 44","161","9","147","98","132; 97")
ids_b<-trimws(unlist(strsplit(b,split=";")))
length(ids_b)
dim(a)
which(a=="39")
length(unique(a$cancer))

write.table(a,"~/fatty-acids/outcome_data/data/cancer_datasets_primary_mr_analyses.txt",sep="\t",col.names=TRUE,row.names=FALSE)
IDS<-disc.tab9$ID
length(IDS)
length(grep(";",IDS))
IDS1<-trimws(unlist(strsplit(IDS[grep(";",IDS)],split=";")))
length(IDS1)
IDS<-as.numeric(trimws(unlist(strsplit(IDS,split=";"))))
length(IDS)
length(unique(c(IDS,ids_b)))


a<-a[a$ID %in% IDS,]
ukb.dat<-a[a$study.abbreviation == "UKB",]
ukb<-a$ID[a$study.abbreviation == "UKB"]

ukb.dat[,c("cases","controls","total")]
names(disc.tab9)
Temp<-disc.tab9[grep("UKB",disc.tab9$study.abbreviation),c("cancer","ID","cases","controls","total","study.abbreviation","site","population","subpopulation","cell","pmid","study","system","original_name","Cancer.Group")]
disc.tab9[grep("Non-hodg",disc.tab9$cancer),c("cancer","ID","cases","controls","total","study.abbreviation","site","population","cell","pmid","study","system","original_name","Cancer.Group","subpopulation")]
dim(disc.tab9)
dim(Temp)
Temp$ukb.total.overlap<-114999/Temp$total
min(Temp$ukb.total.overlap)
max(Temp$ukb.total.overlap)
median(Temp$ukb.total.overlap)
Temp[order(Temp$cases),c("ID","ukb.total.overlap")]

# Brain cancer	748	468373	UKB; FinnGen	138; 30
#138: 606	372016

156	337003/(748+468373)

# Myeloid leukaemia	462	372016	UKB	155
114999/(462+372016)

# Lymphoid leukaemia	958	468317	UKB; FinnGen	150; 44
114999/(958+468317)

# Small bowel cancer	156	337003	UKB	161
114999/(156+337003)

# Liver & bile duct cancer		UKB	147
114999/(350+372016)



# for reproductive cancers divide number participants in PUFA GWAS by 2 (assuming equal number of males and females)
# Ovarian cancer	30869	387356	OCAC; OCAC (EAS); UKB; BJ; FinnGen	118; 120; 158; 18; 51
(114999/2)/(30869+387356)

# Prostate cancer	95512	378951	PRACTICAL; UKB; BJ; FinnGen	128; 159; 20; 53
(114999/2)/(95512+378951)

# Breast cancer	139445	398407	BCAC; UKB; FinnGen	6; 139; 31
(114999/2)/(139445+398407)
for(i in 1:nrow(Temp)){
	b<-Temp[i,]
	IDS<-b$ID
	IDS<-as.numeric(trimws(unlist(strsplit(IDS,split=";"))))
	Pos<-unlist(lapply(1:length(IDS),FUN=function(x)
		which(ukb.dat$ID == IDS[x])))
	case.ukb<-ukb.dat$case[Pos]
	con.ukb<-ukb.dat$control[Pos]
	tot.ukb<-ukb.dat$total[Pos]
	overlap<-round(tot.ukb/b$total,
	Min<-min(tot.ukb,114999)
	Min/b$total
	Min/b$total
	id<-b$ID
	c(overlap,id)

}
IDS[!IDS %in% a$ID]

