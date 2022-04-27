library(MRInstruments)
library(TwoSampleMR)
library(ieugwasr)
# ao<-available_outcomes()
ao<-ieugwasr::gwasinfo()

# IDS<-ao$id[which(ao$consortium == "TAG")]
# IDS<-ao$id[which(ao$subcategory %in% c("Autoimmune / inflammatory","Immune system"))]
# IDS<-ao$id[

# Pos<-unlist(lapply(c("inflammatory bowel disease","Crohn's disease", "Crohn’s disease" ,"Crohns disease" , "ulcerative colitis","asthma","esophagitis","eczema","atopic dermatitis","Chronic obstructive pulmonary disease","asbestosis","silicosis","bronchitis","cystitis","bladder inflammation","gingivitis","lichen planus","lichen sclerosus","pancreatitis","Barett’s esophagus","Baretts esophagus","sialadentis","Sjogren","thyroiditis","skin inflammation","psoriasis","hives"),FUN=function(x) grep(x,ao$trait,ignore.case=TRUE)))

# Identify genetic associations for inflammatory conditions of the bowel, lung, esophagus and skin
 # bowel - inflammatory bowel disease, crohn's disease, ulcerative colitis
 # lung - COPD [bronchitis, emphysema], asthma
 # esophagus -Esophagitis
 # skin - hives, psoriasis , eczema 
Pos<-unlist(lapply(c("inflammatory bowel disease","Crohn's disease", "Crohn’s disease" ,"Crohns disease" , "ulcerative colitis","asthma","esophagitis","eczema","atopic dermatitis","Chronic obstructive pulmonary disease","COPD","bronchitis","Barett’s esophagus","Baretts esophagus","skin inflammation","psoriasis","hives"),FUN=function(x) grep(x,ao$trait,ignore.case=TRUE)))

Pos<-unique(Pos)
unique(ao$trait[ao$consortium=="IIBDGC"])
IDS<-ao$id[Pos]
# 124
ao1<-ao[ao$id %in% IDS,]
# ao1[order(ao1$sample_size),]
# ao1$ncase[is.na(ao1$ncase)]<-ao1$sample_size[is.na(ao1$ncase)]
# ao1$trait[which(ao1$ncase== ao1$sample_size)]
# ao1<-ao1[order(ao1$ncase,decreasing=T),]
# ao1<-ao1[!duplicated(ao1$trait),]
IDS<-ao1$id


# proxies<-read.table("~/MR_FattyAcids/data/rs174546_r2proxies_ldlink_15jan21.txt",sep="\t",head=T,stringsAsFactors=F)
# length(unique(proxies$RS_Number[proxies$R2>=0.8]))

load("~/fatty-acids/outcome_data/data/define_fatty_acid_SNPs_v3.rdata")
SNP.proxies.EUR[which(SNP.proxies.EUR$ProxySNP == "rs174546"),]
proxies<-unique(fa.reg2[which(fa.reg2$SNP=="rs174546"),c("SNP","ProxySNP","R2","pop")])
proxies<-proxies[proxies$pop=="EUR",]

Dat<-associations(variants = proxies$ProxySNP,id = IDS)

Datp<-merge(Dat,proxies,by.x="rsid",by.y="ProxySNP")
Dat1<-Datp[Datp$rsid == "rs174546",]
Dat2<-Datp[Datp$rsid != "rs174546",]
Dat2<-Dat2[order(Dat2$R2,decreasing=TRUE),]
Datp<-rbind(Dat1,Dat2)
Dat3<-Datp[!duplicated(Datp$id),]
Dat.ao1<-merge(Dat3,ao1[,names(ao1)!="trait"],by="id")
Dat.ao1$trait2<-NA
Dat1<-Dat.ao1[grep(":",Dat.ao1$trait),]
Traits<-unlist(strsplit(Dat1$trait,split=":"))
Dat1$trait<-trimws(Traits[seq(2,length(Traits),by=2)])
Dat2<-Dat.ao1[grep(":",Dat.ao1$trait,invert=T),]
Dat1<-rbind(Dat1,Dat2)


Dat1$trait2[grep("Inflammatory bowel disease",Dat1$trait,ignore.case=T)]<-"Inflammatory bowel disease"
UC<-Dat1$trait[grep("Ulcerative colitis",Dat1$trait,ignore.case=T)]
UC<-UC[UC != "colitis/not crohns or ulcerative colitis"]
Dat1$trait2[Dat1$trait %in% UC]<-"Ulcerative colitis"
CD<-Dat1$trait[grep("crohn",Dat1$trait,ignore.case=T)]
CD<-CD[CD != "colitis/not crohns or ulcerative colitis" ]
Dat1$trait2[Dat1$trait %in% CD]<-"Crohn's disease"
Dat1$trait2[Dat1$trait=="colitis/not crohns or ulcerative colitis"]<-"Colitis/not crohns or ulcerative colitis"
COPD<-Dat1$trait[unique(unlist(lapply(c("Chronic obstructive pulmonary disease","COPD"),FUN=function(x) grep(x,Dat1$trait,ignore.case=T))))]
COPD<-COPD[COPD != "COPD-associated co-morbidities"]
Dat1$trait2[Dat1$trait %in% COPD]<-"Chronic obstructive pulmonary disease"


# Dat1$trait2[grep("Emphysema/chronic bronchitis",Dat1$trait,ignore.case=T)]<-"Emphysema/chronic bronchitis"

Dat1$trait2[unique(unlist(lapply(c("Chronic bronchitis/emphysema"  ,"Emphysema/chronic bronchitis"),FUN=function(x) grep(x,Dat1$trait,ignore.case=T))))]<-"Emphysema/chronic bronchitis"

Dat1$trait2[grep("Psoriasis",Dat1$trait,ignore.case=T)]<-"Psoriasis"
Dat1$trait2[Dat1$trait %in% c("Bronchitis","bronchitis","Acute bronchitis","Doctor diagnosed chronic bronchitis","Bronchitis, not specified as acute or chronic","Simple and mucoplurulent chronic bronchitis","Unspecified chronic bronchitis")]<-"Bronchitis"
# Dat1$trait2[Dat1$trait== "K80.1 Calculus of gallbladder with other cholecystitis"]<-"Cholecystitis"
# Dat1$trait2[Dat1$trait== "K80.0 Calculus of gallbladder with acute cholecystitis"]<-"Cholecystitis"
# Dat1$trait2[Dat1$trait=="Cystitis"]<-"Cystitis"
# Dat1$trait2[Dat1$trait=="N30 Cystitis" ]<-"Cystitis"
# Dat1$trait2[grep("K81 Cholecystitis",Dat1$trait,ignore.case=T)]<-"Cholecystitis"
# Dat1$trait2[Dat1$trait=="Chlocystitis" ]<-"Chlocystitis"
# Dat1$trait2[Dat1$trait=="Lichen planus" ]<-"Lichen planus"
# Dat1$trait2[Dat1$trait=="Thyroiditis, ILD-related definition" ]<-"Thyroiditis"
# Dat1$trait2[Dat1$trait=="Thyroiditis" ]<-"Thyroiditis"
# Dat1$trait2[Dat1$trait=="Subacute thyroiditis" ]<-"Thyroiditis"

# Dat1$trait2[Dat1$trait=="K85 Acute pancreatitis" ]<-"Pancreatitis"
# Dat1$trait2[Dat1$trait=="pancreatitis" ]<-"Pancreatitis"
# Dat1$trait2[Dat1$trait=="Acute pancreatitis"    ]<-"Pancreatitis"

# Dat1$trait2[Dat1$trait=="Chronic pancreatitis" ]<-"Pancreatitis"
# Dat1$trait2[Dat1$trait=="Alcohol-induced chronic pancreatitis" ]<-"Pancreatitis"
# Dat1$trait2[Dat1$trait=="Acohol-induced acute pancreatitis" ]<-"Pancreatitis"
# Dat1$trait2[Dat1$trait=="Sicca syndrome [Sjogren]" ]<-"Sjogren syndrome"
# Dat1$trait2[Dat1$trait=="Lichen sclerosus et atrophicus" ]<-"Lichen sclerosus"
# Dat1$trait2[Dat1$trait %in% c("Asthma associated co-morbidites") ]<-Dat1$trait[Dat1$trait %in% c("Asthma associated co-morbidites") ]

# Dat1$trait2[Dat1$trait %in% c("Excluding asthma or asthma-COPD overlap","Asthma-related infections","COPD/asthma related infections","Age asthma diagnosed","Medication related adverse effects (Asthma/COPD)","Asthma/COPD-related acute respiratory infections","Asthma-related acute respiratory infections","Psychiatric comorbidites (Asthma/COPD)","Medication related adverse effects (Asthma/COPD)") ]<-NA
AS<-Dat1$trait[grep("Asthma",Dat1$trait,ignore.case=T)]
AS<-AS[AS != "Allergic disease (asthma, hay fever or eczema)"]
Dat1$trait2[Dat1$trait %in% AS]<-"Asthma"
 
Eso<-Dat1$trait[grep("esophagitis",Dat1$trait,ignore.case=T)]
Eso<-Eso[Eso != "K21.9 Gastro-oesophageal reflux disease without oesophagitis" ]
Dat1$trait2[Dat1$trait %in% Eso]<-"Esophagitis"
EC<-Dat1$trait[unique(unlist(lapply(c("dermatitis","eczema"),FUN=function(x) grep(x,Dat1$trait,ignore.case=T))))]
EC<-EC[!EC %in% c("Allergic disease (asthma, hay fever or eczema)","Hayfever, allergic rhinitis or eczema","Hayfever  allergic rhinitis or eczema")]
Dat1$trait2[Dat1$trait %in% EC]<-"Eczema/dermatitis"



# Dat1$trait2[Dat1$trait=="Allergic disease (asthma, hay fever or eczema)" ]<-"Allergic disease (asthma, hay fever or eczema)"
# ao[which(ao$consortium=="IIBDGC"),c("trait","ncase","sample_size")]
Dat1<-Dat1[order(Dat1$ncase,decreasing=T),]
Dat1<-Dat1[!duplicated(Dat1$trait2),]
Dat1$trait[is.na(Dat1$trait2)]
Dat1<-Dat1[!is.na(Dat1$trait2),]
# Dat1[Dat1$trait == "Inflammatory bowel disease",c("ncase","sample_size","pmid","consortium")]


Dat2<-Dat1[Dat1$consortium=="MRC-IEU",]
Dat3<-Dat1[Dat1$consortium!="MRC-IEU",]
# Dat2[Dat2$trait2 == "Emphysema/chronic bronchitis","beta"]
Dat2<-transform_betas(dat=Dat2,effect="beta",effect.se="se")
Dat1<-rbind(Dat2,Dat3)

Dat1[,c("beta","se","trait2","ncase","ncontrol","pmid","ea","nea")]
Dat1$beta<-Dat1$beta*-1
ea<-Dat1$ea 
nea<-Dat1$nea 
eaf<-Dat1$eaf
Dat1$ea <- nea
Dat1$nea <- ea
Dat1$eaf<-1-eaf

library(ggforestplot)
library(ggplot2)

Dat1$name<-paste0(Dat1$trait2,"\nNo. cases=",Dat1$ncase)
Dat1<-Dat1[order(Dat1$ncase,decreasing=T),]
Dat1$organ <-NA
Dat1$organ[Dat1$trait2 %in% c("Asthma","Emphysema/chronic bronchitis","Bronchitis","Chronic obstructive pulmonary disease")]<-"Lung"
Dat1$organ[Dat1$trait2 %in% c("Inflammatory bowel disease","Crohn's disease","Ulcerative colitis","Esophagitis","Colitis/not crohns or ulcerative colitis")]<-"Digestive tract"
Dat1$organ[Dat1$trait2 %in% c("Psoriasis","Eczema/dermatitis")]<-"Skin"
Dat1<-Dat1[order(Dat1$organ),]
Dat1<-Dat1[Dat1$trait2 !="Colitis/not crohns or ulcerative colitis", ]

png("~/fatty-acids/mr/results/plots/inflammation_lookups_opengwas.png", width = 600, height = 800)
	forestplot(df = Dat1,logodds = TRUE,name=name,
	  estimate=beta, se=se,shape=NULL,
		colour = organ,xlab = "OR (95% CI) per copy of rs174546-C")
dev.off()

write.table(Dat1[,c("trait2","trait","organ","beta","se","p","ncase","ncontrol","ea","nea","eaf","consortium","pmid","id")],"~/fatty-acids/mr/results/inflammation_lookups_opengwas.txt",sep="\t",col.names=T,row.names=F,quote=F)

transform_betas<-function(dat=NULL,effect="lnor",effect.se="se"){
	# formula: log OR = beta / (u(1-u)); where u=ncases/(ncases + ncontrol) REPEAT with SE 	
	beta<-dat[,effect]
	se<-dat[,effect.se]
	u<-dat$ncase/(dat$ncase+dat$ncontrol)
	dat[,effect] <- beta / (u * (1 - u))
	dat[,effect.se]<-se / (u * (1 - u)) 	
	return(dat)
}






IDS<-ao1$id[ao1$category == "Disease"]
IDS<-ao$id[which(ao$subcategory == "TAG")]

Dat<-associations(variants = "rs2524299",id = IDS)
Dat<-associations(variants = "rs174546",id = IDS)
Dat1<-merge(Dat,ao1,by="id")

Dat1<-Dat1[order(Dat1$ncase,decreasing=T),]
Dat1<-Dat1[!duplicated(Dat1$trait.x),]

write.table(Dat1,"~/fatty-acids/mr/results/immune_disease_lookups.txt",sep="\t",col.names=T,row.names=F,quote=F)