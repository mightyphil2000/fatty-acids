setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_instruments"
)

cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/

load("extract_sig_snps_pufas_opengwas.RData") #GRCh37
load("extract_sig_snps_pufas.RData")
 
Pos1<-grep("hg19",Dat$File)
Pos2<-grep("hg19",Dat$File,invert=TRUE)
Temp1<-Dat[Pos1,]
Temp2<-Dat[Pos2,]
Files<-unlist(strsplit(Temp1$File,split="hg19/"))

Files<-Files[seq(from=2,to=length(Files),by=2)]
Temp1$File<-Files

Files<-unlist(strsplit(Temp2$File,split="log/"))
Files<-Files[seq(from=2,to=length(Files),by=2)]
Temp2$File<-Files

Dat<-rbind(Temp1,Temp2)
Dat$id<-Dat$File
Dat$open_gwas<-FALSE
Dat2$open_gwas<-TRUE

Dat3<-plyr::rbind.fill(Dat,Dat2)
ids<-unique(Dat3$id)
length(ids)

clump_list<-NULL
for(i in 1:length(ids)){
for(i in 33:length(ids)){	
	print(ids[i])
	dat1<-Dat3[Dat3$id == ids[i],]
	Clump<-ieugwasr::ld_clump( clump_r2 = 0.001,clump_kb = 10000,
	    dplyr::tibble(rsid=dat1$snp, pval=dat1$p, id=dat1$id),
	    plink_bin = "/data/ph14916/.conda/envs/r4/lib/R/library/genetics.binaRies/bin/plink_Linux",
	    bfile = "/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K"
	)#GRCh37
	clump_list[[i]]<-Clump
}

ls /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bed
ls /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim
ls /projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.fam

clump_dat<-do.call(rbind,clump_list)
Dat4<-merge(Dat3,clump_dat,by.x=c("snp","id"),by.y=c("rsid","id"))
head(Dat4)
trait_dat<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/fatty_acid_GWASanalysis_table2.txt",sep="\t",head=TRUE,stringsAsFactors=FALSE,quote="")
trait_dat$filename<-gsub("\xca","",trait_dat$filename)
# grep("pooled_allchr_qc.txt",trait_dat$filename)
trait_dat$filename<-gsub("pooled_allchr_qc.txt","pooled_allchr_qc1.tab",trait_dat$filename)
Pos<-which(!is.na(trait_dat$filename) & trait_dat$filename !="")
trait_dat<-trait_dat[Pos,]
trait_dat<-trait_dat[,names(trait_dat) != "id"]
trait_dat$filename<-trimws(trait_dat$filename)
Dat5<-merge(Dat4,trait_dat[,c("filename","trait","pmid",	"author","consortium","trait","chain","chain.length","sample_size.analysis","population")],by.x="id",by.y="filename",all.x=TRUE)

Dat5$trait<-Dat5$trait.y
Dat5$trait[is.na(Dat5$trait)]<-Dat5$trait.x[is.na(Dat5$trait)]

get_positions_grch37<-function(){
	ukb_bed<-read.table("/projects/MRC-IEU/users/ph14916/ukb_bed/UKBB_10K.bim",sep="\t",head=FALSE,quote="") #GRCh37

	Dat6<-merge(Dat5,ukb_bed,by.x="snp",by.y="V2",all.x=TRUE)
	Dat6$chr[Dat6$snp=="rs143163547"]<-14
	Dat6$position[Dat6$snp=="rs143163547"]<-54616527
	Dat6$chr[is.na(Dat6$chr)]<-Dat6$V1[is.na(Dat6$chr)]
	Dat6$position[is.na(Dat6$position)]<-Dat6$V4[is.na(Dat6$position)]
	return(Dat6)
}
# Dat6[which(Dat6$trait=="arachidonic acid (20:4n6)"),]

Temp<-Dat6[Dat6$chr == 11,]
Dat7<-Dat6[Dat6$chr != 11,]
# we defined the FADS region in GRch37 coordinates:
lower<-61067099	
upper<-62096790
Temp$FADS<-TRUE
Temp$FADS[which(Temp$position < lower | Temp$position>upper)]<-FALSE
Dat7$FADS<-FALSE
Dat8<-rbind(Dat7,Temp)

Dat8$id[Dat8$trait=="Ratio of docosahexaenoic acid to total fatty acids" ]

Dat8<-harmonise_trait_names()

Dat8$n[grep("met-d",Dat8$id)]<-118466 
Dat8$consortium[grep("met-d",Dat8$id)]<-"UK Biobank"
Dat8$population[grep("met-d",Dat8$id)]<-"European"
Dat8$population[grep("met-c",Dat8$id)]<-"European"
Dat8$consortium[grep("met-c",Dat8$id)]<-"Kettunen"

Dat8$ID<-paste(Dat8$consortium,Dat8$trait)
r1<-with(Dat8, tapply(n, consortium, median))
Dat8$median_n<-NA
for(x in names(r1)){
	print(x)
	Dat8$median_n[Dat8$consortium==x]<-r1[names(r1)==x]
}
Dat8<-Dat8[order(Dat8$median_n,decreasing=TRUE),]
eur<-Dat8[Dat8$population=="European",]
eas<-Dat8[Dat8$population=="East Asian",]

met_dat.eur<-eur[!duplicated(eur$trait),c("consortium","trait","median_n","ID")]
met_dat.eas<-unique(eas[,c("consortium","trait","median_n","ID")])
eur1<-eur[eur$ID %in% met_dat.eur$ID, ]

readme<-"eur contains all snp level pufa data; eur1 contains snp level PUFA data but restricted to the single largest study for each trait (when trait is available in more than one study). Dat contains all snp level data for Europeans and East Asians. eas is snp level data for East Asians (all from the SCHS); met_dat describes single largest study available for each trait"

Dat<-Dat8
save(list=c("eur","eur1","eas","Dat","met_dat.eur","met_dat.eas","readme"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_instruments/instruments.Rdata")
# load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/secondary_pufa_instruments/instruments.Rdata")



harmonise_trait_names<-function(){
	Dat$trait[Dat$trait=="Ratio of omega-6 fatty acids to total fatty acids" ]<-"Omega-6 fatty acids"
	Dat$trait[Dat$trait=="Ratio of docosahexaenoic acid to total fatty acids"]<-"docosahexaenoic acid (22:6n3)"
	Dat$trait[Dat$trait=="Ratio of linoleic acid to total fatty acids"]<- "linoleic acid (18:2n6)"   
	Dat$trait[Dat$trait=="Ratio of omega-3 fatty acids to total fatty acids"  ]<- "Omega-3 fatty acids"   
	Dat$trait[Dat$trait=="22:6, docosahexaenoic acid" ]<-  "docosahexaenoic acid (22:6n3)"
	Dat$trait[Dat$trait=="\" X-12442--5,8-tetradecadienoate\""]<-  "tetradecadienoic acid (14:2n9)"
	Dat$trait[Dat$trait=="18:2, linoleic acid (LA)" ]<-  "linoleic acid (18:2n6)"
	Dat$trait<-paste(toupper(substr(Dat$trait, 1, 1)), substr(Dat$trait, 2, nchar(Dat$trait)), sep="")
	return(Dat)
}

	