load("~/fatty-acids/mr/data/instruments.Rdata")

source("~/fatty-acids/mr/scripts/mr_functions.R")

exp_dat<-format_exposure3(dat=eur1,standardise_beta=TRUE,beta="beta",se="se",pval="pval",effect_allele="effect_allele",other_allele="other_allele",eaf="effect_allele_freq",rsid="snp",ID="ID",exposure="trait",samplesize="n")

exposures<-unique(exp_dat$exposure)
length(unique(exp_dat$SNP))
length(unique(exp_dat$exposure))
excl<-c("Stearidonic acid (18:4n3)","Eicosadienoic acid (20:2n6)",
"Adrenic acid (22:4n6)","Docosapentaenoic acid (22:5n6)","Tetradecadienoic acid (14:2n9)","Tetradecadienoic acid (14:2n9)","Dihomo-linolenic acid (20:3n3 or n6)","Other polyunsaturated fatty acids than 18:2","Ratio of omega-6 fatty acids to omega-3 fatty acids")

# N snps associated with 13 seconday PUFAs (that can be instrumented by variation outside and witihn the FADS region)
length(unique(exp_dat$SNP[!exp_dat$exposure %in% excl & !exp_dat$FADS]))

unique(exp_dat$exposure[!exp_dat$exposure %in% excl])

# tabulate numbers for PUFAs that could not be clearly classified as omega 3 or omega 6 PUFAs
other<-c("Tetradecadienoic acid (14:2n9)","Dihomo-linolenic acid (20:3n3 or n6)","Other polyunsaturated fatty acids than 18:2", "Ratio of omega-6 fatty acids to omega-3 fatty acids")
table(exp_dat$exposure[!exp_dat$exposure %in% other])
unique(exp_dat$exposure[!exp_dat$exposure %in% other])
length(unique(exp_dat$SNP[!exp_dat$exposure %in% other]))
length(unique(exp_dat$SNP[!exp_dat$exposure %in% other & exp_dat$population == "European" ]))
length(unique(Dat$snp[!Dat$trait %in% other ]))
length(unique(Dat$trait[!Dat$trait %in% other ]))
length(unique(eas$snp[!eas$trait %in% other]))
length(unique(eas$trait[!eas$trait %in% other]))

eur1[eur1$trait=="Stearidonic acid (18:4n3)","chr"]

any(duplicated(paste(eur1$trait,eur1$snp)))
length(unique(eur1$snp))
length(unique(eas$snp))

# any(duplicated(exp_dat$exposure[exp_dat$FADS]))

r2_list_inclfads<-NULL
# exposure<-exposures[1]
for(exposure in exposures){
	print(exposure)
	test_dat<-exp_dat[exp_dat$exposure == exposure,]
	r2_list_inclfads[[exposure]]<-sum(r2(dat=test_dat))
}


r2_list_exclfads<-NULL
exp_dat_exclfads<-exp_dat[!exp_dat$FADS,]

exposures<-exp_dat_exclfads$exposure
# exposure<-exposures[1]
r2_list_exclfads<-NULL
for(exposure in exposures){
	print(exposure)
	test_dat<-exp_dat_exclfads[exp_dat_exclfads$exposure == exposure,]
	r2_list_exclfads[[exposure]]<-sum(r2(dat=test_dat))
}


# load("~/fatty-acids/mr/results/res_single_ara_la.Rdata")
# load("~/fatty-acids/mr/results/res_single_ara_la_fixed_effects_v2.Rdata")
load("~/fatty-acids/mr/results/mr_secondary_pufas_inclfads.Rdata")
fads1<-res4
fads1$fads<-TRUE
load("~/fatty-acids/mr/results/mr_secondary_pufas_exclfads.Rdata")
fads2<-res4
fads2$fads<-FALSE
Dat<-rbind(fads1,fads2)

Plot_dat<-format_trait(dat=eur1)
Plot_dat<-Plot_dat[order(Plot_dat$PUFA,Plot_dat$chain.length),]

Dups<-unique(Plot_dat$trait[duplicated(Plot_dat$trait)])
Plot_dat<-Plot_dat[Plot_dat$trait %in% Dups,]

n3<-unique(Plot_dat$trait[Plot_dat$PUFA=="Omega 3"])
n6<-unique(Plot_dat$trait[Plot_dat$PUFA=="Omega 6"])

exposure<-names(r2_list_exclfads) 
r2<-as.numeric(unlist(r2_list_exclfads))
excl_fads<-data.frame((cbind(exposure,r2)))
excl_fads$FADS<-"exclude"
exposure<-names(r2_list_inclfads)
r2<-as.numeric(unlist(r2_list_inclfads))
incl_fads<-data.frame(cbind(exposure,r2))
incl_fads$FADS<-"include"
fads_r2<-rbind(incl_fads,excl_fads)

exposure<-names(table(eur1$trait))
nsnps<-as.numeric(table(eur1$trait))
nsnps_tab1<-data.frame(cbind(exposure,nsnps))
nsnps_tab1$FADS<-"include"
exposure<-names(table(eur1$trait[!eur1$FADS]))
nsnps<-as.numeric(table(eur1$trait[!eur1$FADS]))
nsnps_tab2<-data.frame(cbind(exposure,nsnps))
nsnps_tab2$FADS<-"exclude"
nsnps_tab<-rbind(nsnps_tab1,nsnps_tab2)
fads_r2<-merge(fads_r2,nsnps_tab,by=c("exposure","FADS"),all.x=TRUE)

fads_r2<-format_pufa(dat=fads_r2)
fads_r2<-fads_r2[order(fads_r2$PUFA,fads_r2$chain.length,fads_r2$exposure),]
length(unique(eur1$snp))
length(unique(eas$snp))
table(eas$trait)

Temp<-eur1[eur1$trait %in% c("Stearidonic acid (18:4n3)","Eicosadienoic acid (20:2n6)","Adrenic acid (22:4n6)","Docosapentaenoic acid (22:5n6)"),]
length(unique(Temp$snp))
eur1[eur1$trait == "Tetradecadienoic acid (14:2n9)",]
Temp[,c("chr","FADS")]
met_dat.eur
Temp<-eur[!duplicated(eur$ID),]
length(which(table(Temp$trait)==1))
table(Temp$trait,Temp$consortium)
write.table(fads_r2,"~/fatty-acids/mr/results/secondary_instruments_r2.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

r2<-function(dat=NULL){
	var<-1
	2*dat$b^2*dat$eaf.exposure*(1-dat$eaf.exposure)/var
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
	return(dat)
}


format_pufa<-function(dat=NULL){
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
	return(dat)
}


format_trait<-function(dat=NULL){

	dat$PUFA<-NA
	dat$PUFA[grep("n3",dat$trait)]<-"Omega 3"
	dat$PUFA[dat$trait=="Omega-3 fatty acids"]<-"Omega 3"
	dat$PUFA[grep("n6",dat$trait)]<-"Omega 6"
	dat$PUFA[dat$trait=="Omega-6 fatty acids"]<-"Omega 6"
	dat$PUFA[grep("n3 or n6",dat$trait)]<-"other"
	dat$PUFA[is.na(dat$PUFA)] <- "other"

	dat1<-dat[grep(":",dat$trait),]
	dat1<-dat1[dat1$trait!= "Other polyunsaturated fatty acids than 18:2",]
	dat2<-dat[grep(":",dat$trait,invert=TRUE),]
	dat2<-rbind(dat2,dat[dat$trait== "Other polyunsaturated fatty acids than 18:2",])
	Pos<-gregexpr(":",dat1$trait)
	dat1$chain.length<-as.numeric(substr(dat1$trait,start=unlist(Pos)-2,stop=unlist(Pos)-1))
	dat2$chain.length<-NA
	dat<-rbind(dat1,dat2)
	# dat<-dat[dat$population =="European",]
	# dat<-dat[dat$trait=="Arachidonic acid",]
	dat$b<-as.numeric(dat$b)
	dat$se<-as.numeric(dat$se)
	
	return(dat)
}
