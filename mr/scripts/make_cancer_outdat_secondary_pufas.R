source("~/fatty-acids/mr/scripts/mr_functions.R")

Dir<-"~/fatty-acids/outcome_data/secondary_pufas/"
Files<-dir(Dir)

#Files<-Files[11:14]

Dat_list<-NULL
for(i in 1:length(Files)){
	print(i)
	#i<-c(11,12,13,14)
	#i<-11
	print(Files[i])
	# i<-which(Files== "gli66_csi_inst.txt")
	In_file<-paste0(Dir,Files[i])
	Dat<-read.table(In_file,sep="\t",stringsAsFactors=F,head=T,quote="")	
	# Dat<-read.table("~/fatty-acids/mr/data/csi_inst/cll83_csi_inst.txt",sep="\t",stringsAsFactors=FALSE,head=TRUE,quote="")
	# if(any(Dat$ID == 137)) stop("")
	
	if(!"id" %in% names(Dat)) stop("missing id")
	if(!"rsid" %in% names(Dat)) stop("missing rsid")
	if(any(duplicated(Dat$rsid))) stop("dup SNPs")
	Dat$file.outcome<-Files[i]
	Dat_list[[i]]<-Dat
}

Dat<-do.call(plyr::rbind.fill,Dat_list)

# Dat[Dat$study == "ILCCO" & Dat$rsid =="rs1741" ,]
col_names_keep<-c("file.outcome","outcome","population","pmid","study","ncase","ncontrol","UKbiobank","rsid","effect_allele","other_allele","lnor","lnor_se","eaf","p","info1","HWEp","phet","I2","Q","Direction","effect_allele_confirmed","id","z_score","open_gwas","efo","proxy","study")

# Dat$open_gwas[is.na(Dat$open_gwas)]<-FALSE
Dat2<-Dat[,col_names_keep]
if(any(is.na(Dat2$id))) stop("id missing in some studies")

if(5 %in% Dat2$id) stop("ID 5 should not be present")

outcome_dat<-format_outcomes4(dat=Dat2)
outcome_dat$proxy[is.na(outcome_dat$proxy)]<-FALSE
outcome_dat$study
dim(outcome_dat[outcome_dat$outcome=="Colorectal cancer | 60" ,])

save(outcome_dat,file="~/fatty-acids/mr/data/outcome_dat_secondary_pufas.Rdata")



