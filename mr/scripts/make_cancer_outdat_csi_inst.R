

Dir<-"~/fatty-acids/mr/data/csi_inst/"
Files<-dir(Dir)
Dat_list<-NULL
for(i in 1:length(Files)){
	print(i)
	print(Files[i])
	# i<-which(Files== "gli66_csi_inst.txt")
	In_file<-paste0(Dir,Files[i])
	Dat<-read.table(In_file,sep="\t",stringsAsFactors=F,head=T,quote="")	
	# Dat<-read.table("~/fatty-acids/mr/data/csi_inst/cll83_csi_inst.txt",sep="\t",stringsAsFactors=FALSE,head=TRUE,quote="")
	# if(any(Dat$ID == 137)) stop("")
	if(!"ID" %in% names(Dat)) stop("missing id")
	if(!"rsid" %in% names(Dat)) stop("missing rsid")
	if(any(duplicated(Dat$rsid))) stop("dup SNPs")
	Dat$file.outcome<-Files[i]
	Dat_list[[i]]<-Dat
}

Dat<-do.call(plyr::rbind.fill,Dat_list)
unique(Dat$file.outcome[Dat$ID == 25])
col_names_keep<-c("file.outcome","outcome","population","pmid","study","ncase","ncontrol","UKbiobank","rsid","effect_allele","other_allele","lnor","lnor_se","eaf","p","info1","HWEp","phet","I2","Q","Direction","effect_allele_confirmed","ID","z_score","open_gwas","efo")

# Dat$open_gwas[is.na(Dat$open_gwas)]<-FALSE

Dat2<-Dat[,col_names_keep]
if(any(is.na(Dat2$ID))) stop("id missing in some studies")

Dat2<-Dat2[!Dat2$ID %in% c(5),] #study exclusions
unique(Dat2$file.outcome[Dat2$ID == 137])

save(Dat2,file="~/fatty-acids/mr/data/cancer_csi_inst.Rdata")

head(Csi)
head(Dat2)

