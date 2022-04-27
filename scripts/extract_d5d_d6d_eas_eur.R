source("~/fatty_acids/colocalisation/scripts/extract_snps_functions.R")
Files<-c("GLA_to_LA_adjSNP.tab","GLA_to_LA.tab","AA_to_DGLA.tab" )
Files<-paste0("/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/ratios/",Files)
snp_dat1<-lapply(1:length(Files),FUN=function(i) extract_data2(snplist=c("rs968567","rs174546"),File=Files[i],exact_match=TRUE))

snp_dat1<-do.call(rbind,snp_dat1)


Files<-c("AA_to_DGLA_pooled.tab.gz","GLA_to_LA_pooled.tab.gz")
Files<-paste0("/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_pc/ratios/",Files)

snp_dat2<-lapply(1:length(Files),FUN=function(i) extract_data2(snplist=c("rs968567","rs174546"),File=Files[i],exact_match=TRUE))

snp_dat2<-do.call(rbind,snp_dat2)

snp_dat<-rbind(snp_dat1,snp_dat2)

snp_dat<-format_dat()

format_dat<-function(){
	snp_dat$exposure[snp_dat$file %in% c("GLA_to_LA_pooled.tab.gz","GLA_to_LA.tab")]<-"GLA:LA"
	snp_dat$exposure[snp_dat$file == 	"GLA_to_LA_adjSNP.tab"]<-"GLA:LAadj_rs174546"

	snp_dat$exposure[snp_dat$file %in% c("AA_to_DGLA_pooled.tab.gz","AA_to_DGLA.tab")]<-"AA:DGLA"

	snp_dat$population<-"European"
	snp_dat$population[grep("pooled",snp_dat$file)]<-"East Asian"
	return(snp_dat)
}

snp_dat[snp_dat$SNP == "rs174546",c("exposure","population")]
write.table(snp_dat,"~/d5d_d6d_eas_eur_exposure_dat.txt",sep="\t",col.names=T,row.names=F,quote=F)

scp ph14916@bluecrystalp3.acrc.bris.ac.uk:~/d5d_d6d_eas_eur_exposure_dat.txt . 


