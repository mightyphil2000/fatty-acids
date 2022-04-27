

format_results<-function(study=NULL,gtex=TRUE,drop_palindromic_SNPs=FALSE){
	Dat3<-rbind.fill(ref_dat,Dat2)
	if(!is.null(study)) Dat3$study[is.na(Dat3$study)]<-study
	if(gtex)
	{
		Dat4<-Dat3[!is.na(Dat3$ma_count),]
		Dat5<-Dat3[is.na(Dat3$ma_count),]
		Dat4<-Dat4[Dat4$ma_count>=10,]
		Dat<-rbind(Dat5,Dat4)
	}else{
		Dat<-Dat3
	}
	Dat<-Dat[,c("marker","effect_allele","other_allele","eaf","beta","se","study","trait","maf","trait2")]
	if(drop_palindromic_SNPs)
	{
		Alleles<-paste(Dat$effect_allele,Dat$other_allele,sep="")	
		Dat<-Dat[which(!Alleles %in% c("AT","TA","GC","CG")),]
	}
	return(Dat)
}

prep_data_eas2<-function(Dat=NULL){		
	gwis<-data.frame(Dat[1],stringsAsFactors=F)
	bbj<-data.frame(Dat[2],stringsAsFactors=F)
	# eqtl<-data.frame(Dat[3],stringsAsFactors=F)
	can<-data.frame(Dat[3],stringsAsFactors=F)
	bbj2<-effect_allele_bbj_data(dat=bbj)

	ref_dat<-can
	# gwis2<-gwis[gwis$trait !=  "AA:DGLA lnD5Dpooled",]

	# aa_info<-read.table("~/fatty-acids/colocalisation/data/aa_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	# dgla_info<-read.table("~/fatty-acids/colocalisation/data/dgla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	# gla_info<-read.table("~/fatty-acids/colocalisation/data/gla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	# la_info<-read.table("~/fatty-acids/colocalisation/data/la_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	
	# snps1<-snps_pass(Dat=aa_info)
	# snps2<-snps_pass(Dat=dgla_info)

	# Pos<-ref_dat$marker %in% snps1 & ref_dat$marker %in% snps2	
	# ref_dat<-ref_dat[which(Pos),]		
	
	Dat<-do.call(rbind.fill,list(gwis,bbj2))	
	return(list(Dat,ref_dat))
}

prep_data_eas<-function(Dat=NULL){		
	gwis<-data.frame(Dat[1],stringsAsFactors=F)
	bbj<-data.frame(Dat[2],stringsAsFactors=F)
	# eqtl<-data.frame(Dat[3],stringsAsFactors=F)
	can<-data.frame(Dat[3],stringsAsFactors=F)
	bbj2<-effect_allele_bbj_data(dat=bbj)

	ref_dat<-gwis[gwis$trait ==  "AA:DGLA lnD5Dpooled",]
	gwis2<-gwis[gwis$trait !=  "AA:DGLA lnD5Dpooled",]

	# aa_info<-read.table("~/fatty-acids/colocalisation/data/aa_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	# dgla_info<-read.table("~/fatty-acids/colocalisation/data/dgla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	# gla_info<-read.table("~/fatty-acids/colocalisation/data/gla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	# la_info<-read.table("~/fatty-acids/colocalisation/data/la_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	
	# snps1<-snps_pass(Dat=aa_info)
	# snps2<-snps_pass(Dat=dgla_info)

	# Pos<-ref_dat$marker %in% snps1 & ref_dat$marker %in% snps2	
	# ref_dat<-ref_dat[which(Pos),]		
	
	Dat<-do.call(rbind.fill,list(gwis2,bbj2,can))	
	return(list(Dat,ref_dat))
}


# prep_data2<-function(Dat=NULL){		
# 	gwis<-data.frame(Dat[1],stringsAsFactors=F)
# 	gtex<-data.frame(Dat[2],stringsAsFactors=F)
# 	eqtl<-data.frame(Dat[3],stringsAsFactors=F)
# 	crc<-data.frame(Dat[4],stringsAsFactors=F)
# 	gtex2<-effect_allele_gtex_data(dat=gtex)
# 	ref_dat<-gwis[gwis$trait == "AA:DGLA / D5D",]
# 	gwis2<-gwis[gwis$trait != "AA:DGLA / D5D",]

# 	aa_info<-read.table("~/fatty-acids/colocalisation/data/aa_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
# 	dgla_info<-read.table("~/fatty-acids/colocalisation/data/dgla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
# 	gla_info<-read.table("~/fatty-acids/colocalisation/data/gla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
# 	la_info<-read.table("~/fatty-acids/colocalisation/data/la_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	
# 	snps1<-snps_pass(Dat=aa_info)
# 	snps2<-snps_pass(Dat=dgla_info)

# 	Pos<-ref_dat$marker %in% snps1 & ref_dat$marker %in% snps2	
# 	ref_dat<-ref_dat[which(Pos),]		
	
# 	Dat<-do.call(rbind.fill,list(gwis2,gtex2,eqtl,crc))
	
# 	return(list(Dat,ref_dat))
# }
	
	

prep_data<-function(Dat=NULL){		
	gwis<-data.frame(Dat[1],stringsAsFactors=F)
	gtex<-data.frame(Dat[2],stringsAsFactors=F)
	eqtl<-data.frame(Dat[3],stringsAsFactors=F)
	crc<-data.frame(Dat[4],stringsAsFactors=F)
	gtex2<-effect_allele_gtex_data(dat=gtex)
	ref_dat<-gwis[gwis$trait == "AA:DGLA / D5D",]
	gwis2<-gwis[gwis$trait != "AA:DGLA / D5D",]

	aa_info<-read.table("~/fatty-acids/colocalisation/data/aa_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	dgla_info<-read.table("~/fatty-acids/colocalisation/data/dgla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	gla_info<-read.table("~/fatty-acids/colocalisation/data/gla_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	la_info<-read.table("~/fatty-acids/colocalisation/data/la_info_snps80.txt",sep="\t",head=T,stringsAsFactors=F)
	
	snps1<-snps_pass(Dat=aa_info)
	snps2<-snps_pass(Dat=dgla_info)

	Pos<-ref_dat$marker %in% snps1 & ref_dat$marker %in% snps2	
	ref_dat<-ref_dat[which(Pos),]		
	
	Dat<-do.call(rbind.fill,list(gwis2,gtex2,eqtl,crc))
	
	return(list(Dat,ref_dat))
}
	

 
snps_pass<-function(Dat=NULL){
	Dat1<-Dat[Dat$snp!=".",]
	# Dat1[Dat1$V2=="rs174546",]
	# ref[ref$V2=="rs174546",]
	Dat2<-Dat[Dat$snp==".",]
	chr<-unlist(strsplit(ref$V1,split="chr"))
	ref$chr<-chr[chr!=""]
	Dat3<-merge(Dat2,ref,by.x=c("chr","bp"),by.y=c("chr","V2"))
	snps1<-Dat1$snp
	snps2<-Dat3$V4	
	snps<-c(snps1,snps2)
	return(snps)
}
format_data<-function(Dat=NULL){
	Markers<-data.frame(Dat[3])	
	B.matrix<-data.frame(Dat[5])		
	SE.matrix<-data.frame(Dat[6])
	Trait_names<-unlist(Dat[4])
	Dat_list<-NULL	
	for(i in 1:ncol(B.matrix)){
		print(i)
		Dat_list[[i]]<-do.call(cbind,list(B.matrix[,i],SE.matrix[,i],Trait_names[i]))
	}
	Dat2<-data.frame(do.call(rbind,Dat_list),stringsAsFactors=F)
	rownames(Markers)<-NULL
	Dat<-cbind(Markers,Dat2)
	names(Dat)<-c("SNP","chr","pos_grch37","effect_allele","other_allele","eaf","b","se","trait_detail")
	Dat$effect_allele<-toupper(Dat$effect_allele)
	Dat$other_allele<-toupper(Dat$other_allele)
    Pos<-grep("FADS1",Dat$trait_detail)    
    Dat$trait<-NA
    Dat$trait[Pos]<-"FADS1 expression"
    Pos<-grep("FADS2",Dat$trait_detail)
    Dat$trait[Pos]<-"FADS2 expression"    
    Dat$trait[Dat$trait_detail =="Colorectal cancer GECCO/CORECT/CCFR" ]<-"Colorectal cancer"
    Dat$study<-NA
    Dat$study[Dat$trait_detail =="Colorectal cancer GECCO/CORECT/CCFR" ]<-"GECCO/CORECT/CCFR"
    Dat$study[grep("eQTLGen",Dat$trait_detail)]<-"eQTLGen"
	Dat$study[grep("GTEx",Dat$trait_detail)]<-"GTEx"
	Dat$study[grep("CHARGE",Dat$trait_detail)]<-"CHARGE"   
    Dat$trait[Dat$trait_detail =="AA:DGLA / D5D (CHARGE)"  ]<-"AA:DGLA"
    Dat$trait[Dat$trait_detail =="GLA:LA / D6D (CHARGE)"  ]<-"GLA:LA" 
    Dat1<-Dat[grep("expression",Dat$trait_detail),]   
    tissue<-trimws(unlist(strsplit(Dat1$trait_detail,split="in")))
	Dat1$tissue<-tissue[seq(2,length(tissue),by=4)]   	
	Dat2<-Dat[grep("expression",Dat$trait_detail,invert=T),]   
    Dat2$tissue<-NA
    Dat<-rbind(Dat1,Dat2)
    return(Dat)
}

