# install.packages("metap")
# library(metap)
# source("~/fatty-acids/mr/scripts/etract_snps_function")
# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/
# library(meta)
library(ggforestplot)
library(plyr)
library(ggplot2)

snps_all_tab<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/targetSNPs_table_Dec18.txt",sep="\t",head=T,stringsAsFactors=F)

teqtl<-readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/teqtl_fads.txt") # rs7412 # rs735665 # rs964184 # rs6511720

snplist<-list_of_snps()
# which(snplist=="rs735665")
Res<-extract_snps2(snps=snplist)
res.dat<-do.call(rbind.fill,Res)
save(list=c("res.dat","teqtl","snps_all_tab"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/pufas_replicate.Rdata")
cd ~/fatty-acids/mr/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/pufas_replicate.Rdata .

load("~/fatty-acids/mr/data/pufas_replicate.Rdata")
# res.dat$p[which(res.dat$snp=="rs36016715" )],c("trait")] & res.dat$p<5e-8),]

Dat<-get_meta_data(Dat=res.dat)
# min(Dat$p[which(Dat$snp=="rs36016715" )])
Dat[which(Dat$snp=="rs36016715" ),c("trait","p")]
# Dat2<-Dat[Dat$snp %in% snps_all_tab$SNP,]
Dat2<-format_data(Dat=Dat,type.fa=)

Dat3<-find_proxies(Har=Dat2)
Dat4<-harmonise_alleles(HarR=Dat3)
# Dat4[Dat4$SNP  == "rs174549",c("effect_allele_freq","author","effect_allele","other_allele")]

# Dat4<-Dat3[Dat3$p<0.05/161,]
Dat5<-find_region(Dat=Dat4)
Dat6<-exclude_fads_scd_elovl2(Dat=Dat5)

# Dat7<-find_discovery_study(Dat=Dat6)
Dat8_list<-format_results(Dat=Dat6)
Dat8<-data.frame(Dat8_list[1],stringsAsFactors=F)

Dat9<-data.frame(Dat8_list[2],stringsAsFactors=F)
# unique(Dat9[which(is.na(Dat9$z)),c("SNP","beta","se","trait","b_sd","p","p.fisher","consortium")])

write.table(Dat8,"~/fatty-acids/mr/data/pufasnps_repv3.txt",sep="\t",col.names=T,row.names=F,quote=F)
Dat8<-format_results(Dat=Dat6,keep_sig=FALSE)
write.table(Dat8,"~/fatty-acids/mr/data/pufasnps_repv2_keepnonsig.txt",sep="\t",col.names=T,row.names=F,quote=F)

Dat10<-format_results2(Dat=Dat9)

traits<-unique(Dat10$trait)
snps<-unique(Dat10$SNP)

i<-8

# for(i in 1:length(traits)){
	plot<-make_plot(Dat=Dat10,Snp=snps[i])
	pdf(unlist(plot[2]))
		plot[1]	
	dev.off()
# }


length(unique(Dat8$region))
Dat.t<-find_transeqtls(Dat=Dat3)
Dat5<-find_region(Dat=Dat.t)
Dat6<-exclude_fads_scd_elovl2(Dat=Dat5)
Dat7<-find_discovery_study(Dat=Dat6)
Dat8<-format_results(Dat=Dat7)
write.table(Dat8,"~/fatty-acids/mr/data/pufasnps_rep_transeqtls.txt",sep="\t",col.names=T,row.names=F,quote=F)
						 
Res<-extract_snps2(snps="rs735665")
Res1<-do.call(rbind.fill,Res)
Res1<-Res1[order(Res1$p),]
save(list=c("Res1","teqtl","snps_all_tab"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/teqtl_rs735665_v2.Rdata")

cd ~/fatty-acids/mr/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/teqtl_rs735665_v2.Rdata .
load("~/fatty-acids/mr/data/teqtl_rs735665_v2.Rdata")
# Res1$file
Dat<-get_meta_data(Dat=Res1)
Dat2<-format_data(Dat)
Dat8<-format_results(Dat=Dat2)

write.table(Dat8,"~/fatty-acids/mr/data/pufasnps_rep_transeqtls_rs735665.txt",sep="\t",col.names=T,row.names=F,quote=F)

########################
# pQTLs for fatty acids#
########################

# source("~/fatty-acids/mr/scripts/etract_snps_function")
library(TwoSampleMR)
ao<-available_outcomes()
ao$trait[which(ao$id=="prot-a-1013")]
ao2<-ao[grep("fatty acid",ao$trait,ignore.case=T),]
ids<-unique(ao2$id[grep("protein",ao2$trait,ignore.case=T)])

exposure<-  extract_instruments(
    ids,
       p1 = 5e-08,
       clump = TRUE,
       p2 = 5e-08,
       r2 = 0.001,
       kb = 10000,
       access_token = ieugwasr::check_access_token(),
       force_server = TRUE
     )

)
snps<-exposure$SNP
snps<-c("rs138961357", "rs145534968", "rs149081536", "rs2241883")  
Res<-extract_snps2(snps=snps)
res.dat<-do.call(rbind.fill,Res)
save(list=c("res.dat","snps"),file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/pufas_pqtls.Rdata")
cd ~/fatty-acids/mr/data
scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/pufas_pqtls.Rdata .
load("~/fatty-acids/mr/data/pufas_pqtls.Rdata")
Dat<-get_meta_data(Dat=res.dat)
Dat2<-format_data(Dat=Dat,type.fa=c("saturated","monounsaturated","n3 polyunsaturated","n6 polyunsaturated","n3 or n6 polyunsaturated"))

unique(Dat2[Dat2$p<0.05,c("SNP","p","trait")])

Dat8_list<-format_results(Dat=Dat6)
Dat8<-data.frame(Dat8_list[1],stringsAsFactors=F)

Dat9<-data.frame(Dat8_list[2],stringsAsFactors=F)
# unique(Dat9[which(is.na(Dat9$z)),c("SNP","beta","se","trait","b_sd","p","p.fisher","consortium")])

write.table(Dat8,"~/fatty-acids/mr/data/pufasnps_repv3.txt",sep="\t",col.names=T,row.names=F,quote=F)
Dat8<-format_results(Dat=Dat6,keep_sig=FALSE)
write.table(Dat8,"~/fatty-acids/mr/data/pufasnps_repv2_keepnonsig.txt",sep="\t",col.names=T,row.names=F,quote=F)

Dat10<-format_results2(Dat=Dat9)

traits<-unique(Dat10$trait)
snps<-unique(Dat10$SNP)



make_plot<-function(Dat=NULL,Trait=NULL,Snp=NULL){


	if(!is.null(Snps)){
		Plot.Dat<-Dat[Dat$SNP==Snp,]
		Plot.Dat$Name<-Plot.Dat$trait
		plot_file<-paste("~/fatty-acids/mr/results/plots/",Snp,"_PUFA_replication.pdf",sep="")
		plot_file<-gsub(" ","_",plot_file)
		plot_file<-gsub(":",";",plot_file)
		Title<-Snp
	}

	if(!is.null(Trait)){
		Dat.om36<-Dat[Dat$trait %in% c("omega-6 fatty acids","omega-3 fatty acids"),]
		Plot.Dat<-Dat[Dat$trait==Trait,]
		if(!Trait %in% c("omega-6 fatty acids","omega-3 fatty acids")){	
			Plot.Dat<-rbind(Plot.Dat,Dat.om36)
		}
		plot_file<-paste("~/fatty-acids/mr/results/plots/",Trait,"_PUFA_replication.pdf",sep="")
		plot_file<-gsub(" ","_",plot_file)
		plot_file<-gsub(":",";",plot_file)
		Plot.Dat$Name<-Plot.Dat$SNP
		Plot.Dat$Shape<-Plot.Dat$trait
		Title<-Trait
	}
	# Plot.Dat<-Plot.Dat[Plot.Dat$type.fa == "n6 polyunsaturated",]
	# Plot.Dat<-Plot.Dat[!Plot.Dat$trait %in% c("adrenic acid (22:4n6)","docosapentaenoic acid (22:5n6)")  ,]
	# unique(Plot.Dat$consortium)
	 # Plot.Dat<-Plot.Dat[Plot.Dat$consortium!="CHARGE; FraminghamHOS",]
	
	Plot.Dat$Colour<-Plot.Dat$consortium
	Plot.Dat$Shape
	# Plot.Dat$Shape<-"21L" 
	# Plot.Dat$Shape[Plot.Dat$consortium == "CHARGE; FraminghamHOS"]<-"23L"
	unique(Plot.Dat$trait)
	excl.fas<-c("dihomo-linolenic acid (20:3n3 or n6)","alpha or gamma linolenic acid (18:3n3 or 6)") 
	# excl.fas<-c("omega-3 fatty acids","omega-6 fatty acids","dihomo-linolenic acid (20:3n3 or n6)","alpha or gamma linolenic acid (18:3n3 or 6)")
	Plot.Dat<-Plot.Dat[!Plot.Dat$trait %in% excl.fas,]
	# pdf("plot")
	plot<-forestplot(
		df = Plot.Dat,
		logodds = TRUE,
		name=Name,
		estimate=b_sd,
		se=se_sd,
		# shape=Shape,
		colour = Colour,
		xlab = "") +
		labs(title=Title,size=1)
	# Plot<-make_plot(Plot.Dat=Res1,Colour="exposure",meta_analysis=TRUE)
	# pdf(plot_file)

			
		# dev.off()
		return(list(plot,plot_file))
	# }
}

		# +
  # You may also want to add a manual shape scale to mark meta-analysis with a
  # diamond shape
		  # ggplot2::scale_shape_manual(
	   #  values = c(23L, 21L, 21L, 21L, 21L,21L),
	   #  labels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS","dsf"))
		# +
		# theme(plot.title = element_text(size = text.title))+
		# theme(text = element_text(size=text.names))
		# eval(parse(text=Colour))




format_results2<-function(Dat=NULL,snps.keep=Dat8$SNP[Dat8$replicates]){
	
	snps.keep<-c("rs735665",snps.keep) #	"rs735665" is trans eQTL for fADS1/2 that does not replicate at my threshold
	Dat<-Dat[Dat$SNP %in% snps.keep,]

	Dat$region[Dat$SNP=="rs735665"]<-9
	Dat$file[which(Dat$SNP=="rs735665")]<-9
	Dat$ID<-as.numeric(as.factor(paste(Dat$consortium,Dat$trait)))


	# Exclude duplicate regions retaining the index SNP for each region
	Dat1<-Dat[Dat$SNP == "rs780094",]
	Dat2<-Dat[which(Dat$region != 555),]

	Dat<-rbind(Dat1,Dat2)

	Dat1<-Dat[Dat$SNP == "rs79225634",]
	Dat2<-Dat[which(Dat$region != 1915),]
	Dat<-rbind(Dat1,Dat2)

	Dat1<-Dat[Dat$SNP == "rs4985155",]
	Dat2<-Dat[which(!Dat$region %in% c(4857, 4858)),]
	Dat<-rbind(Dat1,Dat2)

	Dat1<-Dat[Dat$SNP == "rs7412",]
	Dat2<-Dat[which(Dat$region != 5419),]
	Dat<-rbind(Dat1,Dat2)

	Dat1<-Dat[Dat$SNP == "rs143988316",]
	Dat2<-Dat[which(!Dat$region %in% c(5367,5368)),]

	Dat<-rbind(Dat1,Dat2)
		
	ID<-paste(Dat$ID,Dat$region)
	# duplicated(ID)
	length(unique(Dat$SNP))
		unique(Dat$trait)



	# which(Dat$SNP=="rs735665")
	 All.fas<-c("omega-3 fatty acids","docosahexaenoic acid (22:6n3)","omega-6 fatty acids",
	 	"dihomo-linolenic acid (20:3n3 or n6)","linoleic acid (18:2n6)",
	 	"gamma-linolenic acid (18:3n6)","arachidonic acid (20:4n6)","eicosapentaenoic acid (20:5n3)",
	 	"stearidonic acid (18:4n3)","docosapentaenoic acid (22:5n3)",
	 	"alpha-linolenic acid (18:3n3)","alpha or gamma linolenic acid (18:3n3 or 6)",
	 	"adrenic acid (22:4n6)","docosapentaenoic acid (22:5n6)","eicosadienoic acid (20:2n6)",
	 	"dihomo-gamma-linolenic acid (20:3n6)")
	 # excl.fas<-c("omega-3 fatty acids","omega-6 fatty acids","dihomo-linolenic acid (20:3n3 or n6)","alpha or gamma linolenic acid (18:3n3 or 6)")
	 excl.fas<-c("dihomo-linolenic acid (20:3n3 or n6)","alpha or gamma linolenic acid (18:3n3 or 6)") 
	 # Dat<-Dat[!Dat$trait %in% excl.fas,]
 	length(unique(Dat$SNP))
	length(unique(Dat$trait))
		length(unique(Dat$consortium))

	Dat$type.fa[grep("n6",Dat$trait)]<-"n6 polyunsaturated"
	Dat$type.fa[grep("n3",Dat$trait)]<-"n3 polyunsaturated"
	Dat$type.fa[grep("or",Dat$trait)]<-"n3 or n6 polyunsaturated"

# "rs780094" 555
# rs79225634 1915
# rs4985155 4857, 4858
# rs7412 5419

	# id<-paste(Dat$trait,Dat$region)
	# Dat[which(is.na(Dat$trait)),]
	

	# table(Dat$trait)
	Pos<-which(is.na(Dat$b_sd))
	Dat$z[Pos]<-Dat$beta[Pos]/Dat$se[Pos]
	Dat[is.na(Dat$effect_allele_freq),]
	Dat$b_sd[Pos]<-b_sd(z=Dat$z[Pos],maf=Dat$effect_allele_freq[Pos],n=Dat$n[Pos])
	Dat$se_sd[Pos]<-Dat$b_sd[Pos]/Dat$z[Pos]
	Dat$consortium[Dat$consortium=="14 cohorts: EGCUT,ERF,FTC,FR97,COROGENE,GenMets,HBCS,KORA,LLS,NTR,NFBC1966,PredictCVD,PROTE,YFS)"]<-"Kettunen"

	
	Dat1<-Dat[grep(")",Dat$trait),]
	Dat2<-Dat[grep(")",Dat$trait,invert=T),]
	Str<-unlist(strsplit(Dat1$trait,split="\\("))
	Chain<-Str[seq(2,length(Str),by=2)]
	Chain<-unlist(strsplit(Chain,split=")"))
	Dat1$chain<-Chain
	Chain<-unlist(strsplit(Dat1$chain,split=":"))
	Chain.length<-Chain[seq(1,length(Chain),by=2)]
	Dat1$chain.length<-Chain.length
	Dat<-rbind(Dat1,Dat2)
	
	Excl<-c("CHARGE; FraminghamHOS; SCHS pooled; Shin","CHARGE; FraminghamHOS; SCHS pooled", "Kettunen; CHARGE; FraminghamHOS; SCHS pooled; Shin")
	Dat<-Dat[!Dat$consortium %in% Excl,]

	# snps that replicated in meta analysed datasets across CHARGE and FRamingham or across CHARGE, Framingham and SCHS:

	snps.cfs<-c("rs10051260","rs79225634","rs36016715","rs3741298","rs821840","rs72654473","rs7412")
	Dat$snps.cfs<-FALSE
	Dat$snps.cfs[Dat$SNP %in% snps.cfs]<-TRUE
	
	Dat<-Dat[,c("SNP","trait","b_sd","se_sd","p","p.fisher","consortium","chain","chain.length","type.fa","region","gene")]
	return(Dat)
}



find_transeqtls<-function(Dat=NULL){
	teqtl<-teqtl[teqtl != "rs6511720"] #rs6511720 identified by taking 161 fatty acids sNPs and looking these up in eqtlgen. where as other teqtls identified FDR<0.05 through hypothesis free analysis
	Dat.t<-Dat[Dat$SNP %in% teqtl,]
	chr11<-fa.reg[fa.reg$Chr == 11,] 
	unique(as.numeric(chr11$SNP.GRCh38.p12) - 123490689 )
	miss<-teqtl[!teqtl %in% Dat$SNP]
	load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata")
	
	# Temp<-fa.reg[fa.reg$SNP =="rs36016715",]
	# unique(Temp[as.numeric(Temp$pval)<5e-8,c("SNP","pval","trait","author")])

	 # unique(fa.reg$region)[order(unique(fa.reg$region))]
	miss %in% fa.reg$SNP #not present amongst target SNPs 
	miss %in% fa.reg$ProxySNP #not present amongst proxies
	return(Dat.t)
}

# Dat<-Dat6

format_results<-function(Dat=NULL,keep_sig=TRUE){
	Dat1<-Dat[order(Dat$p),]
	Dat.list<-NULL
	snps<-unique(Dat$SNP)	
	# i<-which(snps=="rs1077835")
	# i<-which(snps=="rs139872716")
	options(digits = 3) #round scientic numbers 
	Dat_list<-NULL
	for(i in 1:length(snps)){
		print(i)
		Dat<-Dat1[Dat1$SNP == snps[i],]
		# Studies<-unique(Dat$consortium)
		Traits<-unique(Dat$trait)

		Dat$beta[Dat$beta==0]<-0.00000000001
		Dat$z=Dat$beta/Dat$se
		Dat$b_sd<-b_sd(z=Dat$z,maf=Dat$effect_allele_freq,n=Dat$n)
		Dat$se_sd<-Dat$b_sd/Dat$z
		
		# snp2<-unique(Dat$SNP)
		# snp2[!snp2 %in% snp1 ]
		# j<-which(Traits=="adrenic acid (22:4n6)")
		Dat.fe_list<-NULL
		for(j in 1:length(Traits)){
			Dat.p<-Dat[Dat$trait == Traits[j],]
			Dat.p<-Dat.p[order(Dat.p$consortium),]
			if(length(unique(Dat.p$effect_allele))!=1) stop("effect alleles not consistent across studies")
			# Dat.p[,c("SNP","trait","consortium","p")]
			Dat.p$consortium<-gsub("14 cohorts: EGCUT,ERF,FTC,FR97,COROGENE,GenMets,HBCS,KORA,LLS,NTR,NFBC1966,PredictCVD,PROTE,YFS)","Kettunen",Dat.p$consortium)
			Dat.p$consortium<-gsub("TwinsUK/KORA","Shin",Dat.p$consortium)
			# Dat.temp<-Dat.p[Dat.p$consortium %in% c("CHARGE","FraminghamHOS"),]			 
			Dat.all<-NULL
			if(all(c( "Kettunen" ,"CHARGE","FraminghamHOS","SCHS pooled" ,"Shin") %in% Dat.p$consortium)){
				b<-Dat.p$b_sd
				se<-Dat.p$se_sd
				# Meta<-metagen(TE=b,seTE=se,comb.fixed=T,sm="MD")
				w<-1/se^2
				b.fe<-sum(b*w)/(sum(w))
				se.fe<-sqrt(sum(w)^-1)
				z.fe<-abs(b.fe/se.fe)
				p.fe<-pnorm(z.fe,lower.tail=F)*2
				Fisher<-sumlog(Dat.p$p) # Combine p-values by the sum of logs method, also known as Fisher's method, and sometimes as the chi-square (2) method.		
				p.fisher<-Fisher$p
				Studies<-paste(unique(Dat.p$consortium),collapse="; ")
				Dat.all<-data.frame(matrix(c(b.fe,se.fe,z.fe,p.fe,p.fisher,Studies),nrow=1,ncol=6),stringsAsFactors=F)
			}

			Dat.cf<-NULL
			Dat.temp<-Dat.p[Dat.p$consortium %in% c("CHARGE","FraminghamHOS"),]			 
			if(nrow(Dat.temp)==2){
				b<-Dat.temp$b_sd
				se<-Dat.temp$se_sd
				w<-1/se^2
				b.fe<-sum(b*w)/(sum(w))
				se.fe<-sqrt(sum(w)^-1)
				z.fe<-abs(b.fe/se.fe)
				p.fe<-pnorm(z.fe,lower.tail=F)*2
				Fisher<-sumlog(Dat.temp$p) # Combine p-values by the sum of logs method, also known as Fisher's method, and sometimes as the chi-square (2) method.		
				p.fisher<-Fisher$p
				Studies<-paste(unique(Dat.temp$consortium),collapse="; ")
				Dat.cf<-data.frame(matrix(c(b.fe,se.fe,z.fe,p.fe,p.fisher,Studies),nrow=1,ncol=6),stringsAsFactors=F)
				# Dat.fe<-rbind(Dat.fe,Dat.cf2)
			}			

			Dat.cfs<-NULL
			Dat.temp<-Dat.p[Dat.p$consortium %in% c("CHARGE","FraminghamHOS","SCHS pooled"),]			 
			if(nrow(Dat.temp)==3){
				b<-Dat.temp$b_sd
				se<-Dat.temp$se_sd
				w<-1/se^2
				b.fe<-sum(b*w)/(sum(w))
				se.fe<-sqrt(sum(w)^-1)
				z.fe<-abs(b.fe/se.fe)
				p.fe<-pnorm(z.fe,lower.tail=F)*2
				Fisher<-sumlog(Dat.temp$p) # Combine p-values by the sum of logs method, also known as Fisher's method, and sometimes as the chi-square (2) method.		
				p.fisher<-Fisher$p
				Studies<-paste(unique(Dat.temp$consortium),collapse="; ")
				Dat.cfs<-data.frame(matrix(c(b.fe,se.fe,z.fe,p.fe,p.fisher,Studies),nrow=1,ncol=6),stringsAsFactors=F)
				# Dat.fe<-rbind(Dat.fe,Dat.cf2)
			}	


			Dat.cfks<-NULL
			Dat.temp<-Dat.p[Dat.p$consortium %in% c("CHARGE","FraminghamHOS","Kettunen","Shin"),]			 
			if(nrow(Dat.temp)==4){
				b<-Dat.temp$b_sd
				se<-Dat.temp$se_sd
				w<-1/se^2
				b.fe<-sum(b*w)/(sum(w))
				se.fe<-sqrt(sum(w)^-1)
				z.fe<-abs(b.fe/se.fe)
				p.fe<-pnorm(z.fe,lower.tail=F)*2
				Fisher<-sumlog(Dat.temp$p) # Combine p-values by the sum of logs method, also known as Fisher's method, and sometimes as the chi-square (2) method.		
				p.fisher<-Fisher$p
				Studies<-paste(unique(Dat.temp$consortium),collapse="; ")
				Dat.cfks<-data.frame(matrix(c(b.fe,se.fe,z.fe,p.fe,p.fisher,Studies),nrow=1,ncol=6),stringsAsFactors=F)
				# Dat.fe<-rbind(Dat.fe,Dat.cf2)
			}			

			# no examples where snp only present in ket and shin
			# Dat.ks<-NULL
			# Dat.temp<-Dat.p[Dat.p$consortium %in% c("Kettunen","Shin"),]			 
			# if(nrow(Dat.temp)==3){
			# 	b<-Dat.temp$b_sd
			# 	se<-Dat.temp$se_sd
			# 	w<-1/se^2
			# 	b.fe<-sum(b*w)/(sum(w))
			# 	se.fe<-sqrt(sum(w)^-1)
			# 	z.fe<-abs(b.fe/se.fe)
			# 	p.fe<-pnorm(z.fe,lower.tail=F)*2
			# 	Fisher<-sumlog(Dat.temp$p) # Combine p-values by the sum of logs method, also known as Fisher's method, and sometimes as the chi-square (2) method.		
			# 	p.fisher<-Fisher$p
			# 	Studies<-paste(unique(Dat.temp$consortium),collapse="; ")
			# 	Dat.ks<-data.frame(matrix(c(b.fe,se.fe,z.fe,p.fe,p.fisher,Studies),nrow=1,ncol=6),stringsAsFactors=F)
			# 	# Dat.fe<-rbind(Dat.fe,Dat.cf2)
			# }			
				
			Dat.cfss<-NULL
			Dat.temp<-Dat.p[Dat.p$consortium %in% c("CHARGE","FraminghamHOS","SCHS pooled","Shin"),]			 
			if(nrow(Dat.temp)==4){
				b<-Dat.temp$b_sd
				se<-Dat.temp$se_sd
				w<-1/se^2
				b.fe<-sum(b*w)/(sum(w))
				se.fe<-sqrt(sum(w)^-1)
				z.fe<-abs(b.fe/se.fe)
				p.fe<-pnorm(z.fe,lower.tail=F)*2
				Fisher<-sumlog(Dat.temp$p) # Combine p-values by the sum of logs method, also known as Fisher's method, and sometimes as the chi-square (2) method.		
				p.fisher<-Fisher$p
				Studies<-paste(unique(Dat.temp$consortium),collapse="; ")
				Dat.cfss<-data.frame(matrix(c(b.fe,se.fe,z.fe,p.fe,p.fisher,Studies),nrow=1,ncol=6),stringsAsFactors=F)
				# Dat.fe<-rbind(Dat.fe,Dat.cf2)
				# Meta<-metagen(TE=b,seTE=se,comb.fixed=T,sm="MD")
			}	

			Dat.fe<-do.call(rbind,list(Dat.all,Dat.cf,Dat.cfs,Dat.cfss))
			# Dat.fe$SNP <-snps[i]
			# Dat.fe$trait<-Traits[j]	
			Dat.fe_list[[Traits[j]]]<-Dat.fe
		}	

		# Trait<-names(Dat.fe_list)
		Dat_traits<-NULL
		if(class(Dat.fe_list)=="list"){
			Dat_traits<-do.call(rbind,Dat.fe_list)
		}
		Pos<-regexpr(")",rownames(Dat_traits))
		Dat_traits$trait<-substr(rownames(Dat_traits),1,Pos)
		rownames(Dat_traits)<-NULL		
		Dat_traits$SNP<-snps[i]
		Dat_list[[i]]<-Dat_traits
	}

	Datr<-do.call(rbind,Dat_list)
	names(Datr)<-c("b_sd","se_sd","z","p","p.fisher","consortium","trait","SNP")
	Datr$b_sd<-as.numeric(Datr$b_sd)
	Datr$se_sd<-as.numeric(Datr$se_sd)
	Datr$z<-as.numeric(Datr$z)
	Datr$p<-as.numeric(Datr$p)
	Datr$p.fisher<-as.numeric(Datr$p.fisher)				


	# Dat.r<-do.call(rbind,Dat_list)
	
	Datr2<-rbind.fill(Datr,Dat1)

	if(keep_sig){
		Pos1<-which(Datr2$p<0.05/161)
		Pos2<-which(Datr2$p.fisher<0.05/161) 
		Pos3<-unique(c(Pos1,Pos2))
		Datr3<-Datr2[Pos3,]
	}
	
	Dat.list<-NULL
	snps<-unique(Datr3$SNP)

	for(i in 1:length(snps)){
		print(i)
		Dat.temp<-Datr3[Datr3$SNP == snps[i],]	
		
		# charge
		Pos<-which(Dat.temp$consortium=="CHARGE")
		Dat.temp$charge<-NA
		Dat.temp$charge[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$charge<-paste(Dat.temp$charge[!is.na(Dat.temp$charge)],collapse="; ")
		head(Dat.temp)		
		# fhs
		Pos<-which(Dat.temp$consortium=="FraminghamHOS")
		Dat.temp$fhs<-NA
		Dat.temp$fhs[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$fhs<-paste(Dat.temp$fhs[!is.na(Dat.temp$fhs)],collapse="; ")

		Pos<-which(Dat.temp$consortium=="SCHS pooled")
		Dat.temp$schs<-NA
		Dat.temp$schs[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$schs<-paste(Dat.temp$schs[!is.na(Dat.temp$schs)],collapse="; ")

		# ket
		Pos<-which(Dat.temp$consortium=="14 cohorts: EGCUT,ERF,FTC,FR97,COROGENE,GenMets,HBCS,KORA,LLS,NTR,NFBC1966,PredictCVD,PROTE,YFS)")
		Dat.temp$ket<-NA
		Dat.temp$ket[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$ket<-paste(Dat.temp$ket[!is.na(Dat.temp$ket)],collapse="; ")

		#shin 	
		Pos<-which(Dat.temp$consortium=="TwinsUK/KORA")
		Dat.temp$shin<-NA
		Dat.temp$shin[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$shin<-paste(Dat.temp$shin[!is.na(Dat.temp$shin)],collapse="; ")


		# "Kettunen; CHARGE; FraminghamHOS; SCHS pooled; Shin"
		Pos<-which(Dat.temp$consortium=="Kettunen; CHARGE; FraminghamHOS; SCHS pooled; Shin"  )
		Dat.temp$all<-NA
		Dat.temp$all[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2)," pfisher=",formatC(as.numeric(Dat.temp$p.fisher[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$all<-paste(Dat.temp$all[!is.na(Dat.temp$all)],collapse="; ")


		#"CHARGE; FraminghamHOS; SCHS pooled"
		Pos<-which(Dat.temp$consortium=="CHARGE; FraminghamHOS; SCHS pooled"   )
		Dat.temp$cfs<-NA
		Dat.temp$cfs[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2)," pfisher=",formatC(as.numeric(Dat.temp$p.fisher[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$cfs<-paste(Dat.temp$cfs[!is.na(Dat.temp$cfs)],collapse="; ")

		#"Kett & shin
		Pos<-which(Dat.temp$consortium=="CHARGE; FraminghamHOS; SCHS pooled"   )
		Dat.temp$ks<-NA
		Dat.temp$ks[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2)," pfisher=",formatC(as.numeric(Dat.temp$p.fisher[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$ks<-paste(Dat.temp$ks[!is.na(Dat.temp$ks)],collapse="; ")

		# "CHARGE; FraminghamHOS"  
		Pos<-which(Dat.temp$consortium=="CHARGE; FraminghamHOS"   )
		Dat.temp$cf<-NA
		Dat.temp$cf[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2)," pfisher=",formatC(as.numeric(Dat.temp$p.fisher[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$cf<-paste(Dat.temp$cf[!is.na(Dat.temp$cf)],collapse="; ")


		# CHARGE; FraminghamHOS; SCHS pooled; Shin
		Pos<-which(Dat.temp$consortium=="CHARGE; FraminghamHOS; SCHS pooled; Shin" )
		Dat.temp$cfss<-NA
		Dat.temp$cfss[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2)," pfisher=",formatC(as.numeric(Dat.temp$p.fisher[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$cfss<-paste(Dat.temp$cfss[!is.na(Dat.temp$cfss)],collapse="; ")
		
		# CHARGE; FraminghamHOS; Kettunen and Shin
		Pos<-which(Dat.temp$consortium=="CHARGE; FraminghamHOS; Kettunen; Shin" )
		Dat.temp$cfks<-NA
		Dat.temp$cfks[Pos]<-paste(Dat.temp$trait[Pos]," p=",	formatC(as.numeric(Dat.temp$p[Pos]), format = "e", digits = 2)," pfisher=",formatC(as.numeric(Dat.temp$p.fisher[Pos]), format = "e", digits = 2),sep="")
		Dat.temp$cfks<-paste(Dat.temp$cfks[!is.na(Dat.temp$cfks)],collapse="; ")

		Dat.list[[i]]<-Dat.temp[1,]
	}
	
	Dat.list<-do.call(rbind,Dat.list)

	Dat.list<-Dat.list[,c("SNP","charge","fhs","schs","ket","shin","all","cf","cfs","cfss","cfks")]
	
	charge<-Dat.list$charge!=""
	charge[charge]<-1
	charge[!charge]<-0

	fhs<-Dat.list$fhs!=""
	fhs[fhs]<-1
	fhs[!fhs]<-0
	
	schs<-Dat.list$schs!=""
	schs[schs]<-1
	schs[!schs]<-0

	ket<-Dat.list$ket!=""
	ket[ket]<-1
	ket[!ket]<-0

	shin<-Dat.list$shin!=""
	shin[shin]<-1
	shin[!shin]<-0

	All<-Dat.list$all!=""
	All[All]<-1
	All[!All]<-0

	cf<-Dat.list$cf!=""
	cf[cf]<-1
	cf[!cf]<-0

	cfs<-Dat.list$cfs!=""
	cfs[cfs]<-1
	cfs[!cfs]<-0
	
	cfss<-Dat.list$cfss!=""
	cfss[cfss]<-1
	cfss[!cfss]<-0
	
	cfks<-Dat.list$cfks!=""
	cfks[cfks]<-1
	cfks[!cfks]<-0

	Rep<-data.frame(matrix(c(charge,fhs,schs,ket,shin,All,cf,cfs,cfss),ncol=9,nrow=length(charge),byrow=FALSE),stringsAsFactors=F)
	names(Rep)<-c("charge","fhs","schs","ket","shin","All","cf","cfs","cfss")
	Rep$test<-unlist(lapply(1:nrow(Rep),FUN=function(i)
		sum(Rep[i,])))
	Rep$rep<-"no"
	Rep$rep[Rep$test>1]<-"replicates P<0.05/161"
	
	# set rep to no for test=2 if the 2 studies overlap
	Rep$id<-1:nrow(Rep)
	Rep1<-Rep[Rep$test==1,]
	Rep2<-Rep[Rep$test==2,]
	Rep3<-Rep[Rep$test==3,]
	Rep4<-Rep[Rep$test==4,]
	Rep5<-Rep[Rep$test==5,]
	Rep6<-Rep[Rep$test>5,]
	Pos.excl1<-which(Rep2$ket==1 & Rep2$All==1)
	Pos.excl2<-which(Rep2$shin==1 & Rep2$cfss==1)
	Pos.excl3<-which(Rep2$charge==1 & Rep2$cf==1)
	Pos.excl4<-which(Rep2$fhs==1 & Rep2$cf==1)
	Pos.excl5<-which(Rep2$schs==1 & Rep2$cfs==1)
	Pos.excl<-c(Pos.excl1,Pos.excl2,Pos.excl3,Pos.excl4,Pos.excl5) 
	Rep2$rep[Pos.excl]<-"no"
	Pos.excl1<-which(Rep3$shin==1 & Rep3$All==1 & Rep3$cfss == 1)
	Pos.excl2<-which(Rep3$charge==1 & Rep3$cf == 1 & Rep3$cfs==1)
	Pos.excl<-c(Pos.excl1,Pos.excl2)
	Rep3$rep[Pos.excl]<-"no"
	Pos.excl1<-which(Rep5$charge == 1 & Rep5$All==1 & Rep5$cf == 1 & Rep5$cfs == 1 & Rep5$cfss == 1)
	Rep5$rep[Pos.excl1]<-"no"
	Rep<-do.call(rbind,list(Rep1,Rep2,Rep3,Rep4,Rep5,Rep6))
	Rep<-Rep[order(Rep$id),]
	# Pos.charge<-names(Rep)=="charge"
	# Rep$rep[which(Rep[,Pos.charge]==0)]<-"no" #cannot replicate if not associated with PUFA in CHARGE
	Dat.list<-cbind(Dat.list,Rep$rep)
	names(Dat.list)[names(Dat.list) == "Rep$rep"]<-"rep"
	Dat.list$ARA<-FALSE
	Dat.list$ARA[grep("arachidonic",Dat.list$charge)]<-TRUE
	Dat.list$replicates<-FALSE
	Dat.list$replicates[Dat.list$rep =="replicates P<0.05/161"] <-TRUE
	# Dat.list<-Dat.list[Dat.list$replicates,]
	Dat.list<-Dat.list[,!names(Dat.list) %in% c("rep")]

	# add candidate genes and regions
	load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata")
	meta.dat<-unique(fa.reg[,c("SNP","Chr","SNP.GRCh38.p12","SNP.GRCh37","region","gene")])
	# unique(fa.reg$region[grep("ELOVL2",fa.reg$gene)])
	Dat.m<-merge(Dat.list,meta.dat,by="SNP",all.x=TRUE)
	Datr2<-Datr2[,names(Datr2)!="region"]
	Datr2.m<-merge(Datr2,meta.dat,by="SNP",all.x=TRUE)
	Dat.m<-Dat.m[order(Dat.m$region),]
	Dat.m<-Dat.m[order(Dat.m$replicates,decreasing=T),]
	
	return(list(Dat.m,Datr2.m))
}


find_discovery_study<-function(Dat=NULL){
	load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata")
	pufas<-fa.reg[fa.reg$type %in% c("n6 polyunsaturated"  ,"n3 polyunsaturated","desaturase activity","polyunsaturated","n3 or n6 polyunsaturated"),]
	pufas<-pufas[pufas$SNP %in% Dat$SNP,]
	
	pufas<-pufas[pufas$author %in% c("Shin","Kettunen"),]
	names(pufas)[names(pufas) == "pval"]<-"p"
	pufas$author
	Dat2<-rbind.fill(Dat,unique(pufas[,c("SNP","trait","author","p","consortium")]))
	return(Dat2)
}


exclude_fads_scd_elovl2<-function(Dat=NULL){
	load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata")
	regions.excl<-unlist(lapply(c("FADS","ELOVL2","SCD"),FUN=function(x) unique(fa.reg$region[grep(x,fa.reg$gene)])))
	Dat<-Dat[!Dat$region %in% regions.excl,]
	Dat<-Dat[!Dat$region %in% c(3639,3760,2129),] #also exclude these regions which are close to FADS or ELOVL2 genes. FADS specifically located in region 3761. ELOVL2 in 2149, 2150
	return(Dat)
}

find_region<-function(Dat=NULL){
	load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata")
	fa.reg<-unique(fa.reg[,c("SNP","region")])
	Dat<-merge(Dat,fa.reg,by="SNP",all.x=TRUE)
	# Dat[,c("SNP","region")]
	# Dat[,c("region.x","region.y")]
	# unique(Dat$region) %in% Dat$SNP
	return(Dat)
}


format_data<-function(Dat=NULL,type.fa=c( "n3 polyunsaturated","n6 polyunsaturated","n3 or n6 polyunsaturated")){
	# restrict to PUFAs and get rid of duplicates	
	Dat2<-Dat[Dat$type.fa %in% type.fa,]
	unique(Dat$type.fa)
	Dat2<-Dat2[!is.na(Dat2$beta),]

	Dat2$consortium[Dat2$consortium=="InCHIANTI"  ]<-"CHARGE"
	Dat2<-Dat2[!Dat2$consortium %in%  c("SCHS case","SCHS control" ),]
	Dat2<-Dat2[order(Dat2$n,decreasing=T),]
	id2<-paste(Dat2$snp,Dat2$consortium)
	Dat2<-Dat2[id2 != c("rs12806663 SCHS pooled"),]
	Dat3<-Dat2[Dat2$author == "Tanaka" & Dat2$trait =="eicosadienoic acid (20:2n6)",]
	Dat2<-Dat2[Dat2$author !="Tanaka",]
	Dat2<-rbind(Dat2,Dat3)

	Dat2$id<-paste(Dat2$consortium,Dat2$trait,Dat2$snp)
	any(duplicated(Dat2$id))
	Dat2<-Dat2[!duplicated(Dat2$id),]
	names(Dat2)[names(Dat2) == "snp"]<-"SNP"
	Dat2$other_allele[Dat2$other_allele=="TRUE"]<-"T"
	return(Dat2)
}


find_proxies<-function(Har=NULL,snps=snps_all_tab$SNP) { #proxy lookup for missing SNPs
	# unique(Har1$snp)
	# unique(Har$snp)
	load("~/fatty-acids-mr/instruments/define_fatty_acid_SNPs_v3.rdata")
	Har1<-Har[Har$SNP %in% snps,]
	Har2<-Har[!Har$SNP %in% snps,] #proxy/alias SNPs 
	Hars735665<-Har[Har$SNP=="rs735665",]
	# table(Har2$file)
	L.p<-NULL
	L.a<-NULL
	Files<-unique(Har$file)
	for(i in 1:length(Files)){
		print(i)
		print(Files[i])
		Population<-unique(Har$population[Har$file==Files[i]])
		# Population<-"East Asian"
		Pop<-c("EUR","EAS")
		Pos<-which(c("European","East Asian") == Population )
		
		snps.miss<-snps[!snps %in% Har1$SNP[Har1$file==Files[i]]]
		
		proxies<-fa.reg[fa.reg$SNP %in% snps.miss,]

		proxies<-proxies[proxies$pop == Pop[Pos],] 

		Har3<-Har2[which(Har2$file==Files[i]),]
		names(Har3)[names(Har3)=="SNP"]<-"ProxySNP"
		m.proxy<-merge(Har3,proxies[,c("SNP","ProxySNP","R2","pop","Correlated_Alleles")],by="ProxySNP")
		m.proxy<-m.proxy[!duplicated(m.proxy),] #lots of duplicate proxy SNPs. Duplicates seem to correspond to unique snp-trait associations
		names(Har3)[names(Har3)=="ProxySNP"]<-"alias_rs"
		m.alias<-merge(Har3,proxies[,c("SNP","alias_rs","R2")],by="alias_rs")
		m.alias<-m.alias[!duplicated(m.alias),]
		m.proxy<-m.proxy[order(m.proxy$SNP,m.proxy$info,m.proxy$R2,decreasing=T),]
		m.alias<-m.alias[order(m.alias$SNP,m.alias$info,m.alias$R2,decreasing=T),]
		for(j in 1:length(unique(m.proxy$SNP))){
			# m.proxy$pop
			print(j)
			print(unique(m.proxy$SNP)[j])
			Sign<-sign(m.proxy[m.proxy$SNP==unique(m.proxy$SNP)[j],c("beta")])
			# m.proxy[m.proxy$SNP==unique(m.proxy$SNP)[j],]
			# Info<-paste(m.proxy[m.proxy$SNP==unique(m.proxy$SNP)[j],c("info")],collapse="; ")
			# R2<-paste(m.proxy[m.proxy$SNP==unique(m.proxy$SNP)[j],c("R2")],collapse="; ")
			m.proxy$Sign[m.proxy$SNP==unique(m.proxy$SNP)[j]]<-all(Sign>0) | all(Sign<0)
		}
		L.p[[i]]<-m.proxy
		L.a[[i]]<-m.alias
	}
	# c("SNP","ProxySNP","info","R2","beta.pc","eaf")
	L.p<-do.call(rbind.fill,L.p)
	L.a<-do.call(rbind.fill,L.a) #no aliases
	L.p[L.p$file== "score_c182n6_pooled_allchr_qc1.tab",]
	L.p2<-L.p[!is.na(L.p$info),]
	L.p3<-L.p[is.na(L.p$info),]
	L.p2<-L.p2[order(L.p2$SNP,L.p2$info,L.p2$R2,decreasing=T),]
	L.p2<-L.p2[!duplicated(paste(L.p2$file,L.p2$SNP)),]
	L.p3<-L.p3[order(L.p3$SNP,L.p3$R2,decreasing=T),]
	L.p3<-L.p3[!duplicated(paste(L.p3$file,L.p3$SNP)),]
	proxies<-rbind(L.p3,L.p2)
	# proxies2<-L.p[!duplicated(paste(L.p$file,L.p$SNP)),]

	proxies1<-proxies[!is.na(proxies$info),]
	proxies2<-proxies[is.na(proxies$info),]
	proxies1<-proxies1[proxies1$info == 1,] #restrict proxies to genotyped SNPs. 
	proxies<-rbind(proxies1,proxies2)
	min(proxies$R2) 
	# names(Har1)[names(Har1) == "snp"]<-"SNP"
	Har1.r<-rbind.fill(Har1,proxies)


	# ##############################
	# Fix Proxy SNP effect alleles###
	#################################
	# Recode effect and other allele for Proxy SNPs so that it reflects the effect alleles of the target SNP. I checked dbSNP for several SNPs to make sure the proxy and target alleles were consistent with respect to effect allele frequency and the Correlated_Alleles columns. I did not find any inconsistencies, indicating that I can reliably switch the effect alleles to the target alleles using the correlated_alleles column, i.e. the alleles need to be switched to the other alleles with the correlated alleles column. 

	proxies<-unique(Har1.r$ProxySNP)
	proxies<-proxies[!is.na(proxies)]
	Pos<-which(!is.na(Har1.r$ProxySNP))
	Har1.r$id.snproxy<-NA
	Har1.r$id.snproxy[Pos]<-paste(Har1.r$file[Pos],Har1.r$SNP[Pos],Har1.r$ProxySNP[Pos])

	# need to use a different identifier because alleles for proxy SNPs sometimes different betwee studies. Need to include study in identifier or pmid or file.  
	snproxies<-unique(Har1.r$id.snproxy[!is.na(Har1.r$id.snproxy)])
	# i<-which(snproxies=="score_c226n3_case_allchr_qc.txt rs2727271 rs2727270")
	L.proxy<-NULL
	i<-which(snproxies== "score_c182n6_pooled_allchr_qc1.tab rs11604424 rs66505542")
	for(i in 1:length(snproxies)){
		print(i)
		print(snproxies[i])
		# snproxies[grep("rs66505542",snproxies)]
		ProxyDat<- unique(Har1.r[which(Har1.r$id.snproxy == snproxies[i]),c("id.snproxy","SNP","ProxySNP","effect_allele","other_allele","Correlated_Alleles")])
		
		unique(Har1.r[Har1.r$SNP=="rs2727271",c("effect_allele","other_allele","effect_allele_freq")])
		Correlated_Alleles<-ProxyDat$Correlated_Alleles
		alleles<-unlist(strsplit(Correlated_Alleles,split=","))
		allele1<-alleles[1]
		allele1_1<-unlist(strsplit(allele1,split="="))[1]
		allele1_2<-unlist(strsplit(allele1,split="="))[2]
		allele2<-alleles[2]
		allele2_1<-unlist(strsplit(allele2,split="="))[1]
		allele2_2<-unlist(strsplit(allele2,split="="))[2]
		alleleA<-c(allele1_1,allele2_1)
		alleleB<-c(allele1_2,allele2_2)
		ea<-ProxyDat$effect_allele
		oa<-ProxyDat$other_allele
		
		if(ea %in% alleleA &  oa %in% alleleA){ #ea and oa are in alleleA 
			print(ea %in% alleleA &  oa %in% alleleA)
			ProxyDat$ea_new<-alleleB[match(ea,alleleA)]
			ProxyDat$oa_new<-alleleB[match(oa,alleleA)]
			L.proxy[[i]]<-ProxyDat
		}

		if(ea %in% alleleB &  oa %in% alleleB) {#ea and oa are in alleleB 
			print(ea %in% alleleB &  oa %in% alleleB)
			ProxyDat$ea_new<-alleleA[match(ea,alleleB)]
			ProxyDat$oa_new<-alleleA[match(oa,alleleB)]
			L.proxy[[i]]<-ProxyDat
		}	
	}

	ProxyDat<-do.call(rbind,L.proxy)
	names(ProxyDat)[names(ProxyDat) == "effect_allele"]<-"ea"
	names(ProxyDat)[names(ProxyDat) == "other_allele"]<-"oa"
	ProxyDat<-ProxyDat[,c("id.snproxy","ea_new","oa_new")]
	Har1.r.proxy<-Har1.r[!is.na(Har1.r$id.snproxy),]
	Har1.r.Xproxy<-Har1.r[is.na(Har1.r$id.snproxy),]
	Har1.r.m<-merge(Har1.r.proxy,ProxyDat,by="id.snproxy")
	Har1.r.m$effect_allele<-Har1.r.m$ea_new
	Har1.r.m$other_allele<-Har1.r.m$oa_new
	Har1.r.m<-Har1.r.m[,!names(Har1.r.m) %in% c("ea_new","oa_new")]
	# HarR.schs.m[HarR.schs.m$SNP=="rs2727271",c("SNP","ProxySNP","effect_allele","other_allele","Correlated_Alleles","file")]
	HarR.schs.r<-rbind(Har1.r.Xproxy,Har1.r.m)
	HarR.schs.r2<-rbind.fill(HarR.schs.r,Hars735665)
	# Some proxies missing because of deletion/insertions in correlated alleles
	# Temp<-unique(HarR.schs.r$ProxySNP)
	# Temp<-Temp[!is.na(Temp)] 
	# proxies[!proxies %in% Temp]
	return(HarR.schs.r2)
}



get_meta_data<-function(Dat){
	# had to replace some weird character but only works on epi franklin and not blue crystal 3. 
	# fat<-read.table("~/fatty-acids/mr/data/fatty_acid_GWASanalysis_table2.txt",sep="\t",stringsAsFactors=F,head=T)

	# # test<-fat$filename[1:10]
	# fat$filename<-gsub(".txt",".tab",fat$filename)
	# fat$filename<-gsub("pooled_allchr_qc","pooled_allchr_qc1",fat$filename)
	# fat$filename<-gsub("\xca","",fat$filename)
	# fat$filename<-trimws(fat$filename)
	# fat$filename<-gsub("tab.gz.tab","tab",fat$filename)
	# fat$filename<-gsub("tab.tab","tab",fat$filename)
	# write.table(fat,"~/fatty-acids-mr/instruments/fatty_acid_GWASanalysis_table_cleaned.txt",sep="\t",quote=F,col.names=T,row.names=F)
	res.dat<-Dat
	fat<-read.table("~/fatty-acids-mr/instruments/fatty_acid_GWASanalysis_table_cleaned.txt",sep="\t",stringsAsFactors=F,head=T)
	res.dat$file<-gsub(".txt.tab",".tab",res.dat$file)
	res.dat$file<-gsub("pos.txt.gz.tab","pos.tab",res.dat$file)
	fat$filename[!fat$filename %in% res.dat$file]
	res.dat$file[!res.dat$file %in% fat$filename]
	
	# fat$filename[grep("metal.pos",fat$filename)]
	# res.dat$file[grep("metal.pos",res.dat$file)]

	names(fat)[names(fat) == "type"]<-"type.fa"
	Dat<-merge(res.dat,fat[,c("pmid","author","consortium","trait","type.fa","chain","chain.length","SigSNPs","	sample_size.analysis","population","sex","year","filename")],by.x="file",by.y="filename",all.x=TRUE)
	return(Dat)
}


list_of_snps<-function(){
	miss<-teqtl[!teqtl %in% snps_all_tab$SNP]
	Pos<-unlist(lapply(1:length(miss), FUN=function(x) grep(miss[x],snps_all_tab$ProxySNP)))
	# snps_all_tab[Pos,c("SNP","ProxySNP")]

	snp<-snps_all_tab$SNP
	proxy<-unique(unlist(strsplit(snps_all_tab$ProxySNP,split="; ")))
	proxy<-proxy[proxy!=""]
	extra_teqtl<-miss[1] #this SNP is a trans eQTL for FADS1 but is not associated with PUFAs

	snplist<-unique(c(snp,proxy,extra_teqtl))
	return(snplist)
}


extract_snps2<-function(snps=NULL){
	Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/charge/hg19/imputed"
	Files<-dir(Dir)
	Files<-Files[grep(".tab",Files)]
	Files1<-paste(Dir,"/",Files,sep="")

	Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/dorajoo/beta_pc"
	Files<-dir(Dir)
	Files<-Files[grep(".tab",Files)]
	Files2<-paste(Dir,"/",Files,sep="")

	Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tanaka/hg19/imputed"
	Files<-dir(Dir)
	Files<-Files[grep(".tab",Files)]
	Files3<-paste(Dir,"/",Files,sep="")
	Files<-c(Files1,Files2,Files3)

	Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/tintle/hg19/imputed"
	Files<-dir(Dir)
	Files<-Files[grep(".tab",Files)]
	Files4<-paste(Dir,"/",Files,sep="")

	Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/ket/hg19"
	Files<-dir(Dir)
	Files<-Files[grep(".tab",Files)]
	Files5<-paste(Dir,"/",Files,sep="")

	Dir<-"/newprojects/mrcieu/research/data/MR/fatty_acids/data/shin/hg19/imputed"
	Files<-dir(Dir)
	Files<-Files[grep(".tab",Files)]
	Files6<-paste(Dir,"/",Files,sep="")

	Files<-c(Files1,Files2,Files3,Files4,Files5,Files6)
	Res_list<-NULL
	for(i in 1:length(Files)){
		print(i)
		Res_list[[i]]<-extract_snps(snplist=snps,File=Files[i],exact_match=TRUE,file_sep="\t")
	}
	return(Res_list)
}

b_sd<-function(z,maf,n){
    sqrt(((z ^ 2) / (z ^ 2 + n - 2)) /(2 * maf * (1 - maf)))* sign(z)
}


sumlog <-function(p) {
   keep <- (p > 0) & (p <= 1)
   invalid <- sum(1L * keep) < 2
   if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(chisq = NA_real_, df = NA_integer_,
         p = NA_real_, validp = p[keep])
   } else {
      lnp <- log(p[keep])
      chisq <- (-2) * sum(lnp)
      df <- 2 * length(lnp)
      if(length(lnp) != length(p)) {
         warning("Some studies omitted")
      }
      res <- list(chisq = chisq, df = df,
         p = pchisq(chisq, df, lower.tail = FALSE), validp = p[keep])
    }
   class(res) <- c("sumlog", "metap")
   res
}
print.sumlog <- function(x, ...) {
   cat("chisq = ", x$chisq, " with df = ", x$df, " p = ", x$p, "\n")
   invisible(x)
}


harmonise_alleles<-function(HarR=NULL){
	# Harmonise effect alleles across SNPs so that effect allele is the minor allele in CHARGE. If SNP not present in CHARGE then recode to minor allele in Framingham If not present in Framingham use KEttunen. If not present in Kettunen use TwinsUK/KORA
	
	HarR1<-make_ea_ma(Dat=HarR,study="TwinsUK/KORA")
	HarR2 <-make_ea_ma(Dat=HarR1,study="14 cohorts: EGCUT,ERF,FTC,FR97,COROGENE,GenMets,HBCS,KORA,LLS,NTR,NFBC1966,PredictCVD,PROTE,YFS)")
	HarR3 <-make_ea_ma(Dat=HarR2,study="FraminghamHOS") #function converts effect allele to minor allele in specified study
	HarR4 <-make_ea_ma(Dat=HarR3,study="CHARGE")
	HarR4[which(HarR4$effect_allele_freq>0.5),]
	unique(HarR4$consortium)
	

	# # HarR.schs.r2 <-make_ea_ma(dat=HarR,study="FraminghamHOS") #function converts effect allele to minor allele in specified study
	# HarR.schs.r3 <-make_ea_ma(dat=HarR.schs.r2,study="TwinsUK/KORA")
	# HarR.schs.r4 <-make_ea_ma(dat=HarR.schs.r3,study="Kettunen")
	# HarR.schs.r5 <-make_ea_ma(dat=HarR.schs.r4,study="CHARGE")

	IDs<-unique(HarR4$SNP)
	L<-NULL
	# i<-2
	for(i in 1:length(IDs)){
		print(i)
		print(IDs[i])
		Dat<-HarR4[HarR4$SNP == IDs[i],] 
		print(Dat$consortium[1])
		print(Dat$effect_allele_freq[1])
		ea<-Dat$effect_allele
		oa<-Dat$other_allele
		if(any(Dat$effect_allele!=ea[1] & Dat$effect_allele!=oa[1])){ #SNPs on different strands between studies
			stop("SNPs on different strands")
		}
		Pos<-which(Dat$effect_allele!=ea[1]) #recode effect allele to reflect effect allele in first row, which corresponds to CHARGE consortium with effect allele always the minor allele
		Dat$effect_allele[Pos]<-oa[Pos]
		Dat$other_allele[Pos]<-ea[Pos]
		Eaf<-1-Dat$effect_allele_freq[Pos]
		Dat$effect_allele_freq[Pos]<-Eaf
		Beta<-Dat$beta[Pos]*-1
		Dat$beta[Pos]<-Beta
		# print(paste("number of unique effect alleles=",length(unique(Dat$effect_allele)),sep=""))
		L[[i]]<-Dat
	}
	HarEA<-do.call(rbind.fill,L)
	snps1<-unique(HarEA$SNP)
	if(all(unlist(lapply(1:length(snps1),FUN=function(x)
		length(unique(HarEA$effect_allele[HarEA$SNP == snps1[x]]))!=1)))) stop("effect alleles not same across studies")
	if(all(unlist(lapply(1:length(snps1),FUN=function(x)
		length(unique(HarEA$other_allele[HarEA$SNP == snps1[x]]))!=1)))) stop("other alleles not same across studies")

	# effect allele frequency is also always the minor allele except in the SCHS which is not surprising because it is a study from Chinese in Singapore. This means that for the European origin studies effect alleles are consistent across studies even for palindromic SNPs. For SCHS it is impossible to be certain that effect alleles are the same as European studies for palindromic SNPs. Only by confirming reference strand could we do this. 

	if(any(HarEA$effect_allele_freq>0.5)){
		if(unique(HarEA$consortium[HarEA$effect_allele_freq>0.5]) != "SCHS pooled") stop("minor allele >0.5 in non SCHS studies")
	}

	return(HarEA)
}

make_ea_ma<-function(Dat=NULL,study=NULL){
	dat1<-Dat[Dat$consortium == study,]
	maf<-dat1$effect_allele_freq
	Pos.change<-which(maf>0.5)
	Beta<-dat1$beta[Pos.change]*-1
	dat1$beta[Pos.change]<-Beta
	Eaf<-1-dat1$effect_allele_freq[Pos.change]
	dat1$effect_allele_freq[Pos.change]<-Eaf
	EA<-dat1$effect_allele[Pos.change]
	OA<-dat1$other_allele[Pos.change]
	dat1$effect_allele[Pos.change]<-OA
	dat1$other_allele[Pos.change]<-EA
	dat2<-Dat[Dat$consortium != study,]
	dat3 <- rbind(dat1,dat2) 
	return(dat3)
}
