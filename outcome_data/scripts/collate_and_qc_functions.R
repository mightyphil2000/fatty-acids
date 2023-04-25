
compute_nstudies<-function(dat=NULL){
	Dir<-dat$Direction	
	Dir2<-gsub("\\?","",Dir)
	dat$Nstudies<-nchar(Dir2)
	IDS<-unique(dat$ID)
	for(i in 1:length(IDS)){
		dat$Nstudies_median[dat$ID == IDS[i]]<-median(dat$Nstudies[dat$ID == IDS[i] ])		
	}
	return(dat)
}
	
	
basic_qc<-function(dat=NULL){ 
	# unique(dat$ID[is.na(dat$Other.Allele)])
	# ID <- unique(dat$ID[is.na(dat$mean_info)])
	# ID2 <-unique(dat$ID[!is.na(dat$mean_info)])
	# ID <- unique(dat$ID[is.na(dat$HWEp)])
	# ID2 <-unique(dat$ID[!is.na(dat$HWEp)])
	# ID <- unique(dat$ID[is.na(dat$phet)])
	# ID2 <-unique(dat$ID[!is.na(dat$phet)])
	# which(dat$phet<0.001)
	# which(dat$HWEp<0.001)
	# ID2 <-unique(dat$ID[!is.na(dat$phet)])
	# ID<-ID[!ID %in% ID2]
	# length(ID)
	# head(dat)
	# head(dat[dat$ID == 165,])

	# ID2<-c(1:165,c(967,1499))
	# ID2[!ID2 %in% ID]
	# load("~/fatty-acids/outcome_data/data/harmonised_data.Rdata")

	# dat<-prune_info80(Dat=dat)
	# ID1<-unique(dat1$ID)
	# ID2<-unique(dat2$ID)
	# length(ID1)
	# ID2<-ID2[!ID2 %in% ID1]	
	dat1<-dat[!is.na(dat$eaf),]
	dat2<-dat[is.na(dat$eaf),]

	
	MAF<-dat1$eaf
	MAF[MAF>0.5]<-1-MAF[MAF>0.5]
	dat1$MAC_case<-MAF*dat1$ncase*2
	dat1$MAC_control<-MAF*dat1$ncontrol*2

	id1<-unique(dat1$ID[which(dat1$MAC_case<50)]) #which datasets have counts<50
	id2<-unique(dat1$ID[which(dat1$MAC_case>=50)]) #which datasets have count>=50
	length(id2)-length(id1) #n datasets where all SNPs have mac >=50
	id2[!id2 %in% id1] #dataset where all SNPs have mac>=50
	id1[!id1 %in% id2] #dataset where all SNPs have mac<50

	dat1<-dat1[dat1$MAC_case>=50 & dat1$MAC_control>=50,]
	# dat1<-dat1[which(dat1$eaf > 0.001),]
	dat2$MAC_case<-NA
	dat2$MAC_control<-NA
	dat<-rbind(dat1,dat2)

	
	# Pos1<-which(dat$lnor<1 & dat$lnor > -1 )
	# dat<-dat[Pos1,]
	# unique(dat4$study[Pos])
	dat<-dat[order(dat$ID),]
	dat$z<-dat$lnor/dat$se
	return(dat)
}



collate_dat<-function(postqc=TRUE){
	if(postqc){
		setwd("~/fatty-acids/outcome_data/data/harmonised/")
		
	}
	if(!postqc){
		setwd("~/fatty-acids/outcome_data/data/harmonised_preqc/")
	}

	Files<-dir()
	Dat_list<-NULL
	for(i in 1:length(Files)){
		print(i)
		print(Files[i])
		Dat<-read.table(Files[i],sep="\t",stringsAsFactors=F,head=T,quote="")			
		Dat$file.outcome<-Files[i]
		Dat_list[[i]]<-Dat
	}

	Dat<-do.call(plyr::rbind.fill,Dat_list)

	# unique(Dat[,names(Dat)[grep("file",names(Dat),ignore.case=T)]])
	col_names_keep<-c("file.outcome","outcome","population","pmid","study","ncase","ncontrol","UKbiobank","rsid","Effect.Allele","Other.Allele","lnor","se","eaf","p","info","info1","info2","info3","HWEp","phet","I2","Q","Direction","effect_allele_confirmed","ID","test_statistic","all_summary_stats","summary_set","open_gwas","efo")
	# Dat$open_gwas[is.na(Dat$open_gwas)]<-FALSE
	Dat<-Dat[,col_names_keep]
	
	# set open_gwas to TRUE for all UKB datasets 
	Dat$open_gwas[which(Dat$study == "UKB")]<-TRUE
	Dat$open_gwas[which(Dat$ID == 163)]<-FALSE 
	Dat$open_gwas[which(Dat$ID == 73)]<-FALSE #mistakenly coded as in open gwas
	Dat$open_gwas[which(Dat$ID==128)]<-TRUE #mistakenly coded as FALSE for Open GWAS
	# datasets obtained via correspondence
	unique(Dat$ID[which(Dat$open_gwas)])

	length(unique(Dat$study[which(Dat$open_gwas)]))
	Dat$correspondence<-FALSE
	Dat$correspondence[Dat$ID %in% c(1,2,3,4,5,21,22,23,24,25,26,27,57,59,60,61,62,63,64,65,66,67,68,69,70,71,73,80,81,82,83,84,85,86,87,89,94,95,96,96.2,97,98,99,100,101,102,103,104,105,107,121,123,124,125,128,129,130,132,133,134,165)]<-TRUE

	sort(unique(Dat$pmid[Dat$correspondence]))

	# datasets downloaded from the GWAS catalog
	Dat$gwas_catalog<-FALSE
	Dat$gwas_catalog[Dat$ID %in% c(88,90,91,92,120,163)]<-TRUE

	# GWAS consortia in Open GWAS
	Consortia<-unique(Dat[which(Dat$open_gwas & !Dat$study %in% c("BJ","UKB","FinnGen","UKB")),c("study")])
	Consortia<-Consortia[!Consortia %in% c("GCS","NBS","TCS")]
	Dat$consortium<-FALSE
	Dat$consortium[Dat$study %in% Consortia]<-TRUE
	return(Dat)
}



prune_info80<-function(Dat=NULL){
	Dat$info<-as.numeric(Dat$info)
	Dat$info1<-as.numeric(Dat$info1)
	Dat$info2<-as.numeric(Dat$info2)
	Dat$info3<-as.numeric(Dat$info3)
	
	L<-NULL
	for(i in 1:nrow(Dat)){
		# print(i)
		L[[i]]<-mean(c(Dat$info[i],Dat$info1[i],Dat$info2[i],Dat$info3[i]), na.rm=T)
	}

	Dat$mean_info<-unlist(L)
	
	# drop rows with info score <0.80
	Dat1<-Dat[which(Dat$mean_info>=0.80), ]
	Dat2<-Dat[is.na(Dat$mean_info), ]
	Dat<-rbind(Dat1,Dat2)
	return(Dat)
}
# library(cowplot)
# library(rlang)

# legend<-cowplot::get_legend(Plot)
# Plot2<-Plot
# Plot_list<-Plot1
# dev.off()
make_cow_plot<-function(Plot_list=NULL,Title="",Xlab="",Ylab="",out_file=NULL,return_plot=FALSE,width=1000,height=1000,Title_size=0,Title_axis_size=10,Subtitle="",Subtitle_size=0){
	
	Plot<-cowplot::plot_grid(plotlist=Plot_list)

	if(Title!="") { 
		title <- cowplot::ggdraw() + 
				cowplot::draw_label(
					Title,
					# fontface = 'bold',
					fontface = 'plain',
					x = 0,
					hjust = 0,
					size=Title_size)  +
				ggplot2::theme(
				# add margin on the left of the drawing canvas,
				# so title is aligned with left edge of first plot
					plot.margin = ggplot2::margin(0, 0, 0, 7)
					)

		subtitle <- cowplot::ggdraw() + 
		cowplot::draw_label(
			Subtitle,
			# fontface = 'bold',
			fontface = 'plain',
			x = 0,
			hjust = 0,
			# element = "plot.subtitle",
			size=Subtitle_size)  +
		ggplot2::theme(
		# add margin on the left of the drawing canvas,
		# so title is aligned with left edge of first plot
			plot.margin = ggplot2::margin(0, 0, 0, 7)
			)
		# subtitle <- ggdraw() +
  # 						draw_label_theme("By census tract, 2016",
  #                  theme = theme_georgia(), 
  #                  element = "plot.subtitle",

  #                  x = 0.05, hjust = 0, vjust = 1)

		Plot<-cowplot::plot_grid(title,subtitle, Plot,ncol = 1,rel_heights = c(0.05,0.05, 1))
	}
	y.grob <- grid::textGrob(Ylab, 
	                   gp=grid::gpar(fontface="bold", col="black", fontsize=Title_axis_size), rot=90)

	x.grob <- grid::textGrob(Xlab, 
	                   gp=grid::gpar(fontface="bold", col="black", fontsize=Title_axis_size))

	if(!return_plot){
		png(out_file, width = width, height = height)
			gridExtra::grid.arrange(gridExtra::arrangeGrob(Plot, left = y.grob, bottom = x.grob))
		dev.off()	
	}
	if(return_plot){
		return(Plot)
	}

}

# 	Plot<-cowplot::plot_grid(plotlist=Plot_list)

# 	space to the left of the legend
# plot_grid(Plot, legend, rel_widths = c(3, .4))


make_cow_plot2<-function(Plot_list=NULL,Title="",Xlab="",Ylab="",out_file=NULL,return_plot=FALSE,width=1000,height=1000,Title_size=0,Title_axis_size=10,bycols=TRUE,Legend_z=FALSE,Legend_eaf=FALSE,Tiff=FALSE){

	# Plot<-cowplot::plot_grid(plotlist=Plot_list[[1]])
	if(bycols){
		Plot<-cowplot::plot_grid(plotlist=Plot_list,nrow=3,ncol=2)
	}
	if(!bycols){
		Plot<-cowplot::plot_grid(plotlist=Plot_list)
	}
	if(Title!="") { 
		title <- cowplot::ggdraw() + 
				cowplot::draw_label(
					Title,
					fontface = 'bold',
					# fontface = 'plain',
					x = 0,
					hjust = 0,
					size=Title_size)  +
				ggplot2::theme(
				# add margin on the left of the drawing canvas,
				# so title is aligned with left edge of first plot
					plot.margin = ggplot2::margin(0, 0, 0, 7)
					)

		Plot<-cowplot::plot_grid(title, Plot,ncol = 1,rel_heights = c(0.05, 1))
	}

# Left -> legend around the plot
# basic + theme(legend.position = "bottom")
# # Right -> inside the plot area
# basic + theme(
#     legend.position = c(.95, .95),
#     legend.justification = c("right", "top"),
#     legend.box.just = "right",
#     legend.margin = margin(6, 6, 6, 6)
#     )

	if(Legend_z){ #use BCAC or OCAC study which has best legend
		# i<-which(IDS==6)
		i<-which(IDS==118)
		Dat1<-Dat[Dat$ID == IDS[i],]	
		legend_plot<-make_plot_gwas_catalog_zscores(dat=Dat1,efo=unique(Dat1$efo),Title_size_subplot=12,Title=Title,Ylab="",Xlab="",Title_xaxis_size=0,legend=TRUE)
		Legend <- get_legend(legend_plot + theme(legend.box.margin = margin(0, 0, 0, 12)))# create some
		Plot<-plot_grid(Plot, Legend, rel_widths = c(3, .4))

		# Legend <- get_legend(legend_plot + 
		# 	theme(
		# 	    legend.position = c(.30,1.0), #(right,top)			    
		# 	    legend.justification = c("left", "top"),
		# 	    legend.box.just = "left",
		# 	    legend.margin = margin(0,-2,0,0,unit="cm"),			    
		# 	    ))
		# Plot<-plot_grid(Plot, Legend, rel_widths = c(1, .05))
	}

	if(Legend_eaf){ #use BCAC or OCAC study which has best legend
		# i<-which(IDS==6)
		i<-which(IDS==118)
		Dat1<-Dat[Dat$ID == IDS[i],]	
		legend_plot<-make_plot_gwas_catalog_eaf(dat=Dat1,efo=unique(Dat1$efo),Title_size_subplot=12,Title=Title,Ylab="",Xlab="",Title_xaxis_size=0,legend=TRUE)
		Legend <- get_legend(legend_plot + theme(legend.box.margin = margin(0, 0, 0, 12)))# create some
		Plot<-plot_grid(Plot, Legend, rel_widths = c(3, .4))
	}

	y.grob <- textGrob(Ylab, 
	                   gp=gpar(fontface="bold", col="black", fontsize=Title_axis_size), rot=90)

	x.grob <- textGrob(Xlab, 
	                   gp=gpar(fontface="bold", col="black", fontsize=Title_axis_size))

	if(!return_plot){
		if(!Tiff){
			png(out_file, width = width, height = height,)
				grid.arrange(arrangeGrob(Plot, left = y.grob, bottom = x.grob))
			dev.off()	
		}
		if(Tiff){
			tiff(out_file, width = width, height = height,)
				grid.arrange(arrangeGrob(Plot, left = y.grob, bottom = x.grob))
			dev.off()	
		}
	}
	if(return_plot){
		return(Plot)
	}

}



make_plot_lists<-function(Dat=NULL,refstudy1=NULL,refstudy2=NULL){
	IDS<-unique(Dat$ID)
	Plot_list1<-NULL
	Plot_list2<-NULL
	# East Asians with full GWAS data
	for(i in 1:length(IDS)){
		Dat1<-Dat[Dat$ID==IDS[i],]
		Plot_list1[[i]]<-find_allele_errors_wrapper(dat=Dat1,refstudy=refstudy1)		
		Plot_list2[[i]]<-find_allele_errors_wrapper(dat=Dat1,refstudy=refstudy2)	
	}
	return(list(Plot_list1,Plot_list2))
}


make_plot_lists2<-function(Dat=NULL,refstudy=NULL,subtitle_off=FALSE,Title_subx_axis_size=12,Title_size_subplot=8,Ylab="",Xlab=""){
	IDS<-unique(Dat$ID)
	Plot_list1<-NULL
	for(i in 1:length(IDS)){
		Dat1<-Dat[Dat$ID==IDS[i],]
		Plot_list1[[i]]<-find_allele_errors_wrapper(dat=Dat1,refstudy=refstudy,subtitle_off=subtitle_off,Title_subx_axis_size=Title_subx_axis_size,Title_size_subplot=Title_size_subplot)		
	}
	return(Plot_list1)
}


find_allele_errors_wrapper<-function(dat=NULL,refstudy=NULL,Title_subx_axis_size=Title_subx_axis_size,Title_size_subplot=Title_size_subplot,Title="",Ylab=Ylab,Xlab=Xlab){
	# Pos<-which(c("European","East Asian") == population) 
	# population<-c("EUR","EAS")[Pos]
	# refstudy="CHARGE"
	ref.dat<-load_refdata(refstudy=refstudy)

	Plot<-find_allele_error(target_dat=dat,ref_dat=ref.dat,snp="rsid",eaf="eaf",ref_dat_maf="maf",target_dat_effect_allele = "Effect.Allele",target_dat_other_allele= "Other.Allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",target_study="study",ID="ID",ref_study="study",Title_subx_axis_size=Title_subx_axis_size,Title_size_subplot=Title_size_subplot,Title=Title,Ylab=Ylab,Xlab=Xlab)	
	# out_file<-unique(paste0(refstudy,"_",dat$study,"_",dat$ID,".png"))
	# out_file<-unlist(Dat_list[3])
	# out_file<-gsub("/","_",	out_file)

	# if(population=="superpops"){
	# 	png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/superpops_",out_file,".png"))
	# 		print(Dat_list[1])
	# 	dev.off()	
	# }

	# if(population!="superpops"){
	# 	png(paste0("~/fatty-acids/outcome_data/results/plots/eaf_vs_mafref/fasnps_",out_file,".png"))
	# 		print(Dat_list[1])
	# 	dev.off()	
	# }
	
	# allele_errors<-unlist(Dat_list[2])		
	# write.table(allele_errors,paste0("~/fatty-acids/outcome_data/results/snps_allele_errors_",out_file,".txt"),col.names=T,row.names=F,quote=F,sep="\t")
	# ,Dat_list[2])
	return(Plot)
}
# eaf="eaf";snp="rsid";ref_dat_maf="maf";target_dat_effect_allele="Effect.Allele";target_dat_other_allele="Other.Allele";ref_dat_minor_allele="minor_allele";ref_dat_major_allele="major_allele";outcome="outcome";ID="ID";target_dat_population="population";ref_dat_population="population";target_study="study";ref_study="study";

# should exclude palindromic SNPs which can cause apparent conflicts when target and reference datasets on different strands for some SNPs. drop palindromic SNPs or show in different shape?
make_plot_maf<-function(refstudy=NULL,target_dat=NULL,eaf="eaf",snp="rsid",ref_dat_maf="maf",target_dat_effect_allele="Effect.Allele",target_dat_other_allele="Other.Allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",ID=NULL,target_dat_population="population",ref_dat_population="population",target_study="study",ref_study="study",Title_xaxis_size=10,Title_size=10,Title="",Ylab="",Xlab="",cowplot_title=NULL,Dir="~/fatty-acids/outcome_data/data/"){	

	ref_dat<-load_refdata(refstudy=refstudy,Dir=Dir)

	names(target_dat)[names(target_dat) == target_dat_population]<-"target_dat_population"
	names(ref_dat)[names(ref_dat) == ref_dat_population]<-"ref_dat_population"
	names(target_dat)[names(target_dat) == target_study]<-"target_study"
	names(ref_dat)[names(ref_dat) == ref_study]<-"ref_study"


	if(any(names(ref_dat) %in% c(target_dat_effect_allele,target_dat_other_allele,target_dat_effect_allele))) warning("effect allele, other allele or eaf present in refererence dataset with same name as in target dataset")
	

	dat.m<-merge(ref_dat,target_dat,by=snp)
	Pos<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])
	
	# harmonise study to ref dataset minor allele. need to clean this up as a new function with flip_strand and harmonise_allele functions like in the GWAS catalog functions
	dat.m[,eaf]<-as.numeric(dat.m[,eaf])
	dat.m[,ref_dat_maf]<-as.numeric(dat.m[,ref_dat_maf])
	dat.m[,eaf][Pos]<-1-dat.m[,eaf][Pos]
	EA<-dat.m[,target_dat_effect_allele][Pos]
	OA<-dat.m[,target_dat_other_allele][Pos]
	dat.m[,target_dat_effect_allele][Pos]<-OA
	dat.m[,target_dat_other_allele][Pos]<-EA
	
	# harmonise SNPs on different strands
	Pos1<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])
	Pos2<-which(dat.m[,target_dat_effect_allele] == dat.m[,ref_dat_minor_allele])
	dat.m1<-dat.m[Pos1,]
	dat.m2<-dat.m[Pos2,]

	
	# dat.m1[,c("Effect.Allele","Other.Allele","minor_allele2","major_allele2")]
	Pos<-which(dat.m1[,target_dat_effect_allele] != dat.m1[,"minor_allele2"])
	# dat.m1[Pos,c("Effect.Allele","Other.Allele","minor_allele2","major_allele2")]
	dat.m1[,eaf][Pos]<-1-dat.m1[,eaf][Pos]
	EA<-dat.m1[,target_dat_effect_allele][Pos]
	OA<-dat.m1[,target_dat_other_allele][Pos]
	dat.m1[,target_dat_effect_allele][Pos]<-OA
	dat.m1[,target_dat_other_allele][Pos]<-EA

	dat.m<-rbind(dat.m1,dat.m2)

	# dat.m[,c("Effect.Allele","Other.Allele","minor_allele","major_allele","eaf","maf")]

# Pos<-which(dat.m1$Effect.Allele  != dat.m1$minor_allele2 &  dat.m1$Effect.Allele  != dat.m1$minor_allele)
# dat.m1[Pos,c("Effect.Allele","Other.Allele","minor_allele","minor_allele2")]
	dat.m$alleles<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	dat.m$alleles2<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	
	dat.m.test<-dat.m
	outcome_plot<-outcome
	if(is.null(outcome)){
		outcome_plot<-""
	}

	if(!is.null(outcome) & !is.null(target_study) & !is.null(ID)){
		outfile_name<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID=",dat.m.test[,ID]))
		outfile_name<-gsub("\\|","",outfile_name)
		outfile_name<-gsub(" ","_",outfile_name)
		outfile_name<-gsub("__","_",outfile_name)
		outfile_name<-gsub("=","",outfile_name)
		outfile_name<-gsub("/","_",outfile_name)
		outcome_plot<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		target_study<-unique(paste0(dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		# outcome_plot2<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test[,study]," | ID: ",dat.m.test[,ID]))
	}

	Plot_list<-NULL
	dat.m.test<-dat.m.test[order(dat.m.test$ref_dat_population),]
	Pops<-unique(dat.m.test$ref_dat_population)
	# if(length(Pops)>1){
	for(i in 1:length(Pops)){
		# print(pop)
		# dat1<-dat[dat$ref_dat_population==pop,]			
		dat1<-dat.m.test[dat.m.test$ref_dat_population==Pops[i], ]
		pop2<-c("European","East Asian","African","American","South Asian","Global","European","European")
		Pops2<-c("EUR","EAS","AFR","AMR","SAS","ALL","EUR2","EUR1")
		j<-which(Pops2 %in% Pops[i])
		# if(Xlab==""){
		# Xlab<-paste0(pop2[j]," MAF")
		Xlab<-pop2[j]

		# }
		# Title<-pop2[i]
		# Title<-gsub("ALL","Global pop",Title)
		# if(Title == ""){
		# 	Title<-target_study
		# }

		# if(subtitle_off) Title<-""

		Colour<-rep("black",nrow(dat1))
		Colour[which(dat1[,eaf]>0.5)]<-"blue"
		# Colour[dat1[,target_dat_effect_allele]!=dat1[,ref_dat_minor_allele]]<-"red"	
		
		# fix harmonisation functions above so that efffect allele strand flipped to minor_allele (not harmonised with minor_allele2)
		Colour[dat1$Effect.Allele!=dat1$minor_allele & dat1$Effect.Allele!=dat1$minor_allele2]<-"red"	

		Diff<-abs(dat1[,eaf]-dat1[,ref_dat_maf])
		Colour[which(Diff>0.10)]<-"red"
		Shape<-rep(19,nrow(dat1))
		Shape[which(dat1$alleles %in% c("AT","TA","GC","CG"))]<-1
		dat1$eaf<-dat1[,eaf]
		dat1$maf<-dat1[,ref_dat_maf]

		Title_size1<-Title_size
		Subtitle_size1<-8
		if(length(unique(dat.m.test$ref_dat_population)) > 1){
			Title_size1<-0
			Subtitle_size1<-0
		}

		
		# Temp<-dat1[dat1$eaf > 0.5,c("Effect.Allele","minor_allele","Other.Allele","major_allele","eaf","maf","alleles")]
		# dim(Temp)
		# length(which(Temp$alleles %in% c("CG","GC","TA","AT")))

		# snps<-dat1$rsid[dat1$eaf>0.55]
		# dat1[dat1$eaf>0.55,c("minor_allele","major_allele","Effect.Allele","Other.Allele","maf","eaf")]
		# dat2<-dat[dat$rsid %in% snps, ]
		# dat2[,c("rsid","Effect.Allele","Other.Allele","eaf")]
		# head(ref_dat)
		# dat2.m<-merge(dat2,ref_dat,by="rsid")
		# head(dat2.m[,c("rsid","Effect.Allele","Other.Allele","eaf", "maf", "minor_allele2","major_allele2" )])
		
		Subtitle<-unique(paste0("Reported population: ",dat1$target_dat_population))
		

		Plot<-ggplot2::ggplot(dat1, ggplot2::aes(x=maf, y=eaf)) + ggplot2::geom_point(colour=Colour) +ggplot2::ggtitle(Title) +ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size1,hjust = 0))+ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = Subtitle_size1))

			if(length(Pops)> 1){
				Plot_list[[i]]<-Plot
			}
		}	

		if(length(unique(dat.m.test$ref_dat_population)) > 1){
			if(is.null(cowplot_title)){
				cowplot_title<-target_study
			}	
			Plot<-make_cow_plot(Plot_list=Plot_list,Title=cowplot_title,Xlab="",Ylab="",return_plot=TRUE,Title_size=Title_size,Subtitle=Subtitle,Subtitle_size=8)
		}
		# Title_axis_size
	return(Plot)
}

make_plot_maf2<-function(refstudy=NULL,target_dat=NULL,eaf="eaf",snp="rsid",ref_dat_maf="maf",target_dat_effect_allele="Effect.Allele",target_dat_other_allele="Other.Allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome="outcome",ID=NULL,target_dat_population="population",ref_dat_population="population",target_study="study",ref_study="study",Title_xaxis_size=8,Title_size=10,Title="",Ylab="",Xlab="",cowplot_title=NULL,Dir="~/fatty-acids/outcome_data/data/",legend=TRUE){	

	ref_dat<-load_refdata(refstudy=refstudy,Dir=Dir)

	names(target_dat)[names(target_dat) == target_dat_population]<-"target_dat_population"
	names(ref_dat)[names(ref_dat) == ref_dat_population]<-"ref_dat_population"
	names(target_dat)[names(target_dat) == target_study]<-"target_study"
	names(ref_dat)[names(ref_dat) == ref_study]<-"ref_study"


	if(any(names(ref_dat) %in% c(target_dat_effect_allele,target_dat_other_allele,target_dat_effect_allele))) warning("effect allele, other allele or eaf present in refererence dataset with same name as in target dataset")
	
	if("maf" %in% names(target_dat)) {
		names(target_dat)[names(target_dat) == "maf"]<-"maf_target_dat"
	}

	dat.m<-merge(ref_dat,target_dat,by=snp)

	Pos<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])
	
	# harmonise study to ref dataset minor allele. need to clean this up as a new function with flip_strand and harmonise_allele functions like in the GWAS catalog functions
	dat.m[,eaf]<-as.numeric(dat.m[,eaf])
	dat.m[,ref_dat_maf]<-as.numeric(dat.m[,ref_dat_maf])
	dat.m[,eaf][Pos]<-1-dat.m[,eaf][Pos]
	EA<-dat.m[,target_dat_effect_allele][Pos]
	OA<-dat.m[,target_dat_other_allele][Pos]
	dat.m[,target_dat_effect_allele][Pos]<-OA
	dat.m[,target_dat_other_allele][Pos]<-EA
	
	# harmonise SNPs on different strands
	Pos1<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])
	Pos2<-which(dat.m[,target_dat_effect_allele] == dat.m[,ref_dat_minor_allele])
	dat.m1<-dat.m[Pos1,]
	dat.m2<-dat.m[Pos2,]

	
	# dat.m1[,c("Effect.Allele","Other.Allele","minor_allele2","major_allele2")]
	Pos<-which(dat.m1[,target_dat_effect_allele] != dat.m1[,"minor_allele2"])
	# dat.m1[Pos,c("Effect.Allele","Other.Allele","minor_allele2","major_allele2")]
	dat.m1[,eaf][Pos]<-1-dat.m1[,eaf][Pos]
	EA<-dat.m1[,target_dat_effect_allele][Pos]
	OA<-dat.m1[,target_dat_other_allele][Pos]
	dat.m1[,target_dat_effect_allele][Pos]<-OA
	dat.m1[,target_dat_other_allele][Pos]<-EA

	dat.m<-rbind(dat.m1,dat.m2)

	# dat.m[,c("Effect.Allele","Other.Allele","minor_allele","major_allele","eaf","maf")]

# Pos<-which(dat.m1$Effect.Allele  != dat.m1$minor_allele2 &  dat.m1$Effect.Allele  != dat.m1$minor_allele)
# dat.m1[Pos,c("Effect.Allele","Other.Allele","minor_allele","minor_allele2")]
	dat.m$alleles<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	dat.m$alleles2<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	
	dat.m.test<-dat.m
	outcome_plot<-outcome
	if(is.null(outcome)){
		outcome_plot<-""
	}

	if(!is.null(outcome) & !is.null(target_study) & !is.null(ID)){
		outfile_name<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID=",dat.m.test[,ID]))
		outfile_name<-gsub("\\|","",outfile_name)
		outfile_name<-gsub(" ","_",outfile_name)
		outfile_name<-gsub("__","_",outfile_name)
		outfile_name<-gsub("=","",outfile_name)
		outfile_name<-gsub("/","_",outfile_name)
		outcome_plot<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		target_study<-unique(paste0(dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		# outcome_plot2<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test[,study]," | ID: ",dat.m.test[,ID]))
	}

	Plot_list<-NULL
	dat.m.test<-dat.m.test[order(dat.m.test$ref_dat_population),]
	Pops<-unique(dat.m.test$ref_dat_population)
	# if(length(Pops)>1){
	for(i in 1:length(Pops)){
		# i<-1
		# print(pop)
		# dat1<-dat[dat$ref_dat_population==pop,]			
		dat1<-dat.m.test[dat.m.test$ref_dat_population==Pops[i], ]
		pop2<-c("European","East Asian","African","American","South Asian","Global","European","European")
		Pops2<-c("EUR","EAS","AFR","AMR","SAS","ALL","EUR2","EUR1")
		j<-which(Pops2 %in% Pops[i])
		# if(Xlab==""){
		# Xlab<-paste0(pop2[j]," MAF")
		Xlab<-pop2[j]

		# }
		# Title<-pop2[i]
		# Title<-gsub("ALL","Global pop",Title)
		# if(Title == ""){
		# 	Title<-target_study
		# }

		# if(subtitle_off) Title<-""

		Colour<-rep("black",nrow(dat1))
		Colour[which(dat1[,eaf]>0.5)]<-"blue"
		# Colour[dat1[,target_dat_effect_allele]!=dat1[,ref_dat_minor_allele]]<-"red"	
		
		# fix harmonisation functions above so that efffect allele strand flipped to minor_allele (not harmonised with minor_allele2)
		Colour[dat1$Effect.Allele!=dat1$minor_allele & dat1$Effect.Allele!=dat1$minor_allele2]<-"red"	

		Diff<-abs(dat1[,eaf]-dat1[,ref_dat_maf])
		Colour[which(Diff>0.10)]<-"red"
		Shape<-rep(19,nrow(dat1))
		Shape[which(dat1$alleles %in% c("AT","TA","GC","CG"))]<-1
		dat1$eaf<-dat1[,eaf]
		dat1$maf<-dat1[,ref_dat_maf]

		Title_size1<-Title_size
		Subtitle_size1<-8
		if(length(unique(dat.m.test$ref_dat_population)) > 1){
			Title_size1<-0
			Subtitle_size1<-0
		}

		
		# Temp<-dat1[dat1$eaf > 0.5,c("Effect.Allele","minor_allele","Other.Allele","major_allele","eaf","maf","alleles")]
		# dim(Temp)
		# length(which(Temp$alleles %in% c("CG","GC","TA","AT")))

		# snps<-dat1$rsid[dat1$eaf>0.55]
		# dat1[dat1$eaf>0.55,c("minor_allele","major_allele","Effect.Allele","Other.Allele","maf","eaf")]
		# dat2<-dat[dat$rsid %in% snps, ]
		# dat2[,c("rsid","Effect.Allele","Other.Allele","eaf")]
		# head(ref_dat)
		# dat2.m<-merge(dat2,ref_dat,by="rsid")
		# head(dat2.m[,c("rsid","Effect.Allele","Other.Allele","eaf", "maf", "minor_allele2","major_allele2" )])
		
		Subtitle<-unique(paste0("Reported population: ",dat1$target_dat_population))
		
		Shape2<-Colour
		Shape2[Shape2=="red"]<-1
		Shape2[Shape2=="blue"]<-2
		Shape2[Shape2=="black"]<-3		
		# Shape2<-as.numeric(Shape2)
		Plot<-ggplot2::ggplot(dat1, ggplot2::aes(x=maf, y=eaf)) + ggplot2::geom_point(aes(shape=Shape2)) +ggplot2::ggtitle(Title) +ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size1,hjust = 0))+ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = Subtitle_size1)) + 
			ggplot2::scale_shape_manual(name = "Allele frequency conflict",
		                     labels = c("High","Moderate","None") ,
		                     values = c(3,2,1))+
			ggplot2::theme(legend.title=ggplot2::element_text(size=8))+
			ggplot2::theme(legend.text=ggplot2::element_text(size=8))

		if(!legend){
			Plot<-ggplot2::ggplot(dat1, ggplot2::aes(x=maf, y=eaf)) + ggplot2::geom_point(aes(shape=Shape2)) +ggplot2::ggtitle(Title) +ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size1,hjust = 0))+ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle)+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = Subtitle_size1)) + 
				ggplot2::scale_shape_manual(name = "Allele frequency conflict",
			                     labels = c("High","Moderate","None") ,
			                     values = c(3,2,1))+
				ggplot2::theme(legend.title=ggplot2::element_text(size=8))+
				ggplot2::theme(legend.text=ggplot2::element_text(size=8))+
				  theme(legend.position="none")+
				  	ggplot2::theme(text = element_text(size=8)) 
		}


			
# as.factor(unique(Shape2)))
			# scale_shape_manual(labels = c("T999", "T888"), values = c("blue", "red"))


			# scale_fill_discrete(labels=c("High", "Moderate", "None"))

			if(length(Pops)> 1){
				Plot_list[[i]]<-Plot
			}
		}	

		if(length(unique(dat.m.test$ref_dat_population)) > 1){
			if(is.null(cowplot_title)){
				cowplot_title<-target_study
			}	
			Plot<-make_cow_plot(Plot_list=Plot_list,Title=cowplot_title,Xlab="",Ylab="",return_plot=TRUE,Title_size=Title_size,Subtitle=Subtitle,Subtitle_size=8)
		}
		# Title_axis_size
	return(Plot)
}



find_allele_error<-function(target_dat=NULL,ref_dat=NULL,eaf="eaf",snp="rsid",ref_dat_maf="maf",target_dat_effect_allele="effect_allele",target_dat_other_allele="other_allele",ref_dat_minor_allele="minor_allele",ref_dat_major_allele="major_allele",outcome=NULL,study=NULL,ID=NULL,target_dat_population="population",ref_dat_population="population",target_study=NULL,ref_study=NULL,Title_subx_axis_size=Title_subx_axis_size,Title_size_subplot=Title_size_subplot,Title=Title,Ylab=Ylab,Xlab=Xlab){	

	names(target_dat)[names(target_dat) == target_dat_population]<-"target_dat_population"
	names(ref_dat)[names(ref_dat) == ref_dat_population]<-"ref_dat_population"
	names(target_dat)[names(target_dat) == target_study]<-"target_study"
	names(ref_dat)[names(ref_dat) == ref_study]<-"ref_study"

	if(any(names(ref_dat) %in% c(target_dat_effect_allele,target_dat_other_allele,target_dat_effect_allele))) warning("effect allele, other allele or eaf present in refererence dataset with same name as in target dataset")
	
	dat.m<-merge(ref_dat,target_dat,by=snp)

	Pos<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])
	dat.m[,eaf]<-as.numeric(dat.m[,eaf])
	dat.m[,ref_dat_maf]<-as.numeric(dat.m[,ref_dat_maf])
	dat.m[,eaf][Pos]<-1-dat.m[,eaf][Pos]
	EA<-dat.m[,target_dat_effect_allele][Pos]
	OA<-dat.m[,target_dat_other_allele][Pos]
	dat.m[,target_dat_effect_allele][Pos]<-OA
	dat.m[,target_dat_other_allele][Pos]<-EA
	Pos<-which(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])

	# if(sum(Pos) != 0) warning(paste0("effect allele different to minor allele suggesting either strand conflict or other allele error: ",paste(unique(dat.m[,snp][Pos]),collapse=", ")))

	# dat.m[Pos,c("rsid","Effect.Allele","Other.Allele","eaf","minor_allele","major_allele","MAF" ,"ref_dat_population")]
	# dat.m2<-dat.m[Pos,]
	# Dat.m2<-flip_strand(dat=dat.m2,allele1=ref_dat_minor_allele,allele2=ref_dat_major_allele)

	# Pos2<-which(dat.m2[,target_dat_effect_allele] != dat.m2[,ref_dat_minor_allele])
	# if(sum(Pos2) != 0){
	# 	Pos2<-which(dat.m2[,target_dat_other_allele] != dat.m2[,ref_dat_minor_allele])
	# }
	
	dat.m$alleles<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	dat.m$alleles2<-paste0(dat.m[,ref_dat_minor_allele],dat.m[,ref_dat_major_allele])
	####################################################
	# restrict to maf < 0.40 and non palindromics#######	####################################################
	# Pos<-which(dat.m[,eaf]<=0.40 | dat.m[,eaf]>=0.60 )
	# Pos<-which(dat.m[,ref_dat_maf]<=0.45)
	# dat.m.test<-dat.m[Pos,]
	
	# palindromes<-c("CG","GC","TA","AT")
	# Pos<-which(!alleles %in% palindromes)
	# dat.m.test<-dat.m.test[Pos,]
	
	# if(any(dat.m[,target_dat_effect_allele] != dat.m[,ref_dat_minor_allele])) stop("alt allele different from minor allele in reference population")
	dat.m.test<-dat.m
	outcome_plot<-outcome
	if(is.null(outcome)){
		outcome_plot<-""
	}

	if(!is.null(outcome) & !is.null(target_study) & !is.null(ID)){
		outfile_name<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID=",dat.m.test[,ID]))
		outfile_name<-gsub("\\|","",outfile_name)
		outfile_name<-gsub(" ","_",outfile_name)
		outfile_name<-gsub("__","_",outfile_name)
		outfile_name<-gsub("=","",outfile_name)
		outfile_name<-gsub("/","_",outfile_name)
		outcome_plot<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		target_study<-unique(paste0(dat.m.test$target_study," | ID: ",dat.m.test[,ID]))		
		# outcome_plot2<-unique(paste0(dat.m.test[,outcome]," | ", dat.m.test[,study]," | ID: ",dat.m.test[,ID]))
	}

	# Plot<-make_plot(dat=dat.m.test,target_dat_population = unique(dat.m.test$target_dat_population),eaf=eaf,maf=ref_dat_maf,ea=target_dat_effect_allele,ma=ref_dat_minor_allele,outcome_plot=outcome_plot,target_study=target_study)

	Plot_list<-NULL
	dat.m.test<-dat.m.test[order(dat.m.test$ref_dat_population),]
	Pops<-unique(dat.m.test$ref_dat_population)
	# if(length(Pops)>1){
	for(i in 1:length(Pops)){
		# print(pop)
		# dat1<-dat[dat$ref_dat_population==pop,]			
		dat1<-dat.m.test[dat.m.test$ref_dat_population==Pops[i], ]


		# pop<-dat$ref_dat_population
		# Ylab<-paste0("MAF in ",target_study)
		# Ylab<-""
		# Xlab<-paste0("MAF in ",pop," (",unique(dat1$ref_study),")")
		# Xlab<-paste0("MAF in ",Pops[i])
		# Xlab<-Pops[i]
		# Xlab<-""
		pop2<-c("European","East Asian","African","American","South Asian","Global","European","European")
		Pops2<-c("EUR","EAS","AFR","AMR","SAS","ALL","EUR2","EUR1")
		j<-which(Pops2 %in% Pops[i])
		if(Xlab==""){
			Xlab<-pop2[j]
		}
		# Title<-pop2[i]
		# Title<-gsub("ALL","Global pop",Title)
		if(Title == ""){
			Title<-target_study
		}

		# if(subtitle_off) Title<-""

		Colour<-rep("black",nrow(dat1))
		Colour[which(dat1[,eaf]>0.5)]<-"blue"
		Colour[dat1[,target_dat_effect_allele]!=dat1[,ref_dat_minor_allele]]<-"red"	
		Diff<-abs(dat1[,eaf]-dat1[,ref_dat_maf])
		Colour[which(Diff>0.10)]<-"red"
		Shape<-rep(19,nrow(dat1))
		Shape[which(dat1$alleles %in% c("AT","TA","GC","CG"))]<-1
		dat1$eaf<-dat1[,eaf]
		dat1$maf<-dat1[,ref_dat_maf]
		
		Plot<-ggplot2::ggplot(dat1, ggplot2::aes(x=maf, y=eaf)) + ggplot2::geom_point(colour=Colour) +ggplot2::ggtitle(Title) +theme(plot.title = element_text(size = Title_size_subplot,hjust = 0))+ggplot2::labs(y= Ylab, x =Xlab)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_subx_axis_size)) 

			if(length(Pops)> 1){
				Plot_list[[i]]<-Plot
			}
		}	


		if(length(unique(dat.m.test$ref_dat_population)) > 1){
			Plot<-make_cow_plot(Plot_list=Plot_list,Title=Title,Xlab="",Ylab="",return_plot=TRUE)
		}

# ggplot() +
#   ggtitle("Use theme(plot.title = element_text(hjust = 0.5)) to center") +
#   theme(plot.title = element_text(hjust = 0.5))

# +ggplot2::theme(axis.title=ggplot2::element_text(size=5)) 
# +ggplot2::labs(y= Ylab, x =Xlab) 

	# outfile_name=outfile_name

	# Pos<-NULL
	# snp_mismatch<-0
	# # Diff<-abs(dat.m.test[,eaf]-dat.m.test[,ref_dat_maf])
	# Pos<-which(dat.m.test[,eaf]>=0.5)
	# # Pos<-which(Diff>=0.10)
	# # dat.m.test[Pos,]
	# # af_diff<-abs(dat.m.test[,eaf]-dat.m.test[,ref_dat_maf])
	# # dat.m.test[af_diff>0.15,]
	# if(sum(Pos)!=0){		
	# 	dat.m.test2<-dat.m.test[Pos,]				
	# 	if(any(dat.m.test2[,eaf]<=0.45 | dat.m.test2[,eaf]>=0.55)){
	# 		Pos2<-which(dat.m.test2[,eaf]<=0.45 | dat.m.test2[,eaf]>=0.55)
	# 		warning(paste("minor allele conflict for these SNPs: ",paste(unique(dat.m.test2[,snp][Pos2]),collapse=", ")))
	# 		snp_mismatch<-unique(dat.m.test2[,snp][Pos2])
	# 	}
	# }

	# warning(paste("minor allele conflict for these SNPs: ",paste(unique(dat.m.test[,snp][Pos]),collapse=", ")))
	# snp_mismatch<-unique(dat.m.test2[,snp][Pos2])



# return(list("plot"=Plot,"allele_errors"=snp_mismatch,outfile_name=outfile_name))
	return(Plot)
	# return(list("plot"=Plot,"allele_errors"=snp_mismatch,outfile_name=outfile_name))
}





load_refdata<-function(refstudy=NULL,Dir="~/fatty-acids/outcome_data/data/"){
	if(refstudy=="EUR_1000G"){
		ref.dat<-load_plinkfrq(File<-paste0(Dir,"fatty_acid_snps_eur.frq"),population="EUR1",study=refstudy)	
	}
	if(refstudy=="CHARGE"){
		load(paste0(Dir,"ref_dat_charge_imputed.RData"))
		ref.dat<-load_ref(Dat=Charge,population="EUR2",study=refstudy)			
	}
# ref.dat<-plyr::rbind.fill(ref.dat_eur,ref.dat_charge)
	if(refstudy=="EAS_1000G"){
		ref.dat<-load_plinkfrq(File<-paste0(Dir,"fatty_acid_snps_eas.frq"),population="EAS1",study=refstudy)

	}

	if(refstudy=="SCHS"){
		load(paste0(Dir,"ref_dat_schs.RData"))
		ref.dat<-load_ref(Dat=SCHS,population="EAS2",study=refstudy)	
	}
# ref.dat<-plyr::rbind.fill(ref.dat_eas,ref.dat_schs)
	if(refstudy=="ALL_1000G"){
		load(paste0(Dir,"refdat_1000G_superpops.Rdata"))
		ref.dat<-format_refdat_1000G_superpops(refdat_1000G_superpops=refdat_1000G_superpops,study="1000G")			
		ref.dat<-make_refdat(ref.dat=ref.dat)
	}
	return(ref.dat)
}
# make_refdat


compare_gwas_catalog_hits<-function(dat=NULL,efo_id=NULL,efo=NULL,return_plot=FALSE,Title_size_subplot=10,Title=NULL,Ylab="",Xlab=""){
	Plot1<-make_plot_gwas_catalog_zscores(dat=dat,efo_id=efo_id,efo=efo,Title_size_subplot=Title_size_subplot,Title=Title,Ylab=Ylab,Xlab=Xlab)
	Plot2<-make_plot_gwas_catalog_eaf(dat=dat,efo=efo,Title_size_subplot=Title_size_subplot,Title=Title,Ylab=Ylab,Xlab=Xlab)	
	file_out<-unique(paste0(dat$outcome,"_",dat$study,"_ID_",dat$ID))
	file_out<-gsub("/","",file_out)

	# if(!return_plot){
	# 	png(paste0("~/fatty-acids/outcome_data/results/plots/gwas_catalog_hits/",file_out,".png"),width=800)
	# 		print(Plot)
	# 	dev.off()
	# }
	return(list(zscore=Plot1,eaf=Plot2))
}

make_plot_gwas_catalog_zscores2<-function(dat=NULL,efo_id=NULL,efo=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),legend=TRUE,Title=Title,Title_size_subplot=Title_size_subplot,Ylab=Ylab,Xlab=Xlab,Title_xaxis_size=0){
	# if(!is.null(trait)){
		# gwas_catalog<-gwas_catalog_hits(trait=trait)
	# }
# gwas_catalog_ancestral_group="East Asian"
	if(!is.null(efo)){
		gwas_catalog<-gwas_catalog_hits(efo=efo)
	}

	# efo_id<-CheckSumStats::get_efo(unique(dat$outcome))$efo_id
	# efo_id<-"EFO_0000096"
	if(!is.null(efo_id)){
		gwas_catalog<-gwas_catalog_hits(efo_id=efo_id)
	}

	
	Dat.m<-merge(gwas_catalog,dat,by="rsid")
	Dat.m<-Dat.m[!is.na(Dat.m$Effect.Allele.x),]
	Dat.m<-Dat.m[nchar(Dat.m$Effect.Allele.y)==1,]
	Dat.m<-Dat.m[nchar(Dat.m$Other.Allele)==1,]
	Alleles<-paste0(Dat.m$Effect.Allele.y,Dat.m$Other.Allele)
	Dat.m<-Dat.m[!Alleles %in% c("AT","TA","GC","CG"),]
	Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]		
	Dat.m<-harmonise_effect_allele(Dat=Dat.m)

	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y	
	if(any(Pos)) {
		Dat.m<-flip_strand(Dat=Dat.m)
	}
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y
	if(any(Pos)){
		Dat.m<-harmonise_effect_allele(Dat=Dat.m)
	}	
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y

	if(any(Pos)) {
		stop("effect alleles not full harmonised")	
		# Dat.m[Pos,c("rsid","Effect.Allele.x","Effect.Allele.y","Other.Allele")]
	}

	Dat.m$z.y<-Dat.m$lnor.y/Dat.m$se.y
	Dat.m$z.x<-Dat.m$lnor.x/Dat.m$se.x

	# head(Dat.m[,c("p.x","z.x","p.y","z.y")])
	# max(Dat.m$p.x)
	# dim(Dat.m)
	# Ylab<-""
	# Xlab<-""

	Z_scores<-rep("black",nrow(Dat.m))
	Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"blue"
	Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >= 3.890592 & abs(Dat.m$z.x) >= 3.890592 )]<-"red" # Z score of 3.890592 = 2 sided p value of 0.0001	
	
	# Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >=  4.891638  & abs(Dat.m$z.x) >=  4.891638 )]<-"red"

	ancestry1<-Dat.m$ancestral_group
	
	if(is.null(Title)){
		Title<-paste0(unique(dat$study)," | " ,unique(dat$ID) , " | EFO: ", efo)
	}

	labels_colour<-unique(Z_scores)
	labels_colour[labels_colour == "red"]<-"high"
	labels_colour[labels_colour == "blue"]<-"moderate"
	labels_colour[labels_colour == "black"]<-"none"
	values_colour<-unique(Z_scores)
	Pos<-order(values_colour)
	values_colour<-values_colour[Pos]
	labels_colour<-labels_colour[Pos]

	# unique(ancestry1)[order(unique(ancestry1))]
	# 1:length(unique(ancestry1))
	# ancestry2<-c("European")
	
	# Shape<-ancestry1
	# Shape[Shape=="European"]<-15
	# Shape[Shape=="East Asian"]<-16
	# Shape<-as.numeric(Shape)
	# # Shape<-unique(Shape)[order(unique(Shape))]

	labels_shape<-unique(ancestry1)[order(unique(ancestry1))]
	values_shape<-labels_shape
	values_shape[values_shape == "European"]<-15
	values_shape[values_shape == "East Asian"]<-16
	values_shape<-as.numeric(values_shape)


	values_shape2<-Z_scores
	values_shape2[values_shape2=="red"]<-1
	values_shape2[values_shape2=="blue"]<-2
	values_shape2[values_shape2=="black"]<-3
	labels_shape2<-c("High","Moderate","None")
	Z_scores2<-values_shape2
	values_shape2<-unique(as.numeric(values_shape2))

	# Shape2<-Shape
	# Shape2<-unique(Shape2)[order(unique(Shape2))]

	labels_colour2<-unique(ancestry1)[order(unique(ancestry1))]
	values_colour2<-labels_colour2
	values_colour2[values_colour2 == "European"]<-"black"
	values_colour2[values_colour2 == "East Asian"]<-"gray"
	# values_colour2<-as.numeric(values_colour2)

	Subtitle<-paste0(Dat.m$outcome," | ",Dat.m$population)

		# ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=z.x, y=z.y,colour=Z_scores,shape=ancestry1))
		if(legend){
			Plot<-ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=z.x, y=z.y,colour=ancestry1,shape=Z_scores2)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"),
				)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
			 ggplot2::scale_shape_manual(name = "Effect size conflict",
		                     labels = labels_shape2,
		                     # labels = unique(ancestry1)[order(unique(ancestry1))],
		                     # labels = c("European","East Asian"),
		                     values = values_shape2) + 
		                     # values = 1:length(Shape2)) + 
		 	ggplot2::scale_colour_manual(name="GWAS catalog ancestry",
			              labels=labels_colour2,
			              values=values_colour2)+
		 	ggplot2::theme(legend.title=ggplot2::element_text(size=8))+
			ggplot2::theme(legend.text=ggplot2::element_text(size=8))	

			
		}

	if(!legend){
		Plot<-ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=z.x, y=z.y,colour=Z_scores,shape=ancestry1)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		 ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
	                    labels = labels_shape,	                     
	                     values = values_shape) + 
	 	ggplot2::scale_colour_manual(name="Effect size conflict",
		              labels=labels_colour,
		              values=values_colour)+
	 	ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 		legend.text=ggplot2::element_text(size=8),plot.subtitle = ggplot2::element_text(size = 8),
	 		legend.position = "none")
	}
	 # ggplot2::scale_colour_manual(name="Z score conflict",
  #                     labels=unique(Z_scores)[order(unique(Z_scores))] ,
  #                     values=unique(Z_scores)[order(unique(Z_scores))]) 	 

	  # ggplot2::scale_colour_manual(name="Z score conflict",
                      # labels=c("none", "moderate","high"),
                      # values=c("black","blue", "red")) 	 
  	 
  	
	return(Plot)
}

make_plot_gwas_catalog_zscores<-function(dat=NULL,efo_id=NULL,efo=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),legend=TRUE,Title=Title,Title_size_subplot=Title_size_subplot,Ylab=Ylab,Xlab=Xlab,Title_xaxis_size=0){
	# if(!is.null(trait)){
		# gwas_catalog<-gwas_catalog_hits(trait=trait)
	# }
# gwas_catalog_ancestral_group="East Asian"
	if(!is.null(efo)){
		gwas_catalog<-gwas_catalog_hits(efo=efo)
	}

	# efo_id<-CheckSumStats::get_efo(unique(dat$outcome))$efo_id
	# efo_id<-"EFO_0000096"
	if(!is.null(efo_id)){
		gwas_catalog<-gwas_catalog_hits(efo_id=efo_id)
	}

	
	Dat.m<-merge(gwas_catalog,dat,by="rsid")
	Dat.m<-Dat.m[!is.na(Dat.m$Effect.Allele.x),]
	Dat.m<-Dat.m[nchar(Dat.m$Effect.Allele.y)==1,]
	Dat.m<-Dat.m[nchar(Dat.m$Other.Allele)==1,]
	Alleles<-paste0(Dat.m$Effect.Allele.y,Dat.m$Other.Allele)
	Dat.m<-Dat.m[!Alleles %in% c("AT","TA","GC","CG"),]
	Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]		
	Dat.m<-harmonise_effect_allele(Dat=Dat.m)

	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y	
	if(any(Pos)) {
		Dat.m<-flip_strand(Dat=Dat.m)
	}
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y
	if(any(Pos)){
		Dat.m<-harmonise_effect_allele(Dat=Dat.m)
	}	
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y

	if(any(Pos)) {
		stop("effect alleles not full harmonised")	
		# Dat.m[Pos,c("rsid","Effect.Allele.x","Effect.Allele.y","Other.Allele")]
	}

	Dat.m$z.y<-Dat.m$lnor.y/Dat.m$se.y
	Dat.m$z.x<-Dat.m$lnor.x/Dat.m$se.x

	# head(Dat.m[,c("p.x","z.x","p.y","z.y")])
	# max(Dat.m$p.x)
	# dim(Dat.m)
	# Ylab<-""
	# Xlab<-""

	Z_scores<-rep("black",nrow(Dat.m))
	Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)))]<-"blue"
	Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >= 3.890592 & abs(Dat.m$z.x) >= 3.890592 )]<-"red" # Z score of 3.890592 = 2 sided p value of 0.0001	
	
	# Z_scores[which(sign(Dat.m$z.y) != sign(as.numeric(Dat.m$z.x)) & abs(Dat.m$z.y) >=  4.891638  & abs(Dat.m$z.x) >=  4.891638 )]<-"red"

	ancestry1<-Dat.m$ancestral_group
	
	if(is.null(Title)){
		Title<-paste0(unique(dat$study)," | " ,unique(dat$ID) , " | EFO: ", efo)
	}

	labels_colour<-unique(Z_scores)
	labels_colour[labels_colour == "red"]<-"high"
	labels_colour[labels_colour == "blue"]<-"moderate"
	labels_colour[labels_colour == "black"]<-"none"
	values_colour<-unique(Z_scores)
	Pos<-order(values_colour)
	values_colour<-values_colour[Pos]
	labels_colour<-labels_colour[Pos]

	# unique(ancestry1)[order(unique(ancestry1))]
	# 1:length(unique(ancestry1))
	# ancestry2<-c("European")
	
	# Shape<-ancestry1
	# Shape[Shape=="European"]<-15
	# Shape[Shape=="East Asian"]<-16
	# Shape<-as.numeric(Shape)
	# # Shape<-unique(Shape)[order(unique(Shape))]

	labels_shape<-unique(ancestry1)[order(unique(ancestry1))]
	values_shape<-labels_shape
	values_shape[values_shape == "European"]<-15
	values_shape[values_shape == "East Asian"]<-16
	values_shape<-as.numeric(values_shape)

	# Shape2<-Shape
	# Shape2<-unique(Shape2)[order(unique(Shape2))]
	

	Subtitle<-paste0(Dat.m$outcome," | ",Dat.m$population)

		# ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=z.x, y=z.y,colour=Z_scores,shape=ancestry1))
		if(legend){
			Plot<-ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=z.x, y=z.y,colour=Z_scores,shape=ancestry1)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"),
				)+
			ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
			 ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
		                     labels = labels_shape,
		                     # labels = unique(ancestry1)[order(unique(ancestry1))],
		                     # labels = c("European","East Asian"),
		                     values = values_shape) + 
		                     # values = 1:length(Shape2)) + 
		 	ggplot2::scale_colour_manual(name="Effect size conflict",
			              labels=labels_colour,
			              values=values_colour)+
		 	ggplot2::theme(legend.title=ggplot2::element_text(size=8))+
			ggplot2::theme(legend.text=ggplot2::element_text(size=8))
		}

	if(!legend){
		Plot<-ggplot2::ggplot(Dat.m) + ggplot2::geom_point(ggplot2::aes(x=z.x, y=z.y,colour=Z_scores,shape=ancestry1)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		 ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
	                    labels = labels_shape,	                     
	                     values = values_shape) + 
	 	ggplot2::scale_colour_manual(name="Effect size conflict",
		              labels=labels_colour,
		              values=values_colour)+
	 	ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 		legend.text=ggplot2::element_text(size=8),plot.subtitle = ggplot2::element_text(size = 8),
	 		legend.position = "none")
	}
	 # ggplot2::scale_colour_manual(name="Z score conflict",
  #                     labels=unique(Z_scores)[order(unique(Z_scores))] ,
  #                     values=unique(Z_scores)[order(unique(Z_scores))]) 	 

	  # ggplot2::scale_colour_manual(name="Z score conflict",
                      # labels=c("none", "moderate","high"),
                      # values=c("black","blue", "red")) 	 
  	 
  	
	return(Plot)
}


flip_strand<-function(Dat=NULL){
	Pos<-Dat$Effect.Allele.x!=Dat$Effect.Allele.y	
	strand1<-c("A","T","G","C")
	strand2<-c("T","A","C","G")
	# lnor.y<-Dat$lnor.y[Pos]*-1
	# Dat$lnor.y[Pos]<-lnor.y
	ea<-Dat$Effect.Allele.y[Pos]
	oa<-Dat$Other.Allele[Pos]
	Dat$Effect.Allele.y[Pos]<-strand2[match(ea,strand1)]
	Dat$Other.Allele[Pos]<-strand2[match(oa,strand1)]				
	return(Dat)
}
		
harmonise_effect_allele<-function(Dat=NULL){
	Pos<-Dat$Effect.Allele.x!=Dat$Effect.Allele.y 	
	lnor.y<-Dat$lnor.y[Pos]*-1
	Dat$lnor.y[Pos]<-lnor.y
	oa<-Dat$Effect.Allele.y[Pos]
	ea<-Dat$Other.Allele[Pos]
	Dat$Effect.Allele.y[Pos]<-ea
	Dat$Other.Allele[Pos]<-oa
	eaf<-1-Dat$eaf.y[Pos]
	Dat$eaf.y[Pos]<-eaf	
	return(Dat)
}

make_plot_gwas_catalog_eaf2<-function(dat=NULL,efo_id=NULL,efo=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),Title_size_subplot=Title_size_subplot,Title=Title,Ylab=Ylab,Xlab=Xlab,Title_xaxis_size=0,Legend=TRUE){
	# if(!is.null(trait)){
		# gwas_catalog<-gwas_catalog_hits(trait=trait)
	# }

	if(!is.null(efo)){
		gwas_catalog<-gwas_catalog_hits(efo=efo)
	}

	if(!is.null(efo_id)){
		gwas_catalog<-gwas_catalog_hits(efo_id=efo_id)
	}

	Dat.m<-merge(gwas_catalog,dat,by="rsid")		
	Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]	
	Pos<-which(Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y)
	lnor.y<-Dat.m$lnor.y[Pos]*-1
	Dat.m$lnor.y[Pos]<-lnor.y
	oa<-Dat.m$Effect.Allele.y[Pos]
	ea<-Dat.m$Other.Allele[Pos]
	Dat.m$Effect.Allele.y[Pos]<-ea
	Dat.m$Other.Allele[Pos]<-oa
	eaf<-1-Dat.m$eaf.y[Pos]
	Dat.m$eaf.y[Pos]<-eaf

	Dat.m<-Dat.m[!is.na(Dat.m$Effect.Allele.x),]
	Dat.m<-Dat.m[nchar(Dat.m$Effect.Allele.y)==1,]
	Dat.m<-Dat.m[nchar(Dat.m$Other.Allele)==1,]
	Alleles<-paste0(Dat.m$Effect.Allele.y,Dat.m$Other.Allele)
	Dat.m<-Dat.m[!Alleles %in% c("AT","TA","GC","CG"),]
	Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]		
	Dat.m<-harmonise_effect_allele(Dat=Dat.m)

	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y	
	if(any(Pos)) {
		Dat.m<-flip_strand(Dat=Dat.m)
	}
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y
	if(any(Pos)){
		Dat.m<-harmonise_effect_allele(Dat=Dat.m)
	}	
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y

	if(any(Pos)) {
		stop("effect alleles not full harmonised")	
		# Dat.m[Pos,c("rsid","Effect.Allele.x","Effect.Allele.y","Other.Allele")]
	}

	Dat.m$z.y<-Dat.m$lnor.y/Dat.m$se.y
	Dat.m$z.x<-Dat.m$lnor.x/Dat.m$se.x
	ancestry1<-Dat.m$ancestral_group
	
	if(is.null(Title)){
		Title<-paste0(unique(dat$study)," | " ,unique(dat$ID) , " | EFO: ", efo)
	}	

	Dat.m2<-Dat.m[!is.na(Dat.m$eaf.x),]
	Dat.m2<-Dat.m2[!is.na(Dat.m2$eaf.y),]
	ancestry2<-Dat.m2$ancestral_group
	Pos1<-Dat.m2$eaf.x<0.5 & Dat.m2$eaf.y>0.5 | Dat.m2$eaf.x>0.5 & Dat.m2$eaf.y<0.5
	EAF<-rep("black",nrow(Dat.m2))
	if(any(Pos1)) EAF[Pos1]<-"blue"	 
	Pos2<-Dat.m2$eaf.x<0.40 & Dat.m2$eaf.y>0.60 | Dat.m2$eaf.x>0.60 & Dat.m2$eaf.y<0.40
	if(any(Pos2)) EAF[Pos2]<-"red"	 
	# Dat.m2[Pos2,c("eaf.x","eaf.y","Effect.Allele.x","Effect.Allele.y","Other.Allele","z.x","z.y")]

	# Ylab<-paste0("EAF in ",unique(dat$study)," | ID:" ,unique(dat$ID))
	# Xlab<-paste0("EAF in gwas catalog")
	# Ylab<-""
	# Xlab<-""
	labels_colour<-unique(EAF)
	labels_colour[labels_colour == "red"]<-"high"
	labels_colour[labels_colour == "blue"]<-"moderate"
	labels_colour[labels_colour == "black"]<-"none"
	values_colour<-unique(EAF)
	Pos<-order(values_colour)
	values_colour<-values_colour[Pos]
	labels_colour<-labels_colour[Pos]

	labels_shape<-unique(ancestry2)[order(unique(ancestry2))]
	values_shape<-labels_shape
	values_shape[values_shape == "European"]<-15
	values_shape[values_shape == "East Asian"]<-16
	values_shape<-as.numeric(values_shape)

	values_shape2<-unique(EAF)
	EAF2<-EAF
	EAF2[EAF2 == "red"]<-1
	EAF2[EAF2 == "blue"]<-2
	EAF2[EAF2 == "black"]<-3
	values_shape2[values_shape2 == "red"]<-1
	values_shape2[values_shape2 == "blue"]<-2
	values_shape2[values_shape2 == "black"]<-3
	labels_shape2<-unique(EAF)
	values_shape2<-as.numeric(values_shape2)
	# Pos<-order(labels_shape2)
	# labels_shape2<-labels_shape2[Pos]
	# values_shape2<-values_shape2[Pos]
	labels_shape2[labels_shape2=="red"]<-"high"
	labels_shape2[labels_shape2=="blue"]<-"moderate"
	labels_shape2[labels_shape2=="black"]<-"none"

	labels_colour2<-unique(ancestry2)[order(unique(ancestry2))]
	values_colour2<-labels_colour2
	values_colour2[values_colour2 == "European"]<-"black"
	if("East Asians" %in% unique(ancestry2)) { 
		values_colour2[values_colour2 == "East Asian"]<-"gray"
	}
	
	# ancestry2[ancestry2=="European"]<-"black"
	# ancestry2[ancestry2=="East Asian"]<-"gray"
	# values_colour2<-as.numeric(values_colour2)

	Subtitle<-paste0(Dat.m2$outcome," | ",Dat.m2$population)

	
	# ancestry2<-"European"

	if(Legend){
		Plot<-ggplot2::ggplot(Dat.m2) + ggplot2::geom_point(ggplot2::aes(x=eaf.x, y=eaf.y,colour=ancestry2,shape=EAF2)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
		ggplot2::scale_shape_manual(name = "EAF conflict",
		             labels = labels_shape2,
		             values = values_shape2) + 
		ggplot2::scale_colour_manual(name="GWAS catalog ancestry",
		              labels=labels_colour2,
		              values=values_colour2)+
	 	ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 		legend.text=ggplot2::element_text(size=8))
	 }

	if(!Legend){
		# Dat.m2$eaf.y
		Plot<-ggplot2::ggplot(Dat.m2) + ggplot2::geom_point(ggplot2::aes(x=eaf.x, y=eaf.y,colour=EAF,shape=ancestry2)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
		ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
		             labels = labels_shape,
		             values = values_shape) + 
		ggplot2::scale_colour_manual(name="EAF conflict",
		              labels=labels_colour,
		              values=values_colour)+
	 	ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 		legend.text=ggplot2::element_text(size=8),	
	 		legend.position = "none")
	 }
	

	
	return(Plot)
}

make_plot_gwas_catalog_eaf<-function(dat=NULL,efo_id=NULL,efo=NULL,gwas_catalog_ancestral_group=c("European","East Asian"),Title_size_subplot=Title_size_subplot,Title=Title,Ylab=Ylab,Xlab=Xlab,Title_xaxis_size=0,Legend=TRUE){
	# if(!is.null(trait)){
		# gwas_catalog<-gwas_catalog_hits(trait=trait)
	# }

	if(!is.null(efo)){
		gwas_catalog<-gwas_catalog_hits(efo=efo)
	}

	if(!is.null(efo_id)){
		gwas_catalog<-gwas_catalog_hits(efo_id=efo_id)
	}

	Dat.m<-merge(gwas_catalog,dat,by="rsid")		
	Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]	
	Pos<-which(Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y)
	lnor.y<-Dat.m$lnor.y[Pos]*-1
	Dat.m$lnor.y[Pos]<-lnor.y
	oa<-Dat.m$Effect.Allele.y[Pos]
	ea<-Dat.m$Other.Allele[Pos]
	Dat.m$Effect.Allele.y[Pos]<-ea
	Dat.m$Other.Allele[Pos]<-oa
	eaf<-1-Dat.m$eaf.y[Pos]
	Dat.m$eaf.y[Pos]<-eaf

	Dat.m<-Dat.m[!is.na(Dat.m$Effect.Allele.x),]
	Dat.m<-Dat.m[nchar(Dat.m$Effect.Allele.y)==1,]
	Dat.m<-Dat.m[nchar(Dat.m$Other.Allele)==1,]
	Alleles<-paste0(Dat.m$Effect.Allele.y,Dat.m$Other.Allele)
	Dat.m<-Dat.m[!Alleles %in% c("AT","TA","GC","CG"),]
	Dat.m<-Dat.m[Dat.m$ancestral_group %in% gwas_catalog_ancestral_group,]		
	Dat.m<-harmonise_effect_allele(Dat=Dat.m)

	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y	
	if(any(Pos)) {
		Dat.m<-flip_strand(Dat=Dat.m)
	}
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y
	if(any(Pos)){
		Dat.m<-harmonise_effect_allele(Dat=Dat.m)
	}	
	Pos<-Dat.m$Effect.Allele.x!=Dat.m$Effect.Allele.y

	if(any(Pos)) {
		stop("effect alleles not full harmonised")	
		# Dat.m[Pos,c("rsid","Effect.Allele.x","Effect.Allele.y","Other.Allele")]
	}

	Dat.m$z.y<-Dat.m$lnor.y/Dat.m$se.y
	Dat.m$z.x<-Dat.m$lnor.x/Dat.m$se.x
	ancestry1<-Dat.m$ancestral_group
	
	if(is.null(Title)){
		Title<-paste0(unique(dat$study)," | " ,unique(dat$ID) , " | EFO: ", efo)
	}	

	Dat.m2<-Dat.m[!is.na(Dat.m$eaf.x),]
	Dat.m2<-Dat.m2[!is.na(Dat.m2$eaf.y),]
	ancestry2<-Dat.m2$ancestral_group
	Pos1<-Dat.m2$eaf.x<0.5 & Dat.m2$eaf.y>0.5 | Dat.m2$eaf.x>0.5 & Dat.m2$eaf.y<0.5
	EAF<-rep("black",nrow(Dat.m2))
	if(any(Pos1)) EAF[Pos1]<-"blue"	 
	Pos2<-Dat.m2$eaf.x<0.40 & Dat.m2$eaf.y>0.60 | Dat.m2$eaf.x>0.60 & Dat.m2$eaf.y<0.40
	if(any(Pos2)) EAF[Pos2]<-"red"	 
	# Dat.m2[Pos2,c("eaf.x","eaf.y","Effect.Allele.x","Effect.Allele.y","Other.Allele","z.x","z.y")]

	# Ylab<-paste0("EAF in ",unique(dat$study)," | ID:" ,unique(dat$ID))
	# Xlab<-paste0("EAF in gwas catalog")
	# Ylab<-""
	# Xlab<-""
	labels_colour<-unique(EAF)
	labels_colour[labels_colour == "red"]<-"high"
	labels_colour[labels_colour == "blue"]<-"moderate"
	labels_colour[labels_colour == "black"]<-"none"
	values_colour<-unique(EAF)
	Pos<-order(values_colour)
	values_colour<-values_colour[Pos]
	labels_colour<-labels_colour[Pos]

	labels_shape<-unique(ancestry2)[order(unique(ancestry2))]
	values_shape<-labels_shape
	values_shape[values_shape == "European"]<-15
	values_shape[values_shape == "East Asian"]<-16
	values_shape<-as.numeric(values_shape)

	Subtitle<-paste0(Dat.m2$outcome," | ",Dat.m2$population)

	if(!Legend){
		# Dat.m2$eaf.y
		Plot<-ggplot2::ggplot(Dat.m2) + ggplot2::geom_point(ggplot2::aes(x=eaf.x, y=eaf.y,colour=EAF,shape=ancestry2)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
		ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
		             labels = labels_shape,
		             values = values_shape) + 
		ggplot2::scale_colour_manual(name="EAF conflict",
		              labels=labels_colour,
		              values=values_colour)+
	 	ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 		legend.text=ggplot2::element_text(size=8),	
	 		legend.position = "none")
	 }


	if(Legend){
		Plot<-ggplot2::ggplot(Dat.m2) + ggplot2::geom_point(ggplot2::aes(x=eaf.x, y=eaf.y,colour=EAF,shape=ancestry2)) +ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=Subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size_subplot, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size),plot.subtitle = ggplot2::element_text(size = 8))+
		ggplot2::scale_shape_manual(name = "GWAS catalog ancestry",
		             labels = labels_shape,
		             values = values_shape) + 
		ggplot2::scale_colour_manual(name="EAF conflict",
		              labels=labels_colour,
		              values=values_colour)+
	 	ggplot2::theme(legend.title=ggplot2::element_text(size=8),
	 		legend.text=ggplot2::element_text(size=8))
	 }

	# ggplot2::scale_colour_manual(name="EAF conflict",
	#               labels=c("none","moderate","high"),
	#               values=c("black","blue", "red"))

	
	return(Plot)
}

# Plot4<-Plot		

predict_lnor_and_plot1<-function(dat=NULL,Xlab="",Ylab="",threshold=NULL,maf_filter=FALSE){

	plot_dat<-predict_lnor(lnor=dat$lnor,se=dat$se,n=dat$ncase,p=dat$eaf,cases=dat$ncase,controls=dat$ncontrol)
	for(i in 1:ncol(plot_dat)){
		plot_dat[,i][plot_dat[,i] == "Inf" | plot_dat[,i] == "-Inf" ]<-NA
	}
	plot_dat<-plot_dat[complete.cases(plot_dat),]
	outfile_name<-unique(paste(dat$outcome,dat$study,dat$ID,sep="_"))
	outfile_name<-gsub(" ","_",outfile_name)
	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))
	Plot<-make_plot2(dat=plot_dat[,c("lnor_obs","lnor_pred1","eaf")],outfile_name=outfile_name,outcome_name=outcome_name,Xlab=Xlab,Ylab=Ylab,threshold=threshold,maf_filter=maf_filter)
	return(Plot)
	# outfile_name<-gsub("/","_",outfile_name)

	# png(paste0("~/fatty-acids/outcome_data/results/plots/predicted_lnor/",outfile_name,"_pi_method.png"))
	# 	print(Plot)
	# dev.off()

	
}


predict_lnor_and_plot2<-function(dat=NULL,Xlab="",Ylab="",threshold=NULL,maf_filter=FALSE){
	plot_dat<-predict_lnor(lnor=dat$lnor,se=dat$se,n=dat$ncase,p=dat$eaf,cases=dat$ncase,controls=dat$ncontrol)
	for(i in 1:ncol(plot_dat)){
		plot_dat[,i][plot_dat[,i] == "Inf" | plot_dat[,i] == "-Inf" ]<-NA
	}
	plot_dat<-plot_dat[complete.cases(plot_dat),]
	outfile_name<-unique(paste(dat$outcome,dat$study,dat$ID,sep="_"))
	outfile_name<-gsub(" ","_",outfile_name)
	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))

	Plot<-make_plot2(dat=plot_dat[,c("lnor_obs","lnor_pred2","eaf")],outfile_name=outfile_name,outcome_name=outcome_name,Xlab=Xlab,Ylab=Ylab,threshold=threshold,maf_filter=maf_filter)
	# png(paste0("~/fatty-acids/outcome_data/results/plots/predicted_lnor/",outfile_name,"_method_2.png"))
	# 	print(Plot)
	# dev.off()		
	return(Plot)
}


predict_lnor_and_plot_bias2<-function(dat=NULL,Xlab="",Ylab="",threshold=NULL,maf_filter=FALSE){
	plot_dat<-predict_lnor(lnor=dat$lnor,se=dat$se,n=dat$ncase,p=dat$eaf,cases=dat$ncase,controls=dat$ncontrol)
	for(i in 1:ncol(plot_dat)){
		plot_dat[,i][plot_dat[,i] == "Inf" | plot_dat[,i] == "-Inf" ]<-NA
	}
	plot_dat<-plot_dat[complete.cases(plot_dat),]
	plot_dat$bias<-((plot_dat$lnor_pred2-plot_dat$lnor_obs )/plot_dat$lnor_obs)*100	
	outfile_name<-unique(paste(dat$outcome,dat$study,dat$ID,sep="_"))
	outfile_name<-gsub(" ","_",outfile_name)
	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))

	Plot<-make_plot2(dat=plot_dat[, c("lnor_obs","bias","eaf")],outfile_name=outfile_name,outcome_name=outcome_name,Xlab=Xlab,Ylab=Ylab,linear_regression=FALSE,subtitle=subtitle,threshold=threshold,maf_filter=maf_filter)
	# png(paste0("~/fatty-acids/outcome_data/results/plots/predicted_lnor/",outfile_name,"_method_2.png"))
	# 	print(Plot)
	# dev.off()		
	return(Plot)
}

make_boxplot<-function(dat=NULL){
		p <- ggplot(dat, aes(x=lnor_obs, y=bias)) + 
	  		geom_boxplot()
	  	return(p)
	 }



predict_lnor<-function(lnor=NULL,se=NULL,n=NULL,p=NULL,cases=NULL,controls=NULL){
	t.stat<-lnor/se
	lnor_pred1<-sqrt(t.stat^2 / (t.stat^2 + n) / 2 / p*(1-p) * pi^2/3)
		  # sqrt(t^2 / (t^2 + n) / 2 / p*(1-p) * pi^2/3)
	lnor_pred1<-lnor_pred1*sign(t.stat)	

# 	exposed cases = (1-(1-maf)^2)*ncases
# unexposed cases = (1-maf)^2*ncases
# exposed controls = (1-(1-maf)^2)*ncontrols
# unexposed controls = (1-maf)^2*ncontrols

	# n1<-(p^2+2*p*(1-p))*cases
	n1<- (1-(1-p)^2)*cases #exposed cases
	n2<- (1-p)^2*cases #unexposed cases
	n3<- (1-(1-p)^2)*controls #exposed controls
	n4<- (1-p)^2*controls #unexposed controls
	lnse_pred2<-sqrt(1/n1+1/n2+1/n3+1/n4)	
	lnor_pred2<-t.stat*lnse_pred2
	lnor_pred2<-lnor_pred2/1.4 #method 2 overestimates log odds ratio by 40% on average
	Dat<-data.frame(do.call(cbind,list(lnor_obs=lnor,lnor_pred1=lnor_pred1,lnor_pred2=lnor_pred2,eaf=p)))
	return(Dat)
}

predlnor_model<-function(dat=NULL){
	dat<-predict_lnor(lnor=dat$lnor,se=dat$se,n=dat$ncase,p=dat$eaf,cases=dat$ncase,controls=dat$ncontrol)
	for(i in 1:ncol(dat)){
		dat[,i][dat[,i] == "Inf" | dat[,i] == "-Inf" ]<-NA
	}
	Model<-summary(lm(lnor_pred2~lnor_obs,dat))
	int<-Model$coefficients[1,1]
	slope<-Model$coefficients[2,1]
	return(c(int,slope))
}

make_plot_predlnor<-function(dat=NULL,Xlab="",Ylab="",linear_regression=TRUE,subtitle="",threshold=NULL,maf_filter=FALSE,bias=FALSE,Title_size=0,Title=NULL,Title_xaxis_size=10){

	# outfile_name<-unique(paste(dat$outcome,dat$study,dat$ID,sep="_"))
	# outfile_name<-gsub(" ","_",outfile_name)
	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))

	dat<-predict_lnor(lnor=dat$lnor,se=dat$se,n=dat$ncase,p=dat$eaf,cases=dat$ncase,controls=dat$ncontrol)
	dat$bias<-(dat$lnor_pred2-dat$lnor_obs )/dat$lnor_obs*100	
	# Bias1<-((dat$lnor_pred2-dat$lnor_obs )/dat$lnor_obs)*100	
	# Bias2<-(dat$lnor_pred2-dat$lnor_obs )/dat$lnor_obs*100	

	for(i in 1:ncol(dat)){
		dat[,i][dat[,i] == "Inf" | dat[,i] == "-Inf" ]<-NA
	}

	dat<-dat[complete.cases(dat),]
	
	if(is.null(Title)){
	# if(!is.null(outcome_name)){
		Title<-outcome_name
	}

	if(maf_filter){
		maf<-dat$eaf
		maf[maf>0.5]<-1-maf[maf>0.5]
		dat<-dat[maf>threshold,]
	}
	
	Colour<-rep("black",nrow(dat))
	Colour<-dat$eaf
	Colour[Colour>0.5]<-1-Colour[Colour>0.5]
	Colour[Colour<=0.10]<-"red"
	Colour[Colour>0.10 & Colour<=0.20]<-"orange"
	Colour[Colour>0.20 & Colour<=0.30]<-"green"
	Colour[Colour>0.30 & Colour<=0.40]<-"blue"
	Colour[Colour>0.40 & Colour<=0.50]<-"black"
	Shape<-rep(19,nrow(dat))

	dat$Y<-dat$lnor_pred2
	dat$X<-dat$lnor_obs

	if(bias){
		dat$X<-dat$bias
		Med<-round(summary(dat$bias)[3],1)
		p25<-round(summary(dat$bias)[2],1)
		p75<-round(summary(dat$bias)[5],1)
		Min<-round(summary(dat$bias)[1],1)
		Max<-round(summary(dat$bias)[6],1)
		subtitle<-paste0("Median bias=",Med,"% (IQR:",p25,"%, ",p75,"% | min=",Min,"%, max=",Max,"%)")
	}

	# dat$X<-dat[,1]
	# dat$Y<-dat[,2]
	if(bias){	
		linear_regression<-FALSE
	}

	if(linear_regression){	
		Model<-summary(lm(Y~X,dat))
		int<-Model$coefficients[1,1]
		slope<-Model$coefficients[2,1]
		subtitle<-paste0("intercept=",round(int,3)," | ","slope=",round(slope,3))
	}
		# Title_xaxis_size
	Plot<-ggplot2::ggplot(dat, ggplot2::aes(x=X, y=Y)) + ggplot2::geom_point(colour=Colour,shape=Shape) + ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
	ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
	ggplot2::scale_colour_manual(name = "MAF",
		# labels = c("0.11-0.20","0.21-0.3","0.31-0.4","0.41-0.5"),
        # values = c("orange","yellow","blue","black"))
        labels = c("0-0.10","0.11-0.20","0.21-0.3","0.31-0.4","0.41-0.5"),
        values = c("red","orange","yellow","blue","black"))
	return(Plot)
}


make_plot_predlnor3<-function(dat=NULL,Xlab="",Ylab="",linear_regression=TRUE,subtitle="",threshold=NULL,maf_filter=FALSE,bias=FALSE,Title_size=0,Title=NULL,Title_xaxis_size=10,legend=TRUE,standard_errors=FALSE){

	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))
	dat$bias<-(dat$lnor_sh-dat$lnor )/dat$lnor*100	
	# summary(dat$bias)

	# for(i in 1:ncol(dat)){
	# 	dat[,i][dat[,i] == "Inf" | dat[,i] == "-Inf" ]<-NA
	# }
	
	# dat<-dat[complete.cases(dat),]
	
	if(is.null(Title)){
	# if(!is.null(outcome_name)){
		Title<-outcome_name
	}

	if(maf_filter){
		maf<-dat$eaf
		maf[maf>0.5]<-1-maf[maf>0.5]
		dat<-dat[maf>=threshold,]
	}
	
	MAF<-rep("black",nrow(dat))
	MAF<-dat$eaf
	MAF[MAF>0.5]<-1-MAF[MAF>0.5]
	MAF[MAF<=0.10]<-"0.01-0.10"
	MAF[MAF>0.10 & MAF<=0.20]<-"0.11-0.20"
	MAF[MAF>0.20 & MAF<=0.30]<-"0.21-0.30"
	MAF[MAF>0.30 & MAF<=0.40]<-"0.31-0.40"
	MAF[MAF>0.40 & MAF<=0.50]<-"0.41-0.50"
	Shape<-rep(19,nrow(dat))

	dat$Y<-dat$lnor_sh
	dat$X<-dat$lnor

	if(standard_errors){
		dat$Y<-dat$se_sh
		dat$X<-dat$se
	}

	if(bias){
		dat$X<-dat$bias
		Med<-round(summary(dat$bias)[3],1)
		p25<-round(summary(dat$bias)[2],1)
		p75<-round(summary(dat$bias)[5],1)
		Min<-round(summary(dat$bias)[1],1)
		Max<-round(summary(dat$bias)[6],1)
		subtitle<-paste0("Median bias=",Med,"% (IQR:",p25,"%, ",p75,"% | min=",Min,"%, max=",Max,"%)")
	}

	# dat$X<-dat[,1]
	# dat$Y<-dat[,2]
	if(bias){	
		linear_regression<-FALSE
	}

	Values<-c("red","orange","purple","blue","black")
	Labels<-c("0.01-0.10","0.11-0.20","0.21-0.30","0.31-0.40","0.41-0.50")

	if(linear_regression){	
		Model<-summary(lm(Y~X,dat))
		int<-Model$coefficients[1,1]
		slope<-Model$coefficients[2,1]
		subtitle<-paste0("intercept=",round(int,3)," | ","slope=",round(slope,3))
	}		

	if(!legend){
		Plot<-ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=MAF))+theme(legend.position = "none")  + ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		ggplot2::scale_colour_manual(name = "MAF",
	        labels = Labels,
	        values = Values)
	}


	if(legend){
		Plot<-ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=MAF))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+scale_colour_grey(start=0.8,end=0.2)
		# +
		# 		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		# 		ggplot2::scale_colour_manual(name = "MAF",
		# 	        labels = Labels,
		# 	        values = Values)
	}

	return(Plot)
}


make_plot_predlnor2<-function(dat=NULL,Xlab="",Ylab="",linear_regression=TRUE,subtitle="",threshold=NULL,maf_filter=FALSE,bias=FALSE,Title_size=0,Title=NULL,Title_xaxis_size=10,legend=TRUE,standard_errors=FALSE){

	outcome_name<-unique(paste0(dat$outcome," | " ,dat$study," | ",dat$ID))
	dat$bias<-(dat$lnor_sh-dat$lnor )/dat$lnor*100	
	# summary(dat$bias)

	# for(i in 1:ncol(dat)){
	# 	dat[,i][dat[,i] == "Inf" | dat[,i] == "-Inf" ]<-NA
	# }
	
	# dat<-dat[complete.cases(dat),]
	
	if(is.null(Title)){
	# if(!is.null(outcome_name)){
		Title<-outcome_name
	}

	if(maf_filter){
		maf<-dat$eaf
		maf[maf>0.5]<-1-maf[maf>0.5]
		dat<-dat[maf>=threshold,]
	}
	
	MAF<-rep("black",nrow(dat))
	MAF<-dat$eaf
	MAF[MAF>0.5]<-1-MAF[MAF>0.5]
	MAF[MAF<=0.10]<-"0.01-0.10"
	MAF[MAF>0.10 & MAF<=0.20]<-"0.11-0.20"
	MAF[MAF>0.20 & MAF<=0.30]<-"0.21-0.30"
	MAF[MAF>0.30 & MAF<=0.40]<-"0.31-0.40"
	MAF[MAF>0.40 & MAF<=0.50]<-"0.41-0.50"
	Shape<-rep(19,nrow(dat))

	dat$Y<-dat$lnor_sh
	dat$X<-dat$lnor

	if(standard_errors){
		dat$Y<-dat$se_sh
		dat$X<-dat$se
	}

	if(bias){
		dat$X<-dat$bias
		Med<-round(summary(dat$bias)[3],1)
		p25<-round(summary(dat$bias)[2],1)
		p75<-round(summary(dat$bias)[5],1)
		Min<-round(summary(dat$bias)[1],1)
		Max<-round(summary(dat$bias)[6],1)
		subtitle<-paste0("Median bias=",Med,"% (IQR:",p25,"%, ",p75,"% | min=",Min,"%, max=",Max,"%)")
	}

	# dat$X<-dat[,1]
	# dat$Y<-dat[,2]
	if(bias){	
		linear_regression<-FALSE
	}

	Values<-c("red","orange","purple","blue","black")
	Labels<-c("0.01-0.10","0.11-0.20","0.21-0.30","0.31-0.40","0.41-0.50")

	if(linear_regression){	
		Model<-summary(lm(Y~X,dat))
		int<-Model$coefficients[1,1]
		slope<-Model$coefficients[2,1]
		subtitle<-paste0("intercept=",round(int,3)," | ","slope=",round(slope,3))
	}		

	if(!legend){
		Plot<-ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=MAF))+theme(legend.position = "none")  + ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
		ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
		ggplot2::scale_colour_manual(name = "MAF",
	        labels = Labels,
	        values = Values)
	}
	if(legend){
		Plot<-ggplot2::ggplot(dat) + ggplot2::geom_point(ggplot2::aes(x=X, y=Y,colour=MAF))+ ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size, face = "plain"))+
				ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+
				ggplot2::scale_colour_manual(name = "MAF",
			        labels = Labels,
			        values = Values)
	}

	return(Plot)
}


make_plot2<-function(dat=NULL,outfile_name=NULL,outcome_name=NULL,Xlab="",Ylab="",linear_regression=TRUE,subtitle="",threshold=threshold,maf_filter=maf_filter){
	
	Title<-""
	if(!is.null(outcome_name)){
		Title<-outcome_name
	}
	

	# Colour[which(dat[,eaf]>0.5)]<-"blue"
	# Colour[which(dat$alleles %in% c("AT","TA","GC","CG"))]<-"green"
	# Colour[which(dat[,ea] != dat[,ma])]<-"red"
	# Colour[which(dat$alleles %in% c("AT","TA","GC","CG"))]<-"green"
	
	# Shape[which(dat$alleles %in% c("AT","TA","GC","CG"))]<-1
	# dat$alleles[dat$rsid == "rs7774711"]
	if(maf_filter){
		maf<-dat$eaf
		maf[maf>0.5]<-1-maf[maf>0.5]
		dat<-dat[maf>threshold,]
	}
	
	if(any(names(dat) %in% "bias")){
		Med<-round(summary(dat$bias)[3],1)
		p25<-round(summary(dat$bias)[2],1)
		p75<-round(summary(dat$bias)[5],1)
		Min<-round(summary(dat$bias)[1],1)
		Max<-round(summary(dat$bias)[6],1)
		subtitle<-paste0("Median bias=",Med,"% (IQR:",p25,"%, ",p75,"% | min=",Min,"%, max=",Max,"%)")
	}

	Colour<-rep("black",nrow(dat))
	Colour<-dat$eaf
	Colour[Colour>0.5]<-1-Colour[Colour>0.5]
	Colour[Colour<=0.10]<-"red"
	Colour[Colour>0.10 & Colour<=0.20]<-"orange"
	Colour[Colour>0.20 & Colour<=0.30]<-"green"
	Colour[Colour>0.30 & Colour<=0.40]<-"blue"
	Colour[Colour>0.40 & Colour<=0.50]<-"black"
	Shape<-rep(19,nrow(dat))

	dat$X<-dat[,1]
	dat$Y<-dat[,2]	
	if(linear_regression){	
		Model<-summary(lm(Y~X,dat))
		int<-Model$coefficients[1,1]
		slope<-Model$coefficients[2,1]
		subtitle<-paste0("intercept=",round(int,3)," | ","slope=",round(slope,3))
	}
	# dat$X<-dat$eaf
	# dat$Y<-dat$eaf

	# Plot<-make_boxplot(dat=dat)
	
	Plot<-ggplot2::ggplot(dat, ggplot2::aes(x=X, y=Y)) + ggplot2::geom_point(colour=Colour,shape=Shape) + ggplot2::ggtitle(Title) +ggplot2::labs(y= Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold"))+ggplot2::scale_colour_manual(name = "MAF",
		# labels = c("0.11-0.20","0.21-0.3","0.31-0.4","0.41-0.5"),
        # values = c("orange","yellow","blue","black"))
        labels = c("0-0.10","0.11-0.20","0.21-0.3","0.31-0.4","0.41-0.5"),
        values = c("red","orange","yellow","blue","black"))
	return(Plot)
}


make_plot<-function(dat=NULL,outcome=NULL,eaf=NULL,maf=NULL,ea=NULL,ma=NULL,outcome_plot=outcome_plot,target_dat_population =NULL,reference_study=NULL,target_study=NULL){

	# Plots<-NULL
	# for(pop in unique(dat$ref_dat_population)){
		# print(pop)
		# dat1<-dat[dat$ref_dat_population==pop,]	
		dat1<-dat
		pop<-dat$ref_dat_population
		Ylab<-paste0("MAF in ",target_study)
		Xlab<-paste0("MAF in ",pop," (",unique(dat1$ref_study),")")
		pop2<-c("European","East Asian","African","American","South Asian","Global")
		i<-which(c("EUR","EAS","AFR","AMR","SAS","ALL") %in% pop)
		# Title<-pop2[i]
		# Title<-gsub("ALL","Global pop",Title)
		Title<-""

		Colour<-rep("black",nrow(dat1))
		Colour[which(dat1[,eaf]>0.5)]<-"blue"
		Colour[which(dat1[,eaf]>=0.6)]<-"red"
		Colour[dat1[,ea]!=dat1[,ma]]<-"red"	
		Diff<-abs(dat1[,eaf]-dat1[,maf])
		Colour[which(Diff>0.10)]<-"red"
		Shape<-rep(19,nrow(dat1))
		Shape[which(dat1$alleles %in% c("AT","TA","GC","CG"))]<-1
		dat1$eaf<-dat1[,eaf]
		dat1$maf<-dat1[,maf]
		
		Plots<-ggplot2::ggplot(dat1, ggplot2::aes(x=maf, y=eaf)) + ggplot2::geom_point(colour=Colour) +ggplot2::labs(y= Ylab, x =Xlab) +ggplot2::theme(axis.title=ggplot2::element_text(size=10)) 

		# Plots[[pop]]<-ggplot2::ggplot(dat1, ggplot2::aes(x=maf, y=eaf)) + ggplot2::geom_point(colour=Colour) +ggplot2::labs(y= Ylab, x =Xlab) +ggplot2::theme(axis.title=ggplot2::element_text(size=10)) 
		# face="bold"
		# +ggplot2::ggtitle(Title) + ggplot2::theme(plot.title = ggplot2::element_text(size = 1, face = "bold"))
		# axis.text=element_text(size=12),
	# }	

	# # reported_population<-ao$population[which(ao$id %in% File0 )]
	# Title2<-paste0(outcome_plot," | reported population: ",target_dat_population)

	# 		Plots2<-cowplot::plot_grid(plotlist=Plots)
	# 		title <- cowplot::ggdraw() + 
	# 			cowplot::draw_label(
	# 				Title2,
	# 				fontface = 'bold',
	# 				x = 0,
	# 				hjust = 0,
	# 				size = 12)  +
	# 			ggplot2::theme(
	# 			# add margin on the left of the drawing canvas,
	# 			# so title is aligned with left edge of first plot
	# 				plot.margin = ggplot2::margin(0, 0, 0, 7)
	# 				)

	# 		Plot3<-cowplot::plot_grid(title, Plots2,ncol = 1,
	# 		# rel_heights values control vertical title margins
	# 			rel_heights = c(0.05, 1)
	# 		)
		
		
	# Plot<-ggplot(dat, aes(x=maf, y=eaf)) + geom_point(colour=Colour,shape=Shape) +ggtitle(Title) +labs(y= Ylab, x =Xlab) + theme(plot.title = element_text(size = 10, face = "bold"))
	# return(Plot3)
	return(Plots)
}


predz_vs_obsz<-function(dat=NULL,Title_size=0,Title="",Ylab="",Xlab="",Title_xaxis_size=0){
	dat<-dat[abs(dat$p)<=1,]
	dat$z.p<-qnorm(dat$p/2,lower.tail=F)
	# dat[dat$z.p == "NaN", ]
	# dat$z<-dat$test_statistic
	dat$z.lnor<-abs(dat$lnor/dat$se)
	Colour<-"black"
	Cor<-cor(dat$z.p,dat$z.lnor)
	if(Cor<0.99){
		Colour<-"red"
	}
	subtitle<-paste0("Pearson correlation coefficient=",round(Cor,2))
	# +ggplot2::labs(y= Ylab, x =Xlab,)
	# plot(dat$z.p,dat$z.lnor)
	plot<-ggplot2::ggplot(dat, ggplot2::aes(x=z.p, y=z.lnor)) + ggplot2::geom_point(colour=Colour) + ggplot2::ggtitle(Title) + ggplot2::theme(axis.title=ggplot2::element_text(size=Title_xaxis_size))+ggplot2::labs(y=Ylab, x =Xlab,subtitle=subtitle) + ggplot2::theme(plot.title = ggplot2::element_text(size = Title_size),plot.subtitle = ggplot2::element_text(size = 8))
	
	return(plot)
}

# save(dat,file="~/MR_FattyAcids/data/simulated_lnor_seanharrison.Rdata")

# dat2=Test

pred_lnor_sh<-function(dat2=NULL){
	# load("~/harmonised_data.Rdata")
	# dat<-Dat
	# dat<-Dat[!is.na(Dat$eaf),]
	# Pos<-nrow(dat)/2
	# dat1<-dat[1:Pos,]
	# Pos<-Pos+1
	# dat2<-dat[Pos:nrow(dat),]
	if(!any(names(dat2) == "z")){
		dat2$z<-dat2$lnor/dat2$se
	}
	
	log_or<-NULL
	log_or_se<-NULL
	# ID<-NULL
	# snp<-NULL
	snp_n<-nrow(dat2)
	

	for(i in 1:snp_n)
	# for(i in 1:10)
	{
		# i<-1
		print(paste0("Analysing SNP " ,i, " of ",snp_n))
		# qui {

		#Odds of the outcome	
		n_ncase<-dat2$ncase[i]
		n_total<-dat2$ncase[i]+dat2$ncontrol[i]
		odds <- n_ncase/(n_total-n_ncase)
		# maf<-dat2$eaf[i]
		maf<-dat2$maf[i]
		Z<-dat2$z[i]

		# Ns given the MAF

		# p^2 + 2pq + q^2 = 1

		n0 <- n_total*maf^2
		n1 <- n_total*2*maf*(1-maf)
		n2 <- n_total*(1-maf)^2
		N <- n_total
		z <- Z

		#Simulate values of the log-OR, and estimate the Z score

		# preserve
		# clear
		n<-1:1000000
		if(z >= 0) 
		{
			# gen x = _n*0.000001
			x <- n*0.000001
		}else{
			x <- n*-0.000001
		}

		# gen n = _n 

		p0 <- 1/(1+exp(-(log(odds) - x*(n1+2*n2)/N)))
		p1 <- 1/(1+exp(-(log(odds) - x*(n1+2*n2)/N)-x))
		p2 <- 1/(1+exp(-(log(odds) - x*(n1+2*n2)/N)-2*x))

		a <- n0*p0*(1-p0)+n1*p1*(1-p1)+n2*p2*(1-p2)
		b <- n1*p1*(1-p1)+4*n2*p2*(1-p2)
		c <- (n1*p1*(1-p1)+2*n2*p2*(1-p2))^2

		se <- sqrt(a/(a*b-c))
		y = abs(x/se-z)

		y<- abs(x/sqrt((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)/((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)*(n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 4*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)-((n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 2*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)^2)))-z)

		complete<- 0

		j <- 0

		while(complete == 0) 
		{

			# qui su _y
			# Mean<-as.numeric(summary(y)[4])
			# qui su n if _y == r(min)
			Min<-as.numeric(summary(y))[1] #finds the minimum value of _y
			n<-which(y == Min)
			
			# If the minimum isn’t the last observation
			# if r(mean) < 1000000 
			if(n < 1000000)
			{
				complete <- 1
			}else{# If the minimum is the last observation, it hasn’t reached the actual minimum yet, so increase/decrease the observations
				j <- j+1
				if(z >= 0)
				{
					x <- n*0.000001 + (1000000-100)*0.000001*j
				}else{
					x <- -(n*0.000001 + (1000000-100)*0.000001*j)
				}
			}
		 
			# Updat2e the difference between the estimated and observed Z scores

			y <- abs(x/sqrt((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)/((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)*(n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 4*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)-((n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 2*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)^2)))-z)

			

		}

		#While loop complete, so minimum difference found

		# su _y
		# su _x if _y == r(min)
		# local x = r(mean)

		Min<-as.numeric(summary(y)[1])
		x_local<-summary(x[which(y==Min)])[4]
		
		y_se <- sqrt((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)/((n0*odds*exp(x*(n1+2*n2)/N)/(odds+exp(x*(n1+2*n2)/N))^2 + n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)*(n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 4*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)-((n1*odds*exp(x*(n1+2*n2-N)/N)/(odds+exp(x*(n1+2*n2-N)/N))^2 + 2*n2*odds*exp(x*(n1+2*n2-2*N)/N)/(odds+exp(x*(n1+2*n2-2*N)/N))^2)^2)))

		# su _y
		# su _y_se if _y == r(min)
		Min<-as.numeric(summary(y)[1])
		se_local<-summary(y_se[which(y==Min)])[4]

		# restore
		# End of simulation, bring back the original dat2a and updat2e the results

		# replace _log_or = `x’ in `i’
		# replace _log_or_se = `se’ in `i’

		log_or[[i]] <- x_local
		log_or_se[[i]] <- se_local
		# ID[[i]]<-dat2$ID[i]
		# snp[[i]]<-dat2$rsid[i]		
	}
	return(list(log_or,log_or_se))
}



prep_data_plot_predlnor_sh<-function(ID=NULL){
	load("~/fatty-acids/outcome_data/data/dat1_predor_sh.Rdata")
	load("~/fatty-acids/outcome_data/data/dat2_predor_sh.Rdata")
	Dat<-rbind(dat1,dat2)
	Dat<-Dat[Dat$ID == ID,]

	# 1
	# MAF<-Dat$eaf
	# MAF[MAF>0.5]<-1-MAF[MAF>0.5]
	# Dat<-Dat[MAF>=0.01,]
	MAF<-Dat$eaf
	MAF[MAF>0.5]<-1-MAF[MAF>0.5]
	Dat$MAC_case<-MAF*Dat$ncase*2
	Dat$MAC_control<-MAF*Dat$ncontrol*2
	Dat<-Dat[Dat$MAC_case>=50 & Dat$MAC_control>=50,]
	Dat<-Dat[Dat$lnor_sh <= 1.999 & Dat$lnor_sh>= -1.999,] #lnor_sh ==1.999 is an artifiact. 
	# Model<-summary(lm(lnor_sh~lnor,Dat))
	# Int<-Model$coefficients[1,1]
	# Slope<-Model$coefficients[2,1]
	return(Dat)
}
	
	
load_data_info<-function(){
	load("~/fatty-acids/outcome_data/data/Data_info_lnorsh.Rdata")
	load("~/fatty-acids/outcome_data/data/Datb_info_lnorsh.Rdata")
	load("~/fatty-acids/outcome_data/data/Datc_info_lnorsh.Rdata")
	load("~/fatty-acids/outcome_data/data/Datd_info_lnorsh.Rdata")
	Dat<-do.call(rbind,list(Data,Datb,Datc,Datd))
	# Dat[which(is.na(Dat$bias)),]
	Dat<-Dat[Dat$lnor != "-Inf",]
	Dat<-Dat[Dat$lnor != "Inf",]
	Dat<-Dat[Dat$lnor != 0,]
	Dat$bias<-(Dat$lnor_sh-Dat$lnor )/Dat$lnor*100	
	# Dat[which(Dat$bias == "Inf"),]
	return(Dat)
}


estimate_bias<-function(dat=NULL){
	dat<-dat[which(dat$lnor != "-Inf"),]
	dat<-dat[dat$lnor != "Inf",]
	dat<-dat[dat$lnor != 0,]
	dat$bias<-(dat$lnor_sh-dat$lnor )/dat$lnor*100	
	return(dat)
}


fix_info<-function(dat=NULL){
	Dat2<-load_data_info()
	IDS<-unique(Dat2$ID)
	# remove studies that have information on imputation score quality but which were filtered on score <0.80. 
	Dat<-Dat[!Dat$ID %in% IDS,]
	Dat3<-plyr::rbind.fill(Dat2,Dat)
	Dat<-Dat3
	return(Dat)
}
