
find_snp_positions<-function(SNP=NULL,bed_format=FALSE,Chr=NULL,Position=NULL){
	if(!is.null(SNP)){
		Chr<-ref$V1[which(ref$V2 == SNP)]
		ref1<-ref[ref$V1 == Chr,]
		pos<-which(ref1$V2 == SNP)
		minbp<-ref1$V4[pos] - 500000
		maxbp<-ref1$V4[pos] + 500000
		if(!bed_format){
			snplist<-ref1$V2[ref1$V4> minbp & ref1$V4<maxbp ]
		}
		if(bed_format){
			snplist<-ref1[ref1$V4> minbp & ref1$V4<maxbp,c("V1","V4") ]
			snplist<-paste(snplist$V1,"_",snplist$V4,"_",sep="")
		}
	}
	if(is.null(SNP)){
		ref1<-ref[ref$V1 == Chr,]
		if(length(Position) == 2){
			minbp<-Position[1]
			maxbp<-Position[2]
			maxbp-minbp
			snplist<-ref1$V2[ref1$V4> minbp & ref1$V4<maxbp ]
		}
	}
	return(snplist)
}

extract_data<-function(snplist=NULL,type=NULL,file_dir=NULL,wk_dir=NULL,file_list=NULL,exact_match=TRUE){
	library(plyr)
	# library(dplyr)
	if(is.null(wk_dir))
	{
		setwd("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/colocalisation")
	}else{
		setwd(wk_dir)
	}
	if(type=="raw")	
	{
		if(is.null(file_dir))
		{
			Dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/imputation_harmonisation/notimputed"
		}else
		{
			Dir<-file_dir
		}
	}
	
	if(type=="imputed")
	{
		if(is.null(file_dir))
		{
			Dir<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/imputation_harmonisation"
		}else
		{
			Dir<-file_dir
		}
	}
	
	if(type == "not_fatty_acids"){
		Dir<-file_dir
	}

	if(is.null(file_list)){
		Files<-dir(Dir)[grep("tab",dir(Dir))]	
		Files<-paste(Dir,Files,sep="/")
	}else{
		Files<-paste(Dir,file_list,sep="/")
	}

	if(is.null(snplist))
	{
		snplist2<-"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_server.txt"
	}else
	{
		snplist2<-snplist
	}

	if(!exact_match){
		system(paste("fgrep -f ",snplist2,Files[1],"> temp1.txt"))		
	}
	
	if(exact_match){
		system(paste("fgrep -wf ",snplist2,Files[1],"> temp1.txt"))	
	}
	# system("wc temp1.txt")
	File1<-paste("\"",Files[1],"\"",sep="")
	system(paste("awk -F, '{$(NF+1)=",File1,";}1' OFS=\"\t\" temp1.txt > output.txt",sep=""))
	system(paste("head -1", Files[1],"> filehead.txt"))
	system("awk -F, '{$(NF+1)=\"file\";}1' OFS=\"\t\" filehead.txt > fhtemp.txt")
	system(paste("awk -F, '{$(NF+1)=",File1,";}1' OFS=\"\t\" fhtemp.txt > filehead.txt",sep=""))

	

	for(i in 2:length(Files))
	# for(i in 2:5)
	{
		print(Files[i])
		# system(paste("gunzip",Files[i],sep=" "))
		# File<-gsub("tab.gz","tab",Files[i])
		if(!exact_match){
			system(paste("fgrep -f ",snplist2,Files[i],"> temp1.txt"))	
		}
		if(exact_match){
			system(paste("fgrep -wf ",snplist2,Files[i],"> temp1.txt"))	
		}
		File1<-paste("\"",Files[i],"\"",sep="")
		system(paste("awk -F, '{$(NF+1)=",File1,";}1' OFS=\"\t\" temp1.txt >> output.txt",sep=""))
		system(paste("head -1", Files[i],"> fhtemp.txt"))
		# system(paste("head -1", Files[i],"> filehead.txt"))
		system("awk -F, '{$(NF+1)=\"file\";}1' OFS=\"\t\" fhtemp.txt > fhtemp_2.txt")
		# system("awk -F, '{$(NF+1)=\"file\";}1' OFS=\"\t\" filehead.txt > fhtemp.txt")
		# system("mv fhtemp2.txt fhtemp.txt")
		system(paste("awk -F, '{$(NF+1)=",File1,";}1' OFS=\"\t\" fhtemp_2.txt >> filehead.txt",sep=""))
			# system(paste("gzip",File))
	}


	# system(paste("head -1", Files[1],"> filehead.txt"))
	# system("cat filehead.txt")
	# system("awk -F, '{$(NF+1)=\"file\";}1' OFS=\"\t\" filehead.txt > fhtemp.txt")
	# system("cat fhtemp.txt")

	# system(paste("head -1", Files[i],"> fhtemp.txt"))
	# system("cat fhtemp.txt")	
	# system("awk -F, '{$(NF+1)=\"file\";}1' OFS=\"\t\" fhtemp.txt > fhtemp_2.txt")
	# system("cat fhtemp_2.txt")	

	
	# system("cat output.txt")
	

	L<-NULL 
	L.d<-NULL
# Files<-gsub("tab.gz","tab",Files)
	for(i in 1:length(Files))
	# for(i in 1:5)
	{
		print(i)
		print(Files[i])
		system(paste("grep -w \"",Files[i],"\" output.txt > temp1.txt",sep=""))
		Res1<-readLines("filehead.txt")
		ColNames<-unlist(strsplit(Res1[i],split="\t"))
		Res<-read.table("temp1.txt",sep="\t",head=F,stringsAsFactors=F)
		if(unique(Res[,ncol(Res)]) != ColNames[length(ColNames)]) stop(paste("column headers incorrect! for", Files[i]))
		ColNames<-ColNames[1:length(ColNames)-1]
		names(Res)<-ColNames
		if(class(Res) != "data.frame") stop(paste("not a dataframe"))
		if(all(!duplicated(Res$snp)))
		{
			print("no duplicate SNPs")
		}
		# dup<-Res$snp[duplicated(Res$snp)]
		# Res[Res$snp %in% dup,c("type","beta","se","p")]
		if(any(duplicated(Res$snp))){ #duplicated SNPs only occur in imputed files (kettunen was not imputed by DIST and contains no duplicates and no column called type. aLl duplicates occur in inchianti files)
			print("duplicated")
			Res<-Res[order(Res$type),]
			Dups<-unique(Res$snp[duplicated(Res$snp)])
			print(paste("duplicated SNPs:" ,Dups))
			Res.dup<-Res[Res$snp %in% Dups,]
			Res<-Res[!duplicated(Res$snp),] #drop duplicates. Keep type 0 SNPs. Type 2 SNPs missing beta and se inch/ tanaka files
			L.d[[i]]<-Res.dup
		}
		L[[i]]<-Res
	}

	L<-do.call(rbind.fill,L)
	L$id<-paste(L$snp,L$file)
	return(L)
}
	
#extract_pvalue=c(5e-8,7) #first vector is the pvalue threshold and the second vector is the column containing the p values
extract_data2<-function(snplist=NULL,File=NULL,exact_match=FALSE,file_sep="\t",extract_pvalue=NULL,Head=TRUE,out_dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/"){

	if(!is.null(extract_pvalue)){
		threshold<-extract_pvalue[1]
		column<-extract_pvalue[2]
		sys.cmd<-paste("awk -F \"",file_sep,"\" '$",column, " < ", threshold,"' ",File," > ",out_dir,"output.txt",sep="")
		system(sys.cmd)
		sys.cmd<-paste("head -1 ",File," > ",out_dir,"head.txt",sep="")
		system(sys.cmd)
		sys.cmd<-paste("cat ",out_dir,"head.txt ",out_dir,"output.txt > ",out_dir,"output_head.txt",sep="")
		system(sys.cmd)
		Res<-read.table(out_dir,"output_head.txt",head=T,sep="\t",stringsAsFactors=F)
	}

	if(is.null(extract_pvalue)){
	    if(exact_match){
	    	if(length(snplist)!=1){
	    		write.table(snplist,paste(out_dir,"temp.txt",sep=""),col.names=F,row.names=F,quote=F)
	    		snplist<-paste(out_dir,"temp.txt",sep="")
	    	}
	    	if(sum(grep("gz",File))==1){
	    		system(paste("zcat ",File," | fgrep -wf ",snplist," > ",out_dir,"output.txt",sep=""))	
	    	}else{
	      		system(paste("fgrep -wf ",snplist," ", File," > ",out_dir,"output.txt",sep=""))
	      	}	      	      

	        if(Head){
	        	if(sum(grep("gz",File))==1){
	        		system(paste("zcat ",File," | head -1 > ",out_dir,"filehead2.txt",sep=""))
	        	}else{	        		
		        	system(paste("head -1 ",File," > ",out_dir,"filehead.txt",sep=""))
		        }
	        	system(paste("cat ",out_dir,"filehead.txt ",out_dir,"output.txt > ",out_dir,"output_head.txt",sep=""))	
	        }
	        
	        if(!Head){
	        	sys.cmd<-paste("mv ",out_dir,"output.txt ",out_dir,"output_head.txt",sep="")	        	
	        	system(sys.cmd)
	        }
	    }
	    if(!exact_match){
	        # system(paste("head",snplist))
	        # system(paste("grep 11_61736411",File))
	        system(paste("fgrep -f ",snplist," ", File," > ",out_dir,"output.txt",sep="")) 
	        system(paste("head -1 ",File," > ",out_dir,"filehead.txt",sep=""))
	        system(paste("cat ",out_dir,"filehead.txt ",out_dir,"output.txt > ",out_dir,"output_head.txt",sep=""))
	    }
	    paste(out_dir,"output_head.txt",sep="")
	    if(file_sep=="\t"){
			Res<-read.table(paste(out_dir,"output_head.txt",sep=""),sep="\t",head=Head,stringsAsFactors=F)    	
	    }

	    if(file_sep==" "){
		    Res<-read.table(paste(out_dir,"output_head.txt",sep=""),sep=" ",head=T,stringsAsFactors=F)
		}

	    File_name<-unlist(strsplit(File,split="/"))
	    Res$file<-File_name[length(File_name)]
	}
	return(Res)
}   

#developed for extraction of eQTLs from eQTL gen where all genes and all cist eQTLS are in a single file
# gene="ENSG00000134824"
# snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_coloc.txt"
# File="combined-eQTLs_full.EAF.beta.se.chr.pos.txt"

extract_data3<-function(gene=NULL,snplist=NULL,File=NULL,file_sep="\t",out_dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/"){
	
	Sys.cmd<-paste("grep -w ",gene," ", File," | grep -wf ",snplist, " > ",out_dir,"output.txt", sep="")
	system(Sys.cmd)
 	system(paste("head -1 ",File," > ",out_dir,"filehead.txt",sep=""))
	system(paste("cat ",out_dir,"filehead.txt ",out_dir,"output.txt > ",out_dir,"output_head.txt",sep=""))	   
	Res<-read.table(paste(out_dir,"output_head.txt",sep=""),sep="\t",head=TRUE,stringsAsFactors=F)
	return(Res)
}
	

						   
extract_snps<-function(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt",File=NULL,exact_match=FALSE,file_sep="\t",out_dir="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/",Test.gz=FALSE,fill=FALSE,Comment = "#",Head=TRUE){

	if(length(snplist)==1){
		snplist<-readLines(snplist)
		if(any(duplicated(snplist))) message("duplicate SNPs present in snplist")
	}
	# if(length(snplist)>1){
	write.table(unique(snplist),paste(out_dir,"temp.txt",sep=""),col.names=F,row.names=F,quote=F)
	snplist<-paste(out_dir,"temp.txt",sep="")
	# }
	
	if(Test.gz){
		system(paste("zcat ", File," | fgrep -wf ",snplist," > ", out_dir,"output.txt",sep="")) 
		# system(paste("zgrep -wf ",snplist," ", File," > ",out_dir,"output.txt",sep="")) 
        system(paste("zcat ",File," | head -1 > ",out_dir,"filehead.txt",sep=""))
        system(paste("cat ",out_dir,"filehead.txt ",out_dir,"output.txt > ",out_dir,"output_head.txt",sep=""))        
	}

	if(exact_match & !Test.gz){
        # system(paste("grep -w ",snplist," ", File," > ",out_dir,"output.txt",sep="")) 
        system(paste("fgrep -wf ",snplist," ", File," > ",out_dir,"output.txt",sep="")) 
        system(paste("head -1 ",File," > ",out_dir,"filehead.txt",sep=""))
        system(paste("cat ",out_dir,"filehead.txt ",out_dir,"output.txt > ",out_dir,"output_head.txt",sep="")) 
        # system("cat /projects/MRC-IEU/users/ph14916/fatty_acids_summary/output_head.txt")
    }

    if(!exact_match & !Test.gz){
        # system(paste("head",snplist))
        # system(paste("grep 11_61736411",File))
        # system("cat /projects/MRC-IEU/users/ph14916/fatty_acids_summary/temp.txt")
        system(paste("fgrep -f ",snplist," ", File," > ",out_dir,"output.txt",sep="")) 

        system(paste("head -1 ",File," > ",out_dir,"filehead.txt",sep=""))
        system(paste("cat ",out_dir,"filehead.txt ",out_dir,"output.txt > ",out_dir,"output_head.txt",sep=""))
    }
    
  	if(fill){
	    Res<-read.table(paste(out_dir,"output_head.txt",sep=""),sep=file_sep,head=Head,stringsAsFactors=F,fill=T,comment.char = Comment)
	}else{
		Res<-read.table(paste(out_dir,"output_head.txt",sep=""),sep=file_sep,head=Head,stringsAsFactors=F, comment.char = Comment)
	}
    File_name<-unlist(strsplit(File,split="/"))
    if(nrow(Res)!=0) Res$file<-File_name[length(File_name)]

    if(ncol(Res)==2){
    	warning("Column names look like this. Will try to read file using load_plink function")
    	print(names(Res))
    	Res<-read_plink(File=paste(out_dir,"output_head.txt",sep=""))
    }
    return(Res)
}    

# system(paste("head ", snplist,sep=""))
# system(paste("fgrep 11_61569830 ", File,sep=""))

convert_to_bed<-function(snplist="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/snplist_Europeans_rsidsonly2.txt",read_ref=FALSE ){
    if(read_ref){
        ref<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/data_maf0.01_rs.bim",sep="\t",head=F,stringsAsFactors=F)
    }
    snps<-readLines(snplist)
    ref1<-ref[ref$V2 %in% snps,]
    ref1$bed<-paste(ref1$V1,"_",ref1$V4,sep="")
    return(ref1)
}


read_plink<-function(File=NULL){
	ref<-readLines(File)
	Dat<-NULL
	for(i in 1:length(ref)){
		# print(i)
		A<-unlist(strsplit(ref[i],split=" "))
		A<-A[A!=""]
		Dat[[i]]<-A
	}
	Dat2<-data.frame(do.call(rbind,Dat),stringsAsFactors=F)
	names(Dat2)<-Dat2[1,]
	Dat2<-Dat2[2:nrow(Dat2),]
	return(Dat2)
}