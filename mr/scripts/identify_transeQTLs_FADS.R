setwd("/projects/MRC-IEU/users/ph14916/eQTLGen/trans")

# reg<-read.table("~/fatty-acids-mr/instruments/regional_table_Dec18.txt",sep="\t",head=T,stringsAsFactors=F)
# snps_all<-readLines("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/ALLSNPs_Dec18_europeans.txt")
# snps_all_tab<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/targetSNPs_table_Dec18_europeans.txt",sep="\t",head=T,stringsAsFactors=F)
snps_all_tab<-read.table("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/targetSNPs_table_Dec18.txt",sep="\t",head=T,stringsAsFactors=F)


snps_all<-snps_all_tab$SNP

genes<-c("fads","scd","ELOVL2")
Pos<-unique(unlist(lapply(1:length(genes),FUN=function(x)
		grep(genes[x],snps_all_tab$gene,ignore.case=T,invert=F
			))))
Regions<-unique(snps_all_tab$region[Pos])
snps_all_tab2<-snps_all_tab[!snps_all_tab$region %in% Regions,]
snps_all_tab3<-snps_all_tab2[grep("polyunsaturated",snps_all_tab2$type,ignore.case=T),]


fads_sig<-read.table("fads_sig.txt",sep="\t",head=T,stringsAsFactors=F)
fads_all<-read.table("fads_all.txt",sep="\t",head=T,stringsAsFactors=F)

snps<-fads_sig$SNP
fads_all2<-fads_all[fads_all$SNP %in% snps,]

fads_all2[fads_all2$GeneSymbol == "FADS1",]
fads_all2[fads_all2$GeneSymbol == "FADS2",]
fads_all2[fads_all2$GeneSymbol == "FADS3",]

fads_all2<-fads_all[fads_all$SNP %in% snps_all,]
order(fads_all2$Pvalue)

snps_all_miss<-snps_all[!snps_all %in% unique(fads_all2$SNP)]
proxy<-snps_all_tab[snps_all_tab$SNP %in% snps_all_miss,c("ProxySNP","pop")]
proxy<-proxy[proxy$pop=="EUR",]
proxy<-unique(unlist(strsplit(proxy$ProxySNP,split="; ")))
proxy<-proxy[proxy!=""]

fads_all_proxy<-fads_all[fads_all$SNP %in% proxy,]

P.thresh<-0.05/(161*2)

fads_all3<-fads_all2[fads_all2$Pvalue<P.thresh,]
fads_all4<-fads_all_proxy[fads_all_proxy$Pvalue<P.thresh,]

# fads_all2$GeneSymbol[fads_all2$Pvalue<P.thresh]

snps_eqtl<-unique(c(fads_all3$SNP,fads_all4$SNP))


unique(snps_all_tab3$SNP) #SNPs associated with PUFAs within thr 161 SNP set excluding FADS, ELOVL2 and SCD regions
region1<-snps_all_tab$region[snps_all_tab$SNP %in% snps_eqtl]
snps_eqtl2<-snps_eqtl[!snps_eqtl %in% snps_all_tab$SNP]
Pos<-unique(unlist(lapply(1:length(snps_eqtl2),FUN=function(x)
	grep(snps_eqtl2[x],snps_all_tab$ProxySNP)
	)))
region2<-unique(snps_all_tab$region[Pos])
Regions<-c(region1,region2)
snps_all_tab4<-snps_all_tab3[!snps_all_tab3$region %in% Regions,]
dim(snps_all_tab4)
length(unique(snps_all_tab4$region))
unique(snps_all_tab4[,c("SNP","region","author","pval","gene")])
region_dup<-unique(snps_all_tab4$region[duplicated(snps_all_tab4$region)])
snps_all_tab4<-snps_all_tab4[order(snps_all_tab4$region),]
snps_all_tab4[snps_all_tab4$region %in% region_dup,c("SNP","region","gene")]
Dat_list<-NULL
snps_all_tab4$gene[snps_all_tab4$region %in% region_dup]
for(i in 1:length(region_dup)){
	print(region_dup[i])
	Dat<-snps_all_tab4[snps_all_tab4$region %in% region_dup[i],]
	Dat$gene<-paste(unique(unlist(strsplit(Dat$gene,split="; "))),collapse="; ")
	Dat$SNP<-paste(unique(unlist(strsplit(Dat$SNP,split="; "))),collapse="; ")
	Dat$pval<-paste(unique(unlist(strsplit(Dat$pval,split="; "))),collapse="; ")
	Dat$author<-paste(unique(unlist(strsplit(Dat$author,split="; "))),collapse="; ")
	Dat$trait<-paste(unique(unlist(strsplit(Dat$trait,split="; "))),collapse="; ")
	Dat$type<-paste(unique(unlist(strsplit(Dat$type,split="; "))),collapse="; ")
	Dat_list[[i]]<-Dat
}
Dat<-do.call(rbind,Dat_list)


snps_all_tab5<-snps_all_tab4[!snps_all_tab4$region %in% region_dup,]
Dat2<-rbind(Dat,snps_all_tab5)

rs10740118_region<-snps_all_tab$region[grep("rs10740118",snps_all_tab$ProxySNP)]
Dat2[Dat2$region == rs10740118_region,]
Dat2$gene[Dat2$region ==rs10740118_region]

# region	Chr	SNP.GRCh38.p12	SNP.GRCh37
# gene	trait	type
# author
# start	end
Dat2<-Dat2[!duplicated(Dat2$region),c("region","Chr","start","end","gene","pval","SNP","trait","author")]
# dim(Dat2)
write.table(Dat2,"/projects/MRC-IEU/users/ph14916/fatty_acids_summary/regions_sensitivityhp.txt",sep="\t",col.names=T,row.names=F,,quote=F) #these contains regions previously associated with PUFAs excluding FADS, SCD and ELOVL2 regions and also regions containing a trans eQTL for FADS1 or FADS2

scp ph14916@epi-franklin.epi.bris.ac.uk:/projects/MRC-IEU/users/ph14916/fatty_acids_summary/regions_sensitivityhp.txt .
