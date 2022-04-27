




# harmonise outcome names in data obtained by correspondence
Can_all5<-Can_all 
Can_all5$outcome[which(Can_all5$study == "GAMEON")]<-"Cancer (5 sites)"
Can_all5$outcome<-paste(toupper(substr(Can_all5$outcome, 1, 1)), substr(Can_all5$outcome, 2, nchar(Can_all5$outcome)), sep="")
Can_all5$id<-paste("COR:",as.numeric(as.factor(paste(Can_all5$outcome,Can_all5$study))),sep="")
Pos<-which(is.na(Can_all5$study))
Can_all5$study[Pos]<-Can_all5$id[Pos]
Can_all5$id_outcome<-paste(Can_all5$outcome,Can_all5$id)
Can_all5$outcome[Can_all5$outcome=="B-cell childhood acute lymphoblastic leukemia"]<-"Acute lymphoblastic leukemia"
Can_all5$outcome[Can_all5$outcome=="Childhood acute lymphoblastic leukemia"]<-"Acute lymphoblastic leukemia"
Can_all5$outcome[Can_all5$outcome=="Hepatocellular carcinoma in chronic hepatitis B virus carriers"]<-"Hepatocellular carcinoma"
Can_all5$Effect.Allele<-toupper(Can_all5$Effect.Allele)
Can_all5$Other.Allele<-toupper(Can_all5$Other.Allele)
Can_all5$study[Can_all5$outcome=="Lung cancer adjusted for chip" ]<-"UK Biobank adjusted for chip"
Can_all5$study[Can_all5$outcome=="Lung cancer unadjusted for chip" ]<-"UK Biobank unadjusted for chip"
Can_all5$outcome[Can_all5$outcome=="Lung cancer adjusted for chip" ]<-"Lung cancer"
Can_all5$outcome[Can_all5$outcome=="Lung cancer unadjusted for chip" ]<-"Lung cancer"


Cancer_dat<-rbind.fill(Cancer_mrbase,Can_all5)
# Cancer_dat$id<-paste(Cancer_dat$outcome,Cancer_dat$study)

# unique(paste(Cancer_dat$pmid,Cancer_dat$outcome)[which(is.na(Cancer_dat$eaf))])
# UVeal melaoma missing other allele28781888
# 1 SNP NA for Other Allele in Head and cancer
# mrbase IDs for largest sample set for each cancer site from TRICL,OCAC,BCAC and PRACTICAL
# Temp<-unique(Cancer_dat[!is.na(Cancer_dat$id_mrbase),c("outcome","ncase","id_mrbase","study")])
# Temp[order(Temp$ncase,decreasing=T),]

CRC<-unique(Can_all5$outcome[which(Can_all5$study=="GECCO")])
CRC.drop<-CRC[CRC != "Colorectal cancer"]
Can.drop<-c(CRC.drop,"Oral cancer","Oropharyngeal cancer","Glioblastoma","Non-glioblastoma glioma") #retain single largest analysis for each cancer study / consortium


Can_all6<-Can_all5[!Can_all5$outcome %in% Can.drop,] #restrict to single largest set for GECCO consortium and INHANCE consortium and for Glioma GWAS. 

Cancer_dat2<-rbind.fill(Cancer_mrbase_single,Can_all6) #restrict to single largest cancer at each cancer site in MR-Base and cancers obtained by correspondence. 

# Cancer_dat2$id<-paste(Cancer_dat2$outcome,Cancer_dat2$study)
# Cancer_dat2$study[Cancer_dat2$outcome=="Lung cancer adjusted for chip" ]<-"UK Biobank adjusted for chip"
# Cancer_dat2$study[Cancer_dat2$outcome=="Lung cancer unadjusted for chip" ]<-"UK Biobank unadjusted for chip"
# Cancer_dat2$outcome[Cancer_dat2$outcome=="Lung cancer adjusted for chip" ]<-"Lung cancer"
# Cancer_dat2$outcome[Cancer_dat2$outcome=="Lung cancer unadjusted for chip" ]<-"Lung cancer"

# drop duplicate cancers retaining the duplicate with the most cases
Cancer_dat2<-Cancer_dat2[order(Cancer_dat2$ncase,decreasing=T),] 
Cancer_dat2$outcome[Cancer_dat2$outcome=="Cancer" & Cancer_dat2$study == "UK Biobank"]<-"Cancer (UK Biobank)"
Studies<-unique(Cancer_dat2[,c("ncase","outcome","id")])
Studies.drop<-Studies$id[duplicated(Studies$outcome)]
Cancer_dat3<-Cancer_dat2[!Cancer_dat2$id %in% Studies.drop,]

# All cleaned cancer outcome data. Including all subtypes and duplicate outcomes (from supposedly separate / independent studies)

Cancer_dat$original_outcome<-Cancer_dat$outcome
Cancer_dat$outcome<-Cancer_dat$id_outcome
Cancer_dat<-Cancer_dat[,names(Cancer_dat) != "note"]

# sort(unique(Cancer_dat$outcome))
# sort(unique(Cancer_dat$original_outcome))
save(Cancer_dat,file="~/fatty-acids/data/Cancer_dat.Rdata")
write.table(Cancer_dat,"~/fatty-acids/data/Cancer_dat.txt",sep="\t",col.names=T,row.names=F,quote=F) 

#pruned so that tumour subtypes are removed. Single largest cancer at each site retained. Performed within study, so there are duplicate cancer outcomes between studies. UK Biobank Lung cancer adjusted and unadjusted for chip retained
Cancer_dat2$original_outcome<-Cancer_dat2$outcome
Cancer_dat2$outcome<-Cancer_dat2$id_outcome
Cancer_dat2<-Cancer_dat2[,names(Cancer_dat2) != "note"]


save(Cancer_dat2,file="~/fatty-acids/data/Cancer_dat2.Rdata") 
write.table(Cancer_dat2,"~/fatty-acids/data/Cancer_dat2.txt",sep="\t",col.names=T,row.names=F,quote=F) 

#pruned to have no outcome duplicates and largest duplicate retained. Each cancer outcome is unique
Cancer_dat3$original_outcome<-Cancer_dat3$outcome
Cancer_dat3$outcome<-Cancer_dat3$id_outcome
Cancer_dat3<-Cancer_dat3[,names(Cancer_dat3) != "note"]


save(Cancer_dat3,file="~/fatty-acids/data/Cancer_dat3.Rdata") 
write.table(Cancer_dat3,"~/fatty-acids/data/Cancer_dat3.txt",sep="\t",col.names=T,row.names=F,quote=F) 


