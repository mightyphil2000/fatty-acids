library(ggplot2)
library(viridis)
library(data.table)
library(hrbrthemes)
# source("~/fatty-acids/outcome_data/scripts/harmonise_outcomes_functions.R")
source("~/fatty-acids/outcome_data/scripts/collate_and_qc_functions.R")
source("~/fatty-acids/outcome_data/scripts/descriptive_stats_functions.R")
Dat<-collate_dat(postqc=FALSE) 
# length(unique(Dat$ID))
dat<-read.table("~/MR_FattyAcids/presentations/mrc-ieu-researchers/data_analysis_for_ieuresearchers05052022.txt",head=TRUE,sep="\t",stringsAsFactors=FALSE,skip=1,fill=TRUE,quote="")

dat<-clean_data()

length(unique(dat$study.abbreviation))
dim(dat)

png("~/MR_FattyAcids/presentations/mrc-ieu-researchers/plot_ncases.png",width=800,
	height=500)
	hist(dat$cases)
dev.off()


dat3<-cancer_system()
sum(dat3$N)

# dat3<-cancer_group()
dat3<-cancer_source(miss_effect_allele=TRUE)
sum(dat3$N[grep("Open",dat3$source)])
dat3<-cancer_source2(miss_effect_allele=TRUE)
sum(dat3$N[grep("Open",dat3$source)])

dat3<-cancer_scope(miss_effect_allele=TRUE)
dat3<-cancer_ancestry()
dat3<-cancer_metadata()
dat3<-cancer_issue()
dat3<-cancer_issue2()

dat3<-cancer_consortia(qc_pipeline_fail=TRUE) #exclude cancer datasets that failed the QC pipeline
dim(dat3)
median(dat3$cases)
min(dat3$cases)
max(dat3$cases)
length(which(dat3$cases>10000))
length(which(dat3$cases>1000 & dat3$cases<10000))
length(which(dat3$cases>1000))
length(which(dat3$cases<1000))


cancer_types(Dat=dat,qc_pipeline_fail=TRUE)
cancer_types_cases()

total_cases_independent()


median(dat3$cases)
sum(dat3$cases)
min(dat3$cases)
max(dat3$cases)


dat3<-cancer_consortia2(qc_pipeline_fail=TRUE)


# N studies supplying one dataset
length(which(dat3$N==1))+1 #plus CGEMS Prostate cancer, which was excluded because did not report effect allele, other allele of EAF


dat3<-cancer_site()
sum(dat3$N)
dim(dat3)
dat3$ID<-factor(dat3$site,level=unique(dat3$site),ordered=T)

Plot<-ggplot(dat3, aes(y=N, x=ID)) + 
    geom_bar(position="dodge", stat="identity") +
    # scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=10,vjust= 0.5),axis.text.y = element_text(angle = 90,size=10))+
    theme(legend.position="none")+
    xlab("")+
    ylab("")

png("~/fatty-acids/outcome_data/results/plots/plot_site.png",width=500,
	height=1000)
	Plot
dev.off() 


dat3$ID<-factor(dat3$Study,level=unique(dat3$Study),ordered=T)

Plot<-ggplot(dat3, aes(y=N, x=ID)) + 
    geom_bar(position="dodge", stat="identity") +
    # scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=10,vjust= 0.5),axis.text.y = element_text(angle = 90,size=10))+
    theme(legend.position="none")+
    xlab("")+
    ylab("")

png("~/fatty-acids/outcome_data/results/plots/plot_study2.png",width=500,
	height=1000)
	Plot
dev.off() 

unique(dat$cancer)
unique(dat$site)

dat<-dat[order(dat$cases,decreasing=TRUE),]
dat$site[grep("respiratory",dat$site) ]<-"chest"

table(dat$system)
dat$system[dat$system %in% c("skeletal","endocrine","nervous","urinary")]<-"other"
dat$system[dat$system == "integumentary"]<-"skin"
dat$ID<-1:nrow(dat)
dat$ID<-factor(dat$ID,level=unique(dat$ID),ordered=T)
dat$site<-factor(dat$site,level=unique(dat$site),ordered=T)

Plot<-ggplot(dat, aes(fill=system,y=cases, x=ID)) + 
    geom_bar(position="dodge", stat="identity") +
    # scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=25,vjust= 0.5),axis.text.y = element_text(angle = 90,size=45))+
    theme(legend.position="none")+
    xlab("")+
    ylab("")+
    scale_x_discrete(limits=dat$ID, labels = dat$site)

png("~/fatty-acids/outcome_data/results/plots/plot_datasets.png",width=2500,
	height=4000)
	Plot
dev.off()    
    # scale_x_discrete(breaks=factor(dat$site))
    # ,limits=factor(dat$site)

# write.table(dat3,"~/MR_FattyAcids/presentations/mrc-ieu-researchers/top_consortia_biobanks.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


dat3$ID<-paste0(dat3$site,"/",dat3$study,"")
dat3$ID<-factor(dat3$ID,level=unique(dat3$ID),ordered=T)

dat3[dat3$cases == 2442,]
Plot<-ggplot(dat3, aes(y=cases, x=ID)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=45),axis.text.y = element_text(angle = 90,size=45))+
    xlab("")+
    ylab("")

png("~/fatty-acids/outcome_data/results/plots/plot_study.png",width=2500,
	height=4000)
	Plot
dev.off()

dim(dat3)


dat3$ID<-factor(dat3$QC_issue,level=unique(dat3$QC_issue),ordered=T)

Plot<-ggplot(dat3, aes(y=N, x=ID)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=15,vjust= 0.5))+
    xlab("")+
    ylab("")


png("~/MR_FattyAcids/presentations/mrc-ieu-researchers/plot_issue2.png",width=500,
	height=1000)
	Plot
dev.off()


Plot<-ggplot(dat3, aes(y=N, x=data)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=15))+
    xlab("")+
    ylab("")

png("~/MR_FattyAcids/presentations/mrc-ieu-researchers/plot_data2.png",width=500,
	height=1000)
	Plot
dev.off()


Plot<-ggplot(dat3, aes(y=N, x=population)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45,size=15))+
    xlab("")+
    ylab("")

png("~/MR_FattyAcids/presentations/mrc-ieu-researchers/plot_ancestry2.png",width=1000,
	height=500)
	Plot
dev.off()



Plot<-ggplot(dat3, aes(y=N, x=Summary_data_shared)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 45,size=15))+
    xlab("")+
    ylab("")

png("~/MR_FattyAcids/presentations/mrc-ieu-researchers/plot_shared.png",width=1000,
	height=500)
	Plot
dev.off()

dat3$ID<-factor(dat3$system,level=unique(dat3$system),ordered=T)
Plot2<-ggplot(dat3, aes(y=N, x=ID)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=15,vjust=0.5))+
    xlab("")+
    ylab("")


png("~/MR_FattyAcids/presentations/mrc-ieu-researchers/plot_site.png",width=500,
	height=1000)
	Plot2
dev.off()

dat3$source<-factor(dat3$source,level=unique(dat3$source),ordered=T)

Plot3<-ggplot(dat3, aes(y=N, x=source)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E") +
    ggtitle("") +
    theme_ipsum() +
    theme(axis.text.x = element_text(angle = 90,size=15))+
    xlab("")+
    ylab("")



png("~/MR_FattyAcids/presentations/mrc-ieu-researchers/plot_source2.png",width=500,
	height=1000)
	Plot3
dev.off()

