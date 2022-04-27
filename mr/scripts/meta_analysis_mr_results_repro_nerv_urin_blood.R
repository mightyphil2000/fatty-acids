library(meta)


source("~/fatty-acids/mr/scripts/mr_functions.R")

mr_res<-format_metareg_v2()
sort(mr_res$cancer)
# mr_res<-format_metareg_v2(restrict_to_european_studies=TRUE)

# mr_res<-mr_res[mr_res$population != "European",]
M<-create_correlation_matrix(study=mr_res)
M<-M$matrix

rownames(M)<-mr_res$study.abbreviation
colnames(M)<-mr_res$study.abbreviation
# corr_results_list2<-res$corr_results_list2
mr_res$se_decoupled<-decoupling(s=mr_res$se,C=M)

########################################
# Reproductive, blood, nervous, urinary#
########################################

Pos<-which(mr_res$system=="Blood")
length(Pos)
Model<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se[Pos]))
m1<-extract_meta_results(dat=Model)
m1$system<-"Blood cancers"

Model_decoupled<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se_decoupled[Pos]))
m2<-extract_meta_results(dat=Model_decoupled)
m2$system<-"Blood cancers decoupled standard error"
m_blood<-rbind(m1,m2)

paste0(mr_res$cancer[Pos], " (FAMRC ID:",mr_res$ID[Pos],")")
sum(as.numeric(mr_res$cases[Pos]))
sum(as.numeric(mr_res$controls[Pos]))


# reproductive
Pos<-which(mr_res$system=="Reproductive")
Model<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se[Pos]))
m1<-extract_meta_results(dat=Model)
m1$system<-"Reproductive cancers"

Model_decoupled<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se_decoupled[Pos]))
m2<-extract_meta_results(dat=Model_decoupled)
m2$system<-"Reproductive cancers decoupled standard error"
m_repro<-rbind(m1,m2)

paste0(mr_res$cancer[Pos], " (FAMRC ID:",mr_res$ID[Pos],")")
sum(as.numeric(mr_res$cases[Pos]))
sum(as.numeric(mr_res$controls[Pos]))

# Nervous
Pos<-which(mr_res$system=="Nervous")
Model<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se[Pos]))
m1<-extract_meta_results(dat=Model)
m1$system<-"Nervous system cancers"

Model_decoupled<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se_decoupled[Pos]))
m2<-extract_meta_results(dat=Model_decoupled)
m2$system<-"Nervous system cancers decoupled standard error"
m_nerv<-rbind(m1,m2)
paste0(mr_res$cancer[Pos], " (FAMRC ID:",mr_res$ID[Pos],")")
sum(as.numeric(mr_res$cases[Pos]))
sum(as.numeric(mr_res$controls[Pos]))
# Urinary 
Pos<-which(mr_res$system=="Urinary")
Model<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se[Pos]))
m1<-extract_meta_results(dat=Model)
m1$system<-"Urinary system cancers"

Model_decoupled<-metagen(TE=as.numeric(mr_res$b[Pos]),seTE=as.numeric(mr_res$se_decoupled[Pos]))
m2<-extract_meta_results(dat=Model_decoupled)
m2$system<-"Urinary system cancers decoupled standard error"

paste0(mr_res$cancer[Pos], " (FAMRC ID:",mr_res$ID[Pos],")")

sum(as.numeric(mr_res$cases[Pos]))
sum(as.numeric(mr_res$controls[Pos]))
names(mr_res)
mr_res$system[Pos]

m_uri<-rbind(m1,m2)

m_null_cancers<-do.call(rbind,list(m_blood,m_repro,m_nerv,m_uri))

write.table(m_null_cancers,"~/fatty-acids/mr/results/meta_reg_nullcancers.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

