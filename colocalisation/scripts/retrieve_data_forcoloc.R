# remotes::install_github("mrcieu/gwasvcf")
# conda activate r4
setwd("~")
library(gwasvcf)
library(VariantAnnotation)
library(ieugwasr)
associations()

ao<-gwasinfo(access_token = "ya29.a0AfH6SMCcTwFs44nLf4fa0aV5tbB14PLl_8g71OJKlUNE5xY5aQ9b892LVLJJo_BquTTaKSABLoNcsyPuNucYFLlRbClvE6Ae2deldwkzXEZBE-bTO3LsWZmCVZsALSJ-sD4PeLjbt3Y8OnuxYKN7BrCOOaucSLO0iNclEg")
id1<-ao$id[ao$trait=="Malignant neoplasm of respiratory system and intrathoracic organs"]
id2<-ao$id[ao$trait=="Cancer code, self-reported: basal cell carcinoma"]
id3<-ao$id[ao$trait=="Malignant neoplasm of skin"]
id4<-ao$id[ao$consortium=="TRICL"]
IDS<-c(id1,id2,id3,id4,"finn-a-LUNG_CANCER")
tab1 <- ieugwasr::associations(id=IDS, variants="11:61043499-62159523",access_token = "ya29.a0AfH6SMCcTwFs44nLf4fa0aV5tbB14PLl_8g71OJKlUNE5xY5aQ9b892LVLJJo_BquTTaKSABLoNcsyPuNucYFLlRbClvE6Ae2deldwkzXEZBE-bTO3LsWZmCVZsALSJ-sD4PeLjbt3Y8OnuxYKN7BrCOOaucSLO0iNclEg") 

cancer_dat<-tab1
save(cancer_dat,file="~/cancer_dat.RData")

cd ~/fatty-acids/colocalisation/data/cancer
# scp ph14916@bc4login.acrc.bris.ac.uk:~/cancer_data.Rdata . 

rsync -avzh --stats --progress ph14916@bc4login.acrc.bris.ac.uk:~/cancer_dat.RData .  


e.g. in R
library(ieugwasr)
ao<-gwasinfo()
ao$trait[grep("PGC",ao$consortium)][1]
[1] "Autism Spectrum Disorder"
id<-ao$id[grep("PGC",ao$consortium)]
associations(variants="11:61043499-61053600",id[1])
associations(variants="rs174546",id[1])

https://github.com/mrcieu/ieugwasr