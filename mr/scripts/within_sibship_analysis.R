# devtools::install_github("NightingaleHealth/ggforestplot")
# install.packages("ggplot2")
library(ggforestplot)
library(ggplot2)
# require(data.table)
library(data.table)
# install.packages("sandwich")
# require(sandwich)
library(sandwich)
# install.packages("lmtest")
# require(lmtest)
library(lmtest)

setwd("/mnt/storage/private/mrcieu/users/lh14833/shared")
load("/mnt/storage/home/ph14916/cancer_casecontrol_plus_genetic_data.RData")

# generate d5d exposure based on rs174546. d5d is genetically predicted value of d5d activity
Dat<-d5d_function()
# create new cancer variable that combines cases across top cancer findings in the main MR analysis (colorectal cancer, lung cancer, malignant skin cancer). Malignant skin cancer corresponds to non-melanoma and melanoma skin cancer. The scripts that create the case-control cancer data excluded non-malignant cancers be default 
Dat<-combine_cancers_function(cancers=c("overall_colorectal_cancer","overall_lung_cancer","overall_nm_skin_cancer","overall_melanoma"),cancer_name="cancer_sig")
# Dat<-combine_cancers_function(cancers=c("overall_lung_cancer","overall_colorectal_cancer",cancer_name="cancer_crc_lc")
# smoking_cancers<-c("overall_cervical_cancer", "overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer")
# Dat<-combine_cancers_function(cancers=smoking_cancers,cancer_name="cancer_smk")

# j<-which(names(Dat)=="overall_pan_inclc44_cancer")
# Dat<-format_data(analysis=Dat,j=j)    
# Model<-glm(formula = overall_pan_inclc44_cancer ~ d5d, data =Dat,family="binomial")
# summary(Model)
# table(Dat$overall_pan_exclc44_cancer)
# table(Dat$overall_pan_inclc44_cancer)
# table(Dat$cancer_sig)
# table(Dat$cancer_smk)
# dim(dat.m)
# dim(Fam)
# length(unique(Fam$V1))
# length(unique(dat.m$V1))
# length(which(table(dat.m$V1)==6))
# dat.m[1:10,1:10]
Fam<-read.table("data_filtered.fam",sep=" ",head=FALSE,stringsAsFactors=FALSE)
dat.m<-merge(Fam,Dat,by.x="V2",by.y="geneticID")
# glm (family=binomial)

phenotypes <- c("overall_pan_exclc44_cancer", "overall_pan_inclc44_cancer", "overall_lung_cancer", "overall_colorectal_cancer","overall_melanoma_plus_other_malignant_skin_cancer","cancer_sig")

Out<-within_sibs_function2(phenotypes=phenotypes)

save(Out,file="~/sibling_results.Rdata")
# plot_dat(Dat=Out)

within_sibs_function2<-function(phenotypes=NULL){
    Out<-NULL
    for(i in 1:length(phenotypes)){
        Out[[i]]<-within_sibs_function(analysis=dat.m,phenotype=phenotypes[i])
    }
    Out2<-do.call(rbind,Out)
    return(Out2)
}

within_sibs_function<-function(analysis=NULL,phenotype=NULL){
    #Change exposure to the name of your exposure.. can be a single SNP, a score, a phenotype etc.

    #Output data.frame
    output <- data.frame(PHEN=phenotype, BETA_0=NA, BETA_BF=NA, BETA_WF=NA,
                         SE_BETA_0=NA,  SE_BETA_BF=NA, SE_BETA_WF=NA,
                         P_BETA_0=NA, P_BETA_BF=NA, P_BETA_WF=NA,
                         VCV_0=NA, VCV_0_BF=NA, VCV_0_WF=NA, VCV_BF=NA, VCV_BF_WF=NA, VCV_WF=NA)

    j<-which(names(analysis)==phenotype)
    analysis<-format_data(analysis=analysis,j=j)    
    # analysis<-analysis[!is.na(analysis$overall_pan_exclc44_cancer),]
    # analysis<-analysis[!is.na(analysis$overall_colorectal_cancer),]
    # analysis<-analysis[!is.na(analysis$overall_lung_cancer),]
    # analysis$overall_pan_exclc44_cancer<-analysis$overall_pan_exclc44_cancer-1
    # analysis$overall_colorectal_cancer<-analysis$overall_colorectal_cancer-1
    # analysis$overall_lung_cancer<-analysis$overall_lung_cancer-1
    # analysis$overall_colorectal_cancer<-analysis$overall_lung_cancer-1
    # table(analysis$overall_lung_cancer)
    # Mke a matrix with: [FID PHENOTYPE] [individ - family mean ] [family mean]
    merge <- data.table(FID=analysis$FID, OUTCOME=analysis[,j], EXPOSURE=analysis$d5d, AGE=analysis$f.21022.0.0, SEX=analysis$sex2, FAM_MEAN=ave(analysis$rs174546, analysis$FID, FUN=mean))

     # merge <- data.table(FID=analysis$FID, OUTCOME=as.numeric(unlist(analysis[,j,with=FALSE])), EXPOSURE=analysis$rs174546, AGE=analysis$f.21022.0.0, SEX=analysis$sex2, FAM_MEAN=ave(analysis$rs174546, analysis$FID, FUN=mean))

    #Centre exposure around family-mean
    merge2 <- na.omit(merge[,CENTREDEXPOSURE:=EXPOSURE-FAM_MEAN])


    # Run regression for rs174546 on outcomes
    # fit <- lm(formula = OUTCOME ~ FAM_MEAN + CENTREDEXPOSURE, data = merge2)
    ncase<-length(which(merge2[,OUTCOME]==1))
    ncont<-length(which(merge2[,OUTCOME]==0))
    fit <- glm(formula = OUTCOME ~ FAM_MEAN + CENTREDEXPOSURE, data = merge2,family="binomial")

    #Extract coefficients

    output$BETA_0 <- fit$coefficients[1]
    output$BETA_BF <- fit$coefficients[2]
    output$BETA_WF <- fit$coefficients[3]
     
    # save the variance covariance matrix
    vcv_matrix = vcovCL(fit, cluster=merge2$FID)
    if(  is.na(output$BETA_0) | is.na(output$BETA_BF) | is.na(output$BETA_WF) ) {
        output$VCV_0 <- NA
        output$VCV_0_BF <- NA
        output$VCV_0_WF <- NA
        output$VCV_BF <- NA
        output$VCV_BF_WF <- NA
        output$VCV_WF <- NA
    } else {
        output$VCV_0 <- vcv_matrix[1,1]
        output$VCV_0_BF <- vcv_matrix[1,2]
        output$VCV_0_WF <- vcv_matrix[1,3]
        output$VCV_BF <- vcv_matrix[2,2]
        output$VCV_BF_WF <- vcv_matrix[2,3]
        output$VCV_WF <- vcv_matrix[3,3]
    }

    # save the clustered SE's and corresponding p-values
    test_matrix <- coeftest(fit, vcov.=vcv_matrix)
    if(  is.na(output$BETA_0) | is.na(output$BETA_BF) | is.na(output$BETA_WF) ) {
        output$SE_BETA_0 <- NA
        output$SE_BETA_BF <- NA
        output$SE_BETA_WF <- NA
        output$P_BETA_0 <- NA
        output$P_BETA_BF <- NA
        output$P_BETA_WF <- NA
    } else {
        output$SE_BETA_0 <- test_matrix[1,2] 
        output$SE_BETA_BF <- test_matrix[2,2] 
        output$SE_BETA_WF <- test_matrix[3,2] 
        output$P_BETA_0 <- test_matrix[1,4] 
        output$P_BETA_BF <- test_matrix[2,4] 
        output$P_BETA_WF <- test_matrix[3,4] 
    }
    output$ncase<-ncase
    output$ncont<-ncont
    return(output)
}



d5d_function<-function(){
    # change coding so that effect allele is allele c (the major allele and D5D raising allele )
    Dat$d5d<-NA
    Dat$d5d[Dat$rs174546==2]<-0
    Dat$d5d[Dat$rs174546==1]<-1
    Dat$d5d[Dat$rs174546==0]<-2
    Dat$d5d<-Dat$d5d*0.86
    return(Dat)
    # # effect of C allele on D5D activity is 0.86 SD units
    # Dat$rs174546_3<-Dat$rs174546_2 * 0.86
}

format_data<-function(analysis=NULL,j=NULL){
    names(analysis)[names(analysis)=="V1"]<-"FID"
    analysis$sex2[analysis$sex2=="M"]<-1
    analysis$sex2[analysis$sex2=="F"]<-0
    analysis$sex2<-as.numeric(analysis$sex2)    
    analysis[,j]<-analysis[,j]-1
    return(analysis)
}



sig_can_function<-function(){
    cancers<-c("overall_headneck_cancer","overall_colorectal_cancer","overall_lung_cancer","overall_nm_skin_cancer","overall_larynx_cancer","overall_liver_cell_cancer","overall_acute_lymph_leuk_cancer")
    Pos_case<-unique(unlist(lapply(cancers, FUN=function(x) which(Dat[,x]==2))))
    Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
    Dat$sig_cancers<-NA
    Dat$sig_cancers[Pos_case]<-2
    Dat$sig_cancers[Pos_controls]<-1

    return(Dat) 
}


lc_crc_function<-function(){
    lc_crc_cancers<-c("overall_colorectal_cancer","overall_lung_cancer")
    Pos_case<-unique(unlist(lapply(lc_crc_cancers, FUN=function(x) which(Dat[,x]==2))))
    Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
    Dat$lc_crc_cancers<-NA
    Dat$lc_crc_cancers[Pos_case]<-2
    Dat$lc_crc_cancers[Pos_controls]<-1 
    return(Dat) 
}


combine_cancers_function<-function(cancers=NULL,cancer_name=NULL){
    # lung cancer, colorectal cancer, malignant skin cancer and non-melanoma malignant skin cancer
    # names(Dat)[grep("overall",names(Dat))]
      # "overall_breast_cancer","overall_prostate_cancer"   
      #           "overall_endometrial_cancer" 
      #           "overall_pan_inclc44_cancer"
      #            "overall_pan_exclc44_cancer" 
      #            "overall_melanoma"
      #            "overall_ovarian_cancer"
      #            "overall_nm_skin_cancer"      
      #            "overall_myel_leuk_cancer","overall_lymph_leuk_cancer"
      #            "overall_mult_myel" 
      #            "overall_leuk_cancer"
      #            "overall_haem_cancer"
      #            "overall_brain_cancer"
    
    Pos_case<-unique(unlist(lapply(cancers, FUN=function(x) which(Dat[,x]==2))))
    # Dat$overall_pan_inclc44_cancer
    Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
    Dat[,cancer_name]<-NA
    Dat[Pos_case,cancer_name]<-2
    Dat[Pos_controls,cancer_name]<-1  
    return(Dat) 
}

lc_crc_nm_cancers_function<-function(){
    # names(Dat)[grep("overall",names(Dat))]
      # "overall_breast_cancer","overall_prostate_cancer"   
      #           "overall_endometrial_cancer" 
      #           "overall_pan_inclc44_cancer"
      #            "overall_pan_exclc44_cancer" 
      #            "overall_melanoma"
      #            "overall_ovarian_cancer"
      #            "overall_nm_skin_cancer"      
      #            "overall_myel_leuk_cancer","overall_lymph_leuk_cancer"
      #            "overall_mult_myel" 
      #            "overall_leuk_cancer"
      #            "overall_haem_cancer"
      #            "overall_brain_cancer"
    lc_crc_nm_cancers<-c("overall_colorectal_cancer","overall_lung_cancer","overall_nm_skin_cancer")
    Pos_case<-unique(unlist(lapply(lc_crc_nm_cancers, FUN=function(x) which(Dat[,x]==2))))
    # Dat$overall_pan_inclc44_cancer
    Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
    Dat$lc_crc_nm_cancers<-NA
    Dat$lc_crc_nm_cancers[Pos_case]<-2
    Dat$lc_crc_nm_cancers[Pos_controls]<-1  
    return(Dat) 
}

# smoking_cancers_function<-function(){
#     # names(Dat)[grep("overall",names(Dat))]
#       # "overall_breast_cancer","overall_prostate_cancer"   
#       #           "overall_endometrial_cancer" 
#       #           "overall_pan_inclc44_cancer"
#       #            "overall_pan_exclc44_cancer" 
#       #            "overall_melanoma"
#       #            "overall_ovarian_cancer"
#       #            "overall_nm_skin_cancer"      
#       #            "overall_myel_leuk_cancer","overall_lymph_leuk_cancer"
#       #            "overall_mult_myel" 
#       #            "overall_leuk_cancer"
#       #            "overall_haem_cancer"
#       #            "overall_brain_cancer"
#     smoking_cancers<-c("overall_cervical_cancer", "overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer")

#     smoking_cancers[!smoking_cancers %in% names(Dat)]
#     Pos_case<-unique(unlist(lapply(smoking_cancers, FUN=function(x) which(Dat[,x]==2))))
#     Dat$overall_pan_inclc44_cancer
#     Pos_controls<-which(Dat$overall_pan_inclc44_cancer==1)
#     Dat$smoking_cancers<-NA
#     Dat$smoking_cancers[Pos_case]<-2
#     Dat$smoking_cancers[Pos_controls]<-1
#     Dat$non_smoking_cancers_inclc44<-Dat$overall_pan_inclc44_cancer
#     Dat$non_smoking_cancers_inclc44[Pos_case]<-NA
#     Dat$non_smoking_cancers_exclc44<-Dat$overall_pan_exclc44_cancer
#     Dat$non_smoking_cancers_exclc44[Pos_case]<-NA   
#     return(Dat) 
# }


plot_dat<-function(Dat=NULL,text.names=10,text.title=10,Shape=NULL,colour=NULL){
    Dat$PHEN<-c("Pan-cancer\nexclude non-melanoma SC","Pan-cancer","Lung cancer","Colorectal cancer","Malignant skin cancers","Selected cancers")
    Dat<-Dat[order(Dat$ncase,decreasing=TRUE),]
    Dat$PHEN<-paste0(Dat$PHEN,"\nN. cases=",Dat$ncase)
    p<-forestplot(df = Dat,
            logodds = TRUE,
            name=PHEN,
                  estimate=BETA_WF,
                  se=SE_BETA_WF,
                  shape=Shape,
                  colour = colour,
                   xlab = "")+
            # labs(title=Title.plot,size=1)+
            theme(plot.title = element_text(size = text.title))+
            theme(text = element_text(size=text.names))
    return(p)
}
