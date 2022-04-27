library(multcomp)
# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/
load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/UKBB_cancer_aspirin_smoking_plusgenetic_data.RData")

Dat<-combine_cancers_function(cancers=c("overall_colorectal_cancer","overall_lung_cancer","overall_nm_skin_cancer","overall_melanoma"),cancer_name="cancer_sig")
mylogit <- glm( cancer_sig ~ d5d +age_at_recruitment +sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
Model<-summary(mylogit)
table(Dat$cancer_sig)

save(Model,file="/projects/MRC-IEU/users/ph14916/fatty_acids_summary/mr_crc_lc_msc_d5d.Rdata")




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
    
    Pos_case<-unique(unlist(lapply(cancers, FUN=function(x) which(Dat[,x]==1))))
    # Dat$overall_pan_inclc44_cancer
    Pos_controls<-which(Dat$overall_pan_inclc44_cancer==0)
    Dat[,cancer_name]<-NA
    Dat[Pos_case,cancer_name]<-1
    Dat[Pos_controls,cancer_name]<-0
    return(Dat) 
}
