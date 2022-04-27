library(multcomp)
# cd /projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/
load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/UKBB_cancer_aspirin_smoking_plusgenetic_data.RData")

dim(Dat)
table(Dat$overall_colorectal_cancer)
names(Dat)
Dat
# mylogit <- glm(smoking_cancers ~ GS_d5d*GS_eqtlgen_COX2 +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
# summary(mylogit)
# names(Dat)
# mylogit <- glm(overall_lung_cancer ~ GS_d5d*GS_COX1_Lung +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
# summary(mylogit)
# # smoking analysis
# # overall cancer 

# Dat3<-Dat2[Dat2$smk == 1,]
# Dat4<-Dat2[Dat2$smk == 0,]

# names(Dat)


#  "nsaid_baseline_f_6154"                      
#  "aspirin_baseline_f_6154"                    
#  "nsaid_baseline_f_20003"                     
#  "aspirin_baseline_f_20003"  

# table(Dat$nsaid_baseline_f_6154)
# table(Dat$aspirin_baseline_f_6154)
# table(Dat$nsaid_baseline_f_20003)
# table(Dat$aspirin_baseline_f_20003,Dat$aspirin_baseline_f_6154)

i<-which(Cancers == "overall_prostate_cancer")
Cancers<-names(Dat)[grep("overall",names(Dat))]
Cancers<-Cancers[Cancers !=  "overall_ovarian_cancer"]
Cancers<-Cancers[Cancers !=  "overall_prostate_cancer"]
Cancers<-Cancers[Cancers !=  "overall_cervical_cancer"]
Cancers<-Cancers[Cancers !=  "overall_endometrial_cancer" ]
# Cancers<-c("overall_endometrial_cancer",Cancers)



Cancers<-c("overall_colorectal_cancer" ,"overall_lung_cancer"  , "overall_nm_skin_cancer" , "overall_melanoma", "overall_melanoma_plus_other_malignant_skin_cancer","overall_pan_inclc44_cancer","overall_pan_exclc44_cancer","smoking_cancers", "non_smoking_cancers_inclc44","non_smoking_cancers_exclc44","lc_crc_nm_cancers","lc_crc_cancers","sig_cancers" )

smoking_cancers<-c("overall_cervical_cancer", "overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer")
smoking_cancers1<-c("smoking_cancers","overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer")



Cancers<-smoking_cancers1
# table(Dat$overall_ovarian_cancer)
Cancers<-smoking_cancers1

Cancers<-c("overall_pan_inclc44_cancer","overall_pan_exclc44_cancer")

Interactions<-c("smoking2","aspirin2","GS_csi","GS_cpd","GS_gtex_cox1","GS_gtex_cox2","GS_eqtlgen_COX1","GS_eqtlgen_COX2")

Interactions<-c("GS_cpd","CS_csi")


# Res[Res$se0>1,]

# Res[Res$cancer == "smoking_cancers",]
# 270330+4159 

table(Res$cancer)

Res<-model_interactions_function(Dat_test=Dat,Cancers=Cancers)
Smk1<-model_interactions_function(Dat_test=Datsmk1,Cancers=Cancers) #in ever smokers
Smk0<-model_interactions_function(Dat_test=Datsmk0,Cancers=Cancers) #in never smokers 
Asp1<-model_interactions_function(Dat_test=Datasp1,Cancers=Cancers) #in regular aspirin users
Asp0<-model_interactions_function(Dat_test=Datasp0,Cancers=Cancers) #in non-regular aspirin users 

GS_cpd_smk<-cpd_smoking2_interaction(Dat_test=Dat,Cancers=Cancers)
GS_cpd_smk[,c("cancer","b0","se0", "b1","se1","b_int","se_int")]
head(GS_cpd_smk)

Res_main_effects_sm<-model_main_effects_function(Dat_test=Dat,Cancers=Cancers)
head(Res_main_effects)
Res_main_effects[Res_main_effects$interaction=="GS_csi",c("cancer","interaction","b_int","se_int")]
Res_main_effects_sm[Res_main_effects_sm$interaction=="GS_csi",c("cancer","interaction","b_int","se_int")]
Smk1_main_effects<-model_main_effects_function(Dat_test=Datsmk1,Cancers=Cancers)
Smk1_main_effects[,c("cancer","interaction","b_int","se_int")]
Smk0_main_effects<-model_main_effects_function(Dat_test=Datsmk0,Cancers=Cancers)
Asp1_main_effects<-model_main_effects_function(Dat_test=Datasp1,Cancers=Cancers)
Asp0_main_effects<-model_main_effects_function(Dat_test=Datasp0,Cancers=Cancers)

Res_main_effects_smkc<-model_main_effects_function(Dat_test=Dat,Cancers=smoking_cancers1)
Smk1_main_effects_smkc<-model_main_effects_function(Dat_test=Datsmk1,Cancers=smoking_cancers1)
Smk0_main_effects_smkc<-model_main_effects_function(Dat_test=Datsmk0,Cancers=smoking_cancers1)
# Asp1_main_effects_smkc<-model_main_effects_function(Dat_test=Datasp1,Cancers=smoking_cancers)
# Asp0_main_effects_smkc<-model_main_effects_function(Dat_test=Datasp0,Cancers=smoking_cancers)


Res_d5d<-model_d5d_function(Dat_test=Dat)



plot(Res_main_effects$b_int[Res_main_effects$interaction == "GS_csi"],Res_main_effects$b0[Res_main_effects$interaction == "GS_csi"])

Res_main_effects_smkc[Res_main_effects_smkc$interaction == "GS_csi",]

plot(Res_main_effects$b_int[Res_main_effects$interaction == "GS_csi"],Res_main_effects$b0[Res_main_effects$interaction == "GS_csi"])



Temp<-Res_main_effects[Res_main_effects$interaction=="GS_csi",c("cancer","b0","b_int")]
Temp[order(Temp$b0),]
# Is effect of FADS1/2 stronger for smoking-related compared to non-smoking related cancers? This helps address limitation of potential sample overlap amongst studies in summary data analysis

# "overall_pan_inclc44_cancer","overall_pan_exclc44_cancer","smoking_cancers"
mylogit <- glm( overall_pan_inclc44_cancer ~ d5d +aspirin2 +age_at_recruitment +sex+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
Model<-summary(mylogit)

mylogit <- glm(smoking_cancers ~ d5d*smoking2 +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
Model<-summary(mylogit)

mylogit <- glm(smoking_cancers ~ d5d*GS_csi +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
Model<-summary(mylogit)

mylogit <- glm(smoking_cancers ~ d5d+ GS_cpd +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Datsmk1, family = "binomial")
Model<-summary(mylogit)

Res[Res$cancer =="smoking_cancers", ]
Smk1[Smk1$cancer =="smoking_cancers", ]

Res[Res$cancer =="overall_pan_inclc44_cancer", ]
Smk1[Smk1$cancer =="overall_pan_inclc44_cancer", ]

Res[Res$cancer =="lc_crc_cancers", ]
Smk1[Smk1$cancer =="lc_crc_cancers", ]

Res[Res$cancer =="lc_crc_nm_cancers", ]
Smk1[Smk1$cancer =="lc_crc_nm_cancers", ]



head(Res)
head(Smk1)
head(Smk0)
head(Asp1)
head(Asp0)

unique(Res$interaction)
Cancers2<-c("overall_pan_inclc44_cancer","overall_pan_exclc44_cancer","smoking_cancers", "non_smoking_cancers_inclc44","non_smoking_cancers_exclc44")

Cancers1<-c("overall_colorectal_cancer" ,"overall_lung_cancer"  , "overall_nm_skin_cancer" , "overall_melanoma", "overall_melanoma_plus_other_malignant_skin_cancer","lc_crc_nm_cancers","lc_crc_cancers")

Res[Res$cancer %in% c("smoking_cancers", "non_smoking_cancers_inclc44","non_smoking_cancers_exclc44"),]

Res[Res$interaction == "GS_csi" & Res$cancer %in% Cancers2,]

Smk1[Smk1$interaction == "GS_cpd" & Smk1$cancer %in% Cancers2,c("cancer","b0","se0","p0","b1","se1","p1","b_int","se_int","p_int")]

Smk0[Smk0$interaction == "GS_cpd" & Smk0$cancer %in% Cancers2,c("cancer","b0","se0","p0","b1","se1","p1","b_int","se_int","p_int")]



Res[Res$interaction == "aspirin2" & Res$cancer %in% Cancers2,]

Res[Res$interaction == "GS_gtex_cox1" & Res$cancer %in% Cancers2,]
Res[Res$interaction == "GS_gtex_cox2" & Res$cancer %in% Cancers2,]
Res[Res$interaction == "GS_eqtlgen_COX1" & Res$cancer %in% Cancers2,]
Res[Res$interaction == "GS_eqtlgen_COX2" & Res$cancer %in% Cancers2,]

Res[Res$p_int<0.05,]
Smk1[Smk1$p_int<0.05,]
Smk0[Smk0$p_int<0.05,]
Asp0[Asp0$p_int<0.05,]
Asp1[Asp1$p_int<0.05,]

Smk1[Smk1$interaction == "GS_cpd" & Smk1$cancer %in% Cancers2,c("cancer","b0","se0","p0","b1","se1","p1","b_int","se_int","p_int")]

Smk0[Smk0$interaction == "GS_cpd" & Smk0$cancer %in% Cancers2,c("cancer","b0","se0","p0","b1","se1","p1","b_int","se_int","p_int")]

"GS_COX1_Colon_Transverse"  
"GS_COX2_Colon_Transverse" 
 "GS_COX1_Lung"                                      
"GS_COX2_Lung"


# aspirin_baseline_f_6154
table(Dat$aspirin_baseline_f_20003)
table(Dat$aspirin_baseline_f_6154)


j<-which(Interactions == "GS_eqtlgen_COX1")
i<-which(Cancers=="overall_colorectal_cancer")
Dat_test<-Dat

names(Dat)[grep("overall",names(Dat))]
# Combine the individual cancers with most evidence for associations with D5D [P<0.05]
Dat1<-sig_can_function2(bd=Dat)
# Dat1<-sig_can_function1(bd=Dat)
table(Dat1$sig_cancers2)

table(Dat$smoking_cancers)

mylogit <- glm(smoking_cancers ~ d5d + aspirin2 +age_at_recruitment +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
summary(mylogit)

mylogit <- glm(overall_mult_myel ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")

summary(mylogit)

mylogit <- glm(overall_leuk_cancer ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
summary(mylogit)

mylogit <- glm(overall_larynx_cancer ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
summary(mylogit)


Dat_test<-Datcancer_testcancer_test


Res_main_effects[Res_main_effects$cancer == "smoking_cancers",]
Res_main_effects[Res_main_effects$cancer == "non_smoking_cancers_inclc44",]


Res_main_smk_cancers<-model_main_effects_function(Dat_test=Dat,Cancers=smoking_cancers<-c("overall_cervical_cancer", "overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer"))

Smk1_main_smk_cancers<-model_main_effects_function(Dat_test=Datsmk1,Cancers=smoking_cancers<-c("overall_cervical_cancer", "overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer"))

Smk0_main_smk_cancers<-model_main_effects_function(Dat_test=Datsmk0,Cancers=smoking_cancers<-c("overall_cervical_cancer", "overall_kidney_cancer","overall_pancreatic_cancer", "overall_stomach_cancer","overall_aml_cancer", "overall_colorectal_cancer","overall_lung_cancer","overall_oesoph_cancer","overall_liver_cell_cancer","overall_larynx_cancer","overall_bladder_cancer","overall_headneck_cancer"))


model_interactions_function<-function(Dat_test=NULL,Cancers=Cancers)
{
	Results1<-NULL
	Results2<-NULL
	for(i in 1:length(Cancers))
	{
		print(Cancers[i])
		for(j in 1:length(Interactions))
		{
			print(Interactions[j])
			# Dat_test<-Dat

			Dat_test$cancer_test<- Dat_test[,Cancers[i]]
			Dat_test$interaction_test<- Dat_test[,Interactions[j]]
			Results<-NA
			if(length(table(Dat_test$interaction_test))>1)
			{
				mylogit <- glm(cancer_test ~ d5d*interaction_test + age_at_recruitment +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat_test, family = "binomial")
				Model<-summary(mylogit)
				# mylogit1 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
				# mylogit0 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat0, family = "binomial")	
				# # for(j in c("d5d","d5d:smoking2"){
				b0<-Model$coefficients["d5d","Estimate"] #effect in never smokers
				se0<-Model$coefficients["d5d","Std. Error"]
				z0<-Model$coefficients["d5d","z value"]
				p0<-Model$coefficients["d5d","Pr(>|z|)"]
				
				b_int<-Model$coefficients["d5d:interaction_test","Estimate"] #effect change in ever smokes 
				se_int<-Model$coefficients["d5d:interaction_test","Std. Error"]
				z_int<-Model$coefficients["d5d:interaction_test","z value"]
				p_int<-Model$coefficients["d5d:interaction_test","Pr(>|z|)"]
				
				Model_glht<-summary(glht(mylogit, linfct = c("d5d + d5d:interaction_test = 0")))
				# Model_glht<-summary(glht(mylogit, linfct = c("d5d + d5d:smoking2 = 0")))
				b1<-as.numeric(Model_glht$test$coefficients)
				se1<-as.numeric(Model_glht$test$sigma)
				z1<-as.numeric(Model_glht$test$tstat)
				p1<-as.numeric(Model_glht$test$pvalues[1])
				Cases<-length(which(Dat_test$cancer_test==1))
				Controls<-length(which(Dat_test$cancer_test==0))
				Results1[[j]]<-c(Cancers[i],Interactions[j],Cases,Controls,b0,se0,z0,p0,b1,se1,z1,p1,b_int,se_int,z_int,p_int)
				Results<-data.frame(do.call(rbind,Results1))
				rm(list=c("b0","b1","b_int"))								
			}
			Results2[[i]]<-Results			
			
		}
	}
	Res<-data.frame(do.call(rbind,Results2))
	names(Res)<-c("cancer","interaction","cases","controls","b0","se0","z0","p0","b1","se1","z1","p1","b_int","se_int","z_int","p_int")
	
	for(i in c("b0","se0","z0","p0","b1","se1","z1","p1","b_int","se_int","z_int","p_int")){
		print(i)
		Res[,i]<-round(as.numeric(Res[,i]),3)
	}

	return(Res)
}


model_main_effects_function<-function(Dat_test=NULL,Cancers=Cancers)
{
	Results1<-NULL
	Results2<-NULL
	for(i in 1:length(Cancers))
	{
		print(Cancers[i])
		for(j in 1:length(Interactions))
		{
			print(Interactions[j])
			# Dat_test<-Dat

			Dat_test$cancer_test<- Dat_test[,Cancers[i]]
			Dat_test$interaction_test<- Dat_test[,Interactions[j]]
			Results<-NA
			if(length(table(Dat_test$interaction_test))>1)
			{
				# table(Dat$overall_cervical_cancer)
				# table(Dat$smoking2,Dat$overall_cervical_cancer)
				# mylogit <- glm(overall_colorectal_cancer ~ d5d+smoking2 + age_at_recruitment +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat_test, family = "binomial")
			# table(Dat_test$overall_prostate_cancer,Dat_test$GS_csi)		

				mylogit <- glm(cancer_test ~ d5d+interaction_test + age_at_recruitment +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat_test, family = "binomial")
				Model<-summary(mylogit)
				# mylogit1 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
				# mylogit0 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat0, family = "binomial")	
				# # for(j in c("d5d","d5d:smoking2"){
				b0<-Model$coefficients["d5d","Estimate"] #effect in never smokers
				se0<-Model$coefficients["d5d","Std. Error"]
				z0<-Model$coefficients["d5d","z value"]
				p0<-Model$coefficients["d5d","Pr(>|z|)"]
				
				b_int<-Model$coefficients["interaction_test","Estimate"] #effect change in ever smokes 
				se_int<-Model$coefficients["interaction_test","Std. Error"]
				z_int<-Model$coefficients["interaction_test","z value"]
				p_int<-Model$coefficients["interaction_test","Pr(>|z|)"]
				
				# Model_glht<-summary(glht(mylogit, linfct = c("d5d + d5d:interaction_test = 0")))
				# # Model_glht<-summary(glht(mylogit, linfct = c("d5d + d5d:smoking2 = 0")))
				# b1<-as.numeric(Model_glht$test$coefficients)
				# se1<-as.numeric(Model_glht$test$sigma)
				# z1<-as.numeric(Model_glht$test$tstat)
				# p1<-as.numeric(Model_glht$test$pvalues[1])
				Cases<-length(which(Dat_test$cancer_test==1))
				Controls<-length(which(Dat_test$cancer_test==0))
				Results1[[j]]<-c(Cancers[i],Interactions[j],Cases,Controls,b0,se0,z0,p0,b_int,se_int,z_int,p_int)
				Results<-data.frame(do.call(rbind,Results1))
				rm(list=c("b0","b_int"))								
			}
			Results2[[i]]<-Results			
			
		}
	}
	Res<-data.frame(do.call(rbind,Results2))
	names(Res)<-c("cancer","interaction","cases","controls","b0","se0","z0","p0","b_int","se_int","z_int","p_int")
	
	for(i in c("b0","se0","z0","p0","b_int","se_int","z_int","p_int")){
		print(i)
		Res[,i]<-round(as.numeric(Res[,i]),3)
	}

	return(Res)
}


model_d5d_function<-function(Dat_test=NULL)
{
	Results<-NULL
	for(i in 1:length(Cancers))
	{
		print(Cancers[i])
		Dat_test$cancer_test<- Dat_test[,Cancers[i]]
		mylogit <- glm(cancer_test ~ d5d + age_at_recruitment +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat_test, family = "binomial")
		Model<-summary(mylogit)
		# mylogit1 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
		# mylogit0 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat0, family = "binomial")	
		# # for(j in c("d5d","d5d:smoking2"){
		b0<-Model$coefficients["d5d","Estimate"] #effect in never smokers
		se0<-Model$coefficients["d5d","Std. Error"]
		z0<-Model$coefficients["d5d","z value"]
		p0<-Model$coefficients["d5d","Pr(>|z|)"]
		
		# Model_glht<-summary(glht(mylogit, linfct = c("d5d + d5d:interaction_test = 0")))
		# # Model_glht<-summary(glht(mylogit, linfct = c("d5d + d5d:smoking2 = 0")))
		# b1<-as.numeric(Model_glht$test$coefficients)
		# se1<-as.numeric(Model_glht$test$sigma)
		# z1<-as.numeric(Model_glht$test$tstat)
		# p1<-as.numeric(Model_glht$test$pvalues[1])
		Cases<-length(which(Dat_test$cancer_test==1))
		Controls<-length(which(Dat_test$cancer_test==0))
		Results[[i]]<-c(Cancers[i],Cases,Controls,b0,se0,z0,p0)
		rm(list=c("b0"))								
	}		
	Res<-data.frame(do.call(rbind,Results))
	names(Res)<-c("cancer","cases","controls","b0","se0","z0","p0")

	for(i in c("b0","se0","z0","p0")){
		print(i)
		Res[,i]<-round(as.numeric(Res[,i]),5)
	}

	return(Res)
}
# kidney
# head and neck
# liver 
# lung 
# overall_colorectal_cancer

# names(Dat)[grep("overall",names(Dat))]


# sig_can_function<-function(){
# 	lc_crc_cancers<-c("overall_headneck_cancer","overall_kidney_cancer","overall_colorectal_cancer","overall_lung_cancer","overall_liver_cell_cancer")
# 	Pos_case<-unique(unlist(lapply(lc_crc_cancers, FUN=function(x) which(bd[,x]==2))))
# 	Pos_controls<-which(bd$overall_pan_inclc44_cancer==1)
# 	bd$lc_crc_cancers<-NA
# 	bd$lc_crc_cancers[Pos_case]<-2
# 	bd$lc_crc_cancers[Pos_controls]<-1	
# 	return(bd)	
# }



cpd_smoking2_interaction<-function(Dat_test=NULL,Cancers=Cancers)
{
	Results1<-NULL
	Results2<-NULL
	for(i in 1:length(Cancers))
	{
		print(Cancers[i])
		for(j in 1:length(Interactions))
		{
			print(Interactions[j])
			# Dat_test<-Dat

			Dat_test$cancer_test<- Dat_test[,Cancers[i]]
			Dat_test$interaction_test<- Dat_test[,Interactions[j]]
			Results<-NA
			if(length(table(Dat_test$interaction_test))>1)
			{
				mylogit <- glm(cancer_test ~ GS_cpd*smoking2 + age_at_recruitment +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat_test, family = "binomial")
				Model<-summary(mylogit)


				# mylogit1 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
				# mylogit0 <- glm(cancer_test ~ d5d +sex + pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat0, family = "binomial")	
				# # for(j in c("d5d","d5d:smoking2"){
				b0<-Model$coefficients["GS_cpd","Estimate"] #effect in never smokers
				se0<-Model$coefficients["GS_cpd","Std. Error"]
				z0<-Model$coefficients["GS_cpd","z value"]
				p0<-Model$coefficients["GS_cpd","Pr(>|z|)"]
				
				b_int<-Model$coefficients["GS_cpd:smoking2","Estimate"] #effect change in ever smokes 
				se_int<-Model$coefficients["GS_cpd:smoking2","Std. Error"]
				z_int<-Model$coefficients["GS_cpd:smoking2","z value"]
				p_int<-Model$coefficients["GS_cpd:smoking2","Pr(>|z|)"]
				
				Model_glht<-summary(glht(mylogit, linfct = c("GS_cpd + GS_cpd:smoking2 = 0")))
				# Model_glht<-summary(glht(mylogit, linfct = c("d5d + d5d:smoking2 = 0")))
				b1<-as.numeric(Model_glht$test$coefficients)
				se1<-as.numeric(Model_glht$test$sigma)
				z1<-as.numeric(Model_glht$test$tstat)
				p1<-as.numeric(Model_glht$test$pvalues[1])
				Cases<-length(which(Dat_test$cancer_test==1))
				Controls<-length(which(Dat_test$cancer_test==0))
				Results1[[j]]<-c(Cancers[i],Interactions[j],Cases,Controls,b0,se0,z0,p0,b1,se1,z1,p1,b_int,se_int,z_int,p_int)
				Results<-data.frame(do.call(rbind,Results1))
				rm(list=c("b0","b1","b_int"))								
			}
			Results2[[i]]<-Results						
		}
	}
	Res<-data.frame(do.call(rbind,Results2))
	names(Res)<-c("cancer","interaction","cases","controls","b0","se0","z0","p0","b1","se1","z1","p1","b_int","se_int","z_int","p_int")
	
	for(i in c("b0","se0","z0","p0","b1","se1","z1","p1","b_int","se_int","z_int","p_int")){
		print(i)
		Res[,i]<-round(as.numeric(Res[,i]),3)
	}

	return(Res)
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
    
    Pos_case<-unique(unlist(lapply(cancers, FUN=function(x) which(Dat[,x]==1))))
    # Dat$overall_pan_inclc44_cancer
    Pos_controls<-which(Dat$overall_pan_inclc44_cancer==0)
    Dat[,cancer_name]<-NA
    Dat[Pos_case,cancer_name]<-1
    Dat[Pos_controls,cancer_name]<-0
    return(Dat) 
}
