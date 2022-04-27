load("/projects/MRC-IEU/users/ph14916/fatty_acids_summary/cox_smoking/UKBB_cancer_aspirin_smoking_plusgenetic_data.RData")

class(Dat$overall_pan_inclc44_cancer)
summary(Dat$d5d)
class(Dat$GS_csi)
class(Dat$age_at_recruitment)
Dat$sex_num<-Dat$sex
Dat$sex_num[Dat$sex_num=="Female"]<-"0"
Dat$sex_num[Dat$sex_num=="Male"]<-"1"
Dat$sex_num<-as.numeric(Dat$sex_num)

mylogit <- glm( overall_pan_inclc44_cancer ~ d5d*GS_csi +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
Model<-summary(mylogit)

table(Dat$smoking2) #seems like a lot of ever smokers - surveys indicate 59% UK had never smoked in 2017
184002/(184002+151798)*100 #54% UKB participants have never smoked. since is a middle aged sample maybe we do expect fewer never smokers compared to the general population. on other hand, healthier than general pop, and would expct that to push in other direction. 

# Data showcase UKB 2020: 
# total<-501590
# never<-318238
# ever<-198141+55687
# never/(ever+never) #55%. so seems right

summary(Dat$GS_cpd)
hist(Dat$GS_cpd)
sd(Dat$GS_cpd)
Dat$GS_cpd_sd<-(Dat$GS_cpd-mean(Dat$GS_cpd))/sd(Dat$GS_cpd)
Dat$GS_cpd_sd<-Dat$GS_cpd/sd(Dat$GS_cpd)
summary(Dat$GS_cpd_sd)
sd(Dat$GS_cpd_sd)


Dat1<-Dat[Dat$smoking2==1,] 
Dat0<-Dat[Dat$smoking2==0,]

# mylogit <- glm( overall_pan_inclc44_cancer ~ d5d*GS_cpd*smoking2 +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
# Model<-summary(mylogit)

# ever smokers
mylogit <- glm( overall_pan_inclc44_cancer ~ d5d*GS_cpd +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
# mylogit <- glm( overall_pan_inclc44_cancer ~ d5d*smoking2 +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")

Model1<-summary(mylogit)

# Here is an example manual computation of the slope of r holding m at 30.
slope = b[r] + 30*b[r#m] = .43420626 + 30*(-.00681441) = .22977396

-0.3028495 #CS_cpd = 0

slope1 <-   -0.3028495 + 1*0.0573056 #cpd
slope2 <-   -0.3028495 + 2*0.0573056
slope3 <-   -0.3028495 + 3*0.0573056
slope4 <-   -0.3028495 + 4*0.0573056
slope5 <-   -0.3028495 + 5*0.0573056
slope6 <-   -0.3028495 + 6*0.0573056
slope7 <-   -0.3028495 + 6*0.0573056
slope10 <-   -0.3028495 + 10*0.0573056



 -0.3028495 +0.0573056

mylogit <- glm( overall_pan_exclc44_cancer ~ d5d*GS_cpd +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat1, family = "binomial")
Model2<-summary(mylogit)


mylogit <- glm( overall_pan_inclc44_cancer ~ d5d*GS_cpd +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat0, family = "binomial")
Model<-summary(mylogit)


# Never smokers
mylogit <- glm( overall_pan_inclc44_cancer ~ d5d*GS_cpd +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat0, family = "binomial")
Model<-summary(mylogit)

# csi

mylogit <- glm( overall_pan_inclc44_cancer ~ d5d*GS_csi +age_at_recruitment +sex_num+ pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10, data = Dat, family = "binomial")
Model<-summary(mylogit)


 exp(0.0573056)  0.0302311   1.896   0.0580 .  