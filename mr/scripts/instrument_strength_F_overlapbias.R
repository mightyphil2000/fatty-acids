# order
# exposure	instrument includes FADS region?
# Included	
# Alpha-linolenic acid (18:3n3)	no
# Alpha-linolenic acid (18:3n3)	yes
# Docosahexaenoic acid (22:6n3)	no
# Docosahexaenoic acid (22:6n3)	yes
# Docosapentaenoic acid (22:5n3)	no
# Docosapentaenoic acid (22:5n3)	yes
# Omega-3 fatty acids	no
# Omega-3 fatty acids	yes
# Arachidonic acid (20:4n6)	no
# Arachidonic acid (20:4n6)	yes
# Dihomo-gamma-linolenic acid (20:3n6)	no
# Dihomo-gamma-linolenic acid (20:3n6)	yes
# Gamma-linolenic acid (18:3n6)	no
# Gamma-linolenic acid (18:3n6)	yes
# Linoleic acid (18:2n6)	no
# Linoleic acid (18:2n6)	yes
# Omega-6 fatty acids	no
# Omega-6 fatty acids	yes

nsnps<-c(1,2)
r2<-c(0.005,0.0340)
n<-c(8866,8866)
k<-nsnps
F<-round(r2*(n-1-k)/((1-r2)*k ),2)

nsnps<-c(1,2,23,24,2,3,29,30,1,2,2,3,4,5,35,36,47,48)
r2<-c(0.004,0.038,0.016,0.051,0.025,0.099,0.022,0.089,0.005,0.309,0.038,0.119,0.019,0.053,0.024,0.040,0.046,0.046)
n<-c(8866,8866,114999,114999,7791,7791,114999,114999,8631,8631,7596,7596,7596,7596,114999,114999,114999,114999)
k<-nsnps

F<-round(r2*(n-1-k)/((1-r2)*k ),2)


N<-n
K<-k
rsq<-r2
F[Pos]==expf[Pos]
Pos<-which(expf!=F)
F[F!=expf]
expf = (N-K-1)/K * rsq/(1-rsq) # expf is expected value of the F statistic
# N is sample size
# K is number of genetic variants
# rsq is the expected value of R^2
# (otherwise expf can be specified directly)
olsbias<-log(1.10)
overlap.prop<-c(0.3183857)
bias = olsbias*overlap.prop*(1/expf)
# bias is the bias of the IV estimate under the null
# olsbias is the bias of the OLS estimate (observational
# estimate for binary outcome)
# overlap.prop is the proportion of overlap
# between the samples (between 0 and 1)
var = var_y/(N*var_x*rsq) # var is the variance of the IV estimate
# (continuous outcome)
# var_x is the variance of the risk factor
# var_y is the variance of the outcome
var_x<-1
prop.case<-0.1
var = 1/(N*rsq*var_x*prop.case*(1-prop.case))
# var is the variance of the IV estimate
# (binary outcome)
# prop.case is the proportion of cases (between 0 and 1)
type1err = 2-pnorm(1.96+bias/sqrt(var))-pnorm(1.96-bias/sqrt(var))
# type1err is the estimated Type 1 error rate
# rate under the null for a nominal 5% two-sided
# significance level