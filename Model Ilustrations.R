setwd("")
source("functions.R")
##################################################################################################################################################


##################################################
## 		Gibbs for Linear Regession 		##
##################################################

## Simulating data
set.seed(31415)
N<-200
r_beta<-as.matrix(c(1, 0, 2, 0, 3, 2))
r_sigma2<-1.5
aux<-LiR_base(N, r_beta, r_sigma2)
y=aux$y; Covariates=aux$Covariates
nchain=10000; burnin=2000

## Fitting the model
fit<- gibbs_abms(y=y, Covariates=Covariates, family="LiR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1,
a0=1, b0=1, count.iteration=TRUE )

summary_gibbs(fit, BF=TRUE)	#Summary results

##################################################
## 		Gibbs for Logistic regression		##
##################################################

## Simulating data
rm(fit)
set.seed(31415)
N<-200; mu_cov<-0; sigma_cov<-1
r_beta<-c(1, 0, 2, 0, 3, 2)
aux<-LoR_base(N=N, r_beta=r_beta, ni=rep(1, length(y)))
y=aux$y; Covariates=aux$Covariates
nchain=10000; burnin=2000

## Fit
fit<-gibbs_abms(y, Covariates, family="LoR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)),
count.iteration=TRUE )

summary_gibbs(fit, BF=TRUE)	#Summary results

########################################################
## 	Gibbs for Negative binomial regression		##
########################################################

## Simulating data
rm(fit)
set.seed(31415)
N<-200; mu_cov<-0; sigma_cov<-1
r_beta<-c(0.5, -0.8,  1.0,  0,  0.4, -0.7); p<-length(beta); r_r<-2
aux<-NBR_base(N, r_beta, r_r)
y=aux$y; Covariates=aux$Covariates
nchain=10000; burnin<-2000

## Gibbs
fit<-gibbs_abms(y, Covariates, family="NBR", first_excluded=0, nchain=10000, burnin=2000, tau2=1000, rho=1,
a0=1, b0=1, count.iteration=TRUE )

summary_gibbs(fit, BF=TRUE)	#Summary results

##################################################
## 		Gibbs for Quantile Regession 		##
##################################################

## Simulating data
rm(fit)
set.seed(31415)
N<-200; mu_cov<-0; sigma_cov<-1
r_beta<-c(1, 0, 2, 0, 3, 2)
r_sigma2<-1; r_alpha<-0.5
aux<-QR_base(N, r_beta, r_sigma2, r_alpha=r_alpha)
y=aux$y; Covariates=aux$Covariates
nchain<-10000; burnin<-2000

## Gibbs
fit<-gibbs_abms(y, Covariates, family="QR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, alpha=0.5,
a0=1, b0=1, count.iteration=TRUE )

summary_gibbs(fit, BF=TRUE)	#Summary results


########################################################
## 		Gibbs for Skew Normal Regession 		##
########################################################


## Simulating data
rm(fit)
set.seed(3141)
N<-200
r_beta<-c(1, 0, 2, 0, 3, 2)
r_sigma2<-1.5; r_lambda<-4
aux<-SNR_base(N, r_beta, r_sigma2, r_lambda)
y=aux$y; Covariates=aux$Covariates
nchain=10000; burnin=2000

## Gibbs
d=2; b2=1/2			#delta uniform(-1,1) prior
#d=1/2; b2=(pi^2)/4	#delta Jeffrey prior
fit<-gibbs_abms(y, Covariates, family="SNR", first_excluded=0, nchain=10000, burnin=2000, tau2=1000, rho=1,
a0=1, b0=1, d=2, b2=1/2, count.iteration=TRUE )

summary_gibbs(fit, BF=TRUE)	#Summary results