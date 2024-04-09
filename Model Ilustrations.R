setwd("")
source("functions.R")
##################################################################################################################################################

########################################################
## 			Just run this first			##
########################################################

## Logistic Regression Data generator
gen_base_binomial_reg<- function(beta, Covariates, N, ni=rep(1, N))
{
	p<-length(Covariates[1,]) +1; x<-c()
	for(i in 1:(p -1))
	{
			x<-c(x,Covariates[,i])
	}
	X<-matrix(c(rep(1,N),x), ncol = p, nrow = N, byrow=FALSE)


	beta_matrix<-as.matrix(beta)
	prob<- exp(X%*%beta_matrix)/(1 + exp(X%*%beta_matrix))
  
	y<-vector(length=N)
	for(i in 1:N)
	{
		y[i] <- rbinom(n = 1, size = ni[i], prob = prob[i])
	}
  
	base <- data.frame(y = y, # success
                   n_failure = ni - y, # failure
                   ni = ni,  #total number of trials
                   Covariates = Covariates)  #Predictor
	colnames(base)[-(1:3)]<-colnames(Covariates)
	return(base)
}

## Negative Binomial Regression Data generator
gen_base_NegBinomial_reg<- function(N, beta, r, Covariates)
{
	X<-as.matrix(Covariates)
	X<-cbind(rep(1,N),X)
	colnames(X)<-NULL

	beta_matrix<-as.matrix(beta)
	prob<- 1/(1 + exp(X%*%beta_matrix))
  
	y<-vector(length=N)

	for(i in 1:N)
	{
		y[i]<-rnbinom(1, size=r, prob= prob[i])
	}

	base <- data.frame(y = y, # count
                   Covariates = Covariates)  #Predictor
	colnames(base)[-1]<-colnames(Covariates)
	return(base)
}



##################################################################################################################################################

##################################################
## 		Gibbs for Linear Regession 		##
##################################################

## Simulating data
set.seed(31415)
N<-200
r_beta<-as.matrix(c(1, 0, 2, 0, 3, 2))
r_p<-length(r_beta)
r_sigma2<-1.5
X<-matrix( c(rep(1, N), rnorm((r_p -1)*N)), ncol=r_p )
Xbeta<-X%*%r_beta
y<-rnorm(N, mean=Xbeta , sd=sqrt(r_sigma2))
covariables<-X[,2:(length(r_beta))]; colnames(covariables)<-c("X1", "X2", "X3", "X4", "X5")
nchain=10000; burnin=2000

## Fitting the model
fit <- gibbs_LiR(y, covariables, first_excluded=0, nchain, burnin, tau2=1000, rho=1, a0=1, b0=1,
beta.ini=rep(1,length(covariables[1,]) +1), invsigma2.ini=1, count.iteration=TRUE )

summary_gibbs(fit)	#Summary results




##################################################
## 		Gibbs for Logistic regression		##
##################################################

## Simulating data
rm(fit)
set.seed(31415)
N<-200; mu_cov<-0; sigma_cov<-1
beta<-c(1, 0, 2, 0, 3, 2); p<-length(beta)
aux_cov<-rnorm((p-1)*N, mu_cov, sigma_cov)
covariables<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N)); colnames(covariables)<-c("X1", "X2", "X3", "X4", "X5")
base<-gen_base_binomial_reg(beta, covariables, N)
y<-base$y; ni<-base$ni; covariables<-as.matrix(base[,-(1:3)]); first_excluded<-0
nchain=10000; burnin=2000


## Fit
fit<-gibbs_LoR(y, ni, covariables, first_excluded=0, nchain, burnin, tau2=1000, rho=1,
beta.ini=rep(1,length(covariables[1,]) +1), w.ini=rep(1,length(y)), count.iteration=TRUE )
summary_gibbs(fit)	#Summary results




########################################################
## 	Gibbs for Negative binomial regression		##
########################################################

## Simulating data
rm(fit)
set.seed(31415)
N<-200; mu_cov<-0; sigma_cov<-1
beta<-c(0.5, -0.8,  1.0,  0,  0.4, -0.7); p<-length(beta); r<-2
aux_cov<-rnorm((p-1)*N, mu_cov, sigma_cov)
covariables<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N)); colnames(covariables)<-c("X1", "X2", "X3", "X4", "X5")
base<-gen_base_NegBinomial_reg(N, beta, r, covariables)
y<-base$y; covariables<-as.matrix(base[,-1]); first_excluded<-0; 
nchain=10000; burnin<-2000

## Gibbs
fit<-gibbs_NBR(y, covariables, first_excluded=0, nchain, burnin, tau2=1000, rho=1, a0=1, b0=1,
beta.ini=rep(1,length(covariables[1,]) +1), r.ini=1, w.ini=rep(1,length(y)), l.ini=rep(1,length(y)), count.iteration=TRUE )
summary_gibbs(fit)	#Summary results



##################################################
## 		Gibbs for Quantile Regession 		##
##################################################

## Simulating data
rm(fit)
set.seed(31415)
N<-200; mu_cov<-0; sigma_cov<-1
beta<-c(1, 0, 2, 0, 3, 2); p<-length(beta)
sigma2<-1; alpha<-0.5
aux_cov<-rnorm((p-1)*N, mu_cov, sigma_cov)
X<-matrix( c(rep(1, N), aux_cov), ncol=p )
Xbeta<-X%*%beta

y<-vector(length=N)
w<-rexp(N, rate=1/r_sigma2)
varphi<-(1-2*alpha)/(alpha*(1-alpha))
delta2<-2/(alpha*(1-alpha))
y<- Xbeta +varphi*w +sqrt(r_sigma2*delta2*w)*rnorm(N, mean=0, sd=1)	#ALD data

covariables<-X[,2:(length(beta))]; colnames(covariables)<-c("X1", "X2", "X3", "X4", "X5")
nchain<-10000; burnin<-2000

## Gibbs
fit<-gibbs_QR(y, covariables, first_excluded=0, nchain, burnin, alpha=alpha, tau2=1000, rho=1, a0=1, b0=1,
beta.ini=rep(1,length(covariables[1,]) +1), invsigma2.ini=1, w.ini=rep(1,length(y)), count.iteration=TRUE )
summary_gibbs(fit)	#Summary results


########################################################
## 		Gibbs for Skew Normal Regession 		##
########################################################

## Simulating data
rm(fit)
set.seed(3141)
N<-200
r_beta<-c(1, 0, 2, 0, 3, 2)
r_p<-length(r_beta)
r_sigma2<-1.5; r_lambda<-4
r_delta<-r_lambda/sqrt(r_lambda^2 +1)
X<-matrix( c(rep(1, N), rnorm((r_p -1)*N)), ncol=r_p )
Xbeta<-X%*%r_beta
r_w<-abs(rnorm(N, mean=0, sd=1))
r_delta<-r_lambda/sqrt(r_lambda^2 +1)
y<-Xbeta +sqrt(r_sigma2)*r_delta*r_w +sqrt(r_sigma2*(1 -r_delta^2))*rnorm(N, mean=0, sd=1)	#Simulating from Skew Normal
covariables<-X[,2:(length(r_beta))];  colnames(covariables)<-c("X1", "X2", "X3", "X4", "X5")
nchain=1000; burnin=200

## Gibbs
d=2; b2=1/2			#delta uniform(-1,1) prior
#d=1/2; b2=(pi^2)/4	#delta Jeffrey prior
fit<-gibbs_SNR(y, covariables, first_excluded=0, nchain, burnin, tau2=1000, rho=1, d, b2=1,
beta.ini=rep(1,length(covariables[1,]) +1), kappa.ini=1, invzeta2.ini=1, v.ini=1, w.ini=rep(1,length(y)), count.iteration=TRUE )
summary_gibbs(fit)	#Summary results

