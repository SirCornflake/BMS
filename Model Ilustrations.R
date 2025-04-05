setwd("")
source("function.R")
##################################################################################################################################################

################################################################################
##				Some functions, Just run this					##
################################################################################

## Function that generates "y" vector and "Covariates" data.frame for LiR
LiR_base<-function(N, r_beta, r_sigma2)
{
	r_beta<-as.matrix(r_beta)
	p<-length(r_beta)
	X<-matrix( c(rep(1, N), rnorm((p -1)*N)), ncol=p )
	Xbeta<-X%*%r_beta
	y<-rnorm(N, mean=Xbeta, sd=sqrt(r_sigma2))
	Covariates<-X[,2:(length(r_beta))]
	aa<-c()
	for(i in 1:(p-1)){aa<-c(aa, paste0("X",i))}	
	colnames(Covariates)<-aa
	list(y=y, Covariates=Covariates)
}

## Function that generates "y" vector and "Covariates" data.frame for LoR
LoR_base<-function(N, r_beta, ni=rep(1, N))
{
	gen_base_binomial_reg<- function(beta, Covariates, N, ni)
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
	  
		base <- data.frame(y = y, 
	                   n_failure = ni - y, 
	                   ni = ni,  #
	                   Covariates = Covariates)  
		colnames(base)[-(1:3)]<-colnames(Covariates)
		return(base)
	}
	mu_cov<-0; sigma_cov<-1
	p<-length(r_beta)
	aux_cov<-rnorm((p-1)*N, mu_cov, sigma_cov)
	Covariates<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N))
	aa<-c()
	for(i in 1:(p-1)){aa<-c(aa, paste0("X",i))}	
	colnames(Covariates)<-aa
	base<-gen_base_binomial_reg(N=N, beta=r_beta, Covariates=Covariates, ni=ni)
	y<-base$y; Covariates<-as.matrix(base[,-(1:3)]); 
	list(y=y, Covariates=Covariates, ni=ni)
	
}

## Function that generates "y" vector and "Covariates" data.frame for NBR
NBR_base<-function(N, r_beta, r_r)
{

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
	
		base <- data.frame(y = y,
	                   Covariates = Covariates)  
		colnames(base)[-1]<-colnames(Covariates)
		return(base)
	}

	mu_cov<-0; sigma_cov<-1
	p<-length(r_beta)
	aux_cov<-rnorm((p-1)*N, mu_cov, sigma_cov)
	Covariates<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N))
	aa<-c()
	for(i in 1:(p-1)){aa<-c(aa, paste0("X",i))}	
	colnames(Covariates)<-aa
	base<-gen_base_NegBinomial_reg(N, r_beta, r_r, Covariates=Covariates)
	y<-base$y; Covariates<-as.matrix(base[,-1])
	list(y=y, Covariates=Covariates)
}

## Function that generates "y" vector and "Covariates" data.frame for QR
QR_base<-function(N, r_beta, r_sigma2, r_alpha)
{
	mu_cov<-0; sigma_cov<-1
	p<-length(r_beta)
	aux_cov<-rnorm((p-1)*N, mu_cov, sigma_cov)
	X<-matrix( c(rep(1, N), aux_cov), ncol=p )
	Xbeta<-X%*%r_beta
	
	y<-vector(length=N)
	w<-rexp(N, rate=1/r_sigma2)
	varphi<-(1-2*r_alpha)/(r_alpha*(1-r_alpha))
	delta2<-2/(r_alpha*(1-r_alpha))
	y<- Xbeta +varphi*w +sqrt(r_sigma2*delta2*w)*rnorm(N, mean=0, sd=1)	#ALD data
	y<-as.vector(y)

	Covariates<-X[,2:(length(r_beta))]
	aa<-c()
	for(i in 1:(p-1)){aa<-c(aa, paste0("X",i))}	
	colnames(Covariates)<-aa
	list(y=y, Covariates=Covariates, r_alpha=r_alpha)
}

## Function that generates "y" vector and "Covariates" data.frame for SNR
SNR_base<-function(N, r_beta, r_sigma2, r_lambda)
{
	p<-length(r_beta)
	r_delta<-r_lambda/sqrt(r_lambda^2 +1)
	X<-matrix( c(rep(1, N), rnorm((p -1)*N)), ncol=p )
	Xbeta<-X%*%r_beta
	r_w<-abs(rnorm(N, mean=0, sd=1))
	r_delta<-r_lambda/sqrt(r_lambda^2 +1)
	y<-Xbeta +sqrt(r_sigma2)*r_delta*r_w +sqrt(r_sigma2*(1 -r_delta^2))*rnorm(N, mean=0, sd=1)	#Simulating from Skew Normal
	y<-as.vector(y)
	Covariates<-X[,2:(length(r_beta))]
	aa<-c()
	for(i in 1:(p-1)){aa<-c(aa, paste0("X",i))}	
	colnames(Covariates)<-aa
	list(y=y, Covariates=Covariates)

}

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

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

## Fitting the model with model selection, using the Womack prior
fit<- gibbs_abms(y=y, Covariates=Covariates, family="LiR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1,
a0=1, b0=1, count.iteration=TRUE)
summary_gibbs(fit)	#Summary results

## Fitting the model with model selection, using the Beta-Binomial prior
p_selection<- length(r_beta) -1 	#Because the intercept is excluded from the selection process
fit<- gibbs_abms(y=y, Covariates=Covariates, family="LiR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1,
a0=1, b0=1, WomackPrior=FALSE, a_bb=1, b_bb=p_selection^(1), count.iteration=TRUE)
summary_gibbs(fit)	#Summary results

## Fitting the model under the full model (all five predictors), without model selection
fit<- gibbs_abms(y=y, Covariates=Covariates, family="LiR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1,
a0=1, b0=1, model_fixed=c(1,2,3,4,5), count.iteration=TRUE)
summary_gibbs(fit)


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

## Fitting the model with model selection
fit<-gibbs_abms(y, Covariates, family="LoR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)),
count.iteration=TRUE)

summary_gibbs(fit)	#Summary results


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

## Fitting the model with model selection
fit<-gibbs_abms(y, Covariates, family="NBR", first_excluded=0, nchain=10000, burnin=2000, tau2=1000, rho=1,
a0=1, b0=1, count.iteration=TRUE )

summary_gibbs(fit)	#Summary results


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

## Fitting the model with model selection
fit<-gibbs_abms(y, Covariates, family="QR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, alpha=0.5,
a0=1, b0=1, count.iteration=TRUE )

summary_gibbs(fit)	#Summary results


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

## Fitting the model with model selection
d=2; b2=1/2			#delta uniform(-1,1) prior
#d=1/2; b2=(pi^2)/4	#delta Jeffrey prior
fit<-gibbs_abms(y, Covariates, family="SNR", first_excluded=0, nchain=10000, burnin=2000, tau2=1000, rho=1,
a0=1, b0=1, d=2, b2=1/2, count.iteration=TRUE )

summary_gibbs(fit)	#Summary results