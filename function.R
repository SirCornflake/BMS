#rm(list=ls(all=TRUE))
#install.packages("BayesLogit")	#Polya gamma simulation
#install.packages("GIGrvg") 		#Generalized Inverse Gaussian simulation
#install.packages("truncnorm") 	#Truncated Normal simulation
##########################################################################

print.BMS <- function(m){
    print(m$Default)
}



## Womack prior
womack = function(K, rho)
{
  if(length(rho)!=1){stop(paste("'rho' must be a vector of size 1")) }
  if(rho<=0){stop(paste("'rho' must be >0")) }
  if((round(K)!=K) || K<=0){ stop(paste("'K' must be a positive integer")) }
  out = base::rep( -Inf, K + 1 )
  base::names(out) = 0 : K
  out[ K + 1 ] = 0
  for( k in (K-1):0 )
  {
    j = (k+1):K
    bb = out[j+1] + base::lgamma(j+1) - base::lgamma(k+1) - base::lgamma(j+1-k) + base::log(rho)
    out[k+1] = base::log( base::sum( exp(bb - base::max(bb) )) ) + base::max(bb)
    out = out - base::max(out) -base::log(base::sum( exp(out - base::max(out) ) ) )
  }
  base::exp( out + base::lbeta(c(0:K) + 1, K - c(0:K) + 1 ) + base::log(K+1) )
}
## Simulate from the Chinise restaurant distribution
rCRT<-function(n,b,c)
{
  if((round(n)!=n) || n<=0 || length(n) > 1){ stop(paste("'n' must be a positive integer")) }
  if(any(b<0)){stop(paste("'b' must be a non-negative integer")) }
  if(any(c<=0)){stop(paste("'c' must be >0")) }
  if(1!=length(b) || 1!=length(c)){stop(paste("'b' and 'c' must be of size 1")) }

  x<-base::vector(length=n)
  ifelse(b!=0, aux<-matrix(nrow=n, ncol=b), aux<-matrix(nrow=n, ncol=1) )
  for(j in 1:b)
  {
    p =c/(c +j -1)
    if(p>1){p<-1}
    base::ifelse(b!=0, aux[,j]<-stats::rbinom(n,size=1, prob=p), aux[,j]<-base::rep(0,n))
  }
  x<-apply(aux,1,sum)
  return(x)

}


## Gibbs Sampler for Bayesian variable selection models via a spike-and-slab methodology.
gibbs_abms<-function(y, Covariates, family="LiR", first_excluded=0, nchain=10000, burnin=2000, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, WomackPrior=TRUE, a_bb=1, b_bb=1, count.iteration=TRUE )
{
  if((class(y)!="numeric" && class(y)!="integer")  || is.vector(y)==FALSE){stop(paste("'y' must be a numeric or integer vector")) }
  if(class(Covariates)[1]!="data.frame" && (is.matrix(Covariates)==FALSE)){stop(paste("'Covariates' must be a data.frame or a matrix")) }
  if(length(first_excluded)>1 || first_excluded<0){stop(paste("'first_excluded' must be a non-negative integer")) }
  if(length(nchain)>1 || nchain<=0){stop(paste("'nchain' must be a positive integer")) }
  if(burnin<0 || burnin>nchain){stop(paste("'burnin' must be a non-negative integer that equal less than nchain ")) }
  if(tau2<=0){stop(paste("'tau2' must be a positive integer")) }
  if(tau2<=10) warning("It is recomended for 'tau2' equal at least 1000.")
  if(rho<=0){stop(paste("'rho' must be a positive integer")) }
  if(is.null(model_fixed)==FALSE && is.vector(model_fixed)==FALSE){ stop(paste("'model_fixed' must be a vector of size 'p-1' at most")) }
  if(is.vector(model_fixed)==TRUE && length(model_fixed)>(ncol(Covariates)) ){ stop(paste("'model_fixed' must be a vector of size 'p-1' at most")) }
  if(count.iteration!=TRUE && count.iteration!=FALSE){stop(paste("'count.iteration' must be either TRUE or FALSE")) }

  if(a0<=0){stop(paste("'a0' must be a positive real number")) }
  if(b0<=0){stop(paste("'b0' must be a positive real number")) }
  if(length(ni)!=length(y) || any(ni<0)){stop(paste("'ni' must be a positive integer vector with same size as 'y'")) }
  if(alpha>1 || alpha<0){stop(paste("'alpha' must be between 0 and 1")) }
  if(d<=0){stop(paste("'d' must be a positive real number")) }
  if(b2<=0){stop(paste("'b2' must be a positive real number")) }

  t0<-proc.time()

  ## Beta-Binomial pior2 PDF
  BetaBinomialPrior<-function(p, a, b)
  {
  	aux<-vector(length=p+1)
	for(i in 0:p)
	{
		aux[i+1]<-beta(a=a +i, b=p +b -i)/beta(a=a,b=b)
	}
	return(aux)
  }


  ##############################
  ##	   Initiating	    ##
  ##############################
  N <- length(y)
  p<-length(Covariates[1,]) +1
  if(length(colnames(Covariates))==0)
  {
    aux_prednames<-base::vector(length=ncol(Covariates)); for(i in 1:ncol(Covariates)){aux_prednames[i]<-c(paste0("X",i))}
    colnames(Covariates)<-aux_prednames
  }

  X<-cbind(base::rep(1,N), Covariates)
  X<-as.matrix(X)
  colnames(X)<-NULL
  intercept_first_excluded<- first_excluded +1
  p_selection<-p -intercept_first_excluded

  beta.chain<-matrix(1, nrow=nchain, ncol=p, byrow=TRUE)
  beta<-beta.chain[1,]
  Z.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  t.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  gammaOld<-base::rep(1, p_selection)

  if(family=="LiR")
  {
    invsigma2.chain<-matrix(1, nrow=nchain, ncol=1, byrow=TRUE)
    invsigma2<-invsigma2.chain[1,]
    sigma2.chain<-matrix(0, nrow=nchain, ncol=1, byrow=TRUE)
    Z<-matrix(y)
    t.chain[1,]<-invsigma2
  }

  if(family=="LoR")
  {
    w<-rep(1,length(y))
    kappa <- y - ni/2
    Z<-kappa/w
    t.chain[1,]<-w
  }

  if(family=="NBR")
  {
    y_unique<-sort(unique(y))
    y_unique_freq<-as.integer(table(y))
    r.chain<-matrix(1, nrow=nchain, ncol=1, byrow=TRUE)
    r<-r.chain[1,]
    l<-1
    w<-rep(1,length(y))
    Z<-matrix((y-r)/(2*w))
    t.chain[1,]<-w
  }

  if(family=="QR")
  {
    invsigma2.chain<-matrix(1, nrow=nchain, ncol=1, byrow=TRUE)
    invsigma2<-invsigma2.chain[1,]
    sigma2.chain<-matrix(0, nrow=nchain, ncol=1, byrow=TRUE)
    w<-rep(1,length(y))
    Z.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
    t.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
    delta2<-2/(alpha*(1-alpha))
    xi<-(1 -2*alpha)/(alpha*(1-alpha))
    Z<-matrix(y -xi*w)
    t.chain[1,]<-1/(delta2*w)*invsigma2
  }

  if(family=="SNR")
  {
    lambda.chain<-matrix(1, nrow=nchain, ncol=1, byrow=TRUE)
    sigma2.chain<-matrix(2, nrow=nchain, ncol=1, byrow=TRUE)
    w<-rep(1,length(y))
    lambda<-lambda.chain[1,]
    kappa<-1
    invzeta2<-1
    zeta2<-1
    v<-1
    Z<-matrix(y -kappa*w)
    t.chain[1,]<-invzeta2
  }

  model_chain<-matrix(NA, nrow=(nchain- burnin), ncol=p_selection)

  log_tau2<-base::log(tau2)
  frac_tau2<-1/tau2
  if(WomackPrior==TRUE){ log_gamma_prior<-base::log(womack(p_selection,rho)) }else{
  			    log_gamma_prior<-log( BetaBinomialPrior(p=p_selection, a=a_bb, b=b_bb) )
			  }

  InvDiag_tau2_general<-diag(frac_tau2, p)
  Omega_modified_general<-matrix(base::rep(t.chain[1,],each = p), nrow = p, ncol = N, byrow = F)
  tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
  NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
  tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z

  X_index_aux<-base::vector("list", length = 2)
  X_gamma_aux<-base::vector("list", length = 2)
  tXgamma_Omega_aux<-base::vector("list", length = 2)
  chol_V_aux<-base::vector("list", length = 2)
  V_aux<-base::vector("list", length = 2)
  m_aux<-base::vector("list", length = 2)
  X_index_aux_excluded<-c(1:intercept_first_excluded)
  DAlog_weight_model<-base::vector(length=2)


  ######################
  ## 	 Gibbs	  ##
  ######################
  for(i in 1:nchain)
  {
    if(count.iteration==TRUE){cat("  Iteracion", i, "de", nchain, "\r")}

    if(is.null(model_fixed)==TRUE)
	{
		
    		## Add-delete algorithm
    		aux_gamma<-sample(1:p_selection, size=1)
    		gammaCandidate<-gammaOld
    		gammaCandidate[aux_gamma]<-1 -gammaOld[aux_gamma]

    		for(k in 1:2)
    		{
      		if(k==1){gamma_aux<-gammaOld}else{gamma_aux<-gammaCandidate}
      		q_aux<-base::sum(gamma_aux)
      		X_index_aux[[k]]<-unique(c(X_index_aux_excluded, intercept_first_excluded +which(gamma_aux==1)))
      		X_gamma_aux[[k]]<-X[,X_index_aux[[k]]]
      		tXgamma_Omega_aux[[k]]<- tXgamma_Omega_general[X_index_aux[[k]],]
      		tXgamma_Omega_aux_Z<- tXgamma_Omega_Z_general[X_index_aux[[k]]]
      		InvDiag_tau2_aux<-InvDiag_tau2_general[1:(intercept_first_excluded +q_aux), 1:(intercept_first_excluded +q_aux)]
      		chol_V_aux[[k]]<- chol(NotInverse_V_general[X_index_aux[[k]], X_index_aux[[k]]])
      		V_aux[[k]]<- chol2inv(chol_V_aux[[k]])
      		det_V_aux<-1/prod(diag(chol_V_aux[[k]]))^(2)
      		if(det_V_aux<1e-10){det_V_aux<-1e-10}
      		m_aux[[k]]<-V_aux[[k]]%*%tXgamma_Omega_aux_Z
      		DAlog_weight_model[k]<-as.vector( (-q_aux/2)*log_tau2 +(0.5)*base::log(det_V_aux) +
                                          (0.5)*t(m_aux[[k]])%*%tXgamma_Omega_aux_Z +log_gamma_prior[q_aux +1] )
    		}
    		A<- DAlog_weight_model -base::max(DAlog_weight_model)
    		b<-exp(A)
    		DAgamma_prob<-b/base::sum(b)

    		p_gammaOld<-DAgamma_prob[1]
    		p_gammaCandidate<-DAgamma_prob[2]
    		p_selectionCandidate<-min(p_gammaCandidate/p_gammaOld,1)
    		selection<-sample(c("Candidate","Old"), size=1, prob=c(p_selectionCandidate, 1-p_selectionCandidate))
    		if(selection=="Candidate"){gamma<-gammaCandidate; k<-2}else{gamma<-gammaOld; k<-1}
    		q<-base::sum(gamma)
	}else{
		gamma<-rep(0, p-1)
		gamma[model_fixed]<-1
		k<-1
		X_index_aux[[k]]<-unique(c(X_index_aux_excluded, intercept_first_excluded +which(gamma==1)))
      	X_gamma_aux[[k]]<-X[,X_index_aux[[k]]]
		tXgamma_Omega_aux[[k]]<- tXgamma_Omega_general[X_index_aux[[k]],]
      	tXgamma_Omega_aux_Z<- tXgamma_Omega_Z_general[X_index_aux[[k]]]
      	chol_V_aux[[k]]<- chol(NotInverse_V_general[X_index_aux[[k]], X_index_aux[[k]]])
      	V_aux[[k]]<- chol2inv(chol_V_aux[[k]])
		m_aux[[k]]<-V_aux[[k]]%*%tXgamma_Omega_aux_Z			
	}

    ## Updating beta
    beta_index_0<-setdiff(1:p, X_index_aux[[k]])
    if(length(beta_index_0)==0){beta_index_0<-0}
    beta[X_index_aux[[k]]]<-t(mvtnorm::rmvnorm(n = 1, mean = m_aux[[k]], sigma = V_aux[[k]], method="chol"))
    beta[beta_index_0]<-base::rep(0,length(beta_index_0))
    beta.chain[i,] <- beta
    if(i>burnin){model_chain[i -burnin,]<-gamma}

    ## Updating gamma
    gammaOld<-gamma
    X_beta <- as.matrix(X_gamma_aux[[k]])%*%as.matrix(beta[X_index_aux[[k]]])

    if(family=="LiR")
    {
      ## Updating invsigma2=1/sigma2
      a1<-a0 +N/2
      b1<-b0 +base::sum( (y -X_beta)^2 )/2
      invsigma2<-stats::rgamma( 1, shape= a1, rate= b1 )
      invsigma2.chain[i,]<-invsigma2

      ## Updating gamma
      gammaOld<-gamma

      ## Updating other quantities
      Z<-matrix(y)
      Z.chain[i,]<-Z
      t.chain[i,]<-invsigma2
    }

    if(family=="LoR")
    {
      ## Updating w
      w <- BayesLogit::rpg.devroye(num = N, h = ni, z = as.numeric(X_beta))

      ## Updating other quantities
      Z<-kappa/w
      Z.chain[i,]<-Z
      t.chain[i,]<-w
    }

    if(family=="NBR")
    {
      ## Updating w
      w <- BayesLogit::rpg(num = N, h = y +r, z = as.numeric(X_beta))

      ## Updating l
      for(j in 1:length(y_unique))
      {
        aux_l<-rCRT(y_unique_freq[j], b=y_unique[j], c=r)
        index_l<-which(y==y_unique[j])
        l[index_l]<-aux_l
      }

      ## Updating r
      e_Xbeta<-as.vector(exp(X%*%beta))
      r<-stats::rgamma( 1, shape= a0 +base::sum(l), rate= b0 +base::sum(base::log(1 +e_Xbeta)) )
      r.chain[i,]<-r

      ## Updating other quantities
      Z<-matrix((y-r)/(2*w))
      Z.chain[i,]<-Z
      t.chain[i,]<-w
    }

    if(family=="QR")
    {
      ## Updating w
      chi<-invsigma2*((y -X_beta)^2)/delta2; psi<- invsigma2*(xi^2 +2*delta2)/delta2
      for(j in 1:length(y))
      {
        w[j]<-GIGrvg::rgig(1, lambda=1/2, chi=chi[j], psi=psi)
      }

      ## Updating invsigma2=1/sigma2
      a1<-a0 +3*N/2
      b1<-b0 +base::sum( ((y -X_beta -xi*w)^2 +2*delta2*w^2)/(2*delta2*w) )
      invsigma2<-stats::rgamma( 1, shape= a1, rate= b1 )
      invsigma2.chain[i,]<-invsigma2

      ## Updating other quantities
      Z<-matrix(y -xi*w)
      Z.chain[i,]<-Z
      t.chain[i,]<-1/(delta2*w)*invsigma2
    }

    if(family=="SNR")
    {
      ## Updating w
      mu_w<-(y*kappa -X_beta*kappa)/(kappa^2 +zeta2); sigma2_w<- (zeta2)/(kappa^2 +zeta2 )
      w<-truncnorm::rtruncnorm(n=N, a=0, b=Inf, mean=mu_w, sd=sqrt(sigma2_w))

      ## Updating invzeta2=1/zeta2
      aux_zeta1<-0.5*(N+1) +1
      aux_zeta2<-0.5*( (v*kappa^2/b2) +base::sum( (y -X_beta -kappa*w)^2 ) )
      invzeta2<-stats::rgamma( 1, shape= aux_zeta1, rate= aux_zeta2 )
      zeta2<-1/invzeta2

      ## Updating kappa
      mu_kappa<-base::sum(y*w -X_beta*w)/(v/b2 +base::sum(w^2))
      sigma2_kappa<-(zeta2)/(v/b2 +base::sum(w^2))
      kappa<-stats::rnorm(1, mean=mu_kappa, sd=sqrt(sigma2_kappa))

      ## Updating v
      shape_v<- (d+1)/2
      rate_v<- 0.5*( invzeta2*(kappa^2)/b2 +d )
      v<-stats::rgamma(1, shape=shape_v, rate=rate_v)

      ## Updating other quantities
      Z<-matrix(y -kappa*w)
      Z.chain[i,]<-Z
      t.chain[i,]<-invzeta2

      ## Saving original SN parameters
      lambda<-kappa/sqrt(zeta2)
      lambda.chain[i,]<-lambda
      sigma2<-zeta2 +kappa^2
      sigma2.chain[i,]<-sigma2
    }

    Omega_modified_general<-matrix(base::rep(t.chain[i,],each = p), nrow = p, ncol = N, byrow = F)
    tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
    NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
    tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z
  }

  t1<-proc.time()
  Time<-(t1-t0)[3]

  if(first_excluded!=0)
  {
    prednames<-colnames(Covariates)[-(1:first_excluded)]
  }else{prednames<-colnames(Covariates)}

  beta.chain<-beta.chain[(burnin +1):nchain,]

  if(family=="LiR" || family=="QR")
  {
    sigma2.chain<-1/invsigma2.chain[(burnin +1):nchain,]
    Output<-list(tau2=tau2, y=y, Covariates=Covariates, Z_chain=Z.chain, t_chain=t.chain, beta_chain=beta.chain, sigma2_chain=sigma2.chain, model_chain=model_chain, Default=list(family=family, prednames=prednames, Seconds=Time))
  }

  if(family=="LoR")
  {
    Output<-list(tau2=tau2, Covariates=Covariates, Z_chain=Z.chain, t_chain=t.chain, beta_chain=beta.chain, model_chain=model_chain, Default=list(family=family, prednames=prednames, Seconds=Time))
  }

  if(family=="NBR")
  {
    r.chain<-r.chain[(burnin +1):nchain,]
    Output<-list(tau2=tau2, y=y, Covariates=Covariates, Z_chain=Z.chain, t_chain=t.chain, beta_chain=beta.chain, r_chain=r.chain, model_chain=model_chain, Default=list(family=family, prednames=prednames, Seconds=Time))
  }

  if(family=="SNR")
  {
    lambda.chain<-lambda.chain[(burnin +1):nchain,]
    sigma2.chain<-sigma2.chain[(burnin +1):nchain,]
    Output<- list(tau2=tau2, y=y, Covariates=Covariates, Z_chain=Z.chain, t_chain=t.chain, beta_chain=beta.chain, sigma2_chain=sigma2.chain, lambda_chain=lambda.chain, model_chain=model_chain, Default=list(family=family, prednames=prednames, Seconds=Time))
  }

  class(Output) <- 'abms'
  Output
}


## Summary table for Gibbs sampler
summary_gibbs<-function(fit, BF=FALSE)
{
  if(class(fit)!="abms"){stop(paste("'fit' must be a 'abms' class object")) }
  if(BF!=TRUE && BF!=FALSE){stop(paste("'BF' must be either TRUE or FALSE")) }

  ## Auxiliar function for Summary table for Gibbs sampler
  aux_summary_gibbs<-function(fit)
  {
    ExploredModels<-unique(fit$model_chain)
    TableExploredModels<-data.frame(ExploredModels)
    TableExploredModels[,ncol(TableExploredModels) +1]<-base::rep(0,nrow(TableExploredModels))
    colnames(TableExploredModels)<-c(fit$Default$prednames, "Proportion")

    for(j in 1:nrow(ExploredModels))
    {
      num_selected<-0
      for(i in 1:nrow(fit$model_chain))
      {
        if(all(fit$model_chain[i,]==ExploredModels[j,])==TRUE){num_selected<-num_selected +1}
      }
      TableExploredModels[j,length(fit$Default$prednames)+1]<-num_selected/nrow(fit$model_chain)

    }
    TableExploredModels<-TableExploredModels[order(TableExploredModels$Proportion, decreasing=TRUE),]
  }

  aux_BF<-function(fit,ExploredModels)
  {
    aux_ExploredModels<- ExploredModels[, -ncol(ExploredModels)]
    Indexes<-vector(mode="list", length=nrow(aux_ExploredModels))
    loglik<-vector(mode="list", length=nrow(aux_ExploredModels))
    X_index_aux<-vector(mode="list", length=nrow(aux_ExploredModels))
    X_gamma<-vector(mode="list", length=nrow(aux_ExploredModels))
    InvDiag_tau2_aux<-vector(mode="list", length=nrow(aux_ExploredModels))
    X_index_aux<-vector(mode="list", length=nrow(aux_ExploredModels))
    tau2<-fit$tau2
    Covariates<-fit$Covariates
    p<-ncol(Covariates) +1
    first_excluded<-abs(ncol(Covariates) -length(aux_ExploredModels))
    intercept_first_excluded<- first_excluded +1
    p_selection<-p -intercept_first_excluded

    N<-nrow(Covariates)
    X<-as.matrix(cbind(base::rep(1,N), Covariates))
    InvDiag_tau2_general<-diag(1/tau2, p)
    X_index_aux_excluded<-c(1:intercept_first_excluded)

    for(j in 1:nrow(aux_ExploredModels))
    {
      q_aux<-base::sum(aux_ExploredModels[j,])
      X_gamma[[j]]<-as.matrix(fit$Covariates[,which(aux_ExploredModels[j,]==1)])
      InvDiag_tau2_aux[[j]]<-InvDiag_tau2_general[1:(intercept_first_excluded +q_aux), 1:(intercept_first_excluded +q_aux)]
      X_index_aux[[j]]<-unique(c(X_index_aux_excluded, intercept_first_excluded +which(aux_ExploredModels[j,]==1)))

    }
    for(i in 1:nrow(fit$beta_chain))
    {
      Z<-fit$Z_chain[i,]; t_Omega<-fit$t_chain[i,]
      Omega_modified_general<-matrix(base::rep(t_Omega,each = p), nrow = p, ncol = N, byrow = F)
      tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
      NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
      tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z
      for(j in 1:nrow(aux_ExploredModels))
      {

        tXgamma_Omega_aux_Z<- tXgamma_Omega_Z_general[X_index_aux[[j]]]
        chol_V_aux<- chol(NotInverse_V_general[X_index_aux[[j]], X_index_aux[[j]]])  #Cholesky of "not inverted" V
        V_aux<- chol2inv(chol_V_aux)
        m_aux<-V_aux%*%tXgamma_Omega_aux_Z
        det_V_aux<-1/prod(diag(chol_V_aux))^(2)

        loglik[[j]][i]<-0.5*log(det_V_aux) +(0.5)*t(m_aux)%*%tXgamma_Omega_aux_Z -0.5*p_selection*log(tau2)
      }
    }
    New_ExploredModels<-ExploredModels
    for(j in 1:nrow(aux_ExploredModels))
    {
      #Conditional_BF<-mean(exp(loglik[[1]] -loglik[[j]]))
      Marginal_BF<-mean(loglik[[1]]) -mean(loglik[[j]])
      #New_ExploredModels[j,ncol(ExploredModels) +1]<-Conditional_BF
      New_ExploredModels[j,ncol(ExploredModels) +1]<-Marginal_BF
    }
    colnames(New_ExploredModels)[-(1:ncol(ExploredModels))]<-c("Log_Marginal_BF_Estimator")
    return(New_ExploredModels)
  }

  ExploredModels<-aux_summary_gibbs(fit)
  if(BF==TRUE)
  {
    ExploredModels<-aux_BF(fit,ExploredModels)
    most_model<-ExploredModels[1,-((ncol(ExploredModels) -1):ncol(ExploredModels)) ]
  }else{most_model<-ExploredModels[1,-ncol(ExploredModels)]}

  beta_index<-c()
  for(j in 1:nrow(fit$beta_chain))
  {
    if(prod(most_model==fit$model_chain[j,])==1){beta_index<-c(beta_index,j)}
  }
  beta_chain<-fit$beta_chain[beta_index,]
  if(fit$Default$family=="LiR"){parameters<-cbind(fit$sigma2_chain[beta_index], beta_chain); aux_rownames<-c("sigma2")}
  if(fit$Default$family=="LoR"){parameters<-cbind(beta_chain); aux_rownames<-c()}
  if(fit$Default$family=="NBR"){parameters<-cbind(fit$r_chain[beta_index], beta_chain); aux_rownames<-c("r")}
  if(fit$Default$family=="QR"){parameters<-cbind(fit$sigma2_chain[beta_index], beta_chain); aux_rownames<-c("sigma2")}
  if(fit$Default$family=="SNR"){parameters<-cbind(fit$sigma2_chain[beta_index], fit$lambda_chain[beta_index], beta_chain); aux_rownames<-c("sigma2", "lambda")}

  p<-ncol(beta_chain)
  Mean<-apply(parameters, 2, mean)
  Quantile<-t(apply(parameters, 2, stats::quantile, prob=c(0.025,0.975)))
  SD<-t(apply(parameters, 2, stats::sd))
  Table<-as.data.frame(matrix(c(Mean,Quantile, SD), byrow=FALSE, ncol=4))
  colnames(Table)<-c("Mean", "2.5% quantile", "97.5% quantile", "SD")
  rownames(Table)<-c(aux_rownames, "intercept", fit$Default$prednames)

  list(Mean_IC=Table, Explored_Models=ExploredModels)
}

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


## Function for comparing BF Models. Compute the log marginal BF estimator
#Models is a matrix that need to have two or more models. Function compute the BF between the first models and the others
BF_model<-function(Models, logBF=FALSE, y, Covariates, family="LoR", 
first_excluded=0, nchain=10000, burnin=2000, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
a0=1, b0=1, d=2, b2=1/2, count.iteration=TRUE )
{
	aux_BF<-function(fit, aux_ExploredModels)
	{
    		Indexes<-vector(mode="list", length=nrow(aux_ExploredModels))
    		loglik<-vector(mode="list", length=nrow(aux_ExploredModels))
    		X_index_aux<-vector(mode="list", length=nrow(aux_ExploredModels))
    		X_gamma<-vector(mode="list", length=nrow(aux_ExploredModels))
    		InvDiag_tau2_aux<-vector(mode="list", length=nrow(aux_ExploredModels))
    		X_index_aux<-vector(mode="list", length=nrow(aux_ExploredModels))
    		tau2<-fit$tau2
    		Covariates<-fit$Covariates
    		p<-ncol(Covariates) +1
    		first_excluded<-abs(ncol(Covariates) -ncol(aux_ExploredModels))
    		intercept_first_excluded<- first_excluded +1
    		p_selection<-p -intercept_first_excluded

    		N<-nrow(Covariates)
    		X<-as.matrix(cbind(base::rep(1,N), Covariates))
    		InvDiag_tau2_general<-diag(1/tau2, p)
    		X_index_aux_excluded<-c(1:intercept_first_excluded)

    		for(j in 1:nrow(aux_ExploredModels))
    		{
      		q_aux<-base::sum(aux_ExploredModels[j,])
      		X_gamma[[j]]<-as.matrix(fit$Covariates[,which(aux_ExploredModels[j,]==1)])
      		InvDiag_tau2_aux[[j]]<-InvDiag_tau2_general[1:(intercept_first_excluded +q_aux), 1:(intercept_first_excluded +q_aux)]
      		X_index_aux[[j]]<-unique(c(X_index_aux_excluded, intercept_first_excluded +which(aux_ExploredModels[j,]==1)))

    		}
    		for(i in 1:nrow(fit$beta_chain))
    		{
      		Z<-fit$Z_chain[i,]; t_Omega<-fit$t_chain[i,]
      		Omega_modified_general<-matrix(base::rep(t_Omega,each = p), nrow = p, ncol = N, byrow = F)
      		tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
      		NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
      		tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z
      		for(j in 1:nrow(aux_ExploredModels))
      		{
				tXgamma_Omega_aux_Z<- tXgamma_Omega_Z_general[X_index_aux[[j]]]
				chol_V_aux<- chol(NotInverse_V_general[X_index_aux[[j]], X_index_aux[[j]]])  #Cholesky of "not inverted" V
				V_aux<- chol2inv(chol_V_aux)
				m_aux<-V_aux%*%tXgamma_Omega_aux_Z
				det_V_aux<-1/prod(diag(chol_V_aux))^(2)

				loglik[[j]][i]<-0.5*log(det_V_aux) +(0.5)*t(m_aux)%*%tXgamma_Omega_aux_Z -0.5*p_selection*log(tau2)
			}
		}

		New_ExploredModels<- data.frame(aux_ExploredModels)
    		for(j in 1:nrow(aux_ExploredModels))
    		{
      		Conditional_BF<-mean(exp(loglik[[1]] -loglik[[j]]))
      		Marginal_BF<-exp( mean(loglik[[1]]) -mean(loglik[[j]]) )
      		New_ExploredModels[j,ncol( aux_ExploredModels) +1]<-Conditional_BF
      		New_ExploredModels[j,ncol( aux_ExploredModels) +2]<-Marginal_BF
    		}
    		colnames(New_ExploredModels)[-(1:ncol(aux_ExploredModels))]<-c("Conditional_BF", "Marginal_BF_Estimator")
    		list(loglik=loglik, ExploredModels=New_ExploredModels)
  	}

	## Fitting fixing at model1 and model2
	fit<-list()
	for(i in 1:length(Models))
	{
		fit[[i]]<-gibbs_abms(y, Covariates, family=family, first_excluded=first_excluded, nchain=nchain, burnin=burnin, tau2=tau2, rho=rho, ni=ni, alpha=alpha,
		a0=a0, b0=b0, d=d, b2=b2, model_fixed=Models[[i]], count.iteration=TRUE )
	}

	## Creating aux_ExploredModels table
	aux_ExploredModels<-rbind()
	for(i in 1:length(Models))
	{
		aux_gamma<-rep(0, p-1)
		if( all(Models[[i]]==0) ){ }else{ aux_gamma[Models[[i]]]<-1  }
		aux_ExploredModels<-rbind(aux_ExploredModels, aux_gamma)
	}
	aux_gammaName<-c()
	for(i in 1:length(Models)){aux_gammaName<-c(aux_gammaName, paste0("gamma",i))}
	rownames(aux_ExploredModels)<-aux_gammaName

	## BF where fit1 and fit2 where first use the (w,psi) of model1 and then the (w,psi) of model2, etc.
	BF_fit<-list()
	for(i in 1:length(Models))
	{
		BF_fit[[i]]<-aux_BF(fit[[i]], aux_ExploredModels)	#BF using w and psi of fit[[i]]
	}

	## BF where fit1 and fit2 use their respective (w,psi)
	loglik<-vector(mode="list", length(Models))
	for(i in 1:length(Models))
	{
		loglik[[i]]<-BF_fit[[i]]$loglik[[i]]	#loglik of Models[[i]] with it (w,psi)
	}

	log_BF_marginal<- data.frame(aux_ExploredModels)
    	for(j in 1:nrow(aux_ExploredModels))
    	{
      	Marginal_logBF<-mean(loglik[[1]]) -mean(loglik[[j]])
      	log_BF_marginal[j,ncol( aux_ExploredModels) +1]<-Marginal_logBF

    	}
    	colnames(log_BF_marginal)[-(1:ncol(aux_ExploredModels))]<-c("log_Marginal_BF_Estimator")

    	list("log_BF_marginal"=log_BF_marginal)
 
}


