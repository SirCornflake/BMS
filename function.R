#rm(list=ls(all=TRUE))
#install.packages("BayesLogit")	#Polya gamma simulation
#install.packages("GIGrvg") 		#Generalized Inverse Gaussian simulation
#install.packages("truncnorm") 	#Truncated Normal simulation
##########################################################################

print.BMS <- function(m){
    print(m$Default)
}


## Womack prior
womack = function(K, podds)
{
  out = base::rep( -Inf, K + 1 )
  base::names(out) = 0 : K
  out[ K + 1 ] = 0
  for( k in (K-1):0 )
  {
    j = (k+1):K
    bb = out[j+1] + base::lgamma(j+1) - base::lgamma(k+1) - base::lgamma(j+1-k) + base::log(podds)
    out[k+1] = base::log( base::sum( exp(bb - base::max(bb) )) ) + base::max(bb)
    out = out - base::max(out) -base::log(base::sum( exp(out - base::max(out) ) ) )
  }
  base::exp( out + base::lbeta(c(0:K) + 1, K - c(0:K) + 1 ) + base::log(K+1) )
}

## Simulate from the Chinise restaurant distribution
rCRT<-function(n,b,c)
{
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

## Gibbs for Linear regression
gibbs_LiR<-function(y, Covariates, first_excluded=0, nchain=1000, burnin=0, tau2=1000, rho=1, a0=1, b0=1,
                    beta.ini=base::rep(1,length(Covariates[1,]) +1), invsigma2.ini=1, count.iteration=TRUE )
{
  t0<-proc.time()

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


  beta.chain<-matrix(beta.ini, nrow=nchain, ncol=p, byrow=TRUE)
  beta<-beta.chain[1,]
  invsigma2.chain<-matrix(c(invsigma2.ini), nrow=nchain, ncol=1, byrow=TRUE)
  invsigma2<-invsigma2.chain[1,]
  sigma2.chain<-matrix(0, nrow=nchain, ncol=1, byrow=TRUE)
  Z.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  eta.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)

  gammaOld<-base::rep(1, p_selection)

  Z<-matrix(y)

  model_chain<-matrix(NA, nrow=(nchain- burnin), ncol=p_selection)

  log_tau2<-base::log(tau2)
  frac_tau2<-1/tau2
  log_womack_prior<-base::log(womack(p_selection,rho))
  InvDiag_tau2_general<-diag(frac_tau2, p)
  Omega_modified_general<-matrix(base::rep(invsigma2, each = p), nrow = p, ncol = N, byrow = F)
  tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
  NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general		#"Not invertet" V for all possible model
  tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z			#tXgamma_Omega_general%*%Z for all possible models

  X_index_aux<-base::vector("list", length = 2)			#For not computing columns matrix X given gamma again in parameter simulation
  X_gamma_aux<-base::vector("list", length = 2)			#For not computing matrix X given gamma again in parameter simulation
  tXgamma_Omega_aux<-base::vector("list", length = 2)		#For not computing matrix X^(t)%*%Omega given  this again in parameter simulation
  chol_V_aux<-base::vector("list", length = 2)			#For saving Cholesky decomposition of V^-1 and use it to compute the determinant of V and possibly beta simulation
  V_aux<-base::vector("list", length = 2)				#For not computing matrix V again in beta simulation
  m_aux<-base::vector("list", length = 2)				#For not computing vector m again in beta simulation
  X_index_aux_excluded<-c(1:intercept_first_excluded)		#Here we will save the selected columns (Note that we start with the ones that were excluded of the selection process)
  DAlog_weight_model<-base::vector(length=2)	#Vector containig those two log-weight


  ######################
  ## 	 Gibbs	  ##
  ######################
  for(i in 1:nchain)
  {
    if(count.iteration==TRUE){cat("  Iteracion", i, "de", nchain, "\r")}

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
      chol_V_aux[[k]]<- chol(NotInverse_V_general[X_index_aux[[k]], X_index_aux[[k]]])  #Cholesky of "not inverted" V
      V_aux[[k]]<- chol2inv(chol_V_aux[[k]])		#Inverting "not inverted" V using Cholesky
      det_V_aux<-1/prod(diag(chol_V_aux[[k]]))^(2)	#Saving determinant of V using Cholesky
      if(det_V_aux<1e-10){det_V_aux<-1e-10}
      m_aux[[k]]<-V_aux[[k]]%*%tXgamma_Omega_aux_Z
      DAlog_weight_model[k]<-as.vector( (-q_aux/2)*log_tau2 +(0.5)*base::log(det_V_aux) +
                                          (0.5)*t(m_aux[[k]])%*%tXgamma_Omega_aux_Z +log_womack_prior[q_aux +1] )
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


    ## Updating beta
    beta_index_0<-setdiff(1:p, X_index_aux[[k]])
    if(length(beta_index_0)==0){beta_index_0<-0}
    beta[X_index_aux[[k]]]<-m_aux[[k]] +chol(V_aux[[k]])%*%as.matrix(rnorm(length(X_index_aux[[k]])))
    beta[beta_index_0]<-base::rep(0,length(beta_index_0))			#simulated betas under the selection process that were NOT selected (Dirac measure at 0)
    beta.chain[i,] <- beta						#beta chain for this iteration
    if(i>burnin){model_chain[i -burnin,]<-gamma}			#Counting the selected model after burn-in


    ## Updating invsigma2=1/sigma2
    X_beta <- as.matrix(X_gamma_aux[[k]])%*%as.matrix(beta[X_index_aux[[k]]])
    a1<-a0 +N/2
    b1<-b0 +base::sum( (y -X_beta)^2 )/2
    invsigma2<-rgamma( 1, shape= a1, rate= b1 )
    invsigma2.chain[i,]<-invsigma2

    ## Updating gamma
    gammaOld<-gamma

    ## Updating other quantities
    Z<-matrix(y)
    Z.chain[i,]<-Z
    eta.chain[i,]<-invsigma2
    Omega_modified_general<-matrix(base::rep(invsigma2, each = p), nrow = p, ncol = N, byrow = F)
    tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
    NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
    tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z

  }
  beta.chain<-beta.chain[(burnin +1):nchain,]
  sigma2.chain<-1/invsigma2.chain[(burnin +1):nchain,]

  if(first_excluded!=0)
  {
    prednames<-colnames(Covariates)[-(1:first_excluded)]
  }else{prednames<-colnames(Covariates)}
  t1<-proc.time()
  Time<-(t1-t0)[3]

  Output<-list(tau2=tau2, Covariates=Covariates, Z_chain=Z.chain, eta_chain=eta.chain, beta_chain=beta.chain, sigma2_chain=sigma2.chain, model_chain=model_chain, Default=list(Regression="LiR", prednames=prednames, Seconds=Time))
  class(Output) <- 'BMS'
  Output

}

## Gibbs Logistic Regression
gibbs_LoR<-function(y, ni, Covariates, first_excluded=0, nchain=1000, burnin=0, tau2=100, rho=1,
                    beta.ini=base::rep(1,length(Covariates[1,]) +1), w.ini=base::rep(1,length(y)), count.iteration=TRUE )
{
  t0<-proc.time()
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

  beta.chain<-matrix(beta.ini, nrow=nchain, ncol=p, byrow=TRUE)
  beta<-beta.chain[1,]
  Z.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  eta.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)

  gammaOld<-base::rep(1, p_selection)

  w<-w.ini
  kappa <- y - ni/2
  Z<-kappa/w

  model_chain<-matrix(NA, nrow=(nchain- burnin), ncol=p_selection)

  log_tau2<-base::log(tau2)
  frac_tau2<-1/tau2
  log_womack_prior<-base::log(womack(p_selection,rho))
  InvDiag_tau2_general<-diag(frac_tau2, p)
  Omega_modified_general<-matrix(base::rep(w,each = p), nrow = p, ncol = N, byrow = F)
  tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
  NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general		#"Not invertet" V for all possible model
  tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z			#tXgamma_Omega_general%*%Z for all possible models


  X_index_aux<-base::vector("list", length = 2)			#For not computing columns matrix X given gamma again in parameter simulation
  X_gamma_aux<-base::vector("list", length = 2)			#For not computing matrix X given gamma again in parameter simulation
  tXgamma_Omega_aux<-base::vector("list", length = 2)		#For not computing matrix X^(t)%*%Omega given  this again in parameter simulation
  chol_V_aux<-base::vector("list", length = 2)			#For saving Cholesky decomposition of V^-1 and use it to compute the determinant of V and possibly beta simulation
  V_aux<-base::vector("list", length = 2)				#For not computing matrix V again in beta simulation
  m_aux<-base::vector("list", length = 2)				#For not computing vector m again in beta simulation
  X_index_aux_excluded<-c(1:intercept_first_excluded)		#Here we will save the selected columns (Note that we start with the ones that were excluded of the selection process)
  DAlog_weight_model<-base::vector(length=2)	#Vector containig those two log-weight


  ######################
  ## 	 Gibbs	  ##
  ######################
  for(i in 1:nchain)
  {
    if(count.iteration==TRUE){cat("  Iteracion", i, "de", nchain, "\r")}

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
      chol_V_aux[[k]]<- chol(NotInverse_V_general[X_index_aux[[k]], X_index_aux[[k]]])  #Cholesky of "not inverted" V
      V_aux[[k]]<- chol2inv(chol_V_aux[[k]])		#Inverting "not inverted" V using Cholesky
      det_V_aux<-1/prod(diag(chol_V_aux[[k]]))^(2)	#Saving determinant of V using Cholesky
      if(det_V_aux<1e-10){det_V_aux<-1e-10}
      m_aux[[k]]<-V_aux[[k]]%*%tXgamma_Omega_aux_Z
      DAlog_weight_model[k]<-as.vector( (-q_aux/2)*log_tau2 +(0.5)*base::log(det_V_aux) +
                                          (0.5)*t(m_aux[[k]])%*%tXgamma_Omega_aux_Z +log_womack_prior[q_aux +1] )
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


    ## Updating beta
    beta_index_0<-setdiff(1:p, X_index_aux[[k]])
    if(length(beta_index_0)==0){beta_index_0<-0}
    beta[X_index_aux[[k]]]<-m_aux[[k]] +chol(V_aux[[k]])%*%as.matrix(rnorm(length(X_index_aux[[k]])))
    beta[beta_index_0]<-base::rep(0,length(beta_index_0))			#simulated betas under the selection process that were NOT selected (Dirac measure at 0)
    beta.chain[i,] <- beta						#beta chain for this iteration
    if(i>burnin){model_chain[i -burnin,]<-gamma}			#Counting the selected model after burn-in


    ## Updating w
    #beta_index_1<- intercept_first_excluded +which(gamma==1)
    X_beta <- X_gamma_aux[[k]]	%*%beta[X_index_aux[[k]]]
    w <- BayesLogit::rpg.devroye(num = N, h = ni, z = as.numeric(X_beta))


    ## Updating gamma
    gammaOld<-gamma

    ## Updating other quantities
    Z<-kappa/w
    Z.chain[i,]<-Z
    eta.chain[i,]<-w
    Omega_modified_general<-matrix(base::rep(w,each = p), nrow = p, ncol = N, byrow = F)
    tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
    NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
    tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z
  }
  beta.chain<-beta.chain[(burnin +1):nchain,]

  if(first_excluded!=0)
  {
    prednames<-colnames(Covariates)[-(1:first_excluded)]
  }else{prednames<-colnames(Covariates)}
  t1<-proc.time()
  Time<-(t1-t0)[3]

  Output<-list(tau2=tau2, Covariates=Covariates, Z_chain=Z.chain, eta_chain=eta.chain, beta_chain=beta.chain, model_chain=model_chain, Default=list(Regression="LoR", prednames=prednames, Seconds=Time))
  class(Output) <- 'BMS'
  Output
}


## Gibbs for Negative Binomial regression
gibbs_NBR<-function(y, Covariates, first_excluded=0, nchain=1000, burnin=0, tau2=1000, rho=1, a0=1, b0=1,
                    beta.ini=base::rep(1,length(Covariates[1,]) +1), r.ini=1, w.ini=base::rep(1,length(y)), l.ini=base::rep(1,length(y)), count.iteration=TRUE )
{
  t0<-proc.time()
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

  beta.chain<-matrix(beta.ini, nrow=nchain, ncol=p, byrow=TRUE)
  r.chain<-matrix(c(r.ini), nrow=nchain, ncol=1, byrow=TRUE)
  Z.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  eta.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  beta<-beta.chain[1,]
  r<-r.chain[1,]
  l<-l.ini
  w<-w.ini

  gammaOld<-base::rep(1, p_selection)

  Z<-matrix((y-r)/(2*w))

  model_chain<-matrix(NA, nrow=(nchain- burnin), ncol=p_selection)

  log_tau2<-base::log(tau2)
  frac_tau2<-1/tau2
  log_womack_prior<-base::log(womack(p_selection,rho))
  InvDiag_tau2_general<-diag(frac_tau2, p)
  Omega_modified_general<-matrix(base::rep(w,each = p), nrow = p, ncol = N, byrow = F)
  tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
  NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
  tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z
  y_unique<-sort(unique(y))
  y_unique_freq<-as.integer(table(y))

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
      chol_V_aux[[k]]<- chol(NotInverse_V_general[X_index_aux[[k]], X_index_aux[[k]]])  #Cholesky of "not inverted" V
      V_aux[[k]]<- chol2inv(chol_V_aux[[k]])		#Inverting "not inverted" V using Cholesky
      det_V_aux<-1/prod(diag(chol_V_aux[[k]]))^(2)	#Saving determinant of V using Cholesky
      if(det_V_aux<1e-10){det_V_aux<-1e-10}
      m_aux[[k]]<-V_aux[[k]]%*%tXgamma_Omega_aux_Z
      DAlog_weight_model[k]<-as.vector( (-q_aux/2)*log_tau2 +(0.5)*base::log(det_V_aux) +
                                          (0.5)*t(m_aux[[k]])%*%tXgamma_Omega_aux_Z +log_womack_prior[q_aux +1] )
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



    ## Updating beta
    beta_index_0<-setdiff(1:p, X_index_aux[[k]])
    if(length(beta_index_0)==0){beta_index_0<-0}
    beta[X_index_aux[[k]]]<-m_aux[[k]] +chol(V_aux[[k]])%*%as.matrix(rnorm(length(X_index_aux[[k]])))
    beta[beta_index_0]<-base::rep(0,length(beta_index_0))
    beta.chain[i,] <- beta
    if(i>burnin){model_chain[i -burnin,]<-gamma}

    ## Updating w
    X_beta <- as.matrix(X_gamma_aux[[k]])%*%as.matrix(beta[X_index_aux[[k]]])
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
    r<-rgamma( 1, shape= a0 +base::sum(l), rate= b0 +base::sum(base::log(1 +e_Xbeta)) )
    r.chain[i,]<-r

    ## Updating gamma
    gammaOld<-gamma

    ## Updating other quantities
    Z<-matrix((y-r)/(2*w))
    Z.chain[i,]<-Z
    eta.chain[i,]<-w
    Omega_modified_general<-matrix(base::rep(w,each = p), nrow = p, ncol = N, byrow = F)
    tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
    NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
    tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z


  }
  beta.chain<-beta.chain[(burnin +1):nchain,]

  if(first_excluded!=0)
  {
    prednames<-colnames(Covariates)[-(1:first_excluded)]
  }else{prednames<-colnames(Covariates)}
  t1<-proc.time()
  Time<-(t1-t0)[3]

  Output<-list(tau2=tau2, y=y, Covariates=Covariates, Z_chain=Z.chain, eta_chain=eta.chain, beta_chain=beta.chain, r_chain=r.chain, model_chain=model_chain, Default=list(Regression="NBR", prednames=prednames, Seconds=Time))
  class(Output) <- 'BMS'
  Output
}


gibbs_QR<-function(y, Covariates, first_excluded=0, nchain=1000, burnin=0, alpha=0.5, tau2=1000, rho=1, a0=1, b0=1,
                   beta.ini=base::rep(1,length(Covariates[1,]) +1), invsigma2.ini=1, w.ini=base::rep(1,length(y)), count.iteration=TRUE )
{
	t0<-proc.time()
  	##############################
  	##   Initiating    ##
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

  	beta.chain<-matrix(beta.ini, nrow=nchain, ncol=p, byrow=TRUE)
  	beta<-beta.chain[1,]
  	invsigma2.chain<-matrix(c(invsigma2.ini), nrow=nchain, ncol=1, byrow=TRUE)
  	invsigma2<-invsigma2.chain[1,]
  	sigma2.chain<-matrix(0, nrow=nchain, ncol=1, byrow=TRUE)
  	w<-w.ini
	Z.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
	eta.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)

  	gammaOld<-base::rep(1, p_selection)

  	delta2<-2/(alpha*(1-alpha))
  	xi<-(1 -2*alpha)/(alpha*(1-alpha))
  	Z<-matrix(y -xi*w)

  	model_chain<-matrix(NA, nrow=(nchain- burnin), ncol=p_selection)

  	log_tau2<-base::log(tau2)
  	frac_tau2<-1/tau2
  	log_womack_prior<-base::log(womack(p_selection,rho))
  	InvDiag_tau2_general<-diag(frac_tau2, p)
  	Omega_modified_general<-matrix(base::rep(1/(delta2*w)*invsigma2, each = p), nrow = p, ncol = N, byrow = F)
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
  ##  Gibbs  ##
  ######################
  	for(i in 1:nchain)
  	{
    	if(count.iteration==TRUE){cat("  Iteracion", i, "de", nchain, "\r")}

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
       	(0.5)*t(m_aux[[k]])%*%tXgamma_Omega_aux_Z +log_womack_prior[q_aux +1] )
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


    	## Updating beta
    	beta_index_0<-setdiff(1:p, X_index_aux[[k]])
    	if(length(beta_index_0)==0){beta_index_0<-0}
    	#beta[X_index_aux[[k]]]<-m_aux[[k]] +chol(V_aux[[k]])%*%as.matrix(rnorm(length(X_index_aux[[k]])))
    	beta[X_index_aux[[k]]]<-t(mvtnorm::rmvnorm(n = 1, mean = m_aux[[k]], sigma = V_aux[[k]], method="chol"))
    	beta[beta_index_0]<-base::rep(0,length(beta_index_0))
    	beta.chain[i,] <- beta
    	if(i>burnin){model_chain[i -burnin,]<-gamma}

    	## Updating w
    	X_beta <- as.matrix(X_gamma_aux[[k]])%*%as.matrix(beta[X_index_aux[[k]]])
    	chi<-invsigma2*((y -X_beta)^2)/delta2; psi<- invsigma2*(xi^2 +2*delta2)/delta2
    	for(j in 1:length(y))
    	{
      	w[j]<-GIGrvg::rgig(1, lambda=1/2, chi=chi[j], psi=psi)
    	}

    	## Updating invsigma2=1/sigma2
    	a1<-a0 +3*N/2
    	b1<-b0 +base::sum( ((y -X_beta -xi*w)^2 +2*delta2*w^2)/(2*delta2*w) )
    	invsigma2<-rgamma( 1, shape= a1, rate= b1 )
    	invsigma2.chain[i,]<-invsigma2

    	## Updating  gamma
    	gammaOld<-gamma

    	## Updating other quantities
    	Z<-matrix(y -xi*w)
	Z.chain[i,]<-Z
	eta.chain[i,]<-1/(delta2*w)*invsigma2
    	Omega_modified_general<-matrix(base::rep(1/(delta2*w)*invsigma2, each = p), nrow = p, ncol = N, byrow = F)
    	tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
    	NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
    	tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z
  	}
  	beta.chain<-beta.chain[(burnin +1):nchain,]
  	sigma2.chain<-1/invsigma2.chain[(burnin +1):nchain,]

  	if(first_excluded!=0)
  	{
    	prednames<-colnames(Covariates)[-(1:first_excluded)]
  	}else{prednames<-colnames(Covariates)}
  	t1<-proc.time()
  	Time<-(t1-t0)[3]

  	Output<-list(tau2=tau2, y=y, Covariates=Covariates, Z_chain=Z.chain, eta_chain=eta.chain, beta_chain=beta.chain, sigma2_chain=sigma2.chain, model_chain=model_chain, Default=list(Regression="QR", prednames=prednames, Seconds=Time))
	class(Output) <- 'BMS'
	Output
}



## Gibbs for Skew Normal regression
gibbs_SNR<-function(y, Covariates, first_excluded=0, nchain=1000, burnin=0, tau2=1000, rho=1, d, b2=1,
                    beta.ini=base::rep(1,length(Covariates[1,]) +1), kappa.ini=1, invzeta2.ini=1, v.ini=1, w.ini=base::rep(1,length(y)), count.iteration=TRUE )
{
  t0<-proc.time()
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

  beta.chain<-matrix(beta.ini, nrow=nchain, ncol=p, byrow=TRUE)
  lambda.chain<-matrix(c(kappa.ini*sqrt(invzeta2.ini)), nrow=nchain, ncol=1, byrow=TRUE)
  sigma2.chain<-matrix(c(1/invzeta2.ini +kappa.ini^2), nrow=nchain, ncol=1, byrow=TRUE)
  Z.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  eta.chain<-matrix(NA, nrow=nchain, ncol=N, byrow=TRUE)
  beta<-beta.chain[1,]
  w<-w.ini
  lambda<-lambda.chain[1,]
  kappa<-kappa.ini
  invzeta2<-invzeta2.ini
  zeta2<-1/invzeta2
  v<-v.ini

  gammaOld<-base::rep(1, p_selection)

  Z<-matrix(y -kappa*w)

  model_chain<-matrix(NA, nrow=(nchain- burnin), ncol=p_selection)

  log_tau2<-base::log(tau2)
  frac_tau2<-1/tau2
  log_womack_prior<-base::log(womack(p_selection,rho))
  InvDiag_tau2_general<-diag(frac_tau2, p)
  Omega_modified_general<-matrix(base::rep(invzeta2, each = p), nrow = p, ncol = N, byrow = F)
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
                                          (0.5)*t(m_aux[[k]])%*%tXgamma_Omega_aux_Z +log_womack_prior[q_aux +1] )
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


    ## Updating beta
    beta_index_0<-setdiff(1:p, X_index_aux[[k]])
    if(length(beta_index_0)==0){beta_index_0<-0}
    beta[X_index_aux[[k]]]<-m_aux[[k]] +chol(V_aux[[k]])%*%as.matrix(rnorm(length(X_index_aux[[k]])))
    beta[beta_index_0]<-base::rep(0,length(beta_index_0))
    beta.chain[i,] <- beta
    if(i>burnin){model_chain[i -burnin,]<-gamma}

    ## Updating w
    X_beta <- as.matrix(X_gamma_aux[[k]])%*%as.matrix(beta[X_index_aux[[k]]])
    mu_w<-(y*kappa -X_beta*kappa)/(kappa^2 +zeta2); sigma2_w<- (zeta2)/(kappa^2 +zeta2 )
    w<-truncnorm::rtruncnorm(n=N, a=0, b=Inf, mean=mu_w, sd=sqrt(sigma2_w))

    ## Updating invzeta2=1/zeta2
    aux_zeta1<-0.5*(N+1) +1
    aux_zeta2<-0.5*( (v*kappa^2/b2) +base::sum( (y -X_beta -kappa*w)^2 ) )
    invzeta2<-rgamma( 1, shape= aux_zeta1, rate= aux_zeta2 )
    zeta2<-1/invzeta2

    ## Updating kappa
    mu_kappa<-base::sum(y*w -X_beta*w)/(v/b2 +base::sum(w^2))
    sigma2_kappa<-(zeta2)/(v/b2 +base::sum(w^2))
    kappa<-rnorm(1, mean=mu_kappa, sd=sqrt(sigma2_kappa))

    ## Updating v
    shape_v<- (d+1)/2
    rate_v<- 0.5*( invzeta2*(kappa^2)/b2 +d )
    v<-rgamma(1, shape=shape_v, rate=rate_v)

    ## Updating gamma
    gammaOld<-gamma

    ## Updating other quantities
    Z<-matrix(y -kappa*w)
    Z.chain[i,]<-Z
    eta.chain[i,]<-invzeta2

    Omega_modified_general<-matrix(base::rep(invzeta2, each = p), nrow = p, ncol = N, byrow = F)
    tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
    NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
    tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z

    ## Saving original SN parameters
    lambda<-kappa/sqrt(zeta2)
    lambda.chain[i,]<-lambda
    sigma2<-zeta2 +kappa^2
    sigma2.chain[i,]<-sigma2
  }
  beta.chain<-beta.chain[(burnin +1):nchain,]
  lambda.chain<-lambda.chain[(burnin +1):nchain,]
  sigma2.chain<-sigma2.chain[(burnin +1):nchain,]

  if(first_excluded!=0)
  {
    prednames<-colnames(Covariates)[-(1:first_excluded)]
  }else{prednames<-colnames(Covariates)}
  t1<-proc.time()
  Time<-(t1-t0)[3]

 Output<- list(tau2=tau2, y=y, Covariates=Covariates, Z_chain=Z.chain, eta_chain=eta.chain, beta_chain=beta.chain, sigma2_chain=sigma2.chain, lambda_chain=lambda.chain, model_chain=model_chain, Default=list(Regression="SNR", prednames=prednames, Seconds=Time))
 class(Output) <- 'BMS'
 Output
}



## Summary table for Gibbs sampler
summary_gibbs<-function(fit, BF=FALSE)
{

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
			Z<-fit$Z_chain[i,]; eta<-fit$eta_chain[i,]
			Omega_modified_general<-matrix(base::rep(eta,each = p), nrow = p, ncol = N, byrow = F)
			tXgamma_Omega_general<-t(X)*Omega_modified_general[1:ncol(X),]
			NotInverse_V_general<-tXgamma_Omega_general%*%X + InvDiag_tau2_general
			tXgamma_Omega_Z_general<- tXgamma_Omega_general%*%Z
			for(j in 1:nrow(aux_ExploredModels))
			{
				#tXgamma_Omega_aux<- tXgamma_Omega_general[X_index_aux[[j]],]
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
			Conditional_BF<-mean(exp(loglik[[1]] -loglik[[j]]))		
			Marginal_BF<-exp( mean(loglik[[1]]) -mean(loglik[[j]]) )
			New_ExploredModels[j,ncol(ExploredModels) +1]<-Conditional_BF
			New_ExploredModels[j,ncol(ExploredModels) +2]<-Marginal_BF
		}
		colnames(New_ExploredModels)[-(1:ncol(ExploredModels))]<-c("Conditional_BF", "Marginal_BF_Estimator")
		return(New_ExploredModels)
	}

  ExploredModels<-aux_summary_gibbs(fit)
  if(BF==TRUE)
  { 
	ExploredModels<-aux_BF(fit,ExploredModels) 
	most_model<-ExploredModels[1,-((ncol(ExploredModels) -2):ncol(ExploredModels)) ]
  }else{most_model<-ExploredModels[1,-ncol(ExploredModels)]}

  beta_index<-c()
  for(j in 1:nrow(fit$beta_chain))
  {
    if(prod(most_model==fit$model_chain[j,])==1){beta_index<-c(beta_index,j)}
  }
  beta_chain<-fit$beta_chain[beta_index,]
  if(fit$Default$Regression=="LiR"){parameters<-cbind(fit$sigma2_chain[beta_index], beta_chain); aux_rownames<-c("sigma2")}
  if(fit$Default$Regression=="LoR"){parameters<-cbind(beta_chain); aux_rownames<-c()}
  if(fit$Default$Regression=="NBR"){parameters<-cbind(fit$r_chain[beta_index], beta_chain); aux_rownames<-c("r")}
  if(fit$Default$Regression=="QR"){parameters<-cbind(fit$sigma2_chain[beta_index], beta_chain); aux_rownames<-c("sigma2")}
  if(fit$Default$Regression=="SNR"){parameters<-cbind(fit$sigma2_chain[beta_index], fit$lambda_chain[beta_index], beta_chain); aux_rownames<-c("sigma2", "lambda")}

  p<-ncol(beta_chain)
  Mean<-apply(parameters, 2, mean)
  Quantile<-t(apply(parameters, 2, quantile, prob=c(0.025,0.975)))
  SD<-t(apply(parameters, 2, sd))
  Table<-as.data.frame(matrix(c(Mean,Quantile, SD), byrow=FALSE, ncol=4))
  colnames(Table)<-c("Mean", "2.5% quantile", "97.5% quantile", "SD")
  rownames(Table)<-c(aux_rownames, "intercept", fit$Default$prednames)


  list(Mean_IC=Table, Explored_Models=ExploredModels)
}


## Step Criteria manual (only backward)
MyStepCriteria<-function(Regression, PredictorsNames, base, Criteria="AIC", ... )
{
	N<-nrow(base)
	r_p<-ncol(base)
	p<-r_p
	aux_PredictorsNames<-PredictorsNames
	flag<-0; j<-1
	flag_intercept<-0
	SelectedCriteria<-c()
	while(flag==0)
	{	
		lengthPredictors<-length(aux_PredictorsNames)
		aux_formula<-list(); Candidates_PredictorsNames<-list()
		for(i in 1:lengthPredictors)
		{
			Candidates_PredictorsNames[[i]]<-aux_PredictorsNames[-(lengthPredictors +1 -i)]
			if(Regression=="Binomial"){aux_formula[[i]]<-paste0("y/ni ~ ", paste0(Candidates_PredictorsNames[[i]], collapse = "+"))	}else{
				aux_formula[[i]]<-paste0("y ~ ", paste0(Candidates_PredictorsNames[[i]], collapse = "+"))	}
		}
		p<-p-1
		if(p==1){aux_formula<-paste0("y ~ 1")}	#Si p=1 => s�lo tenemos intercepto
		
		if(Regression=="Normal"){ aux_mods<- lapply(aux_formula, function(frml) glm(frml, data = base, family="gaussian")) }
		if(Regression=="Binomial"){ aux_mods<- lapply(aux_formula, function(frml) glm(frml, data = base, family = "binomial", weights = ni)) }
		if(Regression=="NegBinomial"){ aux_mods<- lapply(aux_formula, function(frml) glm.nb(frml, data=base, init.theta=0.5) ) }
		if(Regression=="Quantile"){ aux_mods<- lapply(aux_formula, function(frml) rq(frml, tau=r_alpha, data=base)) }
		if(Regression=="SkewNormal"){ aux_mods<- lapply(aux_formula, function(frml) selm(frml, data = base, family="SN")) }
		aux_Criteria<-unlist(lapply(aux_mods, Criteria_lm, N=N, p=p, Criteria=Criteria))
		aux_index<-which(min(aux_Criteria)==aux_Criteria)
	
		if(j>=2)
		{
			if(SelectedCriteria[j-1]<aux_Criteria[aux_index]){flag<-1; SelectedPredictors<-aux_PredictorsNames}else{
			SelectedCriteria[j]<-aux_Criteria[aux_index]; j<-j+1; old_aux_mods<-aux_mods; old_aux_index<-aux_index
			aux_PredictorsNames<-Candidates_PredictorsNames[[aux_index]] }
		}else{
			SelectedCriteria[j]<-aux_Criteria[aux_index]; j<-j+1
			aux_PredictorsNames<-Candidates_PredictorsNames[[aux_index]]
			old_aux_mods<-aux_mods; old_aux_index<-aux_index
		}
		if(p==1 && flag==0){flag<-1; flag_intercept<-1}	
	}
	if(flag_intercept==0){ SelectedFormula<-paste0("y ~ ", paste0(SelectedPredictors, collapse = "+"))	 }else{
	SelectedFormula<-aux_formula }		
	list(SelectedFormula=SelectedFormula, fit=old_aux_mods[[old_aux_index]])
}


## Transforming a "binomal data base" into a "bernoulli data base"
DupBaseBinomial<-function(base)
{
	y_IndexGreater1<-which(base$y>1)
	
	y_Greater1<-base$y[y_IndexGreater1]
	duptimes<-rep(1,nrow(base))		
	duptimes[y_IndexGreater1]<-y_Greater1
	
	idx <- rep(1:nrow(base), duptimes)
	dupbase <- base[idx,]		
	dupbase$y[which(dupbase$y>1)]<-1		

	return(dupbase)	
}

## Logistic Regression Data generator
gen_base_binomial_reg<- function(beta, covariables, N, ni=rep(1, N))
{
	p<-length(covariables[1,]) +1; x<-c()
	for(i in 1:(p -1))
	{
			x<-c(x,covariables[,i])
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
                   covariables = covariables)  
	colnames(base)[-(1:3)]<-colnames(covariables)
	return(base)
}


## Negative Binomial Regression Data generator
gen_base_NegBinomial_reg<- function(N, beta, r, covariables)
{
	X<-as.matrix(covariables)
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
                   covariables = covariables)  
	colnames(base)[-1]<-colnames(covariables)
	return(base)
}



