rm(list=ls(all=TRUE))
source("functions.R")

#install.packages("BiocManager")
#BiocManager::install("sparseMatrixStats")

library(mombf)		#Rosell, LA method
library(MASS)		#glm.nb (Negative binomial regression), stepAIC()
library(quantreg)		#For Quantile regression
library(sn)		#Frequentist Regression Skew Normal
library(xtable)
############################################################################################################

##################################################################################################
##						 Functions, Just run this						##
##################################################################################################

mcr <- function(x, drop = FALSE)
  {
    xx <- do.call("paste", c(data.frame(x), sep = "\r"))
    tx <- table(xx)
    mx <- base::names(tx)[which(tx == base::max(tx))[1]]
    x[match(mx, xx), , drop = drop]
  }

Criteria_lm<-function(fit, N, p, Criteria)
{	
	if(Criteria=="AIC"){ result<-2*p -2*logLik(fit)[1] }	
	if(Criteria=="BIC"){ result<-p*log(N) -2*logLik(fit)[1] }
	return(result)
}
#Criteria_lm(fit, N, r_p, Criteria="AIC")

mms_apply<-function(matrix, statistic, ...)
{
	if( length(dim(matrix))==0 ){matrix}else{apply(matrix, 2, statistic) }
}


Summary_SimStudy<-function(aux, ExploredModels=TRUE)
{
	##########################
	##   Explored models    ##
	##########################
	r_p<-ncol(aux$Model_Exploration_OurMethod)
	Regression<-aux$Regression
	
	## Obtaining Explored Model Table from output
	ExploredMOurMethod<-data.frame(aux$Model_Exploration_OurMethod)
	ExploredMLA<-data.frame(aux$Model_Exploration_LA)	#DE AH� CAMBIAR ESTO, PUES YA FUE ARREGLADO EN EL CODIGO
	ExploredMStepAIC<-data.frame(aux$Model_Exploration_StepAIC)
	ExploredMStepBIC<-data.frame(aux$Model_Exploration_StepBIC)

	## Changing predictor colnames for Latex
	aux_namespred<-vector(length=r_p-1)
	for(i in 1:(r_p -1)){ aux_namespred[i]<-paste0("$x_",i,"$") }
	colnames(ExploredMOurMethod)[-r_p]<-aux_namespred
	colnames(ExploredMLA)[-r_p]<-aux_namespred
	colnames(ExploredMStepAIC)[-r_p]<-aux_namespred
	colnames(ExploredMStepBIC)[-r_p]<-aux_namespred
	
	## Sorting rows by frequency
	ExploredMOurMethod<- rbind(ExploredMOurMethod[1,],  ExploredMOurMethod[-1,][order(ExploredMOurMethod[-1,]$Frequency, decreasing=TRUE),] )
	ExploredMLA<- rbind(ExploredMLA[1,],  ExploredMLA[-1,][order(ExploredMLA[-1,]$Frequency, decreasing=TRUE),] )
	ExploredMStepAIC<- rbind(ExploredMStepAIC[1,],  ExploredMStepAIC[-1,][order(ExploredMStepAIC[-1,]$Frequency, decreasing=TRUE),] )
	ExploredMStepBIC<- rbind(ExploredMStepBIC[1,],  ExploredMStepBIC[-1,][order(ExploredMStepBIC[-1,]$Frequency, decreasing=TRUE),] )
	
	PercentagesOurMethod<-ExploredMOurMethod$Frequency
	PercentagesLA<-ExploredMLA$Frequency
	PercentagesStepAIC<-ExploredMStepAIC$Frequency
	PercentagesStepBIC<-ExploredMStepBIC$Frequency
	
	## Changing frequency values to percentages in text, for latex
	aux_EM2<-vector(length=nrow(ExploredMOurMethod))
	for(i in 1:nrow(ExploredMOurMethod)){ aux_EM2[i]<-paste0(ExploredMOurMethod$Frequency[i], "%") }
	ExploredMOurMethod$Frequency<-aux_EM2
	aux_EM2<-vector(length=nrow(ExploredMLA))
	for(i in 1:nrow(ExploredMLA)){ aux_EM2[i]<-paste0(ExploredMLA$Frequency[i], "%") }
	ExploredMLA$Frequency<-aux_EM2
	aux_EM2<-vector(length=nrow(ExploredMStepAIC))
	for(i in 1:nrow(ExploredMStepAIC)){ aux_EM2[i]<-paste0(ExploredMStepAIC$Frequency[i], "%") }
	ExploredMStepAIC$Frequency<-aux_EM2
	aux_EM2<-vector(length=nrow(ExploredMStepBIC))
	for(i in 1:nrow(ExploredMStepBIC)){ aux_EM2[i]<-paste0(ExploredMStepBIC$Frequency[i], "%") }
	ExploredMStepBIC$Frequency<-aux_EM2

	## Adding the column where will go the "True model" text, for latex
	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMOurMethod) -1))
	ExploredMOurMethod<- head( cbind(" "=aux_EM, ExploredMOurMethod), 4 )

	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMLA) -1))
	ExploredMLA<- head( cbind(" "=aux_EM, ExploredMLA), 4 )
	
	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMStepAIC) -1))
	ExploredMStepAIC<- head( cbind(" "=aux_EM, ExploredMStepAIC), 4 )
	
	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMStepBIC) -1))
	ExploredMStepBIC<-head( cbind(" "=aux_EM, ExploredMStepBIC), 4 )
	
	## Latex All in one table
	if(Regression=="Normal" || Regression=="Binomial" || Regression=="Quantile")
	{
		ExploredModel<-list(Our_Method=ExploredMOurMethod, LA=ExploredMLA, Step_AIC=ExploredMStepAIC, Step_BIC=ExploredMStepBIC)
	}else{
		ExploredModel<-list(Our_Method=ExploredMOurMethod, Step_AIC=ExploredMStepAIC, Step_BIC=ExploredMStepBIC)
	}
	
	
	#############################
	## 	  Estimation    	 ##
	#############################
	maxOurMethod<-max(PercentagesOurMethod)
	maxLA<-max(PercentagesLA)
	maxStepAIC<-max(PercentagesStepAIC)
	maxStepBIC<-max(PercentagesStepBIC)
	
	flag_OurMethod<-0; flag_LA<-0; flag_StepAIC<-0; flag_StepBIC<-0
	if(length(which(maxOurMethod==PercentagesOurMethod))==1){flag_OurMethod<-1}
	if(length(which(maxLA==PercentagesLA))==1){flag_LA<-1}
	if(length(which(maxStepAIC==PercentagesStepAIC))==1){flag_StepAIC<-1}
	if(length(which(maxStepBIC==PercentagesStepBIC))==1){flag_StepBIC<-1}
	
	## Estimates
	EstimatesOurMethod<-aux$OurMethod_Posterior_Summary[,c(1,3)]
	rep_blank<-nrow(EstimatesOurMethod)
	if(flag_OurMethod==0){EstimatesOurMethod$mean<-rep("- - -", rep_blank); EstimatesOurMethod$sd<-rep("- - -", rep_blank)}
	
	EstimatesLA<-aux$LA_Posterior_Summary[,c(1,3)] 
	if(flag_LA==0){EstimatesLA$mean<-rep("- - -", rep_blank); EstimatesLA$sd<-rep("- - -", rep_blank)}
	
	if(Regression=="Quantile"){EstimatesStepAIC<-rbind(sigma2=c(NaN,NaN),aux$StepAIC_Summary[,c(1,3)])}else{
	EstimatesStepAIC<-aux$StepAIC_Summary[,c(1,3)]
	}
	if(flag_StepAIC==0){EstimatesStepAIC$mean<-rep("- - -", rep_blank); EstimatesStepAIC$sd<-rep("- - -", rep_blank)}
	 
	if(Regression=="Quantile"){EstimatesStepBIC<-rbind(sigma2=c(NaN,NaN),aux$StepBIC_Summary[,c(1,3)])}else{
	EstimatesStepBIC<-aux$StepBIC_Summary[,c(1,3)] 
	}
	if(flag_StepBIC==0){EstimatesStepBIC$mean<-rep("- - -", rep_blank); EstimatesStepBIC$sd<-rep("- - -", rep_blank)}
	
	aux_names<-vector(length=r_p-1)
	for(i in 0:(r_p -1)){ aux_names[i+1]<-paste0("beta_{",i,"}") }
	if(Regression=="Normal"){parameters<-c(r_sigma2, r_beta); row_names<-c("sigma2", aux_names)}
	if(Regression=="Binomial"){parameters<-c(r_beta); row_names<-c(aux_names)}
	if(Regression=="NegBinomial"){parameters<-c(r_r, r_beta); row_names<-c("r", aux_names)}
	if(Regression=="Quantile"){parameters<-c(r_sigma2, r_beta); row_names<-c("sigma2", aux_names)}
	if(Regression=="SkewNormal"){parameters<-c(r_lambda, r_sigma2, r_beta); row_names<-c("lambda", "sigma2", aux_names)}
	if(Regression=="Normal" || Regression=="Binomial")
	{
		Estimates<-data.frame(parameters, EstimatesOurMethod, EstimatesLA, EstimatesStepAIC, EstimatesStepBIC)
		rownames(Estimates)<-row_names
		colnames(Estimates)[-1]<-c("mean_OurMethod", "sd_OurMethod", "mean_LA", "sd_LA", "mean_StepAIC", "sd_StepAIC" ,"mean_StepBIC", "sd_StepBIC")
	}else{
		Estimates<-data.frame(parameters, EstimatesOurMethod, EstimatesStepAIC, EstimatesStepBIC)
		rownames(Estimates)<-row_names
		colnames(Estimates)[-1]<-c("mean_OurMethod", "sd_OurMethod", "mean_StepAIC", "sd_StepAIC" ,"mean_StepBIC", "sd_StepBIC")
	}
	if(ExploredModels==TRUE){OutPut<-list(ExploredModel=ExploredModel, Estimates=Estimates)}else{
	OutPut<-list(Estimates=Estimates)}
	OutPut
}
	


General_Sim_OurMethod<-function(N, r_beta, R, nchain=1000, burnin=0, Regression="Normal", ...)
{
t0<-proc.time()
	first_excluded=0
	p<-length(r_beta)
	intercept_first_excluded<- first_excluded +1	#Aqu� estoy excluyendo el intercepto del proceso de selecci�n
	p_selection<-p -intercept_first_excluded		#n�mero de coeficientes a testear
	DataList<-list()					#Lista que en cada entrada guarda "base" de cada r�plica, esto es, matriz con respuesta y covariables
	RealModel<-+(r_beta!=0)[-1]	#Real model. Erasing intercept, because it's out of the selection process
	PredictorsIndex1_Real<-which(RealModel==1)	#Which beta!=0
	PredictorsNames<-colnames(data.frame(matrix(ncol=p-1)))	#Predictor's names. In "X1,X2,..." format


	## Preparation for Our method
	SelectedModelsCountsOurMethod_rbind<-rbind( c(RealModel,0) )	#Counting selected models in replicas
	colnames(SelectedModelsCountsOurMethod_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsOurMethod_list<-NULL
	CoefsMeanOurMethod_matrix<-matrix(0,ncol=p, nrow=R) #Matriz de la media a posteriori de beta en cada R�plica 
	CoefsSDOurMethod_matrix<-matrix(0,ncol=p, nrow=R) #Matriz de la sd a posteriori de beta en cada R�plica 
	if(Regression=="Normal"){sigma2_mean_vector<-vector(length=R); sigma2_sd_vector<-vector(length=R) }
	if(Regression=="NegBinomial"){r_mean_vector<-vector(length=R); r_sd_vector<-vector(length=R) }
	if(Regression=="Quantile"){sigma2_mean_vector<-vector(length=R); sigma2_sd_vector<-vector(length=R) }
	if(Regression=="SkewNormal"){lambda_mean_vector<-vector(length=R); lambda_sd_vector<-vector(length=R); sigma2_mean_vector<-vector(length=R); sigma2_sd_vector<-vector(length=R) }	
	time_OurMethod<-0

	## Preparation for Rusell LA method
	SelectedModelsCountsLA_rbind<-rbind( c(RealModel,0) )	#Counting selected models in replicas
	colnames(SelectedModelsCountsLA_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsLA_list<-NULL
	CoefsMeanLA_matrix<-matrix(0,ncol=p, nrow=R) 

	if(Regression=="Normal"){Sigma2MeanLA_vector<-vector(length=R) }
	time_LA<-0
	
	## Preparation for StepAIC method
	glmStepAICCoefsMatrix<-matrix(0,ncol=p, nrow=R)  
	glmStepAICSdMatrix<-matrix(0,ncol=p, nrow=R)  
	SelectedModelsCountsStepAIC_rbind<-rbind( c(RealModel,0) )
	colnames(SelectedModelsCountsStepAIC_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsStepAIC_list<-NULL
	if(Regression=="Normal"){glmStepAIC_sigma2EstVector<-vector(length=R); glmStepAIC_sigma2SdVector<-vector(length=R)}
	if(Regression=="NegBinomial"){glmStepAIC_rEstVector<-vector(length=R); glmStepAIC_rSdVector<-vector(length=R)}
	if(Regression=="SkewNormal"){glmStepAIC_lambdaEstVector<-vector(length=R); glmStepAIC_lambdaSdVector<-vector(length=R); glmStepAIC_sigma2EstVector<-vector(length=R); glmStepAIC_sigma2SdVector<-vector(length=R)}
	time_StepAIC<-0

	## Preparation for StepBIC method
	glmStepBICCoefsMatrix<-matrix(0,ncol=p, nrow=R)  
	glmStepBICSdMatrix<-matrix(0,ncol=p, nrow=R)  
	SelectedModelsCountsStepBIC_rbind<-rbind( c(RealModel,0) )	
	colnames(SelectedModelsCountsStepBIC_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsStepBIC_list<-NULL
	if(Regression=="Normal"){glmStepBIC_sigma2EstVector<-vector(length=R); glmStepBIC_sigma2SdVector<-vector(length=R)}
	if(Regression=="NegBinomial"){glmStepBIC_rEstVector<-vector(length=R); glmStepBIC_rSdVector<-vector(length=R)}
	if(Regression=="SkewNormal"){glmStepBIC_lambdaEstVector<-vector(length=R); glmStepBIC_lambdaSdVector<-vector(length=R); glmStepBIC_sigma2EstVector<-vector(length=R); glmStepBIC_sigma2SdVector<-vector(length=R)}
	time_StepBIC<-0
	
	y<-vector(length=N)
	for(k in 1:R)
	{
		#cat("  Replica", k, " ;", " Iteracion cadena", i, "de", nchain, "\r")
		cat("  Replica", k, "\r")

		## Pre-simulation stuff
		aux_cov<-rnorm((p-1)*N, 0, 1)			

		if(Regression=="Normal")
		{
			X<-matrix( c(rep(1, N), rnorm((p -1)*N)), ncol=p )
			Xbeta<-X%*%r_beta
			y<-rnorm(N, mean=Xbeta, sd=sqrt(r_sigma2))
			covariables<-X[,2:(length(r_beta))]
			base<-data.frame(y, covariables)

			## For LA method
			familyReg<-"normal"	
			X_LA<-as.matrix(X); y_LA<-y
		}

		if(Regression=="Binomial")
		{							
			covariables<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N))
			base<-gen_base_binomial_reg(r_beta, covariables, N, ni=ni)
			base<-base[,-2]
			y<-base$y
			X<-data.frame("Intercept"=rep(1,N), covariables)

			## For LA method
			familyReg<-"binomial"
			NewBase<-DupBaseBinomial(base)
			X_LA<-as.matrix(data.frame(Intercept=rep(1, nrow(NewBase)), X=NewBase[,c(-1,-2)] ))
			y_LA<-NewBase$y

		}

		if(Regression=="NegBinomial")
		{
			Check<-FALSE
			while(Check==FALSE)	
			{
				covariables<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N))		
				base<-gen_base_NegBinomial_reg(N, r_beta, r_r, covariables)
				y<-base$y
			
				t0_StepAIC<-proc.time()
				fit_AIC <- try(MyStepCriteria(Regression, PredictorsNames, base, Criteria="AIC"))
				t1_StepAIC<-proc.time()
				time_StepAIC<-time_StepAIC +(t1_StepAIC-t0_StepAIC)[3]

				t0_StepBIC<-proc.time()
				fit_BIC<-try( MyStepCriteria(Regression, PredictorsNames, base, Criteria="BIC") )	
				t1_StepBIC<-proc.time()
				time_StepBIC<-time_StepBIC +(t1_StepBIC-t0_StepBIC)[3]

				if( grepl("Error",fit_AIC)[1]==FALSE && grepl("Error",fit_BIC)[1]==FALSE ){Check<-TRUE}
			}
			X<-data.frame("Intercept"=rep(1,N),covariables)
			fit_AIC<-fit_AIC$fit
			fit_BIC<-fit_BIC$fit
		}

		if(Regression=="Quantile")
		{
			X<-matrix( c(rep(1, N), aux_cov), ncol=p )
			Xbeta<-X%*%r_beta
			w<-rexp(N, rate=1/r_sigma2)
			varphi<-(1-2*r_alpha)/(r_alpha*(1-r_alpha))
			delta2<-2/(r_alpha*(1-r_alpha))
			y<- Xbeta +varphi*w +sqrt(r_sigma2*delta2*w)*rnorm(N, mean=0, sd=1)
			covariables<-X[,2:(length(r_beta))]
			base<-data.frame(y, covariables)

			## For LA method
			familyReg<-"twopiecelaplace"
			X_LA<-as.matrix(X); y_LA<-y
						
		}

		if(Regression=="SkewNormal")
		{
			X<-matrix( c(rep(1, N), rnorm((p -1)*N)), ncol=p )
			Xbeta<-X%*%r_beta
			r_w<-abs(rnorm(N, mean=0, sd=1))
			r_delta<-r_lambda/sqrt(r_lambda^2 +1)
			y<-Xbeta +sqrt(r_sigma2)*r_delta*r_w +sqrt(r_sigma2*(1 -r_delta^2))*rnorm(N, mean=0, sd=1)
			covariables<-X[,2:(length(r_beta))]
			base<-data.frame(y, covariables)		
		}
		X<-as.matrix(X)
		
		DataList[[k]]<-base



		######################################
		##		Our method		##
		######################################
		t0_OurMethod<-proc.time()
		## Gibbs
		if(Regression=="Normal")
		{
			fit_OurMethod<-gibbs_LiR(y, covariables, first_excluded=0, nchain, burnin, tau2=1000, rho=1, a0=1, b0=1,
			beta.ini=rep(1,length(covariables[1,]) +1), invsigma2.ini=1, count.iteration=FALSE )
		}
		if(Regression=="Binomial")
		{
			fit_OurMethod<-gibbs_LoR(y, ni, covariables, first_excluded, nchain=nchain, burnin=burnin, tau2=1000, rho=1,
			beta.ini=rep(1,length(covariables[1,]) +1), w.ini=rep(1,length(y)), count.iteration=FALSE )
		}

		if(Regression=="NegBinomial")
		{
			fit_OurMethod<-gibbs_NBR(y, covariables, first_excluded=0, nchain, burnin, tau2=1000, rho=1, a0=1, b0=1,
			beta.ini=rep(1,length(covariables[1,]) +1), r.ini=1, w.ini=rep(1,length(y)), l.ini=rep(1,length(y)), count.iteration=FALSE )
		}

		if(Regression=="Quantile")
		{
			fit_OurMethod<-gibbs_QR(y, covariables, first_excluded=0, nchain=nchain, burnin=burnin, alpha=r_alpha, tau2=1000, rho=1, a0=1, b0=1,
			beta.ini=rep(1,length(covariables[1,]) +1), invsigma2.ini=1, w.ini=rep(1,length(y)), count.iteration=FALSE )	
		}
		
		if(Regression=="SkewNormal")
		{
			fit_OurMethod<-gibbs_SNR(y, covariables, first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, d=d, b2=b2,
			beta.ini=rep(1,length(covariables[1,]) +1), kappa.ini=1, invzeta2.ini=1, v.ini=1, w.ini=rep(1,length(y)), count.iteration=FALSE )
		}
		t1_OurMethod<-proc.time()
		time_OurMethod<-time_OurMethod +(t1_OurMethod -t0_OurMethod)[3]

		
		## Adding Selected model to the count
		SelectedModelOurMethod<-as.numeric(mcr(fit_OurMethod$model_chain))
		flag<-0; j<-1
		while(flag==0)
		{
			if( all(SelectedModelsCountsOurMethod_rbind[j,-p]==SelectedModelOurMethod) )
			{
				SelectedModelsCountsOurMethod_rbind[j,p]<- SelectedModelsCountsOurMethod_rbind[j,p]+1
				IndexSelectedModelsOurMethod_list[[j]]<-c(IndexSelectedModelsOurMethod_list[[j]], k)
				flag<-1
			}else{j<-j+1}
			if(j>nrow(SelectedModelsCountsOurMethod_rbind))
			{
				SelectedModelsCountsOurMethod_rbind<-rbind( SelectedModelsCountsOurMethod_rbind, c(SelectedModelOurMethod, 1) )
				IndexSelectedModelsOurMethod_list[[j]]<-k
				flag<-2
			}
		}

				
		## Posterior Mean and SD given the most selected model
		index_selected<-which(colSums(t(fit_OurMethod$model_chain) == SelectedModelOurMethod) == ncol(fit_OurMethod$model_chain))
		coefs_chain<-fit_OurMethod$beta_chain[index_selected,]	
		CoefsMeanOurMethod_matrix[k,]<-apply(coefs_chain, 2, mean)
		CoefsSDOurMethod_matrix[k,]<-apply(coefs_chain, 2, sd)

		if(Regression=="Normal" || Regression=="Quantile" )
		{
			sigma2_mean_vector[k]<-mean(fit_OurMethod$sigma2_chain[index_selected])
			sigma2_sd_vector[k]<-sd(fit_OurMethod$sigma2_chain[index_selected])
		}
		if(Regression=="NegBinomial")
		{
			r_mean_vector[k]<-mean(fit_OurMethod$r_chain[index_selected])
			r_sd_vector[k]<-sd(fit_OurMethod$r_chain[index_selected])
		}
		if(Regression=="SkewNormal")
		{
			lambda_mean_vector[k]<-mean(fit_OurMethod$lambda_chain[index_selected])
			lambda_sd_vector[k]<-sd(fit_OurMethod$lambda_chain[index_selected])
			sigma2_mean_vector[k]<-mean(fit_OurMethod$sigma2_chain[index_selected])
			sigma2_sd_vector[k]<-sd(fit_OurMethod$sigma2_chain[index_selected])
		}


		
		######################################
		## 	 	 Rusell LA 		##
		######################################
		t0_LA<-proc.time()
		if( Regression=="Normal" || Regression=="Binomial" || Regression=="Quantile" )
		{

		prCoef <- momprior(tau=1)  
		prDelta <- modelbbprior(1,1)  

		## Fitting the model
		if(Regression=="Quantile")
		{
			priorSkew<- 2*r_alpha -1	#Fixing quantile
			fit_LA <- modelSelection(y=y_LA, x=X_LA, family=familyReg, priorCoef=prCoef, priorDelta=prDelta, priorSkew=priorSkew, method="Laplace",
			verbose==FALSE, center=FALSE, scale=FALSE, niter=nchain, burnin=burnin, includevars=1)
		}else{	
			fit_LA <- modelSelection(y=y_LA, x=X_LA, family=familyReg, priorCoef=prCoef, priorDelta=prDelta, method="Laplace",
			verbose==FALSE, center=FALSE, scale=FALSE, niter=nchain, burnin=burnin, includevars=1)
		}
		t1_LA<-proc.time()
		time_LA<-time_LA +(t1_LA-t0_LA)[3]

		## Coefficients estimates
		if(Regression=="Normal"){ CoefsMeanLA_matrix[k,]<-as.vector(coef(fit_LA)[c(-1, -(p+2)),1]) }
		if(Regression=="Binomial"){ CoefsMeanLA_matrix[k,]<-as.vector(coef(fit_LA)[-1,1]) }

		if(Regression=="Normal" )
		{
			Sigma2MeanLA_vector[k]<-as.vector(coef(fit_LA)[p+1,1])
		}

		SelectedModelLA<-as.vector(fit_LA$postMode)[-1]		
		flag<-0; j<-1
		while(flag==0)
		{
			if( all(SelectedModelsCountsLA_rbind[j,-p]==SelectedModelLA) )
			{
				SelectedModelsCountsLA_rbind[j,p]<- SelectedModelsCountsLA_rbind[j,p]+1
				IndexSelectedModelsLA_list[[j]]<-c(IndexSelectedModelsLA_list[[j]], k)
				flag<-1
			}else{j<-j+1}
			if(j>nrow(SelectedModelsCountsLA_rbind))
			{
				SelectedModelsCountsLA_rbind<-rbind( SelectedModelsCountsLA_rbind, c(SelectedModelLA, 1) )
				IndexSelectedModelsLA_list[[j]]<-k
				flag<-2
			}
		}


		}



		######################################
		## 	 	 StepAIC glm 		##
		######################################
		
		if(Regression!="NegBinomial"){t0_StepAIC<-proc.time()}
		if(Regression=="Normal")
		{
			fit<-glm(y~., data = base, family="gaussian")	
			fit_AIC<-stepAIC(fit, trace=FALSE)
		}

		if(Regression=="Binomial")
		{
			fit<-glm(y/ni ~., data = base, family = "binomial", weights = ni)			
			fit_AIC<-stepAIC(fit, trace=FALSE)										
		}

		if(Regression=="Quantile")
		{
			fit<-rq(y~., tau=r_alpha, data=base)			
			fit_AIC<-stepAIC(fit, trace=FALSE)			
		}
		if(Regression=="SkewNormal")
		{
			fit_AIC<-MyStepCriteria(Regression, PredictorsNames, base, Criteria="AIC" )$fit
		}
		if(Regression!="NegBinomial"){t1_StepAIC<-proc.time(); time_StepAIC<-time_StepAIC +(t1_StepAIC-t0_StepAIC)[3]}
		

		if(Regression=="SkewNormal")
		{ 
			coef_stepAIC<-coef(fit_AIC, "DP")
			PredNamesStepAIC<-head(names(coef_stepAIC)[c(-1)],-2)	
		}
		if(Regression!="SkewNormal")
		{
			coef_stepAIC<-fit_AIC$coefficients								
			PredNamesStepAIC<-names(fit_AIC$coefficients)[-1]
		}


		## Adding Selected model to the count
		IndexSelectedStepAIC<-match(PredNamesStepAIC, PredictorsNames)
		SelectedModelStepAIC<-rep(0, p-1)
		SelectedModelStepAIC[IndexSelectedStepAIC]<-1

		flag<-0; j<-1
		while(flag==0)
		{
			if( all(SelectedModelsCountsStepAIC_rbind[j,-p]==SelectedModelStepAIC) )
			{
				SelectedModelsCountsStepAIC_rbind[j,p]<- SelectedModelsCountsStepAIC_rbind[j,p]+1
				IndexSelectedModelsStepAIC_list[[j]]<-c(IndexSelectedModelsStepAIC_list[[j]], k)
				flag<-1
			}else{j<-j+1}
			if(j>nrow(SelectedModelsCountsStepAIC_rbind))
			{
				SelectedModelsCountsStepAIC_rbind<-rbind( SelectedModelsCountsStepAIC_rbind, c(SelectedModelStepAIC, 1) )
				IndexSelectedModelsStepAIC_list[[j]]<-k
				flag<-2
			}
		}
			
		aux_betaIndex_1<-which(c(1,SelectedModelStepAIC)!=0)	#Indexes of beta's such that beta_j!=0 (including intercept)	
		

		## Coefficient estimated
		aux_betaStepAIC<-rep(0,p)
		if(Regression!="SkewNormal"){aux_betaStepAIC[aux_betaIndex_1]<-as.vector(coefficients(fit_AIC)) }		
		if(Regression=="SkewNormal"){ aux_betaStepAIC[aux_betaIndex_1]<-as.vector(coefficients(fit_AIC, "DP"))[1:length(aux_betaIndex_1)] }
		glmStepAICCoefsMatrix[k,]<-aux_betaStepAIC			
	
		## Other parameters estimates with their SD (Normal, NegBinomial, Skewnormal)
		if(Regression=="NegBinomial")
		{
			glmStepAIC_rEstVector[k]<-summary(fit_AIC)$theta
			glmStepAIC_rSdVector[k]<-summary(fit_AIC)$SE.theta
		}
		if(Regression=="SkewNormal")
		{
			glmStepAIC_lambdaEstVector[k]<-as.vector(coefficients(fit_AIC, "DP"))[length(aux_betaIndex_1) +2]
			glmStepAIC_lambdaSdVector[k]<-sqrt( vcov(fit_AIC, "DP")[length(aux_betaIndex_1) +2, length(aux_betaIndex_1)+2] )
			glmStepAIC_sigma2EstVector[k]<-as.vector(coefficients(fit_AIC, "DP"))[length(aux_betaIndex_1) +1]
			glmStepAIC_sigma2SdVector[k]<-sqrt( vcov(fit_AIC, "DP")[length(aux_betaIndex_1) +1, length(aux_betaIndex_1)+1] )
		}
		if(Regression=="Normal")
		{
			glmStepAIC_sigma2EstVector[k]<-as.vector(summary(fit_AIC)$dispersion)
			p_glmStepAIC<-length(aux_betaIndex_1)
			glmStepAIC_sigma2SdVector[k]<-sqrt(2*(glmStepAIC_sigma2EstVector[k]^2)/(N-p_glmStepAIC))
		}

		## Coefficient's SD estimates
		if(Regression!="Quantile" && Regression!="SkewNormal")	#For Normal, Binomial, and NegBinomial regression
		{
			aux_SdGlmStepAIC<-rep(0,p)
			aux_SdGlmStepAIC[aux_betaIndex_1]<-sqrt(diag(vcov(fit_AIC)))[1:length(aux_betaIndex_1)]	
			glmStepAICSdMatrix[k,]<-aux_SdGlmStepAIC			
		}
		if(Regression=="Quantile")
		{
			aux_SdCoef<-summary.rq(fit_AIC, se="boot")
			aux_SdGlmStepAIC<-rep(0,p)
			aux_SdGlmStepAIC[aux_betaIndex_1]<-as.vector(aux_SdCoef$coefficients[,2])		
			glmStepAICSdMatrix[k,]<-aux_SdGlmStepAIC			
		}
		if(Regression=="SkewNormal")
		{
			aux_SdGlmStepAIC<-rep(0,p)
			aux_SdGlmStepAIC[aux_betaIndex_1]<-sqrt(diag(vcov(fit_AIC, "DP")))[1:length(aux_betaIndex_1)]	
			glmStepAICSdMatrix[k,]<-aux_SdGlmStepAIC			
		}

	
	
		######################################
		## 	 	 StepBIC glm 		##
		######################################

		if(Regression!="NegBinomial"){t0_StepBIC<-proc.time()}
		if(Regression=="Normal")
		{
			## Step BIC
			fit_BIC<-MyStepCriteria(Regression, PredictorsNames, base, Criteria="BIC" )$fit
		}

		if(Regression=="Binomial")
		{
			## Step BIC
			fit_BIC<-MyStepCriteria(Regression, PredictorsNames, base, Criteria="BIC", ni=ni )$fit								
		}

		if(Regression=="Quantile")
		{
			## Step BIC
			fit_BIC<-MyStepCriteria(Regression, PredictorsNames, base, Criteria="BIC", r_alpha=r_alpha)$fit
		}
		if(Regression=="SkewNormal")
		{
			## Step BIC
			fit_BIC<-MyStepCriteria(Regression, PredictorsNames, base, Criteria="BIC" )$fit
		}
		if(Regression!="NegBinomial"){t1_StepBIC<-proc.time(); time_StepBIC<-time_StepBIC +(t1_StepBIC-t0_StepBIC)[3]}

	
		
		if(Regression=="SkewNormal")
		{ 
			coef_stepBIC<-coef(fit_BIC, "DP")
			PredNamesStepBIC<-head(names(coef_stepBIC)[c(-1)],-2)	#Erasing intercept, omega, and alpha
		}
		if(Regression!="SkewNormal")
		{
			coef_stepBIC<-fit_BIC$coefficients							
			PredNamesStepBIC<-names(fit_BIC$coefficients)[-1]
		}

		## Adding Selected model to the count
		IndexSelectedStepBIC<-match(PredNamesStepBIC, PredictorsNames)
		SelectedModelStepBIC<-rep(0, p-1)
		SelectedModelStepBIC[IndexSelectedStepBIC]<-1

		flag<-0; j<-1
		while(flag==0)
		{
			if( all(SelectedModelsCountsStepBIC_rbind[j,-p]==SelectedModelStepBIC) )
			{
				SelectedModelsCountsStepBIC_rbind[j,p]<- SelectedModelsCountsStepBIC_rbind[j,p]+1
				IndexSelectedModelsStepBIC_list[[j]]<-c(IndexSelectedModelsStepBIC_list[[j]], k)
				flag<-1
			}else{j<-j+1}
			if(j>nrow(SelectedModelsCountsStepBIC_rbind))
			{
				SelectedModelsCountsStepBIC_rbind<-rbind( SelectedModelsCountsStepBIC_rbind, c(SelectedModelStepBIC, 1) )
				IndexSelectedModelsStepBIC_list[[j]]<-k
				flag<-2
			}
		}

		aux_betaIndex_1<-which(c(1,SelectedModelStepBIC)!=0)	#Indexes of beta's such that beta_j!=0 (including intercept)	

		## Coefficient estimated
		aux_betaStepBIC<-rep(0,p)
		if(Regression!="SkewNormal"){aux_betaStepBIC[aux_betaIndex_1]<-as.vector(coefficients(fit_BIC)) }		
		if(Regression=="SkewNormal"){ aux_betaStepBIC[aux_betaIndex_1]<-as.vector(coefficients(fit_BIC, "DP"))[1:length(aux_betaIndex_1)] }
		glmStepBICCoefsMatrix[k,]<-aux_betaStepBIC			
	
		## Other parameters estimates with their SD (Normal, NegBinomial, Skewnormal)
		if(Regression=="NegBinomial")
		{
			glmStepBIC_rEstVector[k]<-summary(fit_BIC)$theta
			glmStepBIC_rSdVector[k]<-summary(fit_BIC)$SE.theta
		}
		if(Regression=="SkewNormal")
		{
			glmStepBIC_lambdaEstVector[k]<-as.vector(coefficients(fit_BIC, "DP"))[length(aux_betaIndex_1) +2]
			glmStepBIC_lambdaSdVector[k]<-sqrt( vcov(fit_BIC, "DP")[length(aux_betaIndex_1) +2, length(aux_betaIndex_1)+2] )
			glmStepBIC_sigma2EstVector[k]<-as.vector(coefficients(fit_BIC, "DP"))[length(aux_betaIndex_1) +1]
			glmStepBIC_sigma2SdVector[k]<-sqrt( vcov(fit_BIC, "DP")[length(aux_betaIndex_1) +1, length(aux_betaIndex_1)+1] )
		}
		if(Regression=="Normal")
		{
			glmStepBIC_sigma2EstVector[k]<-as.vector(summary(fit_BIC)$dispersion)
			p_glmStepBIC<-length(aux_betaIndex_1)
			glmStepBIC_sigma2SdVector[k]<-sqrt(2*(glmStepBIC_sigma2EstVector[k]^2)/(N-p_glmStepBIC))
		}


		## Coefficient's SD estimates
		if(Regression!="Quantile" && Regression!="SkewNormal")	#For Normal, Binomial, and NegBinomial regression
		{
			aux_SdGlmStepBIC<-rep(0,p)
			aux_SdGlmStepBIC[aux_betaIndex_1]<-sqrt(diag(vcov(fit_BIC)))[1:length(aux_betaIndex_1)]		 
			glmStepBICSdMatrix[k,]<-aux_SdGlmStepBIC		
		}
		if(Regression=="Quantile")
		{
			aux_SdCoef<-summary.rq(fit_BIC, se="boot")
			aux_SdGlmStepBIC<-rep(0,p)
			aux_SdGlmStepBIC[aux_betaIndex_1]<-as.vector(aux_SdCoef$coefficients[,2])		
			glmStepBICSdMatrix[k,]<-aux_SdGlmStepBIC			
		}
		if(Regression=="SkewNormal")
		{
			aux_SdGlmStepBIC<-rep(0,p)
			aux_SdGlmStepBIC[aux_betaIndex_1]<-sqrt(diag(vcov(fit_BIC, "DP")))[1:length(aux_betaIndex_1)]	
			glmStepBICSdMatrix[k,]<-aux_SdGlmStepBIC			
		}

		## Uncomment this if you want to save the simulated data
		#text_save<-paste0(Regression,"_","Size",N,"_","p",r_p,"",".Rdata")
		#save.image(file=text_save)
	
	}
	t1<-proc.time()

	## True Model
	TruePredNames<-PredictorsNames[which(RealModel==1)]
	TruePredFormula<-paste0(TruePredNames, collapse = "+")

	## Our method
	MostSelectedModelOurMethod<-as.numeric(which(SelectedModelsCountsOurMethod_rbind[,p] ==max(SelectedModelsCountsOurMethod_rbind[,p])) )
	IndexReplicaMostSelectedOurMethod<-IndexSelectedModelsOurMethod_list[[MostSelectedModelOurMethod]]
	mean_posterior_mean<-mms_apply(CoefsMeanOurMethod_matrix[IndexReplicaMostSelectedOurMethod,], mean)		
	median_posterior_mean<-mms_apply(CoefsMeanOurMethod_matrix[IndexReplicaMostSelectedOurMethod,], median)
	sd_posterior_mean<-mms_apply(CoefsSDOurMethod_matrix[IndexReplicaMostSelectedOurMethod,], mean)		

	SelectedModelOurMethodNames<-colnames(SelectedModelsCountsOurMethod_rbind)[which(SelectedModelsCountsOurMethod_rbind[MostSelectedModelOurMethod, -p]==1)]
	SelectedModelOurMethodFormula<-paste0(SelectedModelOurMethodNames, collapse = "+")
	SelectedModelCountOurMethod_summary<-data.frame(SelectedModel=SelectedModelOurMethodFormula, Count=length(IndexReplicaMostSelectedOurMethod), time=time_OurMethod)	


	## Rusell LA
	MostSelectedModelLA<-as.numeric(which(SelectedModelsCountsLA_rbind[,p] ==max(SelectedModelsCountsLA_rbind[,p])) )
	if( length(MostSelectedModelLA)>1 ){MostSelectedModelLA<-MostSelectedModelLA[1]}
	IndexReplicaMostSelectedLA<-IndexSelectedModelsLA_list[[MostSelectedModelLA]]
	mean_CoefsMeanLA<-mms_apply(CoefsMeanLA_matrix[IndexReplicaMostSelectedLA,], mean)		
	median_CoefsMeanLA<-mms_apply(CoefsMeanLA_matrix[IndexReplicaMostSelectedLA,], median)	
	sd_CoefsMeanLA<-mms_apply(CoefsMeanLA_matrix[IndexReplicaMostSelectedLA,], sd)		

	SelectedModelLANames<-colnames(SelectedModelsCountsLA_rbind)[which(SelectedModelsCountsLA_rbind[MostSelectedModelLA, -(p+1)]==1)]
	SelectedModelLAFormula<-paste0(SelectedModelLANames, collapse = "+")
	SelectedModelCountLA_summary<-data.frame(SelectedModel=SelectedModelLAFormula, Count=length(IndexReplicaMostSelectedLA), time=time_LA)	


	## For stepAIC glm method	
	MostSelectedModelStepAIC<-as.numeric(which(SelectedModelsCountsStepAIC_rbind[,p] ==max(SelectedModelsCountsStepAIC_rbind[,p])) )
	if( length(MostSelectedModelStepAIC)>1 ){MostSelectedModelStepAIC<-MostSelectedModelStepAIC[1]}
	IndexReplicaMostSelectedStepAIC<-IndexSelectedModelsStepAIC_list[[MostSelectedModelStepAIC]]
	glmStepAIC_mean<-mms_apply(glmStepAICCoefsMatrix[IndexReplicaMostSelectedStepAIC,], mean)	
	glmStepAIC_median<-mms_apply(glmStepAICCoefsMatrix[IndexReplicaMostSelectedStepAIC,], median)	
	glmStepAIC_sd_iter<-mms_apply(glmStepAICSdMatrix[IndexReplicaMostSelectedStepAIC,], mean)	
	SelectedModelStepAICNames<-colnames(SelectedModelsCountsStepAIC_rbind)[which(SelectedModelsCountsStepAIC_rbind[MostSelectedModelStepAIC, -p]==1)]
	SelectedModelStepAICFormula<-paste0(SelectedModelStepAICNames, collapse = "+")
	SelectedModelCountStepAIC_summary<-data.frame(SelectedModel=SelectedModelStepAICFormula, Count=length(IndexReplicaMostSelectedStepAIC), time=time_StepAIC)	


	## For stepBIC glm method	
	MostSelectedModelStepBIC<-as.numeric(which(SelectedModelsCountsStepBIC_rbind[,p] ==max(SelectedModelsCountsStepBIC_rbind[,p])) )
	if( length(MostSelectedModelStepBIC)>1 ){MostSelectedModelStepBIC<-MostSelectedModelStepBIC[1]}
	IndexReplicaMostSelectedStepBIC<-IndexSelectedModelsStepBIC_list[[MostSelectedModelStepBIC]]
	glmStepBIC_mean<-mms_apply(glmStepBICCoefsMatrix[IndexReplicaMostSelectedStepBIC,], mean)	
	glmStepBIC_median<-mms_apply(glmStepBICCoefsMatrix[IndexReplicaMostSelectedStepBIC,], median)	
	glmStepBIC_sd_iter<-mms_apply(glmStepBICSdMatrix[IndexReplicaMostSelectedStepBIC,], mean)	
	SelectedModelStepBICNames<-colnames(SelectedModelsCountsStepBIC_rbind)[which(SelectedModelsCountsStepBIC_rbind[MostSelectedModelStepBIC, -p]==1)]
	SelectedModelStepBICFormula<-paste0(SelectedModelStepBICNames, collapse = "+")
	SelectedModelCountStepBIC_summary<-data.frame(SelectedModel=SelectedModelStepBICFormula, Count=length(IndexReplicaMostSelectedStepBIC), time=time_StepBIC)	

	
	if(Regression=="NegBinomial")
	{
		mean_posterior_mean<-c(mean(r_mean_vector[IndexReplicaMostSelectedOurMethod]), mean_posterior_mean)
		median_posterior_mean<-c(median(r_mean_vector[IndexReplicaMostSelectedOurMethod]), median_posterior_mean)
		sd_posterior_mean<-c(mean(r_sd_vector[IndexReplicaMostSelectedOurMethod]), sd_posterior_mean)

		glmStepAIC_mean<-c(mean(glmStepAIC_rEstVector),glmStepAIC_mean)
		glmStepAIC_median<-c(median(glmStepAIC_rEstVector),glmStepAIC_median)
		glmStepAIC_sd_iter<-c(mean(glmStepAIC_rSdVector),glmStepAIC_sd_iter)

		glmStepBIC_mean<-c(mean(glmStepBIC_rEstVector),glmStepBIC_mean)
		glmStepBIC_median<-c(median(glmStepBIC_rEstVector),glmStepBIC_median)
		glmStepBIC_sd_iter<-c(mean(glmStepBIC_rSdVector),glmStepBIC_sd_iter)
	}


	if(Regression=="SkewNormal")
	{
		mean_posterior_mean<-c(mean(lambda_mean_vector[IndexReplicaMostSelectedOurMethod]), mean(sigma2_mean_vector[IndexReplicaMostSelectedOurMethod]), mean_posterior_mean)
		median_posterior_mean<-c(median(lambda_mean_vector[IndexReplicaMostSelectedOurMethod]), median(sigma2_mean_vector[IndexReplicaMostSelectedOurMethod]), median_posterior_mean)
		sd_posterior_mean<-c(mean(lambda_sd_vector[IndexReplicaMostSelectedOurMethod]), mean(sigma2_sd_vector[IndexReplicaMostSelectedOurMethod]), sd_posterior_mean)
		
		glmStepAIC_mean<-c(mean(glmStepAIC_lambdaEstVector), mean(glmStepAIC_sigma2EstVector), glmStepAIC_mean)
		glmStepAIC_median<-c(median(glmStepAIC_lambdaEstVector), median(glmStepAIC_sigma2EstVector), glmStepAIC_median)
		glmStepAIC_sd_iter<-c(mean(glmStepAIC_lambdaSdVector), mean(glmStepAIC_sigma2SdVector), glmStepAIC_sd_iter)

		glmStepBIC_mean<-c(mean(glmStepBIC_lambdaEstVector), mean(glmStepBIC_sigma2EstVector), glmStepBIC_mean)
		glmStepBIC_median<-c(median(glmStepBIC_lambdaEstVector), median(glmStepBIC_sigma2EstVector), glmStepBIC_median)
		glmStepBIC_sd_iter<-c(mean(glmStepBIC_lambdaSdVector), mean(glmStepBIC_sigma2SdVector), glmStepBIC_sd_iter)

	}

	if(Regression=="Normal" || Regression=="Quantile")
	{
		mean_posterior_mean<-c(mean(sigma2_mean_vector[IndexReplicaMostSelectedOurMethod]), mean_posterior_mean)
		median_posterior_mean<-c(median(sigma2_mean_vector[IndexReplicaMostSelectedOurMethod]), median_posterior_mean)
		sd_posterior_mean<-c(mean(sigma2_sd_vector[IndexReplicaMostSelectedOurMethod]), sd_posterior_mean)

		if(Regression=="Normal")	#rq() doesn't provide a sigma2 estimate
		{
			mean_CoefsMeanLA<-c(mean(Sigma2MeanLA_vector[IndexReplicaMostSelectedOurMethod]), mean_CoefsMeanLA)
			median_CoefsMeanLA<-c(median(Sigma2MeanLA_vector[IndexReplicaMostSelectedOurMethod]), median_CoefsMeanLA)
			sd_CoefsMeanLA<-c(sd(Sigma2MeanLA_vector[IndexReplicaMostSelectedOurMethod]), sd_CoefsMeanLA)

			glmStepAIC_mean<-c(mean(glmStepAIC_sigma2EstVector), glmStepAIC_mean)
			glmStepAIC_median<-c(median(glmStepAIC_sigma2EstVector), glmStepAIC_median)
			glmStepAIC_sd_iter<-c(mean(glmStepAIC_sigma2SdVector), glmStepAIC_sd_iter)

			glmStepBIC_mean<-c(mean(glmStepBIC_sigma2EstVector), glmStepBIC_mean)
			glmStepBIC_median<-c(median(glmStepBIC_sigma2EstVector), glmStepBIC_median)
			glmStepBIC_sd_iter<-c(mean(glmStepBIC_sigma2SdVector), glmStepBIC_sd_iter)

		}
	}


	CoefNames<-NULL; for(i in 1:length(r_beta)) { 	CoefNames<-c(CoefNames, paste0("beta", i-1))	}
	if(Regression=="NegBinomial"){CoefNames<-c("r", CoefNames)}
	if(Regression=="SkewNormal"){CoefNames<-c("lambda", "sigma2", CoefNames)}
	if(Regression=="Normal" || Regression=="Quantile"){CoefNames<-c("sigma2", CoefNames)}
	OurMethodResult<-data.frame(mean=mean_posterior_mean, median=median_posterior_mean, sd=sd_posterior_mean)				
	LAResult<-data.frame(mean=mean_CoefsMeanLA, median=median_CoefsMeanLA, sd=sd_CoefsMeanLA)
	StepAICMethodResult<-data.frame(mean=glmStepAIC_mean, median=glmStepAIC_median, sd=glmStepAIC_sd_iter)		
	StepBICMethodResult<-data.frame(mean=glmStepBIC_mean, median=glmStepBIC_median, sd=glmStepBIC_sd_iter)
	rownames(OurMethodResult)<-CoefNames			
	if(Regression=="Normal" || Regression=="Binomial"){rownames(LAResult)<-CoefNames}
	if(Regression=="Quantile"){rownames(StepAICMethodResult)<-CoefNames[-1]; rownames(StepBICMethodResult)<-CoefNames[-1]}else{
					rownames(StepAICMethodResult)<-CoefNames; rownames(StepBICMethodResult)<-CoefNames}

	
	aux<-list(DataList=DataList, coefs_PostMean_matrix_OurMethod=CoefsMeanOurMethod_matrix,	
	True_Model= TruePredFormula, "Selected_Model_Count_OurMethod"=SelectedModelCountOurMethod_summary, 
	"Selected_Model_Count_LA"=SelectedModelCountLA_summary,
	"Selected_Model_Count_StepAIC"=SelectedModelCountStepAIC_summary,
	"Selected_Model_Count_StepBIC"=SelectedModelCountStepBIC_summary,
	Model_Exploration_OurMethod=SelectedModelsCountsOurMethod_rbind, 
	Model_Exploration_LA=SelectedModelsCountsLA_rbind,
	Model_Exploration_StepAIC=SelectedModelsCountsStepAIC_rbind,
	Model_Exploration_StepBIC=SelectedModelsCountsStepBIC_rbind,
	"OurMethod_Posterior_Summary"=OurMethodResult, "LA_Posterior_Summary"=LAResult,
	"StepAIC_Summary"=StepAICMethodResult, "StepBIC_Summary"=StepBICMethodResult, Regression=Regression,
	"Time_Minutes"=(t1 -t0)[3]/60)

	return(aux)

}
############################################################################################################
R<-100	#Number of experiments
nchain<-10000; burnin<-2000	#chain large and burn-in
tau2<-1000; rho<-1	#tau and rho value



############################################
## 		LiR Regression			##
############################################
r_sigma2<-2		#value of the true sigma^2
Regression<-"Normal"

## Scenario 1 (Manuscript, for Table 1 and 2)
set.seed(31415)
N<-100
r_beta<-c(1, 0, 2, 0, 3, 2)
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results

## Scenario 2 (Supplementary material, for Table 1)
set.seed(31415)
N<-200
r_beta<-c( sample(c(1.5,-1,2,-2.5), size=70, replace=TRUE), rep(0,30)  )
#r_beta<-c( sample(c(1.5,-1,2,-2.5), size=100, replace=TRUE) )
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results

## Scenario 3 (Supplementary material, for Table 2)
set.seed(31415)
N<-200
r_beta<-c( rep(c(0.1,0.2,0.3,0.4),5), sample(c(1.5,-1,2,-2.5), size=30, replace=TRUE) )
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results


############################################
## 		LoR Regression			##
############################################
Regression<-"Binomial"

## Scenario 1 (Manuscript, for Table 3 and 4)
set.seed(31415)
N<-200
r_beta<-r_beta<-c(1, 0, 2, 0, 3, 2)
p<-length(r_beta)
ni<-rep(1, N)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results


## Scenario 2 (Supplementary material, for Table 3)
set.seed(31415)
N<-500
r_beta<-c(sample(c(-1,1), 20, replace=TRUE))
p<-length(r_beta)
ni<-rep(1, N)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results


## Scenario 3 (Supplementary material, for Table 4)
set.seed(31415)
N<-1000
r_beta<-c(sample(c(-1,1), 30, replace=TRUE), rep(0,20))
p<-length(r_beta)
ni<-rep(1, N)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results


############################################
## 		NB Regression			##
############################################
Regression<-"NegBinomial"

## Scenario 1 (Manuscript, for Table 5 and 6)
set.seed(31415)
N<-200
r_r<-2; r_beta<-c(0.5, -0.8, 1, 0, 0.4, -0.7)
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results

## Scenario 2 (Supplementary material, for Table 5)
set.seed(31415)
N<-1000
r_r<-2; r_beta<-c(rep(c(0.5, -0.8, 0.4, -0.7),5), rep(0,10))
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results


############################################
## 		QR Regression			##
############################################
r_sigma2<-1
Regression<-"Quantile"

## Scenario 1 (Manuscript, for Table 7 and 8)
N<-200
r_beta<-c(1, 0, 2, 0, 3, 2)
r_p<-length(r_beta)
alphaVector<-c(0.1, 0.5, 0.9)
for(j in 1:length(alphaVector))
{
	set.seed(31415)
	r_alpha<-alphaVector[j]
	DataIter<-j
	if(r_alpha==0.1){aux_0.1<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, DataIter)}
	if(r_alpha==0.5){aux_0.5<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, DataIter)}
	if(r_alpha==0.9){aux_0.9<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, DataIter)}
}
Summary_SimStudy(aux_0.1, ExploredModels=TRUE)	#Summary results
Summary_SimStudy(aux_0.5, ExploredModels=TRUE)	#Summary results
Summary_SimStudy(aux_0.9, ExploredModels=TRUE)	#Summary results

## Scenario 2 (Supplementary material, for Table 6)
N<-1000
r_beta<-c(sample(c(1,-1), size=30, replace=TRUE), rep(0,20))
r_p<-length(r_beta)
alphaVector<-c(0.1, 0.5, 0.9)
for(j in 1:length(alphaVector))
{
	set.seed(31415)
	r_alpha<-alphaVector[j]
	DataIter<-j
	if(r_alpha==0.1){aux_0.1<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, DataIter)}
	if(r_alpha==0.5){aux_0.5<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, DataIter)}
	if(r_alpha==0.9){aux_0.9<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, DataIter)}
}
Summary_SimStudy(aux_0.1, ExploredModels=TRUE)	#Summary results
Summary_SimStudy(aux_0.5, ExploredModels=TRUE)	#Summary results
Summary_SimStudy(aux_0.9, ExploredModels=TRUE)	#Summary results


############################################
## 		SNR  Regression			##
############################################
Regression="SkewNormal"
r_sigma2<-1.5; r_lambda<-4	#true value for sigma^2 and lambda
d=2; b2=1/2			#for uniform(-1,1) prior on delta

## Scenario 1 (Manuscript, for Table 9 and 10)
set.seed(31415)
N<-200
r_beta<-r_beta<-c(1, 0, 2, 0, 3, 2)
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results


## Scenario 2 (Supplementary material, for Table 7)
set.seed(31415)
N<-500
r_beta<-c(rep(1,15), rep(-1,15), rep(0,20))
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results


## Scenario 3 (Supplementary material, for Table 8)
set.seed(31415)
N<-500
r_beta<-c(rep(1,30), rep(-1,30), rep(0,40))
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results

