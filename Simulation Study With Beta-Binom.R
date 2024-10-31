rm(list=ls(all=TRUE))
source("functions.R")


library(mombf)		#Rosell, LA method
library(MASS)		#glm.nb, regresion binomial negativa, steAIC()
library(quantreg)	#Quantile Regression
library(sn)		#Frequentist Regression Skew Normal
library(xtable)
############################################################################################################

##################################################################################################
##						 Functions, Just run this						##
##################################################################################################

## Step Criteria manual (only backward)

## Step Criteria manual (only backward)
MyStepCriteria<-function(Regression, PredictorsNames, base, Criteria="AIC", ... )
{
	Criteria_lm<-function(fit, N, p, Criteria)
	{	
		if(Criteria=="AIC"){ result<-2*p -2*logLik(fit)[1] }	
		if(Criteria=="BIC"){ result<-p*log(N) -2*logLik(fit)[1] }
		return(result)
	}

	N<-nrow(base)
	r_p<-ncol(base)
	p<-r_p -1
	aux_PredictorsNames<-PredictorsNames
	flag<-0; j<-1
	flag_intercept<-0
	SelectedCriteria<-c()

	if(Regression=="Normal"){ p_aux<-p+2; full_mods<- glm(y~., data = base, family="gaussian") }
	if(Regression=="Binomial"){ full_mods<- glm(y~., data = base, family = "binomial", weights = ni) }
	if(Regression=="NegBinomial"){ p_aux<-p+2; full_mods<- glm.nb(y~., data=base, init.theta=0.5) }
	if(Regression=="Quantile"){ full_mods<- rq(y~., tau=r_alpha, data=base) }
	if(Regression=="SkewNormal"){ p_aux<-p+3; full_mods<- selm(y~., data = base, family="SN") }
	Criteria_FullModel<-Criteria_lm(full_mods, N, p=p_aux, Criteria)

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
		if(p==0){aux_formula<-paste0("y ~ 1"); Candidates_PredictorsNames<-"Intercept"}
		
		if(Regression=="Normal"){ p_aux<-p+2; aux_mods<- lapply(aux_formula, function(frml) glm(frml, data = base, family="gaussian")) }
		if(Regression=="Binomial"){ aux_mods<- lapply(aux_formula, function(frml) glm(frml, data = base, family = "binomial", weights = ni)) }
		if(Regression=="NegBinomial"){ p_aux<-p+2; aux_mods<- lapply(aux_formula, function(frml) glm.nb(frml, data=base, init.theta=0.5) ) }
		if(Regression=="Quantile"){ aux_mods<- lapply(aux_formula, function(frml) rq(frml, tau=r_alpha, data=base)) }
		if(Regression=="SkewNormal"){ p_aux<-p+3; aux_mods<- lapply(aux_formula, function(frml) selm(frml, data = base, family="SN")) }
		aux_Criteria<-unlist(lapply(aux_mods, Criteria_lm, N=N, p=p_aux, Criteria=Criteria))
		
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
			if(Criteria_FullModel<SelectedCriteria){ flag<-1; SelectedPredictors<-PredictorsNames; old_aux_index<-1; old_aux_mods<-full_mods  }
		}
		if(p==0 && flag==0 && Candidates_PredictorsNames=="Intercept"){flag<-1; flag_intercept<-1}	
	}
	if(p==0 && flag==1){ SelectedPredictors<-"X1" }
	if(flag_intercept==0){ SelectedFormula<-paste0("y ~ ", paste0(SelectedPredictors, collapse = "+"))	 }else{
	SelectedFormula<-aux_formula }		
	list(SelectedFormula=SelectedFormula, fit=old_aux_mods[[old_aux_index]])
}


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
	r_p<-ncol(aux$Model_Exploration_OurMethod_womack)
	Regression<-aux$Regression
	
	## Obtaining Explored Model Table from output
	ExploredMOurMethod_womack<-data.frame(aux$Model_Exploration_OurMethod_womack)
	ExploredMOurMethod_bb<-data.frame(aux$Model_Exploration_OurMethod_bb)
	ExploredMLA<-data.frame(aux$Model_Exploration_LA)	
	ExploredMStepAIC<-data.frame(aux$Model_Exploration_StepAIC)
	ExploredMStepBIC<-data.frame(aux$Model_Exploration_StepBIC)


	## Changing predictor colnames for Latex
	aux_namespred<-vector(length=r_p-1)
	for(i in 1:(r_p -1)){ aux_namespred[i]<-paste0("$x_",i,"$") }
	colnames(ExploredMOurMethod_womack)[-r_p]<-aux_namespred
	colnames(ExploredMOurMethod_bb)[-r_p]<-aux_namespred
	colnames(ExploredMLA)[-r_p]<-aux_namespred
	colnames(ExploredMStepAIC)[-r_p]<-aux_namespred
	colnames(ExploredMStepBIC)[-r_p]<-aux_namespred
	

	## Sorting rows by frequency
	ExploredMOurMethod_womack<- rbind(ExploredMOurMethod_womack[1,],  ExploredMOurMethod_womack[-1,][order(ExploredMOurMethod_womack[-1,]$Frequency, decreasing=TRUE),] )
	ExploredMOurMethod_bb<- rbind(ExploredMOurMethod_bb[1,],  ExploredMOurMethod_bb[-1,][order(ExploredMOurMethod_bb[-1,]$Frequency, decreasing=TRUE),] )
	ExploredMLA<- rbind(ExploredMLA[1,],  ExploredMLA[-1,][order(ExploredMLA[-1,]$Frequency, decreasing=TRUE),] )
	ExploredMStepAIC<- rbind(ExploredMStepAIC[1,],  ExploredMStepAIC[-1,][order(ExploredMStepAIC[-1,]$Frequency, decreasing=TRUE),] )
	ExploredMStepBIC<- rbind(ExploredMStepBIC[1,],  ExploredMStepBIC[-1,][order(ExploredMStepBIC[-1,]$Frequency, decreasing=TRUE),] )
	
	PercentagesOurMethod_womack<-ExploredMOurMethod_womack$Frequency
	PercentagesOurMethod_bb<-ExploredMOurMethod_bb$Frequency
	PercentagesLA<-ExploredMLA$Frequency
	PercentagesStepAIC<-ExploredMStepAIC$Frequency
	PercentagesStepBIC<-ExploredMStepBIC$Frequency
	

	## Changing frequency values to percentages in text, for latex
	aux_EM2<-vector(length=nrow(ExploredMOurMethod_womack))
	for(i in 1:nrow(ExploredMOurMethod_womack)){ aux_EM2[i]<-paste0(ExploredMOurMethod_womack$Frequency[i], "%") }
	ExploredMOurMethod_womack$Frequency<-aux_EM2

	aux_EM2<-vector(length=nrow(ExploredMOurMethod_bb))
	for(i in 1:nrow(ExploredMOurMethod_bb)){ aux_EM2[i]<-paste0(ExploredMOurMethod_bb$Frequency[i], "%") }
	ExploredMOurMethod_bb$Frequency<-aux_EM2

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
	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMOurMethod_womack) -1))
	ExploredMOurMethod_womack<- head( cbind(" "=aux_EM, ExploredMOurMethod_womack), 4 )

	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMOurMethod_bb) -1))
	ExploredMOurMethod_bb<- head( cbind(" "=aux_EM, ExploredMOurMethod_bb), 4 )

	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMLA) -1))
	ExploredMLA<- head( cbind(" "=aux_EM, ExploredMLA), 4 )
	
	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMStepAIC) -1))
	ExploredMStepAIC<- head( cbind(" "=aux_EM, ExploredMStepAIC), 4 )
	
	aux_EM<-c( "True Model", rep(" ", nrow(ExploredMStepBIC) -1))
	ExploredMStepBIC<-head( cbind(" "=aux_EM, ExploredMStepBIC), 4 )
	

	## Latex All in one table
	if(Regression=="Normal" || Regression=="Binomial" || Regression=="Quantile")
	{
		ExploredModel<-list(Our_Method_womack=ExploredMOurMethod_womack, Our_Method_bb=ExploredMOurMethod_bb, LA=ExploredMLA, Step_AIC=ExploredMStepAIC, Step_BIC=ExploredMStepBIC)
	}else{
		ExploredModel<-list(Our_Method_womack=ExploredMOurMethod_womack, Our_Method_bb=ExploredMOurMethod_bb, Step_AIC=ExploredMStepAIC, Step_BIC=ExploredMStepBIC)
	}
	
	
	#############################
	## 	  Estimation    	   ##
	#############################
	maxOurMethod_womack<-max(PercentagesOurMethod_womack)
	maxOurMethod_bb<-max(PercentagesOurMethod_bb)
	maxLA<-max(PercentagesLA)
	maxStepAIC<-max(PercentagesStepAIC)
	maxStepBIC<-max(PercentagesStepBIC)
	
	flag_OurMethod<-0; flag_LA<-0; flag_StepAIC<-0; flag_StepBIC<-0
	if(length(which(maxOurMethod_womack==PercentagesOurMethod_womack))==1){flag_OurMethod_womack<-1}
	if(length(which(maxOurMethod_bb==PercentagesOurMethod_bb))==1){flag_OurMethod_bb<-1}
	if(length(which(maxLA==PercentagesLA))==1){flag_LA<-1}
	if(length(which(maxStepAIC==PercentagesStepAIC))==1){flag_StepAIC<-1}
	if(length(which(maxStepBIC==PercentagesStepBIC))==1){flag_StepBIC<-1}
	
	## Estimates
	EstimatesOurMethod_womack<-aux$OurMethod_Posterior_Summary_womack[,c(1,3)]
	rep_blank<-nrow(EstimatesOurMethod_womack)
	if(flag_OurMethod_womack==0){EstimatesOurMethod_womack$mean<-rep("- - -", rep_blank); EstimatesOurMethod_womack$sd<-rep("- - -", rep_blank)}

	EstimatesOurMethod_bb<-aux$OurMethod_Posterior_Summary_bb[,c(1,3)]
	rep_blank<-nrow(EstimatesOurMethod_bb)
	if(flag_OurMethod_bb==0){EstimatesOurMethod_bb$mean<-rep("- - -", rep_blank); EstimatesOurMethod_bb$sd<-rep("- - -", rep_blank)}
	
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
		Estimates<-data.frame(parameters, EstimatesOurMethod_womack, EstimatesOurMethod_bb, EstimatesLA, EstimatesStepAIC, EstimatesStepBIC)
		rownames(Estimates)<-row_names
		colnames(Estimates)[-1]<-c("mean_OurMethod_womack", "sd_OurMethod_womack", "mean_OurMethod_bb", "sd_OurMethod_bb", "mean_LA", "sd_LA", "mean_StepAIC", "sd_StepAIC" ,"mean_StepBIC", "sd_StepBIC")
	}else{
		Estimates<-data.frame(parameters, EstimatesOurMethod_womack, EstimatesOurMethod_bb, EstimatesStepAIC, EstimatesStepBIC)
		rownames(Estimates)<-row_names
		colnames(Estimates)[-1]<-c("mean_OurMethod_womack", "sd_OurMethod_womack", "mean_OurMethod_bb", "sd_OurMethod_bb", "mean_StepAIC", "sd_StepAIC" ,"mean_StepBIC", "sd_StepBIC")
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
	intercept_first_excluded<- first_excluded +1	
	p_selection<-p -intercept_first_excluded		
	DataList<-list()					
	RealModel<-+(r_beta!=0)[-1]	
	PredictorsIndex1_Real<-which(RealModel==1)	
	PredictorsNames<-colnames(data.frame(matrix(ncol=p-1)))	


	## Preparation for Our method (Womack prior)
	SelectedModelsCountsOurMethod_rbind<-rbind( c(RealModel,0) )	#Counting selected models in replicas
	colnames(SelectedModelsCountsOurMethod_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsOurMethod_list<-NULL
	CoefsMeanOurMethod_matrix<-matrix(0,ncol=p, nrow=R) 
	CoefsSDOurMethod_matrix<-matrix(0,ncol=p, nrow=R) 
	if(Regression=="Normal"){sigma2_mean_vector<-vector(length=R); sigma2_sd_vector<-vector(length=R) }
	if(Regression=="NegBinomial"){r_mean_vector<-vector(length=R); r_sd_vector<-vector(length=R) }
	if(Regression=="Quantile"){sigma2_mean_vector<-vector(length=R); sigma2_sd_vector<-vector(length=R) }
	if(Regression=="SkewNormal"){lambda_mean_vector<-vector(length=R); lambda_sd_vector<-vector(length=R); sigma2_mean_vector<-vector(length=R); sigma2_sd_vector<-vector(length=R) }	
	time_OurMethod<-0

	## Preparation for Our method (Beta-Binomial Prior)
	SelectedModelsCountsOurMethod_bb_rbind<-rbind( c(RealModel,0) )	#Counting selected models in replicas
	colnames(SelectedModelsCountsOurMethod_bb_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsOurMethod_bb_list<-NULL
	CoefsMeanOurMethod_bb_matrix<-matrix(0,ncol=p, nrow=R) 
	CoefsSDOurMethod_bb_matrix<-matrix(0,ncol=p, nrow=R)
	if(Regression=="Normal"){sigma2_bb_mean_vector<-vector(length=R); sigma2_bb_sd_vector<-vector(length=R) }
	if(Regression=="NegBinomial"){r_bb_mean_vector<-vector(length=R); r_bb_sd_vector<-vector(length=R) }
	if(Regression=="Quantile"){sigma2_bb_mean_vector<-vector(length=R); sigma2_bb_sd_vector<-vector(length=R) }
	if(Regression=="SkewNormal"){lambda_bb_mean_vector<-vector(length=R); lambda_bb_sd_vector<-vector(length=R); sigma2_bb_mean_vector<-vector(length=R); sigma2_bb_sd_vector<-vector(length=R) }	
	time_OurMethod_bb<-0

	## Preparation for Rusell LA method
	SelectedModelsCountsLA_rbind<-rbind( c(RealModel,0) )	#Counting selected models in replicas
	colnames(SelectedModelsCountsLA_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsLA_list<-NULL
	CoefsMeanLA_matrix<-matrix(0,ncol=p, nrow=R) #Matriz de la media a posteriori de beta en cada R�plica 
	#CoefsSdLA_matrix<-matrix(0,ncol=p, nrow=R) #Matriz de la sd a posteriori de beta en cada R�plica 
	if(Regression=="Normal"){Sigma2MeanLA_vector<-vector(length=R) }
	#if(Regression=="Normal"){Sigma2MeanLA_vector<-vector(length=R); Sigma2SdLA_vector<-vector(length=R) }
	time_LA<-0
	
	## Preparation for StepAIC method
	glmStepAICCoefsMatrix<-matrix(0,ncol=p, nrow=R)  
	glmStepAICSdMatrix<-matrix(0,ncol=p, nrow=R)   
	SelectedModelsCountsStepAIC_rbind<-rbind( c(RealModel,0) )	#Counting selected models in replicas
	colnames(SelectedModelsCountsStepAIC_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsStepAIC_list<-NULL
	if(Regression=="Normal"){glmStepAIC_sigma2EstVector<-vector(length=R); glmStepAIC_sigma2SdVector<-vector(length=R)}
	if(Regression=="NegBinomial"){glmStepAIC_rEstVector<-vector(length=R); glmStepAIC_rSdVector<-vector(length=R)}
	if(Regression=="SkewNormal"){glmStepAIC_lambdaEstVector<-vector(length=R); glmStepAIC_lambdaSdVector<-vector(length=R); glmStepAIC_sigma2EstVector<-vector(length=R); glmStepAIC_sigma2SdVector<-vector(length=R)}
	time_StepAIC<-0

	## Preparation for StepBIC method
	glmStepBICCoefsMatrix<-matrix(0,ncol=p, nrow=R)   
	glmStepBICSdMatrix<-matrix(0,ncol=p, nrow=R) 
	SelectedModelsCountsStepBIC_rbind<-rbind( c(RealModel,0) )	#Counting selected models in replicas
	colnames(SelectedModelsCountsStepBIC_rbind)<-c(PredictorsNames, "Frequency")
	IndexSelectedModelsStepBIC_list<-NULL
	if(Regression=="Normal"){glmStepBIC_sigma2EstVector<-vector(length=R); glmStepBIC_sigma2SdVector<-vector(length=R)}
	if(Regression=="NegBinomial"){glmStepBIC_rEstVector<-vector(length=R); glmStepBIC_rSdVector<-vector(length=R)}
	if(Regression=="SkewNormal"){glmStepBIC_lambdaEstVector<-vector(length=R); glmStepBIC_lambdaSdVector<-vector(length=R); glmStepBIC_sigma2EstVector<-vector(length=R); glmStepBIC_sigma2SdVector<-vector(length=R)}
	time_StepBIC<-0

	y<-vector(length=N)
	for(k in 1:R)
	{
		cat("  Replica", k, "\r")

		## Pre-simulation stuff
		aux_cov<-rnorm((p-1)*N, 0, 1)				
		#aux_cov<-runif((p-1)*N, 0,1)


		if(Regression=="Normal")
		{
			#X<-matrix( c(rep(1, N), aux_cov), ncol=p )
			X<-matrix( c(rep(1, N), rnorm((p -1)*N)), ncol=p )
			Xbeta<-X%*%r_beta
			y<-rnorm(N, mean=Xbeta, sd=sqrt(r_sigma2))
			Covariates<-X[,2:(length(r_beta))]
			base<-data.frame(y, Covariates)

			## For LA method
			familyReg<-"normal"	
			X_LA<-as.matrix(X); y_LA<-y
		}

		if(Regression=="Binomial")
		{							
			Covariates<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N))
			base<-gen_base_binomial_reg(r_beta, Covariates, N, ni=ni)
			base<-base[,-2]
			y<-base$y
			X<-data.frame("Intercept"=rep(1,N), Covariates)

			## For LA method
			familyReg<-"binomial"
			NewBase<-DupBaseBinomial(base)
			X_LA<-as.matrix(data.frame(Intercept=rep(1, nrow(NewBase)), X=NewBase[,c(-1,-2)] ))
			y_LA<-NewBase$y

		}

		if(Regression=="NegBinomial")
		{
			Check<-FALSE
			while(Check==FALSE)		#"Arreglando" el problema de glm.nb
			{
				Covariates<-data.frame(matrix(aux_cov, ncol=p-1, nrow=N))		
				base<-gen_base_NegBinomial_reg(N, r_beta, r_r, Covariates)
				y<-base$y
			
				## Adelantando pega del "fit them all glm" por tema del problema del glm.nb
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
			X<-data.frame("Intercept"=rep(1,N),Covariates)
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
			Covariates<-X[,2:(length(r_beta))]
			base<-data.frame(y, Covariates)

			## For LA method
			familyReg<-"twopiecelaplace"
			X_LA<-as.matrix(X); y_LA<-y
						
		}

		if(Regression=="SkewNormal")
		{
			#X<-matrix( c(rep(1, N), aux_cov), ncol=p )
			X<-matrix( c(rep(1, N), rnorm((p -1)*N)), ncol=p )
			Xbeta<-X%*%r_beta
			r_w<-abs(rnorm(N, mean=0, sd=1))
			r_delta<-r_lambda/sqrt(r_lambda^2 +1)
			y<-Xbeta +sqrt(r_sigma2)*r_delta*r_w +sqrt(r_sigma2*(1 -r_delta^2))*rnorm(N, mean=0, sd=1)
			Covariates<-X[,2:(length(r_beta))]
			base<-data.frame(y, Covariates)		
		}
		X<-as.matrix(X)

		DataList[[k]]<-base



		######################################
		##			BMS 		##
		##		(Womack prior)		##
		######################################
		t0_OurMethod<-proc.time()
		## Gibbs
		if(Regression=="Normal")
		{
			fit_OurMethod<-gibbs_abms(y, Covariates, family="LiR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=TRUE, a_bb=1, b_bb=1, count.iteration=FALSE )
		}
		if(Regression=="Binomial")
		{
			fit_OurMethod<-gibbs_abms(y, Covariates, family="LoR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=TRUE, a_bb=1, b_bb=1, count.iteration=FALSE )
		}

		if(Regression=="NegBinomial")
		{
			fit_OurMethod<-gibbs_abms(y, Covariates, family="NBR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=TRUE, a_bb=1, b_bb=1, count.iteration=FALSE )
		}

		if(Regression=="Quantile")
		{
			fit_OurMethod<-gibbs_abms(y, Covariates, family="QR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=TRUE, a_bb=1, b_bb=1, count.iteration=FALSE )	
		}
		
		if(Regression=="SkewNormal")
		{
			fit_OurMethod<-gibbs_abms(y, Covariates, family="SNR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=TRUE, a_bb=1, b_bb=1, count.iteration=FALSE )
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


		##################################################
		##			BMS 			##
		##		(Beta binomial prior)		##
		##################################################
		t0_OurMethod_bb<-proc.time()
		## Gibbs
		if(Regression=="Normal")
		{
			fit_OurMethod_bb<-gibbs_abms(y, Covariates, family="LiR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=FALSE, a_bb=1, b_bb=p_selection^(1.1), count.iteration=FALSE )
		}
		if(Regression=="Binomial")
		{
			fit_OurMethod_bb<-gibbs_abms(y, Covariates, family="LoR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=FALSE, a_bb=1, b_bb=p_selection^(1.1), count.iteration=FALSE )
		}

		if(Regression=="NegBinomial")
		{
			fit_OurMethod_bb<-gibbs_abms(y, Covariates, family="NBR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=FALSE, a_bb=1, b_bb=p_selection^(1.1), count.iteration=FALSE )
		}

		if(Regression=="Quantile")
		{
			fit_OurMethod_bb<-gibbs_abms(y, Covariates, family="QR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=FALSE, a_bb=1, b_bb=p_selection^(1.1), count.iteration=FALSE )	
		}
		
		if(Regression=="SkewNormal")
		{
			fit_OurMethod_bb<-gibbs_abms(y, Covariates, family="SNR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, ni=rep(1, length(y)), alpha=0.5,
                     a0=1, b0=1, d=2, b2=1/2, model_fixed=NULL, womack=FALSE, a_bb=1, b_bb=p_selection^(1.1), count.iteration=FALSE )
		}
		t1_OurMethod_bb<-proc.time()
		time_OurMethod_bb<-time_OurMethod_bb +(t1_OurMethod_bb -t0_OurMethod_bb)[3]

		
		## Adding Selected model to the count
		SelectedModelOurMethod_bb<-as.numeric(mcr(fit_OurMethod_bb$model_chain))
		flag<-0; j<-1
		while(flag==0)
		{
			if( all(SelectedModelsCountsOurMethod_bb_rbind[j,-p]==SelectedModelOurMethod_bb) )
			{
				SelectedModelsCountsOurMethod_bb_rbind[j,p]<- SelectedModelsCountsOurMethod_bb_rbind[j,p]+1
				IndexSelectedModelsOurMethod_bb_list[[j]]<-c(IndexSelectedModelsOurMethod_bb_list[[j]], k)
				flag<-1
			}else{j<-j+1}
			if(j>nrow(SelectedModelsCountsOurMethod_bb_rbind))
			{
				SelectedModelsCountsOurMethod_bb_rbind<-rbind( SelectedModelsCountsOurMethod_bb_rbind, c(SelectedModelOurMethod_bb, 1) )
				IndexSelectedModelsOurMethod_bb_list[[j]]<-k
				flag<-2
			}
		}

		index_selected_bb<-which(colSums(t(fit_OurMethod_bb$model_chain) == SelectedModelOurMethod_bb) == ncol(fit_OurMethod_bb$model_chain))
		coefs_chain_bb<-fit_OurMethod_bb$beta_chain[index_selected_bb,]	
		CoefsMeanOurMethod_bb_matrix[k,]<-apply(coefs_chain_bb, 2, mean)
		CoefsSDOurMethod_bb_matrix[k,]<-apply(coefs_chain_bb, 2, sd)

		if(Regression=="Normal" || Regression=="Quantile" )
		{
			sigma2_bb_mean_vector[k]<-mean(fit_OurMethod_bb$sigma2_chain[index_selected_bb])
			sigma2_bb_sd_vector[k]<-sd(fit_OurMethod_bb$sigma2_chain[index_selected_bb])
		}
		if(Regression=="NegBinomial")
		{
			r_bb_mean_vector[k]<-mean(fit_OurMethod_bb$r_chain[index_selected])
			r_bb_sd_vector[k]<-sd(fit_OurMethod_bb$r_chain[index_selected])
		}
		if(Regression=="SkewNormal")
		{
			lambda_bb_mean_vector[k]<-mean(fit_OurMethod_bb$lambda_chain[index_selected])
			lambda_bb_sd_vector[k]<-sd(fit_OurMethod_bb$lambda_chain[index_selected])
			sigma2_bb_mean_vector[k]<-mean(fit_OurMethod_bb$sigma2_chain[index_selected])
			sigma2_bb_sd_vector[k]<-sd(fit_OurMethod_bb$sigma2_chain[index_selected])
		}


		
		######################################
		## 	 	 Rusell LA 			##
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
			verbose=FALSE, center=FALSE, scale=FALSE, niter=nchain, burnin=burnin, includevars=1)
		}else{	
			fit_LA <- modelSelection(y=y_LA, x=X_LA, family=familyReg, priorCoef=prCoef, priorDelta=prDelta, method="Laplace",
			verbose=FALSE, center=FALSE, scale=FALSE, niter=nchain, burnin=burnin, includevars=1)
		}
		t1_LA<-proc.time()
		time_LA<-time_LA +(t1_LA-t0_LA)[3]

		## Coefficients estimates
		if(Regression=="Normal"){ CoefsMeanLA_matrix[k,]<-as.vector(coef(fit_LA)[c(-1, -(p+2)),1]) }
		if(Regression=="Binomial"){ CoefsMeanLA_matrix[k,]<-as.vector(coef(fit_LA)[-1,1]) }

		if(Regression=="Normal" )
		{
			Sigma2MeanLA_vector[k]<-as.vector(coef(fit_LA)[p+2,1])
		}

		## Explored Models
		#IndexSelectedLA<-as.vector(which(fit_LA$postMode==1))
		#SelectedModelLA<-rep(0, p-1)
		#SelectedModelLA[IndexSelectedLA]<-1
		SelectedModelLA<-as.vector(fit_LA$postMode)[-1]		#Excluding intercept
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

		#SelectedModelsCountsLA_rbind
		#IndexSelectedModelsLA_list

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
			PredNamesStepAIC<-head(names(coef_stepAIC)[c(-1)],-2)	#Erasing intercept, omega, and alpha
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

	## Our method (Womack prior)
	MostSelectedModelOurMethod<-as.numeric(which(SelectedModelsCountsOurMethod_rbind[,p] ==max(SelectedModelsCountsOurMethod_rbind[,p])) )
	IndexReplicaMostSelectedOurMethod<-IndexSelectedModelsOurMethod_list[[MostSelectedModelOurMethod]]
	mean_posterior_mean<-mms_apply(CoefsMeanOurMethod_matrix[IndexReplicaMostSelectedOurMethod,], mean)		
	median_posterior_mean<-mms_apply(CoefsMeanOurMethod_matrix[IndexReplicaMostSelectedOurMethod,], median)
	sd_posterior_mean<-mms_apply(CoefsSDOurMethod_matrix[IndexReplicaMostSelectedOurMethod,], mean)			

	SelectedModelOurMethodNames<-colnames(SelectedModelsCountsOurMethod_rbind)[which(SelectedModelsCountsOurMethod_rbind[MostSelectedModelOurMethod, -p]==1)]
	SelectedModelOurMethodFormula<-paste0(SelectedModelOurMethodNames, collapse = "+")
	SelectedModelCountOurMethod_summary<-data.frame(SelectedModel=SelectedModelOurMethodFormula, Count=length(IndexReplicaMostSelectedOurMethod), time=time_OurMethod)	


	## Our method (beta-binomial prior)
	MostSelectedModelOurMethod_bb<-as.numeric(which(SelectedModelsCountsOurMethod_bb_rbind[,p] ==max(SelectedModelsCountsOurMethod_bb_rbind[,p])) )
	IndexReplicaMostSelectedOurMethod_bb<-IndexSelectedModelsOurMethod_bb_list[[MostSelectedModelOurMethod_bb]]
	mean_posterior_mean_bb<-mms_apply(CoefsMeanOurMethod_bb_matrix[IndexReplicaMostSelectedOurMethod_bb,], mean)		
	median_posterior_mean_bb<-mms_apply(CoefsMeanOurMethod_bb_matrix[IndexReplicaMostSelectedOurMethod_bb,], median)	
	sd_posterior_mean_bb<-mms_apply(CoefsSDOurMethod_bb_matrix[IndexReplicaMostSelectedOurMethod_bb,], mean)			

	SelectedModelOurMethodNames_bb<-colnames(SelectedModelsCountsOurMethod_bb_rbind)[which(SelectedModelsCountsOurMethod_bb_rbind[MostSelectedModelOurMethod_bb, -p]==1)]
	SelectedModelOurMethodFormula_bb<-paste0(SelectedModelOurMethodNames_bb, collapse = "+")
	SelectedModelCountOurMethod_bb_summary<-data.frame(SelectedModel=SelectedModelOurMethodFormula_bb, Count=length(IndexReplicaMostSelectedOurMethod_bb), time=time_OurMethod_bb)	


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

		mean_posterior_mean_bb<-c(mean(r_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), mean_posterior_mean_bb)
		median_posterior_mean_bb<-c(median(r_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), median_posterior_mean_bb)
		sd_posterior_mean_bb<-c(mean(r_bb_sd_vector[IndexReplicaMostSelectedOurMethod_bb]), sd_posterior_mean_bb)

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

		mean_posterior_mean_bb<-c(mean(lambda_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), mean(sigma2_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), mean_posterior_mean_bb)
		median_posterior_mean_bb<-c(median(lambda_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), median(sigma2_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), median_posterior_mean_bb)
		sd_posterior_mean_bb<-c(mean(lambda_bb_sd_vector[IndexReplicaMostSelectedOurMethod_bb]), mean(sigma2_bb_sd_vector[IndexReplicaMostSelectedOurMethod_bb]), sd_posterior_mean_bb)
		
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

		mean_posterior_mean_bb<-c(mean(sigma2_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), mean_posterior_mean_bb)
		median_posterior_mean_bb<-c(median(sigma2_bb_mean_vector[IndexReplicaMostSelectedOurMethod_bb]), median_posterior_mean_bb)
		sd_posterior_mean_bb<-c(mean(sigma2_bb_sd_vector[IndexReplicaMostSelectedOurMethod_bb]), sd_posterior_mean_bb)

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
	OurMethodResult_bb<-data.frame(mean=mean_posterior_mean_bb, median=median_posterior_mean_bb, sd=sd_posterior_mean_bb)		
	LAResult<-data.frame(mean=mean_CoefsMeanLA, median=median_CoefsMeanLA, sd=sd_CoefsMeanLA)
	StepAICMethodResult<-data.frame(mean=glmStepAIC_mean, median=glmStepAIC_median, sd=glmStepAIC_sd_iter)		
	StepBICMethodResult<-data.frame(mean=glmStepBIC_mean, median=glmStepBIC_median, sd=glmStepBIC_sd_iter)
	rownames(OurMethodResult)<-CoefNames		
	rownames(OurMethodResult_bb)<-CoefNames		
	if(Regression=="Normal" || Regression=="Binomial"){rownames(LAResult)<-CoefNames}
	if(Regression=="Quantile"){rownames(StepAICMethodResult)<-CoefNames[-1]; rownames(StepBICMethodResult)<-CoefNames[-1]}else{
					rownames(StepAICMethodResult)<-CoefNames; rownames(StepBICMethodResult)<-CoefNames}

	aux<-list(DataList=DataList, coefs_PostMean_matrix_OurMethod=CoefsMeanOurMethod_matrix,	
	True_Model= TruePredFormula, "Selected_Model_Count_OurMethod_womack"=SelectedModelCountOurMethod_summary,
	"Selected_Model_Count_OurMethod_bb"=SelectedModelCountOurMethod_bb_summary, 
	"Selected_Model_Count_LA"=SelectedModelCountLA_summary,
	"Selected_Model_Count_StepAIC"=SelectedModelCountStepAIC_summary,
	"Selected_Model_Count_StepBIC"=SelectedModelCountStepBIC_summary,
	Model_Exploration_OurMethod_womack=SelectedModelsCountsOurMethod_rbind, 
	Model_Exploration_OurMethod_bb=SelectedModelsCountsOurMethod_bb_rbind,
	Model_Exploration_LA=SelectedModelsCountsLA_rbind,
	Model_Exploration_StepAIC=SelectedModelsCountsStepAIC_rbind,
	Model_Exploration_StepBIC=SelectedModelsCountsStepBIC_rbind,
	"OurMethod_Posterior_Summary_womack"=OurMethodResult, "OurMethod_Posterior_Summary_bb"=OurMethodResult_bb,
	"LA_Posterior_Summary"=LAResult,
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
N<-200
r_beta<-c(1, 0, 2, 0, 3, 2)
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, r_sigma2=r_sigma2)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results

## Scenario 2 (Supplementary material, for Table 1)
set.seed(31415)
N<-200
#r_beta<-c( sample(c(1.5,-1,2,-2.5), size=70, replace=TRUE), rep(0,30)  )
r_beta<-c( sample(c(1.5,-1,2,-2.5), size=100, replace=TRUE) )
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, r_sigma2=r_sigma2)
Summary_SimStudy(aux, ExploredModels=TRUE)	#Summary results

## Scenario 3 (Supplementary material, for Table 2)
set.seed(31415)
N<-200
r_beta<-c( rep(c(0.1,0.2,0.3,0.4),5), sample(c(1.5,-1,2,-2.5), size=30, replace=TRUE) )
r_p<-length(r_beta)
aux<-General_Sim_OurMethod(N, r_beta, R, nchain=nchain, burnin=burnin, Regression=Regression, r_sigma2=r_sigma2)
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

