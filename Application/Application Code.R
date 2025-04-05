rm(list=ls(all=TRUE))

setwd("")
source("function.R")
library(dplyr)
library(correlationfunnel)
library(ggplot2)
library(ggpubr)
##################################################################################################################################################

################################################################################
##				Some functions, Just run this					##
################################################################################

## QR Application: Latex Summary table per quantile
QR_LatexTable<-function(fit, covaribles, r_alpha)
{
	Latex<-data.frame(Predictor=c("Variance", "Intercept", colnames(covaribles)))
	for(i in 1:length(r_alpha))
	{
		aa<-summary_gibbs(fit[[i]])
		Mean<-(aa$Mean_IC)$Mean
		Latex<-cbind(Latex, Mean)
	}
	colnames(Latex)[-1]<-r_alpha
	Latex
}

## QR Application: Credible interval plot
QR_CIPlot<-function(fit, Covariates, r_alpha)
{
	CI<-list()
	DF<-rbind()
	CoefNames<-c("Variance", "Intercept", colnames(Covariates))
	#r_alpha<-as.numeric(colnames(LatexTable)[-1])
	for(i in 1:length(r_alpha))
	{
		aa<-summary_gibbs(fit[[i]])
		Mean<-(aa$Mean_IC)$Mean
		for(j in 1:length(CoefNames))
		{
			DF<-rbind()
			if(i==1){ DF<-rbind(aa$Mean_IC[j,]) }
			if(i>1){ DF<-rbind(CI[[j]], aa$Mean_IC[j,]) }
			CI[[j]]<-DF
		}
	}
	
	for(j in 1:length(CoefNames))
	{
		CI[[j]]<-cbind("alpha"=r_alpha, CI[[j]])
		CI[[j]]<-data.frame(CI[[j]])
		colnames(CI[[j]])[c(3,4)]<-c("quantile.2.5.", "quantile.97.5.")
	}
	names(CI)<-CoefNames
	
	My_Theme = 	
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"),
	plot.title = element_text(size = 38, family="mono"),
	axis.title.x = element_text(size = 30),
	axis.text.x = element_text(size = 14),
	axis.text.y = element_text(size = 14),
	axis.title.y = element_text(size = 30))
	
	CIPlot<-list()
	for(j in 1:length(CoefNames))
	{
		if(j==1 | j==5)
		{
			CIPlot[[j]]<-ggplot(CI[[j]], aes(alpha)) + 
			geom_line(aes(y=Mean), colour="#69b3a2") + 
			geom_point(aes(y=Mean), shape=21, color="black", fill="#69b3a2", size=2) +
			geom_ribbon(aes(ymin=quantile.2.5., ymax=quantile.97.5.), alpha=0.1) +
			labs(x = "", y="Estimate", title=CoefNames[j]) + My_Theme
		}
		if(j==9)
		{
			CIPlot[[j]]<-ggplot(CI[[j]], aes(alpha)) + 
			geom_line(aes(y=Mean), colour="#69b3a2") + 
			geom_point(aes(y=Mean), shape=21, color="black", fill="#69b3a2", size=2) +
			geom_ribbon(aes(ymin=quantile.2.5., ymax=quantile.97.5.), alpha=0.1) +
			labs(x = "Quantile", y="Estimate", title=CoefNames[j]) +
			My_Theme
		}
		if(j==10 |j ==11 | j==12)
		{
			CIPlot[[j]]<-ggplot(CI[[j]], aes(alpha)) + 
			geom_line(aes(y=Mean), colour="#69b3a2") + 
			geom_point(aes(y=Mean), shape=21, color="black", fill="#69b3a2", size=2) +
			geom_ribbon(aes(ymin=quantile.2.5., ymax=quantile.97.5.), alpha=0.1) +
			labs(x = "Quantile", y="", title=CoefNames[j]) +
			My_Theme
		}
		if(j==2 | j==3 | j==4 | j==6 | j==7 |j==8)
		{	
			CIPlot[[j]]<-ggplot(CI[[j]], aes(alpha)) + 
			geom_line(aes(y=Mean), colour="#69b3a2") + 
			geom_point(aes(y=Mean), shape=21, color="black", fill="#69b3a2", size=2) +
			geom_ribbon(aes(ymin=quantile.2.5., ymax=quantile.97.5.), alpha=0.1) +
			labs(x = "", y="", title=CoefNames[j]) +
			My_Theme
		}
	}
	ggarrange(CIPlot[[1]], CIPlot[[2]], CIPlot[[3]], CIPlot[[4]], CIPlot[[5]], CIPlot[[6]], CIPlot[[7]],
	CIPlot[[8]], CIPlot[[9]], CIPlot[[10]], CIPlot[[11]], CIPlot[[12]])
}


## QR Application: Predictor plot
QR_PredPlot<-function(base, Latextable)
{	
	response<-base[,1]; LatexUsed<-Latextable
	if(colnames(base)[1]=="pas"){responseName<-"Systolic pressure (mmHg)"}else{responseName<-"Diastolic pressure (mmHg)"}
	base_used<-base[,-1]
	max_predictor_index<-3
	PredPlot<-list()
	
	My_Theme = 	
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"),
	plot.title = element_text(size = 28),
	axis.title.x = element_text(size = 32, family="mono"),
	axis.text.x = element_text(size = 20),
	axis.text.y = element_text(size = 20),
	axis.title.y = element_text(size = 28))
	
	for(predictor_index in 1:max_predictor_index )
	{
		predictor<-base_used[,predictor_index]
		pred_name<-colnames(base_used)[predictor_index]
		if(predictor_index==1){ gg_object<-ggplot(data=base_used, aes(x=age, y=response)) +geom_point(aes(color =response)) +labs(x="age", title="", y=responseName) +
		scale_color_gradient(low = "grey", high = "#03BAE5", guide="none") }
	
		if(predictor_index==2){ gg_object<-ggplot(data=base_used, aes(x=waist, y=response)) +geom_point(aes(color =response)) +labs(x="waist", y="") +
		scale_color_gradient(low = "grey", high = "#03BAE5", guide="none") }
	
		if(predictor_index==3){ gg_object<-ggplot(data=base_used, aes(x=bmi, y=response)) +geom_point(aes(color =response)) +labs(x="bmi", y="") +
		scale_color_gradient(low = "grey", high = "#03BAE5", name="Pressure")  }
	
		for(alphaIndex in 2:(length(r_alpha) +1) )
		{
			Posterior_coef<-LatexUsed[,alphaIndex]	
			Posterior_coef<-Posterior_coef[c(-1,-(predictor_index +2))] 
			#predFixed<-base_used[1,-predictor_index]				
			predFixed<-round(apply(base_used[,-predictor_index], 2, mean))	
			predFixed<-c(1,as.numeric(predFixed))	
			FixedPredCoef<-as.vector(t(as.matrix(predFixed))%*%as.matrix(Posterior_coef))
			betaPred<-LatexUsed[,alphaIndex][predictor_index +2]	
			if(betaPred==0){color<-"green"}else{color<-"red"}
	
			predictor_label<- seq(min(predictor), max(predictor), by=0.1)
			Pred<-cbind(predictor_label)
			#head(Pred)
		
			Pred<-data.frame(cbind(Pred, pred=betaPred*predictor_label+ FixedPredCoef))
			if(betaPred!=0){ gg_object<- gg_object +geom_line(data=Pred,aes(x=predictor_label,y=pred), size=1.5, color="black")+ My_Theme }
			if(betaPred==0){ gg_object<- gg_object +geom_line(data=Pred,aes(x=predictor_label,y=pred), size=1.5, color="black", linetype = "longdash")+ My_Theme }
		}
		PredPlot[[predictor_index]]<-gg_object
	} 
	ggarrange(PredPlot[[1]], PredPlot[[2]], PredPlot[[3]], ncol=3 )
}
##################################################################################################################################################
set.seed(31415)


##################################################
## 		Loading data-base				##
##################################################

base<-read.csv("ens.csv")
head(base)

## pas base
base_pas<-data.frame(pas=base$pas,base[,c(-1,-2)])
#head(base_pas)

## pad base
base_pad<-data.frame(pad=base$pad,base[,c(-1,-2)])
#head(base_pad)




##################################################
## 		Fitting BMS for pad base		##
##		    for each quantile			##
##################################################
y<-base_pad[,1]
Covariates<-base_pad[,-1]
nchain<-10000; burnin<-2000
r_alpha<-c(0.01, seq(0.05,0.95, 0.10))
fit_Pad<-list()
for(i in 1:length(r_alpha))
{
	cat("  Quantile", r_alpha[i], "\r")
	fit_Pad[[i]]<-gibbs_abms(y=y, Covariates=Covariates, family="QR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, alpha=r_alpha[i],
			  a0=1, b0=1, count.iteration=FALSE )
}

##################################################
## 		Fitting BMS for pas base		##
##		    for each quantile			##
##################################################
y<-base_pas[,1]
Covariates<-base_pas[,-1]
nchain<-10000; burnin<-2000
r_alpha<-c(0.01, seq(0.05,0.95, 0.10))
fit_Pas<-list()
for(i in 1:length(r_alpha))
{
	cat("  Quantile", r_alpha[i], "\r")
	fit_Pas[[i]]<-gibbs_abms(y=y, Covariates=Covariates, family="QR", first_excluded=0, nchain=nchain, burnin=burnin, tau2=1000, rho=1, alpha=r_alpha[i],
			  a0=1, b0=1, count.iteration=FALSE )
}


########################################################
## 		Pas and Pad summary Tables			##
########################################################

## Pad, Manuscript Table 6
( LatexTablePad<-QR_LatexTable(fit_Pad, Covariates, r_alpha) )

## Pas, Manuscript Table 14
( LatexTablePas<-QR_LatexTable(fit_Pas, Covariates, r_alpha) )


########################################################
## 		Pas and Pad Credible Plots			##
########################################################

## Pad, Figure 1 Manuscript
QR_CIPlot(fit_Pad, Covariates, r_alpha)

## Pas, Figure 2, Supplementary material
QR_CIPlot(fit_Pas, Covariates, r_alpha)


########################################################
## 		Pas and Pad Prediction Plots			##
########################################################

## Pad, Figure 1, Supplementary material
QR_PredPlot(base_pad, LatexTablePad)

## Pas, Figure 3, Supplementary material
QR_PredPlot(base_pas, LatexTablePas)


