library(ggplot2)	#for ggplot
library(viridis)	#for ggplot
library(hrbrthemes)	#for ggplot


## Models that will be compared. BF between Models[[1]] and the others
Models<-list()
Models[[1]]<-c(2,4,5)
Models[[2]]<-c(2,4)
Models[[3]]<-2

N_vec<-c(200, 500, 1000)	#For which sample size we want to compute the BF
first_excluded<-0
nchain=10000; burnin=2000


############################################
## 			LiR 				##
############################################
set.seed(31415)

r_beta<-c(1,0,2,0,3,2)
r_sigma2<-2
BF_summary_LiR<-vector(mode="list", length=length(N_vec))
for(i in 1:length(N_vec))
{
	N<-N_vec[i]
	ni<-rep(1, N);
	aux<-LiR_base(N, r_beta, r_sigma2)
	BF_summary_LiR[[i]]<-BF_model(Models=Models, y=aux$y, Covariates=aux$Covariates, family="LiR")
}



############################################
## 			LoR 				##
############################################
set.seed(31415)

r_beta<-c(1,0,2,0,3,2)
BF_summary_LoR<-vector(mode="list", length=length(N_vec))
for(i in 1:length(N_vec))
{
	N<-N_vec[i]
	ni<-rep(1, N);
	aux<-LoR_base(N=N, r_beta=r_beta, ni=ni)
	BF_summary_LoR[[i]]<-BF_model(Models=Models, y=aux$y, Covariates=aux$Covariates, ni=aux$ni, family="LoR")
}


############################################
## 			QR 				##
############################################
set.seed(31415)

r_beta<-c(1,0,2,0,3,2)
r_sigma2<-2
r_alpha<-0.1
BF_summary_QR_01<-vector(mode="list", length=length(N_vec))
BF_summary_QR_05<-vector(mode="list", length=length(N_vec))
BF_summary_QR_09<-vector(mode="list", length=length(N_vec))
for(i in 1:length(N_vec))
{
	N<-N_vec[i]
	aux<-QR_base(N, r_beta, r_sigma2, r_alpha=0.1)
	BF_summary_QR_01[[i]]<-BF_model(Models=Models, y=aux$y, Covariates=aux$Covariates, family="QR", alpha=0.1)

	aux<-QR_base(N, r_beta, r_sigma2, r_alpha=0.5)
	BF_summary_QR_05[[i]]<-BF_model(Models=Models, y=aux$y, Covariates=aux$Covariates, family="QR", alpha=0.5)

	aux<-QR_base(N, r_beta, r_sigma2, r_alpha=0.9)
	BF_summary_QR_09[[i]]<-BF_model(Models=Models, y=aux$y, Covariates=aux$Covariates, family="QR", alpha=0.9)
}


############################################
## 			SNR 				##
############################################
set.seed(31415)

r_beta<-c(1,0,2,0,3,2)
r_sigma2<-1.5; r_lambda<-4
BF_summary_SNR<-vector(mode="list", length=length(N_vec))
for(i in 1:length(N_vec))
{
	N<-N_vec[i]
	ni<-rep(1, N);
	aux<-SNR_base(N, r_beta, r_sigma2, r_lambda)
	BF_summary_SNR[[i]]<-BF_model(Models=Models, y=aux$y, Covariates=aux$Covariates, family="SNR")
}


############################################
## 			NBR 				##
############################################
set.seed(31415)
Models<-list()
Models[[1]]<-c(1,2,4,5)
Models[[2]]<-c(1,5)
Models[[3]]<-c(1)


r_beta<-c(0.5, -0.8,  1.0,  0,  0.4, -0.7)
r_r<-2
BF_summary_NBR<-vector(mode="list", length=length(N_vec))
for(i in 1:length(N_vec))
{
	N<-N_vec[i]
	ni<-rep(1, N);
	aux<-NBR_base(N, r_beta, r_r)
	BF_summary_NBR[[i]]<-BF_model(Models=Models, y=aux$y, Covariates=aux$Covariates, family="NBR")
}


#setwd("C:/Users/asus/Dropbox/BMS, Bayesian Analysis Journal/R/For R package")
setwd("C:/Users/Alumno/Dropbox/BMS, Bayesian Analysis Journal/R/For R package")
load("BF_evolution_AllRegression.Rdata")



############################################
## 			gg plot 			##
############################################


## Generating data base (using only quantile 0.1 for QR)
#for SNR, the sqrt of log BF was computed. Otherwise is to largue
names_N_vec<-c()
for(i in 1:length(N_vec)){ names_N_vec<-c(names_N_vec, paste0("n=",N_vec[i])) }

BF_gamma2<-as.data.frame(matrix(NA, nrow=5*length(N_vec), ncol=3))
BF_gamma3<-as.data.frame(matrix(NA, nrow=5*length(N_vec), ncol=3))
colnames(BF_gamma2)<-c("n", "log_BF", "Regression")
colnames(BF_gamma3)<-c("n", "log_BF", "Regression")
for(i in 1:length(N_vec))
{
	
	BF_gamma2[i,2]<-(BF_summary_LiR[[i]]$Under_respective)$log_Marginal_BF_Estimator[2] 
	BF_gamma3[i,2]<-(BF_summary_LiR[[i]]$Under_respective)$log_Marginal_BF_Estimator[3] 

	BF_gamma2[3 +i,2]<-(BF_summary_LoR[[i]]$Under_respective)$log_Marginal_BF_Estimator[2]
	BF_gamma3[3 +i,2]<-(BF_summary_LoR[[i]]$Under_respective)$log_Marginal_BF_Estimator[3]

	BF_gamma2[6 +i,2]<-(BF_summary_NBR[[i]]$Under_respective)$log_Marginal_BF_Estimator[2] 
	BF_gamma3[6 +i,2]<-(BF_summary_NBR[[i]]$Under_respective)$log_Marginal_BF_Estimator[3] 

	BF_gamma2[9 +i,2]<-(BF_summary_QR_01[[i]]$Under_respective)$log_Marginal_BF_Estimator[2]
	BF_gamma3[9 +i,2]<-(BF_summary_QR_01[[i]]$Under_respective)$log_Marginal_BF_Estimator[3]

	BF_gamma2[12 +i,2]<-sqrt( (BF_summary_SNR[[i]]$Under_respective)$log_Marginal_BF_Estimator[2] )
	BF_gamma3[12 +i,2]<-sqrt( (BF_summary_SNR[[i]]$Under_respective)$log_Marginal_BF_Estimator[3] )
}
BF_gamma3[,1]<-rep(N_vec, 5)
BF_gamma2[,1]<-rep(N_vec, 5)
BF_gamma2[,3]<-c( rep("LiR", 3), rep("LoR", 3), rep("NBR", 3), rep("QR_0.1", 3), rep("SNR", 3) )
BF_gamma3[,3]<-c( rep("LiR", 3), rep("LoR", 3), rep("NBR", 3), rep("QR_0.1", 3), rep("SNR", 3) )



## Computing ggplot, setting theme
My_Theme = 	
theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey"),
plot.title = element_text(size = 38, family="mono"),
axis.title.x = element_text(size = 30),
axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
axis.title.y = element_text(size = 30),
legend.title=element_text(size=20),
legend.text=element_text(size=17) )

## ggplot of BF between model1 (true model) and model 2
ggplot(data=BF_gamma2, aes(x=n, y=log_BF, group=Regression, color=Regression)) +
geom_line() +
geom_point(aes(shape=Regression), size=5) +
#ggtitle("BF between model 1 and model 2") +
ylab("Log Bayes Factor") +
xlab("Sample size") +
My_Theme 

## ggplot of BF between model1 (true model) and model 3
ggplot(data=BF_gamma3, aes(x=n, y=log_BF, group=Regression, color=Regression)) +
geom_line() +
geom_point(aes(shape=Regression), size=5) +
#ggtitle("BF between model 1 and model 2") +
ylab("Log Bayes Factor") +
xlab("Sample size") +
My_Theme 


