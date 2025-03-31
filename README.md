# BMS
The following files are as follows:

These files contain the necessary functions to perform a Bayesian model selection methodology based on the spike-and-slab strategy and an augmentation technique for Linear, Logistic, Negative Binomial, Quantile, and Skew Normal Regression. The model considers a response vector "y" of size "n" and "p" predictors to perform coefficient estimation and asses which ones are relevant to explain the response distribution. These files are the extra examples mentioned in the R package "abms", which is contained in CRAN.

1. "function.R": A R script that contains all the necessary functions to perform estimation and model selection in Linear, Logistic, Negative-Binomial, Quantile, and Skew-Normal regression. This functions are in the R package "abms".

2. "Model Ilustrations.R": A R script that can be used to illustrate all five regression models. “summary_gibbs” function provides the posterior mean of parameters and quantile 2.5% and 97.5%, alongside the explored models with their respective proportion of times that was selected.

3. "Simulation Study.R": A R script where the simulation study's tables and figures of the manuscript and supplementary material can be reproduced

4. "Application Code.R": A R script where the application results can be replicated For the Application reproducibility. It is held in the "Application" folder.

5. "ens.csv": A ".csv" file that contains the data set used in the "Application Code.R" file. It is held in the "Application" folder.

6. “ens_description.pdf”: Data dictionary for the "ens" data-set. It is held in the "Application" folder.
