# BMS
The following files are as follows:

1. "function.R": A R script that contains all the necessary functions to perform estimation and model selection in Linear, Logistic, Negative-Binomial, Quantile, and Skew-Normal regression.

2. "Model Ilustrations.R": A R script that can be used to illustrate all five regression models. “summary_gibbs” function provides the posterior mean of parameters and quantile 2.5% and 97.5%, alongside the explored models with their respective proportion of times that was selected.

3. "Application Code.R": A R script where the application results can be replicated For the Application reproducibility. To do that:
	a) Load de database, 
	b) Fit the models
	c) Run Table estimation, credible interval plot, and prediction plot to replicate the Manuscript tables and plots

4. "ens.csv": A "csv" file that contains the data set used in the "Application Code.R" file

5. “ens_description.pdf”: Data dictionary for the "ens" data-set

6. "Simulation Study Without Bayes Factor.R": A R script that replicates the simulation study WITHOUT the Bayes factor computation

7. "Simulation Study With Bayes Factor.R": A R script that replicates the simulation study WITH the Bayes factor computation
