# BMS
The following files are as follows:

These files contain the necessary functions to perform a Bayesian model selection methodology based on the spike-and-slab strategy and an augmentation technique for Linear, Logistic, Negative Binomial, Quantile, and Skew Normal Regression. The model considers a response vector "y" of size "n" and "p" predictors to perform coefficient estimation and asses which ones are relevant to explain the response distribution.

1. "function.R": A R script containing all the necessary functions to perform estimation and model selection in Linear, Logistic, Negative-Binomial, Quantile, and Skew-Normal regression.

2. "Model Ilustrations.R": A R script that can be used to illustrate all five regression models. “summary_gibbs” function provides the posterior mean of parameters and quantile 2.5% and 97.5%, alongside the explored models with their respective proportion of times that was selected.

3. "Application Code.R": A R script where the application results can be replicated For the Application reproducibility. To do that:
	a) Load the database, 
	b) Fit the models
	c) Run Table estimation, credible interval plot, and prediction plot to replicate the Manuscript tables and plots

4. "ens.csv": A "csv" file that contains the data set used in the "Application Code.R" file

5. “ens_description.pdf”: Data dictionary for the "ens" data-set

6. "Simulation Study With Beta-Binom.R": A R script that performs a simulation study where we also compute our Bayesian model selection method with the Beta-Binomial prior

7. "Simulation Study With Beta-Binom.R": A R script that replicates the simulation study where we only consider the Womack prior for our Bayesian model selection method

8. "Comparing BF between models.R": A R script that computes the Bayes factor between the true model and two nested models for the Linear, Logistic, Negative-Binomial, Quantile, and Skew-Normal regression models
