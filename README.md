# Bayesian-Stats-Project2022
### *Content of the repository*
1) *Import_data.R* : the script contains functions to read data and to prepare the dataset for the analysis.
2) *Auxiliary_Functions.R* : the script contains several functions which performs repetitive tasks (e.g. computation of quantities for the full conditional, evaluation of some functions). In particular you can find the functions implementing the different versions of the Metropolis-Hastings step.
3) *Metropolis_within_Gibbs.R* : in this script the MCMC simulation is performed using a Metropolis-within-Gibbs algorithm. Different choice for the model of the coefficients phi of the warping function can be made inside the script. Depending on the choice of the model for the coefficients, different Metropolis Hasting algorithms are implemented to update their values.
5) *Results_and_Plots.R* : in this script the results produced in 3) are analyzed and the plots of the significant quantities are made.  
