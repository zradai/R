
#######################################
#### Aiding functions for MCMCglmm ####
#######################################

A handful of simple functions to aid handling and interpreting models fitted by MCMCglmm (Hadfield, 2010). Some of the functions are not entirely done yet, so not every kind of models might be used with them (e.g. complex multi-response models, multiple interactions, composite variance structures). Any feedback on the usability and potential errors are very welcome!

In case you use whole or partial codes from this repository, please refer to it!


#################
### Functions ###
#################


NOTE: some functions use the function 'split.direct.sum' from the 'MCMCglmm' package:
https://github.com/cran/MCMCglmm/blob/master/R/split.direct.sum.R


# R2.MCMCglmm

Estimates R-squared for an MCMCglmm model, based on the description given by Shinichi Nakagawa; estimates either a single R^2 value (R2.distribution=FALSE) or the posterior distribution of R^2 (R2.distribution=TRUE).
(https://www.researchgate.net/post/How_can_I_calculate_R2_for_an_Bayesian_MCMC_multilevel_model)

# contr.MCMCglmm

Gets factor level contrasts, i.e. posterior distribution of factor level differences.

# corr.GR

Creates random effects and residual correlation matrices from an MCMCglmm model.

# get.BLUPs

Extracts posterior distributions for random intercepts and slopes; might be used to acquire Best Linear Unbiased Predictors (BLUPs).

# get.dims

Gets dimensions for the prior matrices based on model formula.

# get.prior

Gets prior list for a model formula, that can be used in an MCMCglmm as prior. 

# get.sign

Estimates phylogenetic signal (heritability) from an MCMCglmm model using phylogeny (pedigree) in the random variance structure. 

# pMCMC

Calculates P-value for the posterior distribution of parameter estimates (using the Markov chain), based on the proportion of Markov chain samples crossing zero.

# vcov.GR

Creates random effects and residual co-variance matrices from an MCMCglmm model.

# vcv.names

Gets names of variance elements (used by other functions).

############################

Hadfield, J. D. (2010). MCMC methods for multi-response generalized linear mixed models: the MCMCglmm R package. Journal of Statistical Software, 33(2), 1-22.

