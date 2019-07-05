###############################
#### Simulating phenotypes ####
###############################

A few simple function for simulating phenotypes, and corresponding tools to analyse them. Any feedback on the usability and potential errors are welcome!

In case you use whole or partial codes from this repository, please refer to it!

#################
### Functions ###
#################

# eval.sim

Evaluate phenotypes simulated with 'sim.phen'. Produces a table in which linear regression modelling results are shown, with details of the simulation process generating analysed data.

# examples

Examples how to use 'sim.phen' and 'eval.sim':
 1) two simulations with different population sizes, each of which with different associations between X1 and Y.
 2) 100 simulations: 50 with 10 individuals, 50 with 100 individuals in the population. Association between X1 and Y is random; assessing the evaluations the number of models producing Type I error is shown.

# resamp

Randomly sample from a data table, creating a random subset of the original data. 

# run.sim

Run several simulations, based on the input list. Unlike in 'sim.phen', input list elements (for details see 'sim.phen') can have multiple values (e.g. NI = c(10, 20)), and in the simulations each and every combinations of element values will be used in separate simulations, i.e. the total number of phenotype simulations is equal to the product of list element lengths. 

# sim.phen

Simulate phenotypes, based on specifications in the input list:

 - NI # numbers of individuals in the population
 - NR # numbers of repeated measurements from each individual
 - XA # population-average of intercept
 - B1 # population-average slope
 - X1_mu # population-average of X1
 - X1_V # population-level variance in trait 'X1'
 - VI # random intercept variance
 - VS1 # random slope variance
 - VR # residual variance
 - V.method # method for approximating var-cov matrix in lmer ('statter' or 'kenw')

Each list element can have only one value associated with it. In repeated measurement scenarios (NR > 1) X1-values corresponding to separate individuals might be static (i.e. one X1 value over all repeated measures: dyn_X1=FALSE) or changing (dyn_X1=TRUE). For now, in changing X1 resampling of within-individual X1 values is random.




