#################
### Example 1 ###
#################

input.s<-list(  
  NI=rep(c(10, 100), 1), # numbers of individuals in the population
  NR=c(3), # numbers of repeated measurements from each individual
  XA=0, # population-average of intercept
  B1=rep(c(3.1415, -2.72), 1), # population-average slope
  X1_mu=0, # population-average of X1
  X1_V=c(1), # population-level variance in trait 'X1'
  VI=c(1), # random intercept variance
  VS1=c(0), # random slope variance
  VR=c(1), # residual variance
  V.method=c(0) # method for approximating var-cov matrix in lmer ('statter' or 'kenw')
)

sim1 = run.sim(input.s, dyn_X1=FALSE, rescale=FALSE, fix_Y=FALSE) 
eval.sim(errs, plotting=TRUE, rand.slope=FALSE, averaging=FALSE)

#################
### Example 2 ###
#################

input.s2<-list(  
  NI=rep(c(10, 100), 50), # numbers of individuals in the population
  NR=c(3), # numbers of repeated measurements from each individual
  XA=0, # population-average of intercept
  B1=rep(c(0), 1), # population-average slope
  X1_mu=0, # population-average of X1
  X1_V=c(1), # population-level variance in trait 'X1'
  VI=c(1), # random intercept variance
  VS1=c(0), # random slope variance
  VR=c(2), # residual variance
  V.method=c(0) # method for approximating var-cov matrix in lmer ('statter' or 'kenw')
)

sim2 = run.sim(input.s2, dyn_X1=FALSE, rescale=FALSE, fix_Y=FALSE) 
eval = eval.sim(errs, plotting=TRUE, rand.slope=FALSE, averaging=FALSE)

# check number of simulations resulting in Type I error
require("lattice")
eval.ti = eval[eval$t.value>=2 | eval$t.value<=-2,]
table(eval.ti$NI)
bwplot(B1.estimate ~ as.factor(NI), data = eval.ti)

