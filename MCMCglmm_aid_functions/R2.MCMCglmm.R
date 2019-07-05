# R^2 for MCMCglmm
## so far only for simple random intercept models
## (based on Nakagawa's solution: https://www.researchgate.net/post/How_can_I_calculate_R2_for_an_Bayesian_MCMC_multilevel_model)
R2.MCMCglmm<-function(model, R2.distribution=FALSE, check.MwoR=FALSE){
  
  if(class(model)!="MCMCglmm"){stop("This function requires an object of 'MCMCglmm' class!")}
  
  require("coda")
  
  # in order not to mix solutions of random effects with the fixed effects
  # matrix of fixed effects solutions is indexed by the number of FEs
  Nfixef<-model$Fixed$nfl 
  SOLs<-model$Sol[,1:Nfixef]
  
  if(R2.distribution){
    SampSize<-nrow(SOLs)
    VFix.Alt<-c()
    for(i in 1:SampSize){
      V.i<-var(as.vector(SOLs[i,]%*%t(model$X)))
      VFix.Alt<-c(VFix.Alt, V.i)
    }
    
    Vunits<-grep("units",colnames(model$VCV))
    R2.marg00.post<-VFix.Alt/(VFix.Alt+apply(as.matrix(model$VCV[,Vunits]), 1, sum))
    R2.marg.post<-VFix.Alt/(VFix.Alt+apply(model$VCV, 1, sum))
    R2.cond.post<-(VFix.Alt+apply(as.matrix(model$VCV[,-Vunits]), 1, sum))/(VFix.Alt+apply(model$VCV, 1, sum))
    
    R2.marg00<-posterior.mode(as.mcmc(R2.marg00.post))
    R2.marg<-posterior.mode(as.mcmc(R2.marg.post))
    R2.cond<-posterior.mode(as.mcmc(R2.cond.post))
    
    R2m00HPD.l<-HPDinterval(as.mcmc(R2.marg00.post))[1]
    R2m00HPD.u<-HPDinterval(as.mcmc(R2.marg00.post))[2]
    R2mHPD.l<-HPDinterval(as.mcmc(R2.marg.post))[1]
    R2mHPD.u<-HPDinterval(as.mcmc(R2.marg.post))[2]
    R2cHPD.l<-HPDinterval(as.mcmc(R2.cond.post))[1]
    R2cHPD.u<-HPDinterval(as.mcmc(R2.cond.post))[2]
  } else {
    FixMeans<-0
    for(c in 2:length(colnames(SOLs))){
      FixMeans<-FixMeans+mean(SOLs[,c])*model$X[,c]
    }
    
    VFix<-var(as.vector(apply(SOLs, 2, mean)%*%t(model$X)))
    Vunits<-grep("units",colnames(model$VCV))
    R2.marg00<-VFix/(VFix+sum(apply(as.matrix(model$VCV[,Vunits]), 2, mean))) # using only fixed effect variances
    R2.marg<-VFix/(VFix+sum(apply(model$VCV, 2, mean))) # random effect variances are included in the denominator
    R2.cond<-(VFix+sum(apply(as.matrix(model$VCV[,-Vunits]), 2, mean)))/(VFix+sum(apply(model$VCV, 2, mean)))
    
    R2m00HPD.l<-NA
    R2m00HPD.u<-NA
    R2mHPD.l<-NA
    R2mHPD.u<-NA
    R2cHPD.l<-NA
    R2cHPD.u<-NA
  }
  
  if(check.MwoR){
    R2s<-data.frame(type=c("Marginal.woR", "Marginal", "Conditional"), 
                    R.squared=c(R2.marg00, R2.marg, R2.cond), 
                    HPD.lower=c(R2m00HPD.l, R2mHPD.l, R2cHPD.l), 
                    HPD.upper=c(R2m00HPD.u, R2mHPD.u, R2cHPD.u))
  } else {
    R2s<-data.frame(type=c("Marginal", "Conditional"), 
                    R.squared=c(R2.marg, R2.cond), 
                    HPD.lower=c(R2mHPD.l, R2cHPD.l), 
                    HPD.upper=c(R2mHPD.u, R2cHPD.u))
  }
  
  return(R2s)
  
}
