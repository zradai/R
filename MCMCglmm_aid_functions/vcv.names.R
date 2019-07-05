# acquire response and random variable, names for 'corr.GR'
vcv.names<-function(model){
  
  if(class(model)!="MCMCglmm"){stop("This function requires an object of 'MCMCglmm' class!")}
  
  responses<-all.vars(as.formula(model$Fixed)[[2]])
  randoms<-all.vars(as.formula(model$Random))
  randoms<-randoms[randoms!="trait"]
  
  vcv.n<-list(responses=responses, randoms=randoms)
  return(vcv.n)
  
}
