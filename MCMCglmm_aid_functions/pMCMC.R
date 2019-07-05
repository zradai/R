# get a pMCMC value for any posterior distributions
pMCMC<-function(MCMC.chains){ ## one chain as a vector, or matrix (or data frame) of multiple chains
  
  cn.numb<-FALSE
  
  if(is.null(colnames(MCMC.chains))){
    if(is.null(names(MCMC.chains))){
      #warning("No names are given to the different chains! Order of chains will be used instead of names.")
      cn.numb<-TRUE
    }
  } else {
    chain.names<-colnames(MCMC.chains)
  }
  
  if(is.null(ncol(MCMC.chains))){
    MCMC.chains<-matrix(as.vector(MCMC.chains), ncol=1)
    n.chains<-ncol(MCMC.chains)
  } else {
    n.chains<-ncol(MCMC.chains)
  }
  
  if(cn.numb){
    chain.names<-c(1:n.chains)
  }
  
  
  p.mc<-c()
  for(c in 1:n.chains){
    samp.c<-length(MCMC.chains[,c])
    vals.c<-MCMC.chains[,c]
    below.zero<-length(vals.c[vals.c<0])
    above.zero<-length(vals.c[vals.c>0])
    bzf.c<-below.zero/samp.c
    azf.c<-above.zero/samp.c
    
    if(bzf.c<=azf.c){
      p.mc.c<-2*bzf.c
    } else {
      p.mc.c<-2*azf.c
    }
    p.mc<-c(p.mc, p.mc.c)
  }
  
  all.pMCMC<-data.frame(pMCMC=p.mc, variables=chain.names)
  return(all.pMCMC)
  
}
