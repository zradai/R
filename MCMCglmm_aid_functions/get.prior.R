# prepare a non-informative prior based on the output of 'get.dims', OR a formula + data frame
get.prior<-function(dims, random=NULL, data=NULL, type=c("NI", "FI", "NII"), use.data=FALSE){ 
  # NI - non-informative
  # FI - flat-improper
  # NII - non-informative-improper
  
  if(is.null(type) | length(type)>1){type<-"NI"}
  
  if(is.formula(dims)){
    if(is.null(data)){stop("A data frame is needed when formula is specified instead of list!")}
    formula<-as.formula(dims)
    dims<-get.dims(fixed=formula, random=random, data=data)
  }
  
  
  # non-informative
  if(type=="NI"){
    
    pr<-list(B=list(V=diag(dims$Dimensions$Fixed)*1e5, 
                    mu=rep(0, dims$Dimensions$Fixed)),
             R=list(V=diag(dims$Dimensions$Residual),
                    nu=(dims$Dimensions$Residual-1)+0.002))
    
    if(dims$Dimensions$Random>0){
      for(i in 1:dims$Number.of.elements$Random){
        numb.rand.i<-as.character(paste("G",i,sep=""))
        pr$G[[numb.rand.i]]<-list(V=diag(dims$Dimensions$Random),
                                  nu=(dims$Dimensions$Random-1)+0.002)
      }
    }
  }
  
  # flat improper (is sad to be non informative for the mean, but could be informative for the variance! 
  # posterior distribution is equivalent to the likelihood)
  if(type=="FI"){
    pr<-list(B=list(V=diag(dims$Dimensions$Fixed)*1e5, 
                    mu=rep(0, dims$Dimensions$Fixed)),
             R=list(V=diag(dims$Dimensions$Residual),
                    nu=0))
    
    if(dims$Dimensions$Random>0){
      for(i in 1:dims$Number.of.elements$Random){
        numb.rand.i<-as.character(paste("G",i,sep=""))
        pr$G[[numb.rand.i]]<-list(V=diag(dims$Dimensions$Random),
                                  nu=0)
      }
    }
  }
  
  # non-informative improper (marginal posterior mode should coincide with the REML estimator)
  if(type=="NII"){
    pr<-list(B=list(V=diag(dims$Dimensions$Fixed)*1e5, 
                    mu=rep(0, dims$Dimensions$Fixed)),
             R=list(V=diag(dims$Dimensions$Residual)*1e-16,
                    nu=-2))
    
    if(dims$Dimensions$Random>0){
      for(i in 1:dims$Number.of.elements$Random){
        numb.rand.i<-as.character(paste("G",i,sep=""))
        pr$G[[numb.rand.i]]<-list(V=diag(dims$Dimensions$Random)*1e-16,
                                  nu=-2)
      }
    }
  }
  
  return(pr)
  
}
