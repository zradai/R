# get dimensions for the prior matrices
get.dims<-function(fixed, random=NULL, data){
  
  require("stats")
  
  fixed<-as.formula(fixed)
  if(!is.null(random)){random<-as.formula(random)}
  if(!is.null(rcov)){rcov<-as.formula(rcov)}
  
  resps<-all.vars((as.formula(fixed))[[2]])
  preds<-all.vars((as.formula(fixed))[[3]])
  rands<-all.vars(as.formula(random))
  
  N.response<-length(resps)
  
  if(N.response==1){
    mm.d<-model.matrix(fixed, data=data)
    B.dim<-dim(mm.d)[2]
  }
  
  if(N.response>1){
    
    Y<-c()
    traits<-c()
    for(r in resps){
      Y<-c(Y, data[,r])
      traits<-c(traits, rep(r, nrow(data)))
    }
    
    ndf<-data
    for(b in 1:(length(resps)-1) ){
      ndf<-rbind(ndf, data) 
    }
    
    ndf$Y<-Y
    ndf$trait<-traits
    fixed.new<-as.formula(paste("Y ~", as.character(fixed[3])))
    mm.d<-model.matrix(fixed.new, data=ndf)
    B.dim<-dim(mm.d)[2]
    
  }
  
  if(is.null(random)){rndm.dim<-0} else {rndm.dim<-N.response}
  dims<-list(Dimensions=list(Fixed=B.dim, 
                             Random=rndm.dim, 
                             Residual=N.response),
             Number.of.elements=list(Response=N.response,
                                     Predictor=length(preds[preds!="trait"]),
                                     Random=length(rands[rands!="trait"]) ))
  
  return(dims)
  
}
