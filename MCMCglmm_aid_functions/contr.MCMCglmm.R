# get contrasts from an MCMCglmm model 
## one factor predictor usable only so far
## also, categorical interactions are not yet implemented
contr.MCMCglmm<-function(model, data){
  
  if(class(model)!="MCMCglmm"){stop("This function requires an object of 'MCMCglmm' class!")}
  
  require("coda")
  
  resps<-all.vars((as.formula(model$Fixed$formula))[[2]])
  preds<-all.vars((as.formula(model$Fixed$formula))[[3]])
  preds<-preds[preds!="trait"]
  fctrs<-c()
  for(p in preds){
    if(is.factor(data[,p]) | is.character(data[,p])){
      fctrs<-c(fctrs, p)
    }
  }
  
  lvls<-list()
  for(f in 1:length(fctrs)){
    fr<-as.character(fctrs[f])
    lvls[[fr]]<-levels(as.factor(data[,fr]))
  }
  
  if(colnames(model$Sol)[1]=="(Intercept)"){intercept<-TRUE}else{intercept<-FALSE}
  
  if(intercept){
    # get contrasts from a model fitted with intercept
    lvl.chains<-list()
    for(r in resps){
      lvl.chains.r<-data.frame()
      for(f in names(lvls)){
        for(l in lvls[[f]]){
          intrcpt.f<-as.character(lvls[[f]][1])
          if(length(resps)>1){lvl.nm<-paste("trait",r,":",f,l,sep="")} else {lvl.nm<-paste(f,l,sep="")}
          if(l==intrcpt.f){
            lvl.chains.r[1:nrow(model$Sol),lvl.nm]<-model$Sol[,"(Intercept)"]
          } else {
            lvl.chains.r[1:nrow(model$Sol),lvl.nm]<-model$Sol[,"(Intercept)"]+model$Sol[,lvl.nm]
          }
        }
      }
      lvl.chains[[r]]<-lvl.chains.r
    }
    
  } else {
    
    # get contrasts from a model fitted without intercept
    lvl.chains<-list()
    for(r in resps){
      lvl.chains.r<-data.frame()
      for(f in names(lvls)){
        for(l in lvls[[f]]){
          if(length(resps)>1){lvl.nm<-paste("trait",r,":",f,l,sep="")} else {lvl.nm<-paste(f,l,sep="")}
          lvl.chains.r[1:nrow(model$Sol),lvl.nm]<-model$Sol[,lvl.nm]
        }
      }
      lvl.chains[[r]]<-lvl.chains.r
    }
    
  }
  
  all.contrasts<-data.frame(response=0, contrast=0, post.mode=0, pMCMC=0, HPD.lower=0, HPD.upper=0)
  c<-1
  for(r in resps){
    rcols<-ncol(lvl.chains[[r]])
    for(i in 1:rcols){
      for(j in 1:rcols){
        
        if(i!=j){
          contr.ij<-lvl.chains[[r]][,i]-lvl.chains[[r]][,j]
          post.m.cij<-posterior.mode(as.mcmc(contr.ij))
          hpd.l<-as.numeric(HPDinterval(as.mcmc(contr.ij))[1,1])
          hpd.u<-as.numeric(HPDinterval(as.mcmc(contr.ij))[1,2])
          pdf.ij<-pMCMC(contr.ij)
          p.ij<-as.numeric(pdf.ij$pMCMC)
          
          all.contrasts[c,"response"]<-as.character(r)
          all.contrasts[c,"contrast"]<-as.character(paste(colnames(lvl.chains[[r]])[i], colnames(lvl.chains[[r]])[j], sep=" - "))
          all.contrasts[c,"post.mode"]<-post.m.cij
          all.contrasts[c,"pMCMC"]<-p.ij
          all.contrasts[c,"HPD.lower"]<-hpd.l
          all.contrasts[c,"HPD.upper"]<-hpd.u
          c<-c+1
        }
        
      }
    }
  }
  
  return(all.contrasts)
  
}
