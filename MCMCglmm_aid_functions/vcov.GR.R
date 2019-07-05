
##############################################################
###----- Get variance-covariance matrices from models -----###
##############################################################

vcov.GR<-function(model){
  
  if(class(model)!="MCMCglmm"){stop("This function requires an object of 'MCMCglmm' class!")}
  
  #--- random elements ---
  
  resps<-all.vars(model$Fixed$formula[[2]])
  preds<-all.vars(model$Fixed$formula[[3]])
  preds<-preds[preds!="trait"]
  rand.elements<-split.direct.sum(as.character(model$Random$formula[2]))
  rndM.list<-list()
  for(e in rand.elements){
    rnd.form.e<-as.formula(paste("~",e,sep=""))
    rnd.vars.e<-all.vars(rnd.form.e)
    rnd.vars.e<-rnd.vars.e[rnd.vars.e!="units" & rnd.vars.e!="trait"]
    rnd.only.e<-setdiff(rnd.vars.e, preds)
    pos.rnd<-grep(rnd.only.e, colnames(model$VCV))
    vcv.rnd<-posterior.mode(model$VCV)[pos.rnd]
    
    #--- get vcv component names ---
    vcv.nms.e<-names(vcv.rnd)
    nms.Me<-cutoff.names(vcv.nms.e, rnd.only.e)
    nms.Me<-cutoff.trait(nms.Me[[1]])
    
    #--- build matrices ---
    if(substr(e,1,3)=="idh"){
      if(length(vcv.rnd)==1){
        rndM.e<-as.matrix(vcv.rnd)
        colnames(rndM.e)<-"(Intercept)"
        rownames(rndM.e)<-"(Intercept)"
        }else{
        rndM.e<-diag(vcv.rnd)
        colnames(rndM.e)<-nms.Me
        rownames(rndM.e)<-nms.Me
      }
    }
    
    if(substr(e,1,3)=="us("){
      rndM.e<-matrix(vcv.rnd, ncol=sqrt(length(vcv.rnd)))
      
      if(any(grepl("Intercept", nms.Me))){
        nms.to.use<-nms.Me[1:sqrt(length(nms.Me))]
        ntu.list<-strsplit(nms.to.use, "\\(Intercept)")
        ntu.vector<-c("(Intercept)")
        
        if(length(ntu.list)<2){
          colnames(rndM.e)<-"(Intercept)"
          rownames(rndM.e)<-"(Intercept)"
        }else{
          for(n in 2:length(ntu.list)){
            ntu.vector<-c(ntu.vector, substr(ntu.list[[n]], 1, (nchar(ntu.list[[n]])-1)))
          }
          colnames(rndM.e)<-ntu.vector
          rownames(rndM.e)<-ntu.vector
        }
      }else{
        nms.to.use<-nms.Me[1:sqrt(length(nms.Me))]
        ntu.vector<-c()
        for(n in nms.to.use){
          split.nm<-strsplit(n, "\\:")[[1]]
          ntu.vector<-c(ntu.vector, paste(split.nm[1:(length(split.nm)/2)], collapse=":"))
        }
        colnames(rndM.e)<-ntu.vector
        rownames(rndM.e)<-ntu.vector
      }
      
    }
    
    rndM.list[[rnd.only.e]]<-rndM.e
    
  }
  
  #--- residual structure ---
  
  resid.elements<-split.direct.sum(as.character(model$Residual$formula[2]))
  vcv.all<-colnames(model$VCV)
  resid.nms<-vcv.all[grep("units", vcv.all)]
  vcv.rsd<-posterior.mode(model$VCV[,resid.nms])
  resid.nms<-cutoff.trait(resid.nms)
  resid.nms<-cutoff.units(resid.nms)
  
  for(r in resid.elements){
    if(substr(r,1,3)=="idh"){
      rsdM.r<-diag(vcv.rsd)
      colnames(rsdM.r)<-resid.nms
      rownames(rsdM.r)<-resid.nms
    }
    
    if(substr(r,1,3)=="us("){
      rsdM.r<-matrix(vcv.rsd, ncol=sqrt(length(vcv.rsd)))
      rs.to.use<-resid.nms[1:sqrt(length(resid.nms))]
      rtu.vector<-c()
      for(n in rs.to.use){
        split.nm<-strsplit(n, "\\:")[[1]]
        rtu.vector<-c(rtu.vector, paste(split.nm[1:(length(split.nm)/2)], collapse=":"))
      }
      colnames(rsdM.r)<-rtu.vector
      rownames(rsdM.r)<-rtu.vector
    }
  }
  
  GR.matrices<-list(Random.VCV=rndM.list, 
                    Residual.VCV=rsdM.r)
  return(GR.matrices)
  
}


