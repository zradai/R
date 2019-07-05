
######################################################
###----- Get correlation matrices from models -----###
######################################################

corr.GR<-function(model, p.mcmc=FALSE, type="matrix", verbose=TRUE){
  
  # p.mcmc -- only used when 'type="matrix"'; logical, whether pMCMC values should be saved in separate matrices
  # type=c("matrices", "dataframe")
  
  if(class(model)!="MCMCglmm"){stop("This function requires an object of 'MCMCglmm' class!")}
  if(type!="matrix" & type!="dataframe"){stop("Type of output should be either 'matrix' or 'dataframe'!")}
  if(p.mcmc){
    p.mcmc.list.rnd<-list()
    p.mcmc.list.rs<-list()
    }
  
  RandF<-substr(split.direct.sum(as.character(model$Random$formula[2])), 1, 2)
  ResF<-substr(split.direct.sum(as.character(model$Residual$formula[2])), 1, 2)
  if(!any(grepl("us", c(RandF, ResF)))){stop("At least one variance structure should be fitted with the 'us()' function in the model!")}
  
  #--- random elements ---
  
  resps<-all.vars(model$Fixed$formula[[2]])
  preds<-all.vars(model$Fixed$formula[[3]])
  preds<-preds[preds!="trait"]
  rand.elements<-split.direct.sum(as.character(model$Random$formula[2]))
  rndDF.list<-list()
  rndM.list<-list()
  for(e in rand.elements){
    rnd.form.e<-as.formula(paste("~",e,sep=""))
    rnd.vars.e<-all.vars(rnd.form.e)
    rnd.vars.e<-rnd.vars.e[rnd.vars.e!="units" & rnd.vars.e!="trait"]
    rnd.only.e<-setdiff(rnd.vars.e, preds)
    
    if(length(grep(pattern = rnd.only.e, x = colnames(model$VCV)))<4){
      random.corr<-FALSE
      if(verbose){warning("Less than 4 elements in the random covariance matrix of '", rnd.only.e, "' - excluding from analyses!")}
      } else {
        
        random.corr<-TRUE
        
        pos.rnd<-grep(rnd.only.e, colnames(model$VCV))
        vcv.rnd<-model$VCV[,pos.rnd]
        
        #--- get vcv component names ---
        vcv.nms.e<-colnames(vcv.rnd)
        vcv.nms.e2<-gsub(paste(".",rnd.only.e, sep=""), "", vcv.nms.e)
        nms.Me<-cutoff.names(vcv.nms.e, rnd.only.e)
        nms.Me<-cutoff.trait(nms.Me[[1]])
        
        #--- get correlations ---
        corrs.rnd.e<-data.frame()
        for(n in vcv.nms.e2){
          
          vn.n<-strsplit(n, "\\:")[[1]]
          
          if(any(grepl("(Intercept)",vn.n))){
            if(grep("(Intercept)",vn.n)[1]==1){
              vn.n.1<-"(Intercept)"
              vn.n.2<-gsub(pattern = "\\(Intercept):", replacement = "", x = n)
            }else{
              vn.n.1<-gsub(pattern = "\\:\\(Intercept)", replacement = "", x = n)
              vn.n.2<-"(Intercept)"
            }
          }else{
            vn.n.1<-paste(vn.n[1:(length(vn.n)/2)], collapse=":")
            vn.n.2<-paste(vn.n[(1+length(vn.n)/2):length(vn.n)], collapse=":")
          }
          if(vn.n.1==vn.n.2){
            next
          }else{
            cov.name<-paste(paste(vn.n, collapse=":"), rnd.only.e, sep=".")
            vv1<-paste(paste(rep(vn.n.1, 2), collapse=":"), rnd.only.e, sep=".")
            vv2<-paste(paste(rep(vn.n.2, 2), collapse=":"), rnd.only.e, sep=".")
            
            cov.n<-vcv.rnd[,cov.name]
            sqrt.vv.n<-sqrt(vcv.rnd[,vv1]*vcv.rnd[,vv2])
            corr.post<-cov.n/sqrt.vv.n
            rho.n<-posterior.mode(corr.post)
            hpd.l<-HPDinterval(corr.post)[1]
            hpd.u<-HPDinterval(corr.post)[2]
            p.rho<-pMCMC(corr.post)[1,1]
            
            cr<-nrow(corrs.rnd.e)+1
            corrs.rnd.e[cr,"Element.1"]<-gsub(pattern = "trait", replacement = "", vn.n.1)
            corrs.rnd.e[cr,"Element.2"]<-gsub(pattern = "trait", replacement = "", vn.n.2)
            corrs.rnd.e[cr,"rho"]<-rho.n
            corrs.rnd.e[cr,"HPD.lower"]<-hpd.l
            corrs.rnd.e[cr,"HPD.upper"]<-hpd.u
            corrs.rnd.e[cr,"pMCMC"]<-p.rho
          }
        }
        
        rndDF.list[[rnd.only.e]]<-corrs.rnd.e
        
        #--- build matrices ---
        rndM.e<-matrix(rep(1, ncol(vcv.rnd)), ncol=sqrt(ncol(vcv.rnd)))
        
        if(any(grepl("Intercept", nms.Me))){
          
          ntu.vector<-c("(Intercept)", rndDF.list[[rnd.only.e]][1:(sqrt(ncol(vcv.rnd))-1),"Element.1"])
          ntu.vector<-gsub(pattern = "trait", replacement = "", ntu.vector)
          
        }else{
          
          ntu.vector<-c(rndDF.list[[rnd.only.e]][1,"Element.2"], rndDF.list[[rnd.only.e]][1:(sqrt(ncol(vcv.rnd))-1),"Element.1"])
          ntu.vector<-gsub(pattern = "trait", replacement = "", ntu.vector)
          
        }
        
        colnames(rndM.e)<-ntu.vector
        rownames(rndM.e)<-ntu.vector
        
        rhoDF.e<-rndDF.list[[rnd.only.e]]
        for(i in ntu.vector){
          for(j in ntu.vector){
            if(i==j){next}else{
              rho.ij<-rhoDF.e$rho[rhoDF.e$Element.1==i & rhoDF.e$Element.2==j]
              rndM.e[i,j]<-rho.ij
              rndM.e[j,i]<-rho.ij
            }
          }
        }
        
        rndM.list[[rnd.only.e]]<-rndM.e
        
        if(p.mcmc){
          p.mcmc.Me<-matrix(rep(1, ncol(vcv.rnd)), ncol=sqrt(ncol(vcv.rnd)))
          colnames(p.mcmc.Me)<-ntu.vector
          rownames(p.mcmc.Me)<-ntu.vector
          for(r in ntu.vector){
            for(p in ntu.vector){
              if(r==p){next}else{
                P.rp<-rhoDF.e$pMCMC[rhoDF.e$Element.1==r & rhoDF.e$Element.2==p]
                p.mcmc.Me[r,p]<-P.rp
                p.mcmc.Me[p,r]<-P.rp
              }
            }
          }
          
          p.mcmc.list.rnd[[rnd.only.e]]<-p.mcmc.Me
          
        }
        
      }
  }
  
  
  #--- residual structure ---
  
  resid.elements<-split.direct.sum(as.character(model$Residual$formula[2]))
  resid.pos<-grep("units", colnames(model$VCV))
  
  if(length(resid.pos)<4){
    resid.corr<-FALSE
    if(verbose){warning("Less than 4 elements in the residual covariance matrix of - excluding from analyses!")}
  } else {
    
    resid.corr<-TRUE
    
    resid.nms<-colnames(model$VCV)[resid.pos]
    vcv.rsd<-model$VCV[,resid.nms]
    resid.nms<-gsub(pattern = ".units", replacement = "", x = resid.nms)
    
    rsdDF.list<-list()
    rsdM.list<-list()
    for(r in resid.elements){
      
      #--- get correlations ---
      corrs.rsd.e<-data.frame()
      for(d in resid.nms){
        
        rsn.d<-strsplit(d, "\\:")[[1]]
        rsn.d.1<-paste(rsn.d[1:(length(rsn.d)/2)], collapse=":")
        rsn.d.2<-paste(rsn.d[(1+length(rsn.d)/2):length(rsn.d)], collapse=":")
        
        if(rsn.d.1==rsn.d.2){
          next
        }else{
          rscov.name<-paste(d, "units", sep=".")
          rsvv1<-paste(paste(rep(rsn.d.1, 2), collapse=":"), "units", sep=".")
          rsvv2<-paste(paste(rep(rsn.d.2, 2), collapse=":"), "units", sep=".")
          
          cov.d<-vcv.rsd[,rscov.name]
          sqrt.rsvv.d<-sqrt(vcv.rsd[,rsvv1]*vcv.rsd[,rsvv2])
          corr.post.rs<-cov.d/sqrt.rsvv.d
          rs.rho.d<-posterior.mode(corr.post.rs)
          rs.hpd.l<-HPDinterval(corr.post.rs)[1]
          rs.hpd.u<-HPDinterval(corr.post.rs)[2]
          rs.p.rho<-pMCMC(corr.post.rs)[1,1]
          
          cr<-nrow(corrs.rsd.e)+1
          corrs.rsd.e[cr,"Element.1"]<-gsub(pattern = "trait", replacement = "", rsn.d.1)
          corrs.rsd.e[cr,"Element.2"]<-gsub(pattern = "trait", replacement = "", rsn.d.2)
          corrs.rsd.e[cr,"rho"]<-rs.rho.d
          corrs.rsd.e[cr,"HPD.lower"]<-rs.hpd.l
          corrs.rsd.e[cr,"HPD.upper"]<-rs.hpd.u
          corrs.rsd.e[cr,"pMCMC"]<-rs.p.rho
        }
      }
      
      rsdDF.list[["units"]]<-corrs.rsd.e
      
      #--- build matrices ---
      rsdM.r<-matrix(rep(1, ncol(vcv.rsd)), ncol=sqrt(ncol(vcv.rsd)))
      
      rs.ntu.vector<-c(rsdDF.list[["units"]][1,"Element.2"], rsdDF.list[["units"]][1:(sqrt(ncol(vcv.rsd))-1),"Element.1"])
      
      colnames(rsdM.r)<-rs.ntu.vector
      rownames(rsdM.r)<-rs.ntu.vector
      
      rhoDF.r<-rsdDF.list[["units"]]
      for(k in rs.ntu.vector){
        for(l in rs.ntu.vector){
          if(k==l){next}else{
            rho.kl<-rhoDF.r$rho[rhoDF.r$Element.1==k & rhoDF.r$Element.2==l]
            rsdM.r[k,l]<-rho.kl
            rsdM.r[l,k]<-rho.kl
          }
        }
      }
      
      rsdM.list[["units"]]<-rsdM.r
      
      if(p.mcmc){
        p.mcmc.Mr<-matrix(rep(1, ncol(vcv.rsd)), ncol=sqrt(ncol(vcv.rsd)))
        colnames(p.mcmc.Mr)<-rs.ntu.vector
        rownames(p.mcmc.Mr)<-rs.ntu.vector
        for(s in rs.ntu.vector){
          for(t in rs.ntu.vector){
            if(s==t){next}else{
              P.st<-rhoDF.r$pMCMC[rhoDF.r$Element.1==s & rhoDF.r$Element.2==t]
              p.mcmc.Mr[s,t]<-P.st
              p.mcmc.Mr[t,s]<-P.st
            }
          }
        }
        
        p.mcmc.list.rs[["units"]]<-p.mcmc.Mr
        
      }
      
    }
    
  }
  
  
  
  CORR<-list()
  
  if(type=="matrix"){
    if(random.corr){
      CORR[["Random.rho"]]<-rndM.list
      if(p.mcmc){CORR[["Random.pMCMC"]]<-p.mcmc.list.rnd}
      }
    
    if(resid.corr){
      CORR[["Residual.rho"]]<-rsdM.list
      if(p.mcmc){CORR[["Residual.pMCMC"]]<-p.mcmc.list.rs}
      }
  }
  
  if(type=="dataframe"){
    if(random.corr){CORR[["Random"]]<-rndDF.list}
    if(resid.corr){CORR[["Residual"]]<-rsdDF.list}
  }
  
  if(!random.corr & !resid.corr){stop("Either the random or the residual covariance matrix must have at least 4 elements!")}else{
    return(CORR)
  }
  
}


