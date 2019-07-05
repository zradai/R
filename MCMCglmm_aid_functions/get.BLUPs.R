
# get posterior distributions for random intercepts and/or slopes
get.BLUPs<-function(model, post.mode=TRUE){
  
  ## NOTE: in multiresponse models, all random factors should be given either as 'idh(trait):random' or 'us(trait):random'!
  ## because acquiring posterior modes of random levels is done by calling the random variable's name as 'traitRESP.random'
  
  if(class(model)!="MCMCglmm"){stop("This function requires an object of 'MCMCglmm' class!")}
  if(is.null(model$Random$formula)){stop("This function needs a model with random effects!")}
  if(ncol(model$Sol)==model$Fixed$nfl){stop("Posterior distribution of random effects must be saved during model fitting - make sure to set 'pr=TRUE'!")}
  
  rands<-all.vars(as.formula(model$Random$formula))
  fixs<-all.vars(as.formula(model$Fixed$formula))
  resps<-all.vars(as.formula(model$Fixed$formula)[[2]])
  preds<-all.vars(as.formula(model$Fixed$formula)[[3]])
  
  rands.only<-setdiff(rands, fixs)
  rands.only<-rands.only[rands.only!="trait"]
  Rslope.vars<-intersect(rands, fixs)
  Rslope.vars<-Rslope.vars[Rslope.vars!="trait"]
  
  INTS<-data.frame()
  SLPS<-data.frame()
  
  if(length(resps)>1){
    for(y in resps){
      for(r in rands.only){
        
        Rm<-model$Sol[,(model$Fixed$nfl+1):ncol(model$Sol)]
        Rm<-Rm[,grep(r,colnames(Rm))]
        
        if(length(Rslope.vars)>0){
          if(length(grep("(Intercept)",colnames(Rm)))>0){
            R.INTS<-posterior.mode(Rm[,grep("(Intercept)",colnames(Rm))])
            INTS[(nrow(INTS)+1):(length(R.INTS)+nrow(INTS)), "Random.intercept"]<-R.INTS
            INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.level"]<-substr(names(R.INTS), 13, nchar(names(R.INTS)))
            INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.effect"]<-as.character(r)
            INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Response"]<-as.character(y)
          }else{
            R.INTS<-"empty"
          }
          for(s in Rslope.vars){
            RS.name<-paste("trait",y,":",s, sep="")
            if(length(grep(RS.name, colnames(Rm)))>0){
              R.SLPS<-posterior.mode(Rm[,grep(RS.name, colnames(Rm))])
              SLPS[(nrow(SLPS)+1):(length(R.SLPS)+nrow(SLPS)), "Random.slope"]<-R.SLPS
              SLPS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.level"]<-substr(names(R.SLPS), (nchar(RS.name)+nchar(s)), nchar(names(R.SLPS)))
              SLPS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.effect"]<-as.character(r)
              SLPS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Fix.effect"]<-as.character(s)
              SLPS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Response"]<-as.character(y)
              if(R.INTS[1]=="empty"){
                INTS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.intercept"]<-NA
                INTS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.level"]<-substr(names(R.SLPS), (nchar(RS.name)+nchar(s)), nchar(names(R.SLPS)))
                INTS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.effect"]<-as.character(r)
                INTS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Response"]<-as.character(y)
              }
            } else {
              rint.name<-paste("trait",y,".",r, sep="")
              if(length(grep(rint.name, colnames(Rm)))>0){
                R.INTS<-posterior.mode(Rm[,grep(rint.name,colnames(Rm))])
                INTS[(nrow(INTS)+1):(length(R.INTS)+nrow(INTS)), "Random.intercept"]<-R.INTS
                INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.level"]<-substr(names(R.INTS), (nchar(rint.name)-nchar(r))+1, nchar(names(R.INTS)))
                INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.effect"]<-as.character(r)
                INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Response"]<-as.character(y)
              }
            }
          }
        }
        
        if(length(grep("(Intercept)",colnames(Rm)))==0 & length(Rslope.vars)==0){ 
          rint.name<-paste("trait",y,".",r, sep="")
          R.INTS<-posterior.mode(Rm[,grep(rint.name,colnames(Rm))])
          INTS[(nrow(INTS)+1):(length(R.INTS)+nrow(INTS)), "Random.intercept"]<-R.INTS
          INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.level"]<-substr(names(R.INTS), (nchar(rint.name)-nchar(r))+1, nchar(names(R.INTS)))
          INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.effect"]<-as.character(r)
          INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Response"]<-as.character(y)
        }
        
      }
    }
    
  }else{
    
    for(r in rands.only){
      
      Rm<-model$Sol[,(model$Fixed$nfl+1):ncol(model$Sol)]
      Rm<-Rm[,grep(r,colnames(Rm))]
      
      if(length(Rslope.vars)>0){
        if(length(grep("(Intercept)",colnames(Rm)))>0){
          R.INTS<-posterior.mode(Rm[,grep("(Intercept)",colnames(Rm))])
          INTS[(nrow(INTS)+1):(length(R.INTS)+nrow(INTS)), "Random.intercept"]<-R.INTS
          INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.level"]<-substr(names(R.INTS), 13, nchar(names(R.INTS)))
          INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.effect"]<-as.character(r)
        }else{
          R.INTS<-"empty"
        }
        for(s in Rslope.vars){
          RS.name<-paste(s,r,sep=".")
          if(length(grep(RS.name, colnames(Rm)))>0){
            R.SLPS<-posterior.mode(Rm[,grep(RS.name, colnames(Rm))])
            SLPS[(nrow(SLPS)+1):(length(R.SLPS)+nrow(SLPS)), "Random.slope"]<-R.SLPS
            SLPS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.level"]<-substr(names(R.SLPS), (nchar(s)+2), nchar(names(R.SLPS)))
            SLPS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.effect"]<-as.character(r)
            SLPS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Fix.effect"]<-as.character(s)
            if(R.INTS[1]=="empty"){
              INTS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.intercept"]<-NA
              INTS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.level"]<-substr(names(R.SLPS), (nchar(s)+2), nchar(names(R.SLPS)))
              INTS[(nrow(SLPS)-length(R.SLPS)+1):nrow(SLPS), "Random.effect"]<-as.character(r)
              }
          } else {
            if(length(grep(r, colnames(Rm)))>0){
              R.INTS<-posterior.mode(Rm[,grep(r,colnames(Rm))])
              INTS[(nrow(INTS)+1):(length(R.INTS)+nrow(INTS)), "Random.intercept"]<-R.INTS
              INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.level"]<-substr(names(R.INTS), (nchar(r)+2), nchar(names(R.INTS)))
              INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.effect"]<-as.character(r)
            }
          }
        }
      }
      
      if(length(grep("(Intercept)",colnames(Rm)))==0 & length(Rslope.vars)==0){ # get R-intercepts from models with only R-intercepts
        R.INTS<-posterior.mode(Rm[,grep(r,colnames(Rm))])
        INTS[(nrow(INTS)+1):(length(R.INTS)+nrow(INTS)), "Random.intercept"]<-R.INTS
        INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.level"]<-substr(names(R.INTS), 1, nchar(names(R.INTS)))
        INTS[(nrow(INTS)-length(R.INTS)+1):nrow(INTS), "Random.effect"]<-as.character(r)
      }
    }
  }
  
  blups<-list()
  if(nrow(INTS)>0 & nrow(na.omit(INTS))>0){blups[["Random.intercepts"]]<-INTS}
  if(nrow(SLPS)>0 & nrow(na.omit(SLPS))>0){blups[["Random.slopes"]]<-SLPS}
  
  return(blups)
  
}
