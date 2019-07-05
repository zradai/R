# get phylogenetic signal, or heritability, from an MCMCglmm with either a phylogen.tree or pedigree
gen.sign<-function(model, organism){
  if(class(model)!="MCMCglmm"){stop("This function requires an object of 'MCMCglmm' class!")}
  
  resps<-all.vars(as.formula(model$Fixed$formula)[[2]])
  rands<-all.vars(as.formula(model$Random$formula))
  rands<-rands[rands!="trait"]
  
  post.sign<-data.frame()
  rw<-0
  for(y in resps){
    for(o in 1:length(organism)){
      rw<-rw+1
      org.y<-organism[o]
      ptrn<-paste("trait",y,":trait",y,".",org.y, sep="")
      org.pos<-grep(ptrn, colnames(model$VCV))
      all.rands.pos<-c()
      if(length(org.pos)==0){
        ptrn<-org.y
        org.pos<-grep(ptrn, colnames(model$VCV))
        for(r in c(rands, "units")){
          pt.yor<-as.character(r)
          pos.yor<-grep(pt.yor, colnames(model$VCV))
          all.rands.pos<-c(all.rands.pos, pos.yor)
        }
      }else{
        for(r in c(rands, "units")){
          pt.yor<-paste("trait",y,":trait",y,".",r, sep="")
          pos.yor<-grep(pt.yor, colnames(model$VCV))
          all.rands.pos<-c(all.rands.pos, pos.yor)
        }
      }
      
      Vorg.y<-model$VCV[,org.pos]
      Vtotal.y<-apply(model$VCV[,all.rands.pos], 1, sum)
      
      lambda.PD<-Vorg.y/Vtotal.y
      post.sign[rw,"Response"]<-as.character(y)
      post.sign[rw,"Organism"]<-as.character(org.y)
      post.sign[rw,"lambda.postMode"]<-posterior.mode(lambda.PD)
      post.sign[rw,"lambda.HPD.l"]<-HPDinterval(lambda.PD)[1]
      post.sign[rw,"lambda.HPD.u"]<-HPDinterval(lambda.PD)[2]
    }
  }
  
  return(post.sign)
  
}
