
###################################################
###----- Cut off designated variable names -----###
###------ from names from VCV components -------###
###################################################


cutoff.names<-function(name.list, cutoff.name){
  all.cutoff<-list()
  for(n in cutoff.name){
    cutoff.nms<-strsplit(name.list, paste(".", n, sep=""))
    co.nms.2<-c()
    for(c in 1:length(cutoff.nms)){
      co.nms.2<-c(co.nms.2, cutoff.nms[[c]])
    }
    all.cutoff[[n]]<-co.nms.2
  }
  return(all.cutoff)
}

#############################################################


################################################
###----- Cut off 'trait' from VCV names -----###
################################################

cutoff.trait<-function(vcv.names){
  new.names<-c()
  for(n in 1:length(vcv.names)){
    nn.n<-gsub("trait", "", vcv.names[n])
    new.names<-c(new.names, nn.n)
  }
  return(new.names)
}



#############################################################


################################################
###----- Cut off 'units' from VCV names -----###
################################################

cutoff.units<-function(vcv.names){
  new.names<-c()
  for(n in 1:length(vcv.names)){
    nn.n<-gsub(".units", "", vcv.names[n])
    new.names<-c(new.names, nn.n)
  }
  return(new.names)
}



