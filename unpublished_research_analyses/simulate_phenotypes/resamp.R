
########################################################
##### Resample from the output-table of 'sim.phen' #####
########################################################

resamp<-function(table, n, samp.reps=FALSE){
  
  if(samp.reps){ # stratified sampling, using 'n' repeatitions from the specified datatable
    res<-sample(unique(table$Repeat), n, replace=FALSE)
    table.res<-table
    for(r in 1:n){
      res.r<-res[r]
      table.res<-table.res[table.res$Repeat!=res.r,]
    }
    
  } else {
    
    res<-sample(unique(table$Individual), n, replace=FALSE)
    res.rows<-c()
    for(i in 1:nrow(table)){
      for(r in res){
        if(table[i,"Individual"]==r){res.rows<-c(res.rows, i)}
      }
    }
    
    table.res<-table[res.rows,]
    
  }
  
  
  return(table.res)
}
