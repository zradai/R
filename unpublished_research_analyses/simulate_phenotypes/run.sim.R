###########################
##### Run simulations #####
###########################

run.sim<-function(input.s, dyn_X1=FALSE, rescale=FALSE, fix_Y=FALSE){
  
  tt<-data.frame(Individual=0, Repeat=0, Y=0, X1=0, B1=0, e=0, I=0, S1=0, V.R=0, V.I=0, V.S1=0, Method.error=0, SIM=0)
  tt<-tt[-1,]
  
  c<-1  
  c.max<-(
    length(input.s$NI)*
      length(input.s$NR)*
      length(input.s$XA)*
      length(input.s$B1)*
      length(input.s$X1_mu)*
      length(input.s$X1_V)*
      length(input.s$VI)*
      length(input.s$VS1)*
      length(input.s$VR)*
      length(input.s$V.method)
  )
  for(ni in 1:length(input.s$NI)){
    for(nr in 1:length(input.s$NR)){
      for(xa in 1:length(input.s$XA)){
        for(b1 in 1:length(input.s$B1)){
          for(x1m in 1:length(input.s$X1_mu)){
            for(x1v in 1:length(input.s$X1_V)){
              for(vi in 1:length(input.s$VI)){
                for(vs1 in 1:length(input.s$VS1)){
                  for(vr in 1:length(input.s$VR)){
                    for(vm in 1:length(input.s$V.method)){
                      
                      cat("Progress:",c,"/",c.max, paste("(",round((c/c.max)*100,2),"%)", sep=""), "\n" )
                      
                      input<-list(  
                        NI=input.s$NI[ni],
                        NR=input.s$NR[nr],
                        XA=input.s$XA[xa], 
                        B1=input.s$B1[b1],
                        X1_mu=input.s$X1_mu[x1m],
                        X1_V=input.s$X1_V[x1v],
                        VI=input.s$VI[vi],
                        VS1=input.s$VS1[vs1],
                        VR=input.s$VR[vr],
                        V.method=input.s$V.method[vm]
                      )
                      
                      #if(dyn_X1){tt.s<-sim.phen(input, dyn_X1=TRUE)} else {tt.s<-sim.phen(input, dyn_X1=FALSE)}
                      tt.s<-sim.phen(input, dyn_X1=dyn_X1, rescale=rescale, fix_Y=fix_Y)
                      
                      tt.s$SIM<-c
                      c<-c+1
                      tt<-rbind(tt, tt.s)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  cat("Simulations ready!\n")
  return(tt)
  
}
