##############################
#### EVALUATE SIMULATIONS ####
##############################

eval.sim<-function(table, 
                   plotting=TRUE, 
                   rand.slope=FALSE, 
                   resampling=FALSE, 
                   col.per.ind=TRUE, 
                   ind.lines=FALSE, 
                   var.cov.method="statter", # var-cov matrix estimation function (statter or kenw)
                   averaging=FALSE){ 
  # col.per.ind = colors on plots set unique to each individual
  # ind.lines = put regression lines individually on plot (i.e. for all individuals)
  
  require("lmerTest")
  
  cat("Please wait...\n")
  
  if(resampling){
    cat("Please specify how many samples should the simulation use:\n")
    rs.n<-readline()
    table<-resamp(table,rs.n)
    cat(" sampling...\n")
  }
  
  ests<-data.frame(Simulation=0, NI=0, Reps=0, B1.estimate=0, p.value=0, t.value=0, SE=0, c.R.squared=0, rand.Int.var=0, rand.Slope.var=0, resid.var=0, V.method=0)
  if(averaging){
    #est.avg<-data.frame(Simulation=0, NI=0, Reps=0, B1.estimate=0, p.value=0, t.value=0, SE=0, c.R.squared=0, rand.Int.var=0, rand.Slope.var=0, resid.var=0)
    cat(" averaging...\n")
    est.avg<-ests
    t.avg<-table[-c(1:nrow(table)),c("Individual", "Y", "B1", "I", "SIM")]
    tr<-1
    
    for(ss in unique(table$SIM)){
      for(ti in unique(table$Individual[table$SIM==ss])){
        t.avg[tr,"Individual"]<-ti
        t.avg[tr,"Y"]<-mean(table$Y[table$Individual==ti & table$SIM==ss])
        t.avg[tr,"X1"]<-mean(table$X1[table$Individual==ti & table$SIM==ss])
        t.avg[tr,"B1"]<-unique(table$B1[table$Individual==ti & table$SIM==ss])
        t.avg[tr,"S1"]<-unique(table$S1[table$Individual==ti & table$SIM==ss])
        t.avg[tr,"I"]<-unique(table$I[table$Individual==ti & table$SIM==ss])
        t.avg[tr,"SIM"]<-ss
        tr<-tr+1 
        
      }
    }
    
  }
  
  if(plotting){
    if(length(unique(table$SIM))>6){warning("More than 6 simulations: plotting is automatically turned off. Would You like to plot the results anyway? (y/n)")
      plotting.rdln<-readline()
    }
    if(length(unique(table$SIM))>1){
      par(mfrow=c(ceiling(max(table$SIM)/2), 2))
    }
    plotting.rdln<-"yes"
  } else {
    plotting.rdln<-"no"
  }
  
  for(s in 1:max(table$SIM)){
    
    cat("Analysing",paste(s,".", sep=""),"simulation out of",max(table$SIM),"\n")
    
    if(length(unique(table$Repeat[table$SIM==s]))==1){
      m.s<-lm(Y~X1, data=table[table$SIM==s,])
      ests[s,"Simulation"]<-s
      b.est<-round(coef(summary(m.s))[2], 3)
      ests[s,"B1.estimate"]<-b.est
      if(round(coef(summary(m.s))[8], 3)==0){
        p.val<-"P < 0.001"
      }else{p.val<-paste("P: ",round(coef(summary(m.s))[10], 3),sep="")}
      ests[s,"p.value"]<-p.val
      ests[s,"t.value"]<-round(coef(summary(m.s))[6], 2)
      ests[s,"SE"]<-round(coef(summary(m.s))[4], 4)
      ests[s,"resid.var"]<-as.numeric(var(m.s$resid))
      ests[s,"c.R.squared"]<-as.numeric(r.squaredGLMM(m.s)[2])
      ests[s,"NI"]<-max(table$Individual[table$SIM==s])
      ests[s,"Reps"]<-max(table$Repeat[table$SIM==s])
      ests[s,"original.B1"]<-unique(table$B1[table$SIM==s])
    } else {
      if(!rand.slope){
        m.s<-lmer(Y~X1+(1|Individual), data=table[table$SIM==s,])
        ests[s,"Simulation"]<-s
        b.est<-round(coef(summary(m.s))[2], 3)
        ests[s,"B1.estimate"]<-b.est
        if(round(coef(summary(m.s))[10], 3)==0){
          p.val<-"P < 0.001"
        }else{p.val<-paste("P: ",round(coef(summary(m.s))[10], 3),sep="")}
        ests[s,"p.value"]<-p.val
        ests[s,"t.value"]<-round(coef(summary(m.s))[8], 2)
        ests[s,"SE"]<-round(coef(summary(m.s))[4], 4)
        ests[s,"rand.Int.var"]<-as.numeric(VarCorr(m.s)[[1]][1])
        ests[s,"resid.var"]<-as.numeric(attr(VarCorr(m.s), "sc"))^2
        ests[s,"c.R.squared"]<-as.numeric(r.squaredGLMM(m.s)[2])
        ests[s,"NI"]<-max(table$Individual[table$SIM==s])
        ests[s,"Reps"]<-max(table$Repeat[table$SIM==s])
        ests[s,"original.B1"]<-unique(table$B1[table$SIM==s])
        
      } else {
        
        m.s<-lmer(Y~X1+(X1|Individual), data=table[table$SIM==s,])
        
        if(var.cov.method=="statter"){s.ms<-summary(m.s)}     
        if(var.cov.method=="kenw"){s.ms<-summary(m.s, ddf="kenw")}
        
        ests[s,"Simulation"]<-s
        b.est<-round(coef(s.ms)[2], 3)
        ests[s,"B1.estimate"]<-b.est
        if(round(coef(s.ms)[10], 3)==0){
          p.val<-"P < 0.001"
        }else{p.val<-paste("P: ",round(coef(s.ms)[10], 3),sep="")}
        ests[s,"p.value"]<-p.val
        ests[s,"t.value"]<-round(coef(s.ms)[8], 2)
        ests[s,"SE"]<-round(coef(s.ms)[4], 4)
        ests[s,"rand.Int.var"]<-as.numeric(VarCorr(m.s)[[1]][1])
        ests[s,"rand.Slope.var"]<-as.numeric(VarCorr(m.s)[[1]][4])
        ests[s,"resid.var"]<-as.numeric(attr(VarCorr(m.s), "sc"))^2
        ests[s,"c.R.squared"]<-as.numeric(r.squaredGLMM(m.s)[2])
        ests[s,"NI"]<-max(table$Individual[table$SIM==s])
        ests[s,"Reps"]<-max(table$Repeat[table$SIM==s])
        ests[s,"original.B1"]<-unique(table$B1[table$SIM==s])
      }
    }
    
    
    
    if(plotting & plotting.rdln!="n"){
      
      xlm<-c(min(table$X1), max(max(table$X1)))
      ylm<-c(min(table$Y), max(max(table$Y)))
      
      with(table[table$SIM==s,], plot(X1, Y, main=paste("Simulation: ",s,sep=""), xlim=xlm, ylim=ylm))
      
      if(col.per.ind){
        clz<-rainbow(length(unique(table$Individual[table$SIM==s])))
        ind.colz<-c()
        for(c in 1:length(clz)){
          for(i in 1:max(table$Repeat)){       # max. number of cycles = number of measurements per ind.s
            ind.colz<-c(ind.colz,clz[c])
          }
        }
        with(table[table$SIM==s,], points(X1, Y, col=ind.colz))
      }
      
      txt.y<-max(table$Y)-abs(max(table$Y)*0.1)
      txt.x<-min(table$X1)+abs(min(table$X1)*0.1)
      text(txt.x, txt.y, labels=paste("B_est: ",b.est,", ", p.val, sep=""), pos=4)
      abline(a=fixef(m.s)[1], b=fixef(m.s)[2], lwd=3, lty=2)
      
      if(ind.lines){
        coef.int<-coef(m.s)$Individual[,1]
        coef.slp<-coef(m.s)$Individual[,2]
        if(col.per.ind){
          for(i in 1:length(clz)){
            abline(a=coef.int[i], b=coef.slp[i], col=clz[i], lwd=2, lty=4)
          }
        } else {
          for(i in 1:length(unique(table$Individual))){
            abline(a=coef.int[i], b=coef.slp[i], lwd=2, lty=4)
          }
        }
      } 
      
    }
    
    
    if(averaging){
      
      m.s.a<-lm(Y~X1, data=t.avg[t.avg$SIM==s,])
      
      est.avg[s,"Simulation"]<-s
      b.est<-round(coef(summary(m.s.a))[2], 3)
      est.avg[s,"B1.estimate"]<-b.est
      if(round(coef(summary(m.s.a))[8], 3)==0){
        p.val<-"P < 0.001"
      }else{p.val<-paste("P: ",round(coef(summary(m.s.a))[10], 3),sep="")}
      est.avg[s,"p.value"]<-p.val
      est.avg[s,"t.value"]<-round(coef(summary(m.s.a))[6], 2)
      est.avg[s,"SE"]<-round(coef(summary(m.s.a))[4], 4)
      est.avg[s,"rand.Int.var"]<-0
      est.avg[s,"resid.var"]<-(summary(m.s.a)$sigma)**2
      est.avg[s,"c.R.squared"]<-as.numeric(r.squaredGLMM(m.s.a)[2])
      est.avg[s,"NI"]<-max(table$Individual[table$SIM==s])
      est.avg[s,"Reps"]<-max(table$Repeat[table$SIM==s])
      est.avg[s,"original.B1"]<-unique(table$B1[table$SIM==s])
      
      est.avg$sim.type<-"avg"
      
    }
    
    
    
  }
  
  ests$sim.type<-"lmm"
  
  if(averaging){
    ests<-rbind(ests, est.avg)
    ests<-ests[order(ests$sim.type),]
  }
  
  return(ests)
  
  par(mfrow=c(1,1))
  
}
