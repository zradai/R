
###############################
##### Simulate phenotypes #####
###############################

sim.phen<-function(input, dyn_X1=FALSE, rescale=FALSE, fix_Y=FALSE){ # simulate phenotypes
  
  if(is.null(input$XA)){input$XA<-0}
  if(is.null(input$VI)){input$VI<-0}
  if(is.null(input$VS1)){input$VS1<-0}
 
  NI<-input$NI
  NR<-input$NR
  XA<-input$XA
  B1<-input$B1
  X1_mu<-input$X1_mu
  X1_V<-input$X1_V
  VI<-input$VI
  VS1<-input$VS1
  VR<-input$VR
  
  tab<-data.frame(Individual=0, Repeat=0, Y=0, X1=0, B1=0, e=0, I=0, S1=0, V.R=0, V.I=0, V.S1=0, V.method=0, Method.error=0)
  
  if(dyn_X1==TRUE){ # if dynamic X1 is set true, X1 values are sampled again during the repeated measures
                      # if false, X1 values remain the same within individual over all the repeated measures
    cr<-1
    for(i in 1:NI){
      
      I.i<-rnorm(1,0,sqrt(VI)) # get individual-specific diversion from population mean (random intercept of 'i' individual)
      S1.i<-rnorm(1,0,sqrt(VS1)) # get individual-specific diversion from population mean response (random slope of 'i' individual)
      
      # get phenotype values for given individual, 'NR' times
      for(j in 1:NR){
        
        X1.ij<-rnorm(1,X1_mu, sqrt(X1_V))
        e.ij<-rnorm(1,0,sqrt(VR))
        y.ij<-(XA+I.i)+(B1+S1.i)*X1.ij+e.ij
        
        tab[cr,"Individual"]<-i
        tab[cr,"Repeat"]<-j
        tab[cr,"Y"]<-y.ij
        tab[cr,"X1"]<-X1.ij
        tab[cr,"B1"]<-B1
        tab[cr,"e"]<-e.ij
        tab[cr,"I"]<-I.i
        tab[cr,"S1"]<-S1.i
        tab[cr,"V.R"]<-VR
        tab[cr,"V.I"]<-VI
        tab[cr,"V.S1"]<-VS1
        cr<-cr+1
      }
      
    }
  } else {
    X1.pop<-rnorm(NI, X1_mu, sqrt(X1_V))
    cr<-1
    for(i in 1:NI){
      
      I.i<-rnorm(1,0,sqrt(VI))
      S1.i<-rnorm(1,0,sqrt(VS1))
      
      for(j in 1:NR){
        
        X1.ij<-X1.pop[i]
        e.ij<-rnorm(1,0,sqrt(VR))
        y.ij<-(XA+I.i)+(B1+S1.i)*X1.ij+e.ij
        
        tab[cr,"Individual"]<-i
        tab[cr,"Repeat"]<-j
        tab[cr,"Y"]<-y.ij
        tab[cr,"X1"]<-X1.ij
        tab[cr,"B1"]<-B1
        tab[cr,"e"]<-e.ij
        tab[cr,"I"]<-I.i
        tab[cr,"S1"]<-S1.i
        tab[cr,"V.R"]<-VR
        tab[cr,"V.I"]<-VI
        tab[cr,"V.S1"]<-VS1
        cr<-cr+1
      }
      
    }
  }
  
  if(fix_Y){
    for(i in unique(tab$Individual)){
      tab$Y[tab$Individual==i]<-mean(tab$Y[tab$Individual==i])
    }
  }
  
  if(rescale){
    tab$Y<-(tab$Y-mean(tab$Y))/(sd(tab$Y))
    tab$X1<-(tab$X1-mean(tab$X1))/(sd(tab$X1))
  }
  
  tab$SIM<-1
  
  return(tab)
  
}
