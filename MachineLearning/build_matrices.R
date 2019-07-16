# building Theta matrices
build.matrices <- function(L = 3, s = 2, s.input=2, s.output=1, random.initialization=TRUE){
  # 's.output' is the number of units in the last (output) layer
  # 'L' is the total number of layers (input + hidden + output)
  # 's' is either an integer, or a vector of integers, giving the number of
  # hidden units in hidden layers (WITHOUT BIAS UNITS!!) -- 
  # the order of numbers in the vector maps to the hidden layers 
  # (e.g. if L = 5, then number of hidden layers is 3; s = c(5, 4, 3), 
  # i.e. in the 1st hidden layer s=5, in the 2nd s=4, in the 3rd s=3 )
  
  if(length(s)==1 & L>3){
    s<-rep(s, (L-2))
  }
  s<-c(s.input, s, s.output)
  
  # building L-1 Theta matrices
  TT <- list()
  for(t in 1:(L-1)){
    dims.t<-c(s[(t+1)], (s[t]+1)) # row × col
    if(random.initialization){
      
      TT[[t]]<-matrix(runif(n = prod(dims.t), min = -1, max = 1), nrow=dims.t[1], ncol=dims.t[2] )
      if(length(unique(TT[[t]]))!=prod(dims.t)){
        while (length(unique(TT[[t]]))!=prod(dims.t)) {TT[[t]]<-matrix(runif(n = prod(dims.t), min = -1, max = 1), nrow=dims.t[1], ncol=dims.t[2] )}
      }
      
    }else{
      TT[[t]]<-matrix(rep(0, prod(dims.t)), nrow=dims.t[1], ncol=dims.t[2] )
    }
  }
  return(TT)
}
