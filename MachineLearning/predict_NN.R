### Prediction ###

predict_NN = function(Thetas, x, actFUN=AF_identity, linkFUN=LF_Gauss){
  
  L = length(Thetas)+1
  #s = unique( sapply(X = Thetas, FUN = ncol)[2:(length(Thetas)-1)] )
  s.input = ncol(Thetas[[1]])-1
  s.output = nrow(Thetas[[length(Thetas)]])
  
  aL_s = list()
  aL_s[[1]] = x
  for(l in 2:L){
    n.s = nrow(Thetas[[(l-1)]])
    M.l = matrix(rep(0, n.s*nrow(x)), ncol = n.s)
    for(s in 1:n.s){
      
      if(l < L){
        M.l[,s] = actFUN(x = cbind(1,aL_s[[(l-1)]]), theta = Thetas[[(l-1)]][s,], deriv = FALSE)
      }
      if(l == L){
        M.l[,s] = linkFUN(x = cbind(1,aL_s[[(l-1)]]), theta = Thetas[[(l-1)]][s,])
      }
    }
    aL_s[[l]] = M.l
  }
  
  return(aL_s[[L]])
  
}
