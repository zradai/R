### Cost function ###

cost_function = function(x, y, tt, ll=0, actFUN=AF_identity, linkFUN=LF_Gauss, logistic_cost=FALSE){
  
  L = length(tt) + 1
  m = nrow(x)
  aL_s = list()
  aL_s[[1]] = x
  for(l in 2:L){
    n.s = nrow(tt[[(l-1)]])
    M.l = matrix(rep(0, n.s*m), ncol = n.s)
    for(s in 1:n.s){
      
      if(l < L){
        M.l[,s] = actFUN(x = cbind(1, aL_s[[(l-1)]]), theta = tt[[(l-1)]][s,])
      }
      if(l==L){
        M.l[,s] = linkFUN(x = cbind(1, aL_s[[(l-1)]]), theta = tt[[(l-1)]][s,])
      }
    }
    aL_s[[l]] = M.l
  }
  
  tt.noBias = list()
  for(t in 1:length(tt)){
    tt.t = tt[[t]]
    tt.t[,1] = 0
    tt.noBias[[t]] = tt.t
  }
  
  reg.tt.noBias = lapply(X = tt.noBias, FUN = function(x){return(x^2)})
  
  if(logistic_cost){
    J = -(1/m)*sum( y*log(aL_s[[L]]) + (1-y)*log(1 - aL_s[[L]]) ) + # (-1/m)*sum( y*log(c( aL_s[[L]] )) + (1-y)*log(1-c( aL_s[[L]] )) )
      (ll/(2*m))*sum(sapply(X = reg.tt.noBias, FUN = sum))
  }else{
    #J = (1/m)*sum( (aL_s[[L]] - y)^2 ) + 
    J = (1/m)*sum( ( linkFUN(x = aL_s[[L]], mean = FALSE) - linkFUN(x = y, mean = FALSE))^2 ) + 
      (ll/(2*m))*sum(sapply(X = reg.tt.noBias, FUN = sum))
  }
  
  return(J)
  
}
