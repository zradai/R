
train_NN = function(x, y, tt, actFUN=AF_identity, linkFUN = LF_Gauss, algorithm="GD", 
                    aa=0.001, ll=0, n_iter=1000, 
                    logistic_cost=FALSE, verbose=TRUE, R2=FALSE, CrossVal=FALSE, runtime=FALSE){
  if(runtime){run_start = Sys.time()}
  costs = c()
  L = length(tt) + 1
  m = nrow(x)
  n = ncol(x)
  
  for(e in 1:n_iter){
    
    if(verbose & algorithm=="GD"){cat("Progress: ", e, " / ", n_iter, "\n", sep="")}
    
    # cost function
    J.e = cost_function(x = x, y = y, tt = tt, ll = ll, actFUN = actFUN, linkFUN = linkFUN, logistic_cost = logistic_cost)
    costs= c(costs, J.e)
    
    D.list = list()
    for(d in 1:(L-1)){
      D.list[[d]] = 0
    }
    
    # Training algorithm (for now only gradient descent)
    for(i in 1:nrow(x)){
      
      # Forward propagation
      aL_s = list()
      aL_s[[1]] = cbind(1,t(x[i,]))
      for(l in 2:L){
        n.s = nrow(tt[[(l-1)]])
        M.l = matrix(rep(0, n.s), ncol = n.s)
        for(s in 1:n.s){
          if(l < L){
            M.l[,s] = actFUN(x = aL_s[[(l-1)]], theta = tt[[(l-1)]][s,], deriv = FALSE)
          }
          if(l == L){
            M.l[,s] = linkFUN(x = aL_s[[(l-1)]], theta = tt[[(l-1)]][s,])
          }
        }
        if(l==L){aL_s[[l]] = M.l}else{aL_s[[l]] = cbind(1, M.l)}
      }
      
      # Back propagation; order of error terms in the list is d_L, d_L-1, ..., d_2
      # (i.e. first list element is the error term for the last layer)
      d.list = list()
      backwrd.a = c(L:1)
      d.list[[1]] = aL_s[[L]] - y[i]
      for(d in 2:(L-1)){
        ai = backwrd.a[d]
        z = aL_s[[ai]][-1]
        if(nrow(tt[[ai]])==1){
          d.list[[d]] = (c(as.matrix(tt[[ai]][,-1]) %*% d.list[[(d-1)]])) * actFUN(x = z, deriv = TRUE, simple=TRUE)
        }else{
          d.list[[d]] = (c(t(tt[[ai]][,-1]) %*% d.list[[(d-1)]])) * actFUN(x = z, deriv = TRUE, simple=TRUE)
        }
      }
      
      # Accumulate errors
      for(D in 1:length(D.list)){
        ai = backwrd.a[D]-1
        D.list[[D]] = D.list[[D]] + d.list[[D]] %*% aL_s[[ai]]
      }
    }
    
    # Gradients
    grads = list()
    for(t in 1:length(D.list)){
      ai = backwrd.a[t]-1
      D.t = D.list[[t]]
      
      grad.noReg = (1/m)*D.t
      if(nrow(grad.noReg)==1){
        regM = c(0, (lambda/m)*tt[[ ai ]][,-1])
      }else{
        regM = cbind(0, (lambda/m)*tt[[ ai ]][,-1])
      }
      grads[[ai]] = grad.noReg + regM
    }
    
    for(g in 1:length(grads)){
      tt[[g]] = tt[[g]] - aa*grads[[g]]
    }
    
  }
  
  if(runtime){
    run_stop = Sys.time()
    cat("Runtime", run_stop-run_start, "\n")
  }
  
  # return results
  NN_results = list(
    Thetas = tt,
    Costs = costs
  )
  return(NN_results)
}
