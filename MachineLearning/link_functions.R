### Link functions ###

LF_Gauss = function(x, theta, mean=TRUE){
  if(mean){
    fx = x %*% theta
  }else{
    fx = x
  }
  return(fx)
}

LF_gamma = function(x, theta, mean=TRUE){
  if(mean){
    fx = -(x %*% theta)^(-1)
  }else{
    fx = -x^(-1)
  }
  return(fx)
}

LF_invGauss = function(x, theta, mean=TRUE){
  if(mean){
    fx = (x %*% theta)^(-1/2)
  }else{
    fx = x^-(2)
  }
  return(fx)
}

LF_pois = function(x, theta, mean=TRUE, logplus=0){
  if(mean){
    fx = exp(x %*% theta)
  }else{
    fx = log(x+logplus)
  }
  return(fx)
}

LF_binom = function(x, theta, mean=TRUE){
  if(mean){
    fx = 1/(1+exp(-(x %*% theta)))
  }else{
    fx = log(x/(1-x))
  }
  return(fx)
}

