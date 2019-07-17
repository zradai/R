### Activation functions, with derivative functions as well ###

AF_identity = function(x, theta, deriv=FALSE, simple=TRUE){
  theta = as.matrix(theta)
  if(!deriv){
    gz = x %*% theta
    return(gz)
  }else{
    if(is.null(nrow(x))){m = length(x)}else{m = nrow(x)}
    gpz = rep(1, m)
    return(gpz)
  }
}

AF_logistic = function(x, theta, deriv=FALSE, simple=TRUE){
  theta = as.matrix(theta)
  if(!deriv){
    z = x %*% theta
    gz = 1/(1+exp(-z))
    return(gz)
  }else{
    if(simple){
      gpz = x * (1 - x)
    }else{
      z = x %*% theta
      gz = 1/(1+exp(-z))
      gpz = gz * (1 - gz)
    }
    
    return(gpz)
  }
}

AF_tanh = function(x, theta, deriv=FALSE, simple=TRUE){
  theta = as.matrix(theta)
  if(!deriv){
    z = x %*% theta
    gz = tanh(z)
    return(gz)
  }else{
    if(simple){
      gpz = 1 - x^2
    }else{
      z = x %*% theta
      gz = tanh(z)
      gpz = 1 - gz^2
    }
    return(gpz)
  }
}

AF_relu = function(x, theta, deriv=FALSE, simple=TRUE){
  theta = as.matrix(theta)
  if(!deriv){
    z = x %*% theta
    z[z<0] = 0
    gz = z
    return(gz)
  }else{
    if(simple){
      z = x
      z[z>=0] = 1
      z[z<0] = 0
      gpz=z
    }else{
      z = x %*% theta
      z[z>=0] = 1
      z[z<0] = 0
      gpz=z
    }
    return(gpz)
  }
}

AF_softplus = function(x, theta, deriv=FALSE, simple=TRUE){
  theta = as.matrix(theta)
  if(!deriv){
    z = x %*% theta
    gz = log(1+exp(z))
    return(gz)
  }else{
    if(simple){
      gpz = 1/(1+exp(-x))
    }else{
      z = x %*% theta
      gpz = 1/(1+exp(-z))
    }
    return(gpz)
  }
}

AF_sinus = function(x, theta, deriv=FALSE, simple=TRUE){
  theta = as.matrix(theta)
  if(!deriv){
    z = x %*% theta
    gz = sin(z)
    return(gz)
  }else{
    if(simple){
      gpz = cos(x)
    }else{
      z = x %*% theta
      gpz = cos(z)
    }
    return(gpz)
  }
}

AF_gauss = function(x, theta, deriv=FALSE, simple=TRUE){
  theta = as.matrix(theta)
  if(!deriv){
    z = x %*% theta
    gz = exp(-(z^2))
    return(gz)
  }else{
    if(simple){
      gpz = -2 * x * exp(-(x^2))
    }else{
      z = x %*% theta
      gpz = -2 * z * exp(-(z^2))
    }
    return(gpz)
  }
}
