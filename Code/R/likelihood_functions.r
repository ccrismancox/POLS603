###### Normal #####
l.normal <- function(theta, y, X, sum=TRUE){
  ## Normal negative log-likelihood
  ## inputs: theta: a guess at regression parameters with log(s2) last
  ##         y: dependent variable
  ##         X: independent variables
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood
  s2 <- exp(theta[length(theta)])
  beta <- theta[-length(theta)]
  Xb <- drop(X %*% beta)

  LL <- -1/2 *log(2*base::pi) - 1/2 * log(s2) -1/(2*s2) *(y-Xb)^2
  if(sum){
    return(- sum(LL))
  }else{
    return(LL)
  }
}
r.normal <- function(theta, y, X, sum=TRUE){
  ## Poisson negative score
  ## inputs: theta: a guess at regression parameters with log(s2) last
  ##         y: dependent variable
  ##         X: independent variables
  ##         sum: a flag for whether to return the 
  ##              summed negative score or 
  ##              the Jacobian matrix
  ## returns: (negative) score or Jacobian
  s2 <- exp(theta[length(theta)])
  beta <- theta[-length(theta)]
  Xb <- drop(X %*% beta)
  J <- X*(y-Xb)/s2
  if(sum){
    return(-colSums(J))
  }else{
    return(J)
  }
}

#### Poisson #####
l.poisson <- function(beta, y,X, sum=TRUE){
  ## Poisson negative log-likelihood
  ## inputs: b: a guess at regression parameters
  ##         y: dependent variable
  ##         X: independent variables
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood
  XB <- X%*%beta
  lambda <- exp(XB)
  LL <- y*XB - lambda
  if(sum){
    return(-sum(LL))
  }else{
    ## this is here because the main reason for 
    ## having the vector form is for Vuong tests
    return(LL-lgamma(y+1)) 
  }
}

r.poisson <- function(beta, y,X, sum=TRUE){
  ## Poisson negative score
  ## inputs: b: a guess at regression parameters
  ##         y: dependent variable
  ##         X: independent variables
  ##         sum: a flag for whether to return the 
  ##              summed negative score or 
  ##              the Jacobian matrix
  ## returns: (negative) score or Jacobian
  XB <- drop(X%*%beta)
  lambda <- exp(XB)
  score <- X*(y-lambda)
  if(sum){
    return(-colSums(score))
  }else{
    return(score)
  }
}

h.poisson <- function(beta, y,X){
  ## Poisson negative Hessian
  ## inputs: b: a guess at regression parameters
  ##         y: dependent variable
  ##         X: independent variables
  ## returns: (negative) Hessian
  XB <- drop(X%*%beta)
  lambda <- exp(XB)
  bread <- diag(lambda)
  H <- t(X) %*% bread %*% X
  return(H)
}

#### Negative binomial ####
l.negbin <- function(b,y,X, sum=TRUE){
  a <- b[length(b)] #log(r)
  r <- exp(a)
  b <- b[-length(b)]
  
  Xb <- drop(X%*%beta)
  
  LL <- a*r + y*XB- (y+r)*log(r+exp(Xb)) - 
    lgamma(y+1) +lgamma(y+r) - lgamma(r)
  if(sum){
    return(-LL)
  }else{
    return(LL)
  }
}

r.negbin <- function(b,y,X, sum=TRUE){
  a <- b[length(b)] #log(r)
  r <- exp(a)
  b <- b[-length(b)]
  
  Xb <- drop(X%*%beta)
  l <- exp(Xb)
  
  Db <- X *(y*r - l*r)/(r+l)
  Da <- r *( a - y/(r+l) - log(r + l) + digamma(y+r)-digamma(r) + 1 - r/(r+l))
  J <- cbind(Db, Da)
  if(sum){
    return(-colSums(J))
  }else{
    return(J)
  }
}



