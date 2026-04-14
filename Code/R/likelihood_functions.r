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


#### Binomial: Logit and probit #####
l.binomial <- function(b,Z, link=c("logit", "probit"), sum=TRUE){
  ## binomial negative log-likelihood
  ## inputs: b: guess at regression parameters
  ##         Z: the data, (2*y-1)*X
  ##         link: the link function either 
  ##               logit (default) or probit
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood  
  G <- match.arg(G)
  G <- switch(G, #assign G based on link
              "logit"=plogis,
              "probit"=pnorm)
  lp <- G(Z%*%b, log.p=TRUE)
  if(sum){
    return(-sum(lp))
  }else{
    return(lp)
  }
}
r.binomial <- function(b,Z, G=c("logit", "probit"), sum=TRUE){
  ## binomial first derivatives
  ## inputs: b: guess at regression parameters
  ##         Z: the data, (2*y-1)*X
  ##         link: the link function either 
  ##               logit (default) or probit
  ##         sum: a flag to return either the 
  ##           negative score (default) or Jacobian
  ## returns: either the negative score or Jacobian
  link <- match.arg(link)
  mu <- drop(Z %*% b)      
  
  if(link=="logit"){
    Jac <- Z*plogis(mu,lower.tail = FALSE)
  }else{
    p <- pnorm(mu)
    ## If it gets to the point where it tries to divide by 0, 
    ## stop it by replacing it with the smallest number that
    ## the computer knows is not 0
    p[p < .Machine$double.eps] <- .Machine$double.eps
    Jac <- Z*dnorm(mu)/p    
  }
  if(sum){
    return(-colSums(Jac))
  }else{
    return(Jac)
  }
}
H.binomial <- function(b,Z, link=c("logit", "probit")){
  ## binomial Hessian
  ## inputs: b: guess at regression parameters
  ##         Z: the data, (2*y-1)*X
  ##         link: the link function either 
  ##               logit (default) or probit
  ## returns: negative Hessian
  link <- match.arg(link)
  mu <- drop(Z %*% b)      
  
  if(link=="logit"){
    omega <- sqrt(dlogis(mu))
  }else{
    p <- pnorm(mu)
    ## If it gets to the point where it tries to divide by 0, 
    ## stop it by replacing it with the smallest number that
    ## the computer knows is not 0
    p[p < .Machine$double.eps] <- .Machine$double.eps
    ratio <- dnorm(mu)/p
    omega <- sqrt(ratio*(mu + ratio))
  }
  Zstar <- Z*omega
  Hess <- -t(Zstar) %*% Zstar
  
  return(Hess)
}

#### penalized logits ####
l.logit.pml <- function(b,Z, penalty=c("Jeffreys", "Cauchy", "logF")){
  ## bias-reduced negative log-likelihood
  ## inputs: b: guess at regression parameters
  ##         Z: the data, (2*y-1)*X
  ##         penalty: what penalty to use
  ## returns: (negative) log-likelihood 
  penalty <- match.arg(penalty)
  if(penalty=="Jeffreys"){
    ## note that determinant calculates the logged det by default so slightly 
    ## better than log(det(H))
    p <- determinant(H.binomial(b,Z))$mod[1]/2
  }else{
    if(penalty=="Cauchy"){
      p <- sum(c(dcauchy(b[1],scale=10,log=TRUE),
                 dcauchy(b[-1],scale=2.5,log=TRUE)))
    }else{
      p <- sum(b[-1]/2-log(exp(b[-1])+1))
    }
  }
  
  
  return(-sum(plogis(Z%*%b, log=TRUE))-p)
  
}

r.logit.pml <- function(b,Z, penalty=c("Jeffreys", "Cauchy", "logF")){
  ## bias-reduced negative negative score
  ## inputs: b: guess at regression parameters
  ##         Z: the data, (2*y-1)*X
  ##         penalty: what penalty to use
  ## returns: (negative) log-likelihood 
  penalty <- match.arg(penalty)
  if(penalty=="Jeffreys"){
    w <- sqrt(dlogis(drop(Z%*%b)))
    Ztilde <- Z*w
    h <- diag(Ztilde %*% solve(t(Ztilde)  %*% Ztilde) %*% t(Ztilde))
    p <- plogis(drop(Z%*%b ))
    
    return(-colSums(Z*((1-p)*(1+h/2)-p*(h/2))))
  }else{
    if(penalty=="Cauchy"){
      score <- colSums(Z*(1-plogis(drop(Z%*%b))))
      penalty <- c(-2*b[1]/(b[1]^2+100), 
                   -0.32*b[-1]/(0.16*b[-1]^2+1))
      return(-score-penalty)
    }else{
      score <- colSums(Z*(1-plogis(drop(Z%*%b))))
      penalty <- c(0,(1/2-plogis(b[-1])))
      return(-score-penalty)
    }
  }
  
  
  return(-sum(plogis(Z%*%b, log=TRUE))-p)
  ## note that determinant calculates the logged det by default so slightly 
  ## better than log(det(H))
}




### Duration ####
l.duration <- function(theta, y, X,uncen, 
                       dist=c("exponential", "Weibull"),
                       type=c("aft", "ph"),
                       sum=TRUE){
  ## Exponential and weibull duration models
  ## inputs: theta: guess at regression parameters with log(alpha) last 
  ##                if using the Weibull
  ##         y: dependent variable
  ##         X: independent variables
  ##         uncen: dummy variable indicating if the duration is
  ##                uncensored
  ##         dist: Which distribution to use 
  ##         type: AFT (target E[y|X]) or PH (target h(t|X))
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood 
  dist=match.arg(dist)
  type <- match.arg(type)
  
  if(dist=="Weibull"){
    ln.a <- theta[length(theta)]
    beta <- theta[-length(theta)]
  }else{
    ln.a <- 0
    beta <- theta
  }
  a <- exp(ln.a)
  XB <- drop(X%*%beta)
  aft <- ifelse(type=="aft", -a,1)
  lp <- aft * XB
  lam <- exp(lp)
  lS <- -lam*(y^a)
  lhaz <- log(a) + lp + (a-1)*log(y)
  ll <- uncen*lhaz + lS
  if(sum){
    return(-sum(ll))
  }else{
    return(ll)
  }
}
r.duration <-function(theta, y, X,uncen, 
                      dist=c("exponential", "Weibull"),
                      type=c("aft", "ph"),
                      sum=TRUE){
  dist=match.arg(dist)
  if(dist=="Weibull"){
    ln.a <- theta[length(theta)]
    beta <- theta[-length(theta)]
  }else{
    ln.a <- 0
    beta <- theta
  }  
  a <- exp(ln.a)
  XB <- drop(X%*%beta)
  type <- match.arg(type)
  aft <- ifelse(type=="aft", -a,1)
  lp <- aft * XB
  lam <- exp(lp)
  
  
  Dbeta <- X* (aft*uncen-aft*lam*(y^a))
  if(type=="aft"){
    Da <- uncen*(1-a*XB+a*log(y)) -a*lam*log(y)*(y^a) + a*XB*lam*(y^a)
  }else{
    Da <- uncen*(1+a*log(y)) -a*lam*log(y)*(y^a) 
    
  }
  if(dist=="Weibull"){
    J <- cbind(Dbeta, Da)
  }else{
    J <- Dbeta
  }
  if(sum){
    return(-colSums(J))
  }else{
    return(J)
  }
}
