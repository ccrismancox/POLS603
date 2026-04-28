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


#### IV probits ####

dbivnorm <- function(x1, x2, rho = 0){ 
  ## bivariate normal pdf. Needed for the bivariate probit score
  pi <- base::pi
  denom <- 1 - rho^2
  pdf <- 1/(2*pi * sqrt(denom)) * exp( -1/(2*denom)  * (x1^2 + x2^2 - 2 *rho * x1 *x2))
  return(pdf)
}


l.ivprobit.cont.treat <-  function(b, regr, Y, sum=TRUE){
  ## IVprobit with continuous treatment
  ## inputs: b: guess at regression parameters (gamma, beta, rho, sigma_u)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the treatment equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the treatment variable, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood  
  
  Z=regr[[1]] #treatment equation
  X=regr[[2]] #outcome including the treatment
  treatment <- Y[,1] #endogenous variable
  y <- Y[,2] #outcome
  eta <- b[length(b)]
  rho <-  tanh(eta) 
  rho <- ifelse(rho>= 1 | rho <= -1, sign(rho) * (1- .Machine$double.eps), rho)
  lns <- b[length(b)-1]
  sigma <- exp(lns)
  
  gamma <- b[1:ncol(Z)]
  beta <- b[(ncol(Z)+1):(length(b)-2)]
  
  ZG <- drop(Z %*% gamma)
  norm <- (treatment - ZG)/sigma
  m <- drop(X %*% beta + rho * norm)/sqrt(1-rho^2)
  m <- (2*y-1)*m
  
  LL <- pnorm(m, log.p=TRUE) + dnorm(norm, log = TRUE) - lns
  
  if(sum==TRUE){
    return(-sum(LL))
  }else{
    return(LL)
  }
}




r.ivprobit.cont.treat <-  function(b, regr, Y, sum=TRUE){
  ## IVprobit with continuous treatment
  ## inputs: b: guess at regression parameters (gamma, beta, rho, sigma_u)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the treatment equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the treatment variable, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative score or Jacobian
  ## returns: (negative) score or Jacobian
  
  Z=regr[[1]] #treatment equation
  X=regr[[2]] #outcome including the treatment
  treatment <- Y[,1] #endogenous variable
  y <- Y[,2] #outcome
  eta <- b[length(b)]
  rho <-  tanh(eta) 
  rho <- ifelse(rho>= 1 | rho <= -1, sign(rho) * (1- .Machine$double.eps), rho)
  lns <- b[length(b)-1]
  sigma <- exp(lns)
  
  gamma <- b[1:ncol(Z)]
  beta <- b[(ncol(Z)+1):(length(b)-2)]
  
  ZG <- drop(Z %*% gamma)
  norm <- (treatment - ZG)/sigma
  q <- sqrt(1-rho^2)
  m <- drop((X %*% beta + rho * norm)/q)
  m <- (2*y-1)*m
  p <- pnorm(m)
  d <- dnorm(m)
  
  Dgamma <- Z *(norm*q*p - rho*(2*y-1)*d)/(q*sigma*p)
  Dbeta <- X*((2*y-1)*d)/(q*p)
  
  Dsigma <- (-norm*rho*(2*y-1)*d + q*(norm^2-1)*p)/(q*sigma*p)
  Dsigma <- Dsigma*exp(lns)
  
  Drho <- (m*rho + norm*q*(2*y-1))*d/(q^2*p)
  Drho <- Drho*(1-tanh(eta)^2)
  
  
  
  J <- cbind(Dgamma, Dbeta, Dsigma, Drho)
  if(sum){
    return(-colSums(J))
  }else{
    J
  }
  
}


## binary treatment
l.ivprobit.binary <- function(b, regr, Y,sum=TRUE){
  ## IVprobit with binary treatment
  ## inputs: b: guess at regression parameters (gamma, beta, rho)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the treatment equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the treatment variable, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood  
  ## NOTE: REQUIRES THE PACKAGE pbivnorm
  
  X.treatment=regr[[1]] 
  X.outcome=regr[[2]]
  treatment <- Y[,1]
  y <- Y[,2]
  eta <- b[length(b)]
  rho <-  tanh(eta) *(2*y-1)*(2*treatment-1)
  rho <- ifelse(rho>= 1 | rho <= -1, sign(rho) * (1- .Machine$double.eps), rho)
  
  
  beta.treatment <- b[1:ncol(X.treatment)]
  beta.outcome <- b[(ncol(X.treatment)+1):(length(b)-1)]
  
  mu.treatment <- drop((2*treatment-1)*(X.treatment %*% beta.treatment))
  mu.outcome <- drop((2*y-1)*(X.outcome %*% beta.outcome))
  
  LL <-  pbivnorm(mu.treatment, mu.outcome, rho) 
  
  LL[LL<.Machine$double.eps] <- .Machine$double.eps ## avoid 0 likelihood
  LL <- log(LL) 
  if(sum){
    return(-sum(LL))
  }else{
    return(LL)
  }
}




r.ivprobit.binary <- function(b, regr, Y,sum=TRUE){
  ## IVprobit with binary treatment
  ## inputs: b: guess at regression parameters (gamma, beta, rho)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the treatment equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the treatment variable, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative score or Jacobian
  ## returns: (negative) score or Jacobian
  ## NOTE: REQUIRES THE PACKAGE pbivnorm and the dbivnorm function (above)
  ## Reminder:
  ## the derivative of \Phi_2(x, y, rho) w.r.t x is 
  ## phi(x) * \Phi((y-x*rho)/sqrt(1-rho^2)) 
  
  X.treatment=regr[[1]] 
  X.outcome=regr[[2]]
  treatment <- Y[,1]
  y <- Y[,2]
  eta <- b[length(b)]
  rho <-  tanh(eta) *(2*y-1)*(2*treatment-1)
  rho <- ifelse(rho>=1 | rho <= -1, sign(rho) * (1- .Machine$double.eps), rho)
  
  beta.treatment <- b[1:ncol(X.treatment)]
  beta.outcome <- b[(ncol(X.treatment)+1):(length(b)-1)]
  
  mu.treatment <- drop((2*treatment-1)*(X.treatment %*% beta.treatment))
  mu.outcome <- drop((2*y-1)*(X.outcome %*% beta.outcome))
  
  p <- pbivnorm(mu.treatment, mu.outcome, rho)
  p[p<=0] <- .Machine$double.eps
  
  Db1 <- X.treatment*drop((2*treatment-1)*(dnorm(mu.treatment)* pnorm((mu.outcome-mu.treatment*rho)/sqrt(1-rho^2))))/p
  Db2 <- X.outcome*drop((2*y-1)*(dnorm(mu.outcome)*pnorm((mu.treatment-mu.outcome*rho)/sqrt(1-rho^2))))/p
  Deta <- (dbivnorm(mu.treatment, mu.outcome,rho) * (1-tanh(eta)^2) *(2*y-1)*(2*treatment-1))/p
  if(sum){
    return(-colSums(cbind(Db1, Db2,Deta)))
  }else{
    return(cbind(Db1, Db2,Deta))
    
  }
}




##### Selection ####

l.heckit <- function(theta, Y, regr, sum=TRUE){
  ## Selection model with continuous outcome
  ## inputs: b: guess at regression parameters (gamma, beta, sigma^2, rho)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the treatment equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the treatment variable, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or the vector of
  ##              individual log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood   
  rho <- tanh(theta[length(theta)])
  s2 <-  exp(theta[length(theta)-1])
  Z <- regr[[1]]
  X <- regr[[2]]
  s <- Y[,1]
  y <- Y[,2]
  gamma <- theta[1:ncol(Z)]
  beta <- theta[(ncol(Z)+1):(ncol(Z)+ncol(X))]
  e <- drop(y-X %*% beta)
  ZG <- drop(Z %*%gamma)
  
  mu <- (rho*e/sqrt(s2) + ZG)/sqrt(1-rho^2)
  LL <- ifelse(s==1,
               -log(s2)/2 + pnorm(mu, log.p=TRUE) - 
                 log(pi)/2 - log(2)/2 -
                 e^2 / (2*s2),
               pnorm(ZG, lower=FALSE, log.p=TRUE))
  if(sum){
    return(-sum(LL))
  }else{
    return(LL)
  }
}


r.heckit <- function(theta, Y, regr, sum=TRUE){
  ## Selection model with continuous outcome
  ## inputs: b: guess at regression parameters (gamma, beta, sigma^2, rho)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the treatment equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the treatment variable, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative score or Jacobian
  ## returns: (negative) score or Jacobian
  rho <- tanh(theta[length(theta)]) #tanh(eta)
  s2 <-  exp(theta[length(theta)-1]) #exp(tau)
  Z <- regr[[1]]
  X <- regr[[2]]
  s <- Y[,1]
  y <- Y[,2]
  gamma <- theta[1:ncol(Z)]
  beta <- theta[(ncol(Z)+1):(ncol(Z)+ncol(X))]
  e <- drop(y-X %*% beta)
  ZG <- drop(Z %*%gamma)
  mu <- (rho*e/sqrt(s2) + ZG)/sqrt(1-rho^2)
  
  phat <- pnorm(mu)
  phat[phat < .Machine$double.eps] <- .Machine$double.eps
  p0 <- pnorm(mu,lower=FALSE)
  p0[p0 < .Machine$double.eps] <- .Machine$double.eps
  Dbeta <- X*ifelse(s==1,
                    e/s2 - rho*dnorm(mu)/(sqrt(s2*(1-rho^2))*phat),
                    0)
  Dgamma <- Z*ifelse(s==1,
                     dnorm(mu)/(sqrt(1-rho^2)*phat),
                     -dnorm(mu)/(p0))
  
  
  Dtau <- ifelse(s==1,
                 ( e^2 / (2*s2) - 1/2 - (rho*e*dnorm(mu))/(2*sqrt(s2*(1-rho^2))*phat)),
                 0)
  Deta <- ifelse(s==1,
                 (e+ZG*rho*sqrt(s2))*dnorm(mu)/(sqrt(s2*(1-rho^2))*phat),
                 0)
  J <- cbind(Dgamma,Dbeta, Dtau, Deta)
  
  if(sum){
    return(-colSums(J))
  }else{
    return(J)
  }
}


l.biprobit.selection <- function(b, regr, y,sum=TRUE){
  ## Sample selection bivariate probit
  ## inputs: b: guess at regression parameters (gamma, beta, rho)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the selection equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the selection indicator, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood  
  ## NOTE: REQUIRES THE PACKAGE pbivnorm
  
  
  Z.selection=regr[[1]] 
  X.outcome=regr[[2]]
  selection <- y[,1]
  y <- y[,2]
  eta <- b[length(b)]
  rho <-  tanh(eta) 
  
  
  beta <- b[1:ncol(X.selection)]
  gamma <- b[(ncol(X.selection)+1):(length(b)-1)]
  
  mu.selection <- drop((X.selection %*% beta))
  mu.outcome <- (2*y-1)*drop((X.outcome %*% gamma))
  rho.star <- (2*y-1)*rho
  
  LL <-  ifelse(selection==0,
                pnorm(-mu.selection),
                pbivnorm(mu.outcome, mu.selection, rho.star))
  
  LL[LL<.Machine$double.eps] <- .Machine$double.eps
  LL <- log(LL)
  if(sum){
    return(-sum(LL))
  }else{
    return(LL)
  }
}


r.biprobit.gr <- function(b, regr, y,sum=TRUE){
  ## Sample selection bivariate probit
  ## inputs: b: guess at regression parameters (gamma, beta, rho)
  ##         regr:  A length 2 list where 
  ##                item 1 is matrix for the selection equation (Z)
  ##                item 2 is the outcome regressors (X)
  ##          Y: A two column matrix. Column 1 is the selection indicator, 
  ##             2 is the outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ##         sum: a flag for whether to return the 
  ##              summed negative score or Jacobian
  ## returns: (negative) score or Jacobian
  ## NOTE: REQUIRES THE PACKAGE pbivnorm and dbivnorm function (defined above)
  Z=regr[[1]] 
  X=regr[[2]]
  d <- y[,1]
  y <- y[,2]
  eta <- b[length(b)]
  rho <-  tanh(eta) 
  
  
  gamma <- b[1:ncol(Z)]
  beta <- b[(ncol(Z)+1):(length(b)-1)]
  
  pi <- drop((Z %*% gamma))
  mu <- (2*y-1)*drop((X %*% beta))
  r <- (2*y-1)*rho
  p <- pbivnorm(pi, mu, r)
  p[p<.Machine$double.eps] <- .Machine$double.eps
  
  xi.pi <- (mu-pi*r)/sqrt(1-r^2)
  xi.mu  <- (pi-mu*r)/sqrt(1-r^2)
  
  
  Dgamma <- Z*ifelse(d==1,
                     pnorm(xi.pi)/p,
                     -1/pnorm(pi, lower=FALSE)
  )*dnorm(pi)
  Dbeta <- X*(2*y-1)*d*pnorm(xi.mu)*dnorm(mu)/p
  Drho <- d*(2*y-1)*dbivnorm(mu, pi, r)/p
  Deta <- Drho*(1-rho^2)
  if(sum){
    return(-colSums(cbind(Dgamma, Dbeta,Deta)))
  }else{
    return(cbind(Dgamma, Dbeta,Deta))
    
  }
}


#### Panel logit ####
#### See clogit_functions.R for panel logit functions 

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
  ##              summed negative score or Jacobian
  ## returns: (negative) score or Jacobian
  
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




makeR <- function(y,uncen,X){
  ## Helper function to prepare data for Cox PH 
  ## inputs: y: outcome variable
  ##         uncen: an indicator for uncensored observations
  ##         X: N by k matrix of regressors 
  ## outputs: A list with 3 elements 
  ##          R: a list of which observations are in the 
  ##             in the risk set for each failure time
  ##          d: a vector denoting the number of failures at each time
  ##          Xd: a list containing the rows of X for each exit time
  idx=1:length(y)
  d <- c(table(y[uncen==1]))
  fails <- sort(unique(y[y%in% as.numeric(names(d))]))
  
  #create the risk set for each observed failure by finding the observations
  #still in play at each observed failure time
  R <- sapply(fails, \(t){idx[y >= t ]})
  
  ## Extract the X's for each observed failure
  Xd <- sapply(fails, \(t){X[idx[y == t & uncen==1], ,drop=FALSE]})
  
  ## To get a Jacobian we need to flip it so we know which risk sets
  ## each observation is in
  Robs <- sapply(y, \(t){idx[y >=t ]})
  
  
  return(list(R=R, d=d, Xd=Xd, Robs=Robs, uncen=uncen))
}


l.cox <- function(b, X, Robj, sum=TRUE){
  ## Cox PH with Breslow's method for ties
  ## inputs: b: guess at regression parameters 
  ##         X: N by k matrix of regressors 
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood  
  ## NOTES: Requires an Robj, supplied by makeR (above)
  R <- Robj$R
  d <- Robj$d
  Xd <- Robj$Xd
  
  XB <- drop(X%*%b)
  Xbd <- sapply(Xd, \(x){sum(x%*%b)})
  XBrisk <- d*sapply(R, \(x){logSumExp(XB[x])})
  ll <- Xbd - XBrisk
  if(sum){
    return(-sum(ll))
  }else{
    return(ll)
  }
}



r.cox <- function(b, X, Robj,sum=TRUE){
  ## Cox PH with Breslow's method for ties
  ## inputs: b: guess at regression parameters 
  ##         X: N by k matrix of regressors 
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood  
  ## NOTES: Requires an Robj, supplied by makeR (above)
  R <- Robj$R
  d <- Robj$d
  Xd <- Robj$Xd

  
  XB <- drop(X%*%b)
  Db1 <- t(sapply(Xd, \(x){colSums(x)}))
  XBrisk <- d*t(sapply(R, \(x){colSums(X[x,]*exp(XB[x]))/sum(exp(XB[x]))}))
  Dbeta <- Db1 - XBrisk
  if(sum){
    return(-colSums(Dbeta))
  }else{
    ### Note this is not a proper Jacobian because the 
    ## likelihood is a partial likelihood
    ## This is however what you would use in place
    ## of the Jacobian in a sandwich estimator
    
    ## Taken from Lin and Wei (1989, JASA)
    uncen <- Robj$uncen
    N <- length(Robj$uncen)
    R0 <- sapply(Robs,\(r){sum(exp(XB[r]))})
    R1 <- t(sapply(Robs,\(r){colSums(exp(XB[r])*X[r,,drop=FALSE])}))
    J1 <- X-(R1/R0) 
    
    J <- uncen*J1- sapply(Robs,\(r){sum((uncen[r]*exp(XB[r]))/(N*R0[r]))})*J1 
    return(J)
  }
 
  
}




#### Ordered logit/probit #####

l.ordered <- function(theta,y, X, link=c("logit", "probit")){
  ## Ordered logit 
  ## inputs: theta: guess at regression parameters (beta, tau_1, delta)
  ##         X: N by k matrix of regressors 
  ##         y: outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood  
  ## NOTES: y is numeric here and has to start at 0, so should be
  ## 0, 1, 2
  link <- match.arg(link)
  G <- switch(link, #assign G based on link
              "logit"=plogis,
              "probit"=pnorm)
  beta <- theta[1:ncol(X)]
  t0  <- theta[(ncol(X)+1):length(theta)]
  K <- length(unique(y))
  t0[-1] <-exp(t0[-1])
  tau <- c(-Inf, cumsum(t0), Inf)
  XB <- drop(X %*% beta)
  
  take <- sparseMatrix(i=rep(1:nrow(X),2), 
                       j=c(y+1,y+2), x=1)
  cutoffs <- as.vector(tau*t(take))
  cutoffs <- t(matrix(cutoffs[!is.na(cutoffs) & ! cutoffs==0], nrow=2))
  p.hat <- G(cutoffs[,2]-XB) - G(cutoffs[,1]-XB)
  lp <- log(p.hat)
  return(- sum(lp))
}
r.ordered <- function(theta,y, X,link=c("logit", "probit"), sum=TRUE){
  ## Ordered logit 
  ## inputs: theta: guess at regression parameters (beta, tau_1, delta)
  ##         X: N by k matrix of regressors 
  ##         y: outcome variable
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ##         sum: a flag for whether to return the 
  ##              summed negative score or Jacobian
  ## returns: (negative) score or Jacobian
  ## NOTES: y is numeric here and has to start at 0, so should be
  ## 0, 1, 2  link <- match.arg(link)
  G <- switch(link, #assign G based on link
              "logit"=plogis,
              "probit"=pnorm)
  g <- switch(link, #assign g based on link
              "logit"=dlogis,
              "probit"=dnorm)
  
  beta <- theta[1:ncol(X)]
  t0  <- theta[(ncol(X)+1):length(theta)]
  t0[-1] <-exp(t0[-1])
  K <- length(unique(y))
  tau <- c(-Inf, cumsum(t0), Inf)
  XB <- drop(X %*% beta)
  take <- sparseMatrix(i=rep(1:nrow(X),2), 
                       j=c(y+1,y+2), x=1)
  cutoffs <- as.vector(tau*t(take))
  cutoffs <- t(matrix(cutoffs[!is.na(cutoffs) & ! cutoffs==0], nrow=2))
  
  if(link=="logit"){
    Dbeta <- X*(G(cutoffs[,1]-XB) +G(cutoffs[,2]-XB)-1)
    Dt <- -G(cutoffs[,1]-XB) -G(cutoffs[,2]-XB)+1
  }else{
    Dbeta <- X*((g(cutoffs[,2]-XB) -g(cutoffs[,1]-XB))/(G(cutoffs[,1]-XB) -G(cutoffs[,2]-XB)))
    Dt <- (g(cutoffs[,2]-XB) -g(cutoffs[,1]-XB))/(G(cutoffs[,2]-XB) -G(cutoffs[,1]-XB))
  }
  
  for(k in 1:(K-2)){
    Dtk <-(y>=k)*((g(cutoffs[,2]-XB)-(y>k)*g(cutoffs[,1]-XB))/(G(cutoffs[,2]-XB)-G(cutoffs[,1]-XB)))
    Dtk <- Dtk*t0[k+1]
    Dt <- cbind(Dt,Dtk)
  }
  J <- cbind(Dbeta, Dt)
  if(sum){
    return(-colSums(J))
  }else{
    return(J)
  }
  
}


#### Multinomial Logit ####
l.mnl <- function(theta, y, X, Z, W, sum=TRUE){
  ## Multinomial logit 
  ## inputs: theta: guess at the parameters (beta_1, ..., beta_{J-1}, gamma)
  ##         X: N by k matrix of regressors that do not vary with choice
  ##         Z: N by l*J matrix of regressors that do vary with choice. 
  ##             Note that the columns of Z should be 
  ##             [Z_{10}, ... Z_{1J-1}, Z_{20},..,Z_{2J-1},...,Z_{l0}, ... Z_{lJ-1}]
  ##         y:  NxJ matrix of binary outcomes
  ##         W: unused in this function, but used for the score
  ##         sum: a flag for whether to return the 
  ##              summed negative log-likelihood or 
  ##              the vector of log-likelihoods
  ## returns: (negative) log-likelihood or 
  ##          vector of log-likelihood
  
  
  J <- ncol(y) # number of observed choices
  ngamma <- ncol(Z)/J #number of gammas
  gamma <- theta[(length(theta)-ngamma+1):length(theta)]
  
  ## normalize beta[0]=0
  beta <- cbind(0, matrix( theta[1:(length(theta)-ngamma)], ncol=J-1 ))
  
  V <- X%*% beta + # individual-level payoffs
    matrix(matrix(Z, ncol=ngamma) %*% gamma, ncol=J) #  action-specific payoffs
  
  ## logged choice probs
  vmax <- rowMaxs(V) #from matrixStats
  lp <- V- log(rowSums(exp(V-vmax)))-vmax
  
  ## negative likelihood
  if(sum){
    return(-sum(lp*y))
  }else{
    return(rowSums(lp*y))
  }
  
}



split.matrix <- function(x, f){
  lapply(split(x = seq_len(nrow(x)), f = f),
         function(ind) x[ind, , drop = FALSE])
}

r.mnl <- function(theta, y,X,Z, W, sum=TRUE){
  ## Multinomial logit 
  ## inputs: theta: guess at the parameters (beta_1, ..., beta_{J-1}, gamma)
  ##         X: N by k matrix of regressors that do not vary with choice
  ##         Z: N by l*J matrix of regressors that do vary with choice. 
  ##             Note that the columns of Z should be 
  ##             [Z_{10}, ... Z_{1J-1}, Z_{20},..,Z_{2J-1},...,Z_{l0}, ... Z_{lJ-1}]
  ##         y:  NxJ matrix of binary outcomes
  ##         W: a (possibly sparse) matrix that combines X and Z as 
  ##            described in the notes
  ##         sum: a flag for whether to return the 
  ##              summed negative score or Jacobian
  ## returns: (negative) score or Jacobian
  ## NOTES: relies on split.matrix (above) 
  
  J <- ncol(y) # number of observed choices
  N <- nrow(y)
  ngamma <- ncol(Z)/J #number of gammas
  gamma <- theta[(length(theta)-ngamma+1):length(theta)]
  
  ## normalize beta[0]=0
  beta <- cbind(0, matrix( theta[1:(length(theta)-ngamma)], ncol=J-1 ))
  
  V <- X%*% beta + # individual-level payoffs
    matrix(matrix(Z, ncol=ngamma) %*% gamma, ncol=J) #  action-specific payoffs
  
  
  P <- exp(V)/(rowSums(exp(V)))
  
  if(sum){
    Jac <- matrix(colSums(W* (matrix(1, 1, ncol(W)/J)) %x%( y-P)), ncol=ncol(W)/J)[,-c(1:ncol(X))]
    return(-colSums(Jac))
  }else{
    J1 <- W* (matrix(1, 1, ncol(W)/J)) %x%( y-P) #setup the same matrix
    Dbeta <- J1[,(ncol(X)*J+1):((ncol(X)*J*J))] #discard gamma and beta0
    keepDbeta <- apply(Dbeta,2,\(x){all(x==0)}) # drop the all 0s
    Dbeta <- Dbeta[,!keepDbeta] #drop the all 0s
    Dgamma <- J1[,(ncol(X)*J*J+1):ncol(J1)]
    Dgamma <- split.matrix(t(Dgamma), rep(1:ngamma,each=J))
    Dgamma <- sapply(Dgamma, colSums)
    Jac <- cbind(Dbeta, Dgamma)
    return(Jac)
  }
}
