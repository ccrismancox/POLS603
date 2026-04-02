library(Rcpp)
library(matrixStats)
sourceCpp('../cpp/clogit.cpp')
split.matrix <- function(x, f){
  lapply(split(x = seq_len(nrow(x)), f = f),
         function(ind) x[ind, , drop = FALSE])
  
}
l.clogit <- function(beta, X, y, id){
  N <- length(unique(id))
  Ti <- table(id)
  XB <-  drop(X%*%beta)
  eXB <-  split(exp(XB), id)
  
  k <- tapply(y, id, sum)
  lf <- f(Ti,k,eXB)
  
  ll <- sum(XB*y) -sum(lf)
  return(-ll)
}
r.clogit <- function(beta, X, y, id, sum=TRUE){
  N <- length(unique(id))
  T <- table(id)
  XB0 <- drop(X %*% beta)
  XB <- split(XB0, id)
  eXB <- split(exp(XB0), id)
  k <- tapply(y, id, sum)
  p <- ncol(X)
  Xin <- split.matrix(X, id)
  D.lf <- D_f(T, k, eXB, Xin,p)
  Xy <-  do.call(rbind,by(X*y, id, colSums))
  out <- Xy-t(D.lf)
  if(sum){
    return(-colSums(out))
  }else{
    return(out)
  }
}
vcovCL.clogit <- function(mod){
  par <- mod$coefficients
  X <- as.matrix(mod$x)
  y <- 1*(as.character(mod$y)=="1")
  id <- mod$strata
  Ji <- r.clogit(par, X, y, id,sum=FALSE)
  Meat <- crossprod(Ji)         
  bread <- vcov(mod)
  V <- bread %*% Meat %*% bread 
  rownames(V) <- colnames(V) <- names(par)
  return(V)  
}




l.relogit <- function(theta, X, y, id, GH, MC, sum=TRUE){
  s <- exp(theta[length(theta)])
  beta <- theta[-length(theta)]
  ZB <- drop(X%*%beta )*(y*2-1)
  
  if(missing(GH)){
    a <- MC/sqrt(2)
    w <- rep(1/length(a), length(a))
    GH <- FALSE
  }else{
    a <- GH$x
    w <- GH$w
    GH <- TRUE
  }
  sepZB <- split(ZB, factor(id))
  sepY <- split(y, factor(id))
  
  Pi <- mapply(x=sepZB, y=sepY,
                \(x,y){
                  f <- \(u){return(logit.f(u=u, ZB=x, Y=y, s=s, GH=GH))}
                  sum(w* f(a))
                })
  if(sum){
    return(-sum(log(Pi)))
  }else{
    return(-log(Pi))
  }
}

## Auxillary function for the likelihood to integrate
logit.f <- function(u, ZB,Y,s, GH=TRUE){
  w <- ifelse(GH,1/sqrt(pi),1)
  u <- matrix(rep(u, each=length(ZB)), nrow=length(ZB))*(Y*2-1)*sqrt(2)
  return(w *colProds(plogis(ZB+u*s)))
}
logit.d <- function(u, ZB,Y,s,Z){
  R <- length(u)
  Ti <- length(ZB)
  u <- matrix(rep(u, each=Ti), nrow=Ti)
  u <- u*(Y*2-1)*sqrt(2)
  W <-  Z %x% matrix(1, nrow=R, ncol=1)
  W <- cbind(W, c(t(u)))
  Wb <-  ZB %x% matrix(1, nrow=R, ncol=1)
  Wb <- drop(Wb +c(t(u*s)))
  
  d <- plogis(Wb,lower=FALSE)*W
  r <- rep(1:R, Ti)
  return(t(sapply(split.matrix(d, r), colSums)))
}
r.relogit<- function(theta, X, y, id, GH, MC, sum=TRUE){
  
  s <- exp(theta[length(theta)])
  beta <- theta[-length(theta)]
  Z <- (y*2-1)*X
  ZB <- drop(Z%*%beta )
  if(missing(GH)){
    a <- MC/sqrt(2)
    w <- rep(1/length(a), length(a))
    GH <- FALSE
  }else{
    a <- GH$x
    w <- GH$w
    GH <- TRUE
  }
  
  sepZB <- split(ZB, factor(id))
  sepY <- split(y, factor(id))
  sepZ <- split.matrix(Z, factor(id))
  
  J <- t(mapply(x=sepZB, y=sepY,z=sepZ,
               \(x,y,z){
                 f <- \(u){return(logit.f(u=u, ZB=x, Y=y, s=s, GH=GH))}
                 d <- \(u){return(logit.d(u=u, ZB=x, Y=y, s=s,Z=z))}
                 Pir <- (w * f(a))
                 Dir <- d(a)
                 colSums(Pir*Dir)/sum(Pir)
               }) )
  J[,ncol(J)] <- J[,ncol(J)]*s
  
  if(sum){
    return(-colSums(J))
  }else{
    return(J)
  }
}
