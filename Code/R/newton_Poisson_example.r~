#### Poisson regression with Newton's Method #####
rm(list=ls())
set.seed(1)

#### Generate data #### 
N <- 1000
X1 <- rnorm(N)
X2 <- runif(N, -2,2)
beta <- c(-2, 1, -1)
X <- cbind(1, X1, X2)
y <- rpois(N, lambda = exp(X%*%beta))
summary(y)
hist(y,freq = FALSE,breaks = "scott")


#### fuctions #####
## negative likelihood
l <- function(beta, y,X){ 
  XB <- X%*%beta
  lambda <- exp(XB)
  LL <- y*XB - lambda
  return(-sum(LL))
}

## negative score
r <- function(beta, y,X){
  XB <- drop(X%*%beta)
  lambda <- exp(XB)
  score <- X*(y-lambda)
  return(-colSums(score))
}

## negative Hessian
h <- function(beta, y,X){
  XB <- drop(X%*%beta)
  lambda <- exp(XB)
  Xtilde <- X*sqrt(lambda)
  H <- t(Xtilde) %*%Xtilde
  return(H)
}
##########################


### Newton steps ###
eps <- 1e-5
K <- 50
k <- 0
b <- rep(0, 3) ## initial guess
while(max(abs(r(b,y,X))) > eps){
  b <- b - solve(h(b,y,X)) %*% r(b,y,X)
  k <- k + 1
  b
  if(k>=K){
    break
  }
}
print(b) #estimates 
print(r(b,y,X)) #FOC

#### Standard errors #### 
#Jacobian, same as is the score, but no -colSums
J <-function(beta, y,X){ 
  XB <- drop(X%*%beta)
  lambda <- exp(XB)
  score <- X*(y-lambda)
  return(score)
}

Vh <- solve(h(b, y=y,X=X))
Vopg <- solve(t(J(b, y,X)) %*% J(b, y,X))


Vrobust <- Vh %*% solve(Vopg) %*% Vh
SE.hess <- sqrt(diag(Vh))
SE.opg <- sqrt(diag(Vopg))
SE.robust <- sqrt(diag(Vrobust))
round(cbind(SE.hess, SE.opg, SE.robust),3)





##### BFGS example ###### 
x0 <- rep(0,3) #starting values
names(x0) <- paste0("beta", 0:2) #names 

#### fit the model
fit <- optim(x0, #initial guess
             fn = l, #negative log-likelihood function
             gr=r, #gradient/score function
             y=y, #additional arguments for fn and gr
             X=X, #additional arguments for fn and gr
             method="BFGS", #optimization algorithm
             hessian=TRUE) #return the last guess at the hessian Q^{-1}
fit
h(fit$par, y,X) ## pretty close match to BFGS guess 


## est, SEs, and common z tests 
beta.hat <- fit$par
SE <- sqrt(diag(solve(fit$hessian)))
round(cbind(Est = beta.hat, 
            SE = SE,
            z.stat= beta.hat/SE,
            p.val = 2*pnorm(abs( beta.hat/SE), lower=FALSE)),
      4)


#### fitting with glm #####
df1 <- data.frame(y=y, X1=X1, X2=X2)
fit1 <- glm(y~X1+X2, 
            data=df1,
            family=poisson)
summary(fit1)



## with robust standard errors
library(sandwich)
library(lmtest)
VR <- vcovHC(fit1, type="HC0")
coeftest(fit1,VR )



###### ON Numeric derivatives ####
library(numDeriv)
rbind(grad(l,x0, y=y,X=X), #finite differences
      r(x0,y,X)) #coded derivative

system.time(replicate(1000,
                      grad(l,x0, y=y,X=X,
                           method="simple",
                           method.args=list(eps=1e-5))))
system.time(replicate(1000,r(x0,y,X)))

powers <- seq(-21,0, by=1)
err <- matrix(0, 22, 3)
true <- r(c(0,0,0), y,X)
for(i in 1:22){
    h <- 1*10^(powers[i])
    FD <- grad(l,x0, y=y,X=X,
               method="simple",
               method.args=list(eps=h))
    err[i, ] <- c(abs(FD-true))
}
plot(log(abs(err[,1]),base=10)~powers,
     col="red", type="l",
     ylab="Absolute error, powers of ten",
     xlab="Step size, powers of ten",
     ylim=c(-5.7,3))
lines(log(abs(err[,2]),base=10)~powers, col="blue")
lines(log(abs(err[,3]),base=10)~powers, col="darkgreen")
legend("bottomleft",expression(italic(D[beta[0]]),
                               italic(D[beta[1]]),
                               italic(D[beta[2]])),
       col=c("red", "blue", "darkgreen"), lty=1)





##########Overdispersion ####################
set.seed(1)

## decision process
X <- cbind(1, rnorm(5000, mean=5), rnorm(5000, mean=5))
beta <- c(2, .5, .3)
y <- exp(drop(X%*%beta)) + rnorm(5000,sd=2*sqrt(exp(X%*%beta)))
summary(y)

m2 <- glm(y~X[,-1], family=quasipoisson)
summary(m2)

fit <- optim(x0, #initial guess
             fn = l, #negative log-likelihood function
             gr=r, #gradient/score function
             y=y, #additional arguments for fn and gr
             X=X, #additional arguments for fn and gr
             method="BFGS", #optimization algorithm
             hessian=TRUE) #return the last guess at the hessian Q^{-1}
beta.hat <- fit$par
SE <- sqrt(diag(solve(fit$hessian)))

## Using just the Hessian 
round(cbind(Est = beta.hat, 
            SE = SE,
            z.stat= beta.hat/SE,
            p.val = 2*pnorm(abs( beta.hat/SE), lower=FALSE)),
      4)


## Using the glm version
u <- (y - exp(X %*% fit$par))/sqrt(exp(X %*% fit$par))
s2.hat <- mean(u^2)
SE <- sqrt(diag(s2.hat*solve(fit$hessian)))
round(cbind(Est = beta.hat, 
            SE = SE,
            z.stat= beta.hat/SE,
            p.val = 2*pnorm(abs( beta.hat/SE), lower=FALSE)),
      5)


## fully robust
Vh <- solve(fit$hessian)
Vopg <- solve(t(J(beta.hat, y,X)) %*% J(beta.hat, y,X))
SE <- sqrt(diag(Vh %*% solve(Vopg) %*% Vh))
round(cbind(Est = beta.hat, 
            SE = SE,
            z.stat= beta.hat/SE,
            p.val = 2*pnorm(abs( beta.hat/SE), lower=FALSE)),
      5)


