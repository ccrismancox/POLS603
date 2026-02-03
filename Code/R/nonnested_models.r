rm(list=ls())

##### generate data from true model ####
set.seed(1)
N <- 1000
dat <-data.frame(X1=rnorm(N, mean=20, sd=4),
                 X2=runif(N),
                 X3=rlnorm(N, 5, .25),
                 X4=rpois(N, lambda = 15))
dat$X1 <- dat$X1 + .1*dat$X4
dat$X2 <- dat$X2 + .1*dat$X4
beta <- c(.5, 1, .24)
s2.i <- exp(with(dat, cbind(1,X1, X2, X3)) %*% c(-5,.2, 1, 0.01))
y <- with(dat,cbind(1, X3, X4)) %*% beta + rnorm(N,sd=sqrt(s2.i))
dat <- cbind(y,dat)

####################################################################





#### look at the data ####
summary(dat)
hist(y)

#### fit the canddiate models #### 
m1 <- lm(y~X1+X2+X3, data=dat, x=TRUE, y=TRUE)
m2 <- lm(y~X3+X4, data=dat, x=TRUE, y=TRUE)
m3 <- lm(log(y)~X1+X2+X3, data=dat, x=TRUE, y=TRUE)
m4 <- lm(log(y)~X3+X4, data=dat, x=TRUE, y=TRUE)

summary(m1)
summary(m2)
summary(m3)
summary(m4)



#### compare on information criteria #####
k1 <- length(m1$coef)+1
k2 <- length(m2$coef)+1


AIC1 <- -2*logLik(m1) + 2*k1
AIC2 <- -2*logLik(m2) + 2*k2
rbind(c(AIC(m1), AIC(m2)),
      c(AIC1, AIC2))

AIC3 <- -2*logLik(m3) + 2*k1
AIC4 <- -2*logLik(m4) + 2*k2
rbind(c(AIC(m3), AIC(m4)),
      c(AIC3, AIC4))



N <- nrow(dat)
BIC1 <- -2*logLik(m1) + log(N)*k1
BIC2 <- -2*logLik(m2) + log(N)*k2
rbind(c(BIC(m1), BIC(m2)),
      c(BIC1, BIC2))


BIC3 <- -2*logLik(m3) + log(N)*k1
BIC4 <- -2*logLik(m4) + log(N)*k2
rbind(c(BIC(m3), BIC(m4)),
      c(BIC3, BIC4))

### We can't currently compare log to non-log unless we
### use the log-normal likelihood for the log-normal models.
### As in

ln.aic3 <- -2*sum(dlnorm(dat$y, 
                         meanlog = m3$fitted.values, 
                         sdlog = sqrt(mean(m3$residuals^2)),
                         log=TRUE))+2*k1
ln.aic4 <- -2*sum(dlnorm(dat$y, 
                         meanlog = m4$fitted.values, 
                         sdlog = sqrt(mean(m4$residuals^2)),
                         log=TRUE))+2*k1
## These are comparable
c(AIC(m1), AIC(m2),
  ln.aic3, ln.aic4)



###### Vuong test #####
## we need the likelihood ratio at each point
normal.ll <- function(theta, y,X){
  s2 <- theta[length(theta)]
  beta <- theta[-length(theta)]
  mu <- X %*%  beta
  return(-log(2*pi*s2)/2-(y-mu)^2 / (2*s2))
}

log.normal.ll <- function(theta, ln.y,X){
  s2 <- theta[length(theta)]
  beta <- theta[-length(theta)]
  mu <- X %*%  beta
  return(-log(2*pi*s2)/2-ln.y-(ln.y-mu)^2 / (2*s2))
} 

s1.hat <- mean(m1$residuals^2)
s2.hat <- mean(m2$residuals^2)
s3.hat <- mean(m3$residuals^2)
s4.hat <- mean(m4$residuals^2)



LR1 <- normal.ll(c(m1$coef,s1.hat), y=m1$y, X=m1$x)
LR2 <- normal.ll(c(m2$coef,s2.hat), y=m2$y, X=m2$x)

LR3 <- log.normal.ll(c(m3$coef,s3.hat), ln.y=m3$y, X=m3$x)
LR4<- log.normal.ll(c(m4$coef,s4.hat), ln.y=m4$y, X=m4$x)

## 1 versus 2
LR <- LR1-LR2
omega2 <- mean(LR^2) - mean(LR)^2
LR <- sum(LR)
V <- LR/(sqrt(N*omega2))
V

pnorm(V) ## reject the null in favor of alt. that 2 is better 


## 3 versus 4
LR <- LR3 - LR4
omega2 <- mean(LR^2) - mean(LR)^2
LR <- sum(LR)
V <- LR/(sqrt(N*omega2))
pnorm(V) ## reject the null in favor of alt. that 2 is better 

## So far we have evidence that theory 2 is better than theory 1

## 1 versus 3
LR <- LR1 - LR3
omega2 <- mean(LR^2) - mean(LR)^2
LR <- sum(LR)
V <- LR/(sqrt(N*omega2))
pnorm(V) ## fail to reject the null of equal fit when the alt. is that 2 is better
pnorm(V, lower=F) ## reject the null in favor of alt. that 1 is better 


