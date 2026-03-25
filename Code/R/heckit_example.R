library(sandwich)
library(lmtest)
library(GJRM)
library(MASS)
rm(list=ls())
mroz87 <- read.csv("datasets/Mroz87.csv")

## convert to log wages, but leave as 0 for unselected by convention 
mroz87$lnwage <- ifelse(mroz87$lfp==1, log(mroz87$wage), 0)

## regular OLS
ols <- lm(lnwage ~ exper + I(exper^2) + educ, 
          subset=lfp==1,
          data=mroz87)
coeftest(ols, vcovHC)

## what do we see 
# On average, a 1 year increase in education increases 
# hourly wages by about 11%, holding experience fixed

## On average, a 1 year increase in experience increases 
## hourly wages by about [0.04 -0.002 (years of experience)]*100 %
## hold education fixed.

## For polynomials we may want to plot the fitted values and marginals
## we'll keep it rough for now
par(mfrow=c(1,2))


## fitted values
ed.bar <- mean(mroz87[mroz87$lfp==1,]$educ)
## log normal expected value
Xstar <- cbind(1, 0:30,(0:30)^2, ed.bar)
Ey <- exp(Xstar %*% ols$coefficients + mean(ols$residuals^2)/2)
Db <- drop(Ey)*Xstar
Vey <- Db %*% vcovHC(ols) %*% t(Db)
se.Ey <- sqrt(diag(Vey))
lo <- Ey - 1.96*se.Ey
hi <- Ey +1.96*se.Ey
plot(Ey, x=0:30, xlab="Years of experience", 
     ylab="Expected hourly wages", 
     type="l", ylim=c(2,5.6))
lines(lo, x=0:30, lty="dotted",col="blue")
lines(hi, x=0:30, lty="dotted",col="blue")



## marginals
mar <-100*( ols$coefficients["exper"] + 2*(0:30)*ols$coefficients["I(exper^2)"])
v <- 10000 *(vcovHC(ols)["exper","exper"] + 
               (2*(0:30))^2 * vcovHC(ols)["I(exper^2)","I(exper^2)"] +
               2*2*(0:30) * vcovHC(ols)["exper","I(exper^2)"])
lo <- mar-1.96*sqrt(v)
hi <- mar+1.96*sqrt(v)
plot(mar, x=0:30, xlab="Years of experience", 
     ylab="% increase in hourly wages", type="l", ylim=c(-3,7.2))
lines(lo, x=0:30, lty="dotted",col="blue")
lines(hi, x=0:30, lty="dotted",col="blue")
abline(h=0,lty="dashed")

## When experience is low, the effect of increasing it leads
## to larger changes in wages. This levels off around 20ish years
## deminishing returns


## The heckman two-step
step1 <- glm(lfp ~ age + I( age^2 ) + kids5 + educ,
             x=TRUE, y=TRUE,
             family=binomial("probit"), data=mroz87)
V1 <- vcovHC(step1)
coeftest(step1, V1)

## education and kids have opposing effects here 
## Age we don't see anything on the surface, but note that
linearHypothesis(step1, c("age", "I(age^2)"), vcov=V1)
## They are jointly significant, we reject the null that they're both 0.


mroz87$lam.hat <- dnorm(step1$linear.predictors)/step1$fitted.values
step2 <- lm(lnwage ~ exper + I( exper^2 ) + educ +lam.hat,
            x=TRUE, y=TRUE,
            subset=lfp==1,
            data=mroz87)      


bigX <- step2$x
Z <- model.matrix(~ age + I( age^2 )  + kids5 + educ, 
                  data=mroz87[mroz87$lfp==1,])

delta <- with(mroz87[mroz87$lfp==1,], lam.hat*(lam.hat+Z%*%step1$coef))
s2e <- mean(step2$residuals^2) + mean(delta)*step2$coefficients["lam.hat"]^2
rho.hat <- step2$coefficients["lam.hat"]/sqrt(s2e)
DELTA <- diag(1-rho.hat^2 *drop(delta))
A <- t(bigX) %*% DELTA %*% Z
bread <- solve(t(bigX) %*% bigX)
pb <- t(bigX) %*% (diag(nrow(Z))-rho.hat^2*DELTA)%*% bigX
jam <- rho.hat^2 * (A %*% V1 %*% t(A))
V.heckit <- s2e * bread %*% (pb+jam) %*% bread
coeftest(step2, V.heckit)


## We don't actually find much evidence of selection bias
## Note that $\beta_\lambda$ is small and not distinguishable from 0
## estimates are largely unchanged. Effect of experience will be
## nearly identical
## What about the effect of education?

ame.ed <- step2$coefficients['educ']- 
  step1$coefficients['educ']*step2$coefficients['lam.hat']*mean(delta)
ame.ed
# about the same, lol

## What about the FIML
## changes from before
## model is now BSS for bivariate sample selection
## the selection equation must go first in the formulas
fiml <- gjrm(list(lfp ~ log(huswage) + educ 
                  + age +I(age^2)
                  +kids5, 
                  lnwage ~ educ + exper + I( exper^2 )),
             data=mroz87,
             model="BSS",
             margins=c("probit", "N"))

## another difference, for robust SE rather than cluster, 
## just give it all the rows as the identifier
fiml <- adjCov(fiml, id=1:nrow(mroz87))
summary(fiml)

## very similar again. One thing to note here: sigma2 in the summary
## is sigma of the second equation not sigma^2. Confusing!
gamma.hat <- fiml$coefficients[1:ncol(fiml$X1)]
beta.hat <- fiml$coefficients[(ncol(fiml$X1)+1):(ncol(fiml$X1)+ncol(fiml$X2))]
ZG.fiml <- drop(fiml$X1 %*% gamma.hat)
lam.fiml <- dnorm(ZG.fiml)/pnorm(ZG.fiml)
delta.fiml <- lam.fiml*(lam.fiml+ZG.fiml)
ame.ed2 <- beta.hat["educ"] - gamma.hat["educ"]*fiml$sigma*fiml$theta*mean(delta.fiml)
ame.ed2

## standard error: parametric bootstrap
#first is gamma_educ, second is beta_educ
which(names(fiml$coefficients)=="educ") 
Bmat <- mvrnorm(1000, mu=fiml$coefficients, Sigma = fiml$Vb)
ZG.fiml <- fiml$X1 %*% t(Bmat[,1:ncol(fiml$X1)])
lam.fiml <- dnorm(ZG.fiml)/pnorm(ZG.fiml)
delta.fiml <- lam.fiml*(lam.fiml+ZG.fiml)
ame.sim <- Bmat[,8] - Bmat[,3]*exp(Bmat[,"sigma.star"])*tanh(Bmat[,"theta.star"])*colMeans(delta.fiml)
sd(ame.sim)

## nearly identical to the OLS st err as we may expect
## given how everything else suggests no real problem 
## here.