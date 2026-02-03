rm(list=ls())

#' Suppose that $y_1, \ldots, y_N$ are iid observations from a N(1,1).
#' We will demonstrate how to test hypothesis relating to $\mu$ and $\sigma^2$ 
#' using both likelihood ratio tests and Wald tests.
#' 
#' First example
#' $$H_0: \mu=0.5$$
#' against the alternative that 
#' $$H_A: \mu\neq 0.5.$$
#' 
#' For the LR test this means we have the log-likelihood for full model.
#' The full model imposes no restrictions so we just use our formulas for the MLE 
#' in the these two cases
#' \begin{align*}
#' \hat{\mu} &= \bar{y}\\
#' \hat{\sigma^2} &= \frac{1}{N}\sum_{i=1}^N(y_i-\bar{y})^2.
#' \end{align*}
#' 
#' We also need the restricted model. 
#' Here we fit the model as if the null hypothesis is true,
#' so we fix $\mu=0.5$ and estimate 
#' $$\hat{\sigma^2}_R = \frac{1}{N}\sum_{i=1}^N(y_i-0.5)^2.$$
#' 
#' We then use these two parameters vectors to find the full and restricted log likelihoods.

### Generate data ###
set.seed(1)
y <- rnorm(10, 1,1)

mu.mle <- mean(y) ## MLE of mu
s2.mle <- mean((y-mu.mle)^2) ## unrestricted MLE of sigma^2

s2.mle.R <- mean((y-0.5)^2) ## restricted MLE of sigma^2

## LR test
L.unrestricted <- sum(dnorm(y, mu.mle, sqrt(s2.mle), log=TRUE)) #at hat(theta)
L.unrestricted

L.restricted <- sum(dnorm(y, 0.5, sqrt(s2.mle.R), log=TRUE)) #at hat(theta)_R
L.restricted

## test statistic with 1 degree of freedom (1 restriction in the null hypothesis)
LR.test <- 2*(L.unrestricted-L.restricted)
LR.test
pchisq(LR.test, lower=FALSE, df=1)


#' Now for the Wald test, we don't need the restricted model, but we do need 
#' the variance-covariance matrix of the MLEs. 
#' Recall we have three ways to estimate this, 
#' 1. Expected value of the negative Hessian
#' 2. Negative observed Hessian
#' 3. Outer product of gradients (OPG)


## vcov 1: (-E(H(\theta|y)))^{-1}
E.hessian <- matrix(c(-10/s2.mle, 0,0,-10/(2*s2.mle^2)), 2,2)
V1 <- solve(-E.hessian)

## vcov 2: (-(H(\hat\theta|y)))^{-1}

Obs.hessian <- matrix(c(-10/s2.mle, (10*mu.mle-sum(y))/(s2.mle^2),
                    (10*mu.mle-sum(y))/(s2.mle^2), -10/(2*s2.mle^2)), 2,2)
V2 <- solve(-Obs.hessian)

## Vcov 3: OGP
Jac <- cbind((y-mu.mle)/s2.mle,
             -1/(2*s2.mle)+(y-mu.mle)^2 / (2*s2.mle^2))
## In Jac, the rows are score at each observation.
OPG <- t(Jac) %*% Jac
V3 <- solve(OPG)


list(V1, V2, V3)

wald.stat <- (mu.mle -0.5)/sqrt(V1[1,1])
wald.stat
pnorm(abs(wald.stat), lower=FALSE)*2








#' 
#' Second example
#' $$H_0: \mu=0 \& \sigma^2=1$$
#' against the alternative that 
#' $$H_A: \mu\neq 0 \text{ OR } \sigma^2 \neq 1$$
#' 
#' The unrestricted model doesn't change and now the restricted parameter vector
#' includes both $mu$ and $sigma^2$ so we don't need to refit the model
#' 
## LR test
L.unrestricted <- sum(dnorm(y, mu.mle, sqrt(s2.mle), log=TRUE)) #at hat(theta)
L.unrestricted

L.restricted <- sum(dnorm(y, 0, 1, log=TRUE)) #at hat(theta)_R; no refit required
L.restricted

LR.test <- 2*(L.unrestricted-L.restricted)
LR.test
pchisq(LR.test, lower=FALSE, df=2)


## Wald setup
theta.hat <- c(mu.mle, s2.mle)
A <- diag(2)
b <- c(0,1)

## with the Observed Hessian
wald <- t(A%*% theta.hat -b) %*% solve(A %*% V2 %*% A) %*% (A%*% theta.hat -b) 
wald
pchisq(wald, lower=FALSE, df=2)

## With the OPG
wald <- t(A%*% theta.hat -b) %*% solve(A %*% V3 %*% A) %*% (A%*% theta.hat -b) 
wald
pchisq(wald, lower=FALSE, df=2)







c.theta <- exp(mu.mle + s2.mle/2)-1
Dc.theta <- c(exp(mu.mle + s2.mle/2), exp(mu.mle+s2.mle/2)/2)
Wald <- c.theta %*% solve(Dc.theta %*% V2 %*% Dc.theta) %*% c.theta
pchisq(Wald, lower=FALSE, df=1)
