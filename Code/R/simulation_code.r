l <- function(beta, y,X){
  XB <- X%*%beta
  lambda <- exp(XB)
  LL <- y*XB - lambda
  return(-sum(LL))
}
r <- function(beta, y,X){
  XB <- drop(X%*%beta)
  lambda <- exp(XB)
  score <- X*(y-lambda)
  return(-colSums(score))
}
J <-function(beta, y,X){
  XB <- drop(X%*%beta)
  lambda <- exp(XB)
  score <- X*(y-lambda)
  return(score)
}

#### initial skeleton, we'll fill it from here ####
# set.seed(1)
# B <- 250 ## number of simulations
# N <- 500 ## sample size
# beta <- c(1,1) #parameters
# ## Create an empty matrix for storing results
# output <- matrix(0,
#                  nrow=B, # number of simulations
#                  ncol=8) # number of things to save
# for(b in 1:B){
#   ## Generate data
#   X <- cbind(1,rnorm(N))
#   y <- rlnorm(N, X %*%beta, 2)
#   #### DO WHAT WE NEED TO DO TO GET RESULTS####
#   ## save output of simulation b
#   output[b,] <- c(fit1$par[2], #Poisson est of beta1
#                   diag(Vh)[2], #hessian var of beta1
#                   diag(Vopg)[2], #OPG var of beta1
#                   diag(Vrobust)[2], #Robust var of beta1
#                   AME.true, #log normal AME
#                   AME.est, # estimated AME
#                   V.ame, #estimated variance of the AME
#                   time2[3]/time1[3]) #time differences
# }


#### final Monte Carlo code ####
set.seed(1)
B <- 250 ## number of simulations
N <- 500
beta <- c(1,1)
trueAME <- exp(7/2)


## Create an empty matrix for storing results
output <- matrix(0, 
                 nrow=B, # number of simulations 
                 ncol=7) # number of things to save
for(b in 1:B){
  X <- cbind(1,rnorm(N))
  y <- rlnorm(N, X %*%beta, 2)
  
  #Fit models and check times (Q1 and 5)
  time1 <- system.time({
    fit1 <- optim(rep(0,2),
                  fn=l, gr=r,
                  y=y, X=X,
                  method="BFGS",
                  hessian=TRUE)
  })
  time2 <- system.time({ #no score
    fit2 <- optim(rep(0,2),
                  fn=l, 
                  y=y, X=X,
                  method="BFGS",
                  hessian=TRUE)
  })
  
  ## find the variances (Q2)
  Vh <- solve(fit1$hessian)
  OPG <- crossprod(J(fit1$par, y=y,X=X))
  Vopg <- solve(OPG)
  V.robust <- Vh %*% OPG %*% Vh
  
  ## find the AME (Q3)
  AME.est <- fit1$par[2] * mean(exp(X %*% fit1$par))
  
  ## Find the variance of the estimated AME (Q4)
  Xstar <- X
  Xstar[,2] <- Xstar[,2] + 1/fit2$par[2]
  D.ame.beta <- colMeans(Xstar*fit1$par[2]*exp(drop(X%*%fit1$par)))
  V.ame <- D.ame.beta %*% V.robust %*% D.ame.beta
  
  ## save output of simulation b
  output[b,] <- c(fit1$par[2], #Poisson est of beta1
                  diag(Vh)[2], #hessian var of beta1
                  diag(Vopg)[2], #OPG var of beta1             
                  diag(V.robust)[2], #Robust var of beta1 
                  AME.est, # estimated AME
                  V.ame, #estimated variance of the AME
                  time2[3]/time1[3]) #time differences
}



mean(output[,1]) ## A little bit less but not bad

var(output[,1]) ## observed variance 
round(colMeans(output[,2:4]),3) ##variance estimates

trueAME
mean(output[,5]) ## slightly under estimates on average

## out of curiousity, do we do well on percentage 
## change estimates?
c(mean(exp(output[,1])-1), exp(1)-1) #yes

var(output[,5]) ## observed variance 
mean(output[,6]) ##variance estimates

## under estimates our uncertainty


mean(output[,7]) ##nearly 3 times slower with score on average




########### bootstrapping methods ##########
#### Bootstrapping examples ####
## basics
library(readstata13)
library(matrixStats)
library(MASS)
library(dplyr)
## econometric
library(lmtest)
library(sandwich)
## tables and figures
library(xtable)

rm(list=ls())


protests <- read.dta13("datasets/remittances_and_protests1.dta")
protests <- subset(protests, sample==1)


f1 <- banks_protest_count ~ remit*dict + 
  l1gdp + #GDP pc capita
  l1pop+ #population
  l1nbr5 + #neighboring protests
  l12gr+ #Economic growth
  l1migr+ #Net migration
  elec3+ #election
  factor(cowcode)-1

## Poisson: ML-Dummy variable
poi <- glm(f1, data=protests, family=poisson, y=TRUE, x=TRUE)


##  parametric bootstrap
B <- 100 ##number of bootstraps. should be larger, but this is for an example
boot.out1 <- matrix(0, ncol=2, nrow=B)
X0 <- X1 <- poi$x
X0[,"dict"] <- X0[,"remit:dict"] <- 0
X1[,"dict"] <- 1
X1[, "remit:dict"] <- X1[,"remit"]
for(b in 1:B){
  y.sim <- rpois(length(poi$y), lambda = poi$fitted.values)
  boot.mod <-  glm(update(f1, y.sim~.),
                   data=protests, family=poisson)
  poissonAME <-  c(boot.mod$coefficients['remit']*
                     mean(exp(X0 %*% boot.mod$coefficients)),
                   (boot.mod$coefficients['remit'] + 
                      boot.mod$coefficients["remit:dict"])*
                     mean(exp(X1 %*% boot.mod$coefficients)))
  
  boot.out1[b,] <-poissonAME
  
}
se.parametric <- colSds(boot.out1)

## non-parametric bootstrap
boot.out2 <- matrix(0, ncol=4, nrow=B)

dim(protests)
dim(poi$x)
for(b in 1:B){
  ## ordinary bootstrap
  id1 <-  sample(1:nrow(protests), replace=TRUE)
  boot.mod <-  glm(f1,data=protests[id1,],
                   family=poisson,x=TRUE)
  X0 <- X1 <- boot.mod$x
  X0[,"dict"] <- X0[,"remit:dict"] <- 0
  X1[,"dict"] <- 1
  X1[, "remit:dict"] <- X1[,"remit"]
  poissonAME1 <-  c(boot.mod$coefficients['remit']*
                      mean(exp(X0 %*% boot.mod$coefficients)),
                    (boot.mod$coefficients['remit'] + 
                       poi$coefficients["remit:dict"])*
                      mean(exp(X1 %*% boot.mod$coefficients)))
  
  ## clustered
  id2 <-  sample(unique(protests$cowcode), replace=TRUE)
  boot.df <- sapply(id2, \(x){subset(protests, cowcode==x)},
                    simplify=FALSE)
  boot.df <- lapply(1:length(boot.df), 
                    \(x){
                      boot.df[[x]]$cowcode <- x
                      return(boot.df[[x]])
                    })
  boot.df <- do.call(rbind, boot.df)
  boot.mod <-  glm(f1,data=boot.df,
                   family=poisson,x=TRUE)
  X0 <- X1 <- boot.mod$x
  X0[,"dict"] <- X0[,"remit:dict"] <- 0
  X1[,"dict"] <- 1
  X1[, "remit:dict"] <- X1[,"remit"]
  poissonAME2 <-  c(boot.mod$coefficients['remit']*
                      mean(exp(X0 %*% boot.mod$coefficients)),
                    (boot.mod$coefficients['remit'] + 
                       boot.mod$coefficients["remit:dict"])*
                      mean(exp(X1 %*% boot.mod$coefficients)))
  poissonAME2 
  
  boot.out2[b,] <- c(poissonAME1,poissonAME2)
  
}
se.nonparametric <- colSds(boot.out2)

## Clarify: parametric bootstrap
## Here we'll refit the model to drop the all-zero cases because
## we're relying on the estimated variance. 
## For the fixed effects on the all zero units these will also be excessively large (these parameters are unidentified)
protests2 <- protests %>% 
  mutate(sum.y= sum(banks_protest_count), .by=cowcode) %>%
  filter(sum.y > 0)
poi <- update(poi, data=protests2)

Bmat1 <- mvrnorm(B, poi$coef,vcov(poi))
Bmat2 <- mvrnorm(B, poi$coef,vcovCL(poi, protests2$cowcode))
X0 <- X1 <- poi$x
X0[,"dict"] <- X0[,"remit:dict"] <- 0
X1[,"dict"] <- 1
X1[, "remit:dict"] <- X1[,"remit"]
poissonAME1 <-  cbind(Bmat1[,'remit']*
                        colMeans(exp(X0 %*%t(Bmat1))),
                      (Bmat1[,'remit'] + Bmat1[,"remit:dict"])*
                        colMeans(exp(X1 %*%t(Bmat1))))
poissonAME2 <-  cbind(Bmat2[,'remit']*
                        colMeans(exp(X0 %*%t(Bmat2))),
                      (Bmat2[,'remit'] + Bmat2[,"remit:dict"])*
                        colMeans(exp(X1 %*%t(Bmat2))))
se.clarify <- c(colSds(poissonAME1), colSds(poissonAME2))


compMat <- rbind(parametric=c(se.parametric, NA, NA),
                 non.parametric=se.nonparametric,
                 clarify=se.clarify)
colnames(compMat) <- c("Auto", "Demo", "Auto-cluster", "Demo-cluster")
print(compMat)