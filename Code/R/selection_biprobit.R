library(GJRM)
library(pbivnorm)
library(readstata13)
library(sandwich)
library(lmtest)
library(ggplot2)
library(MASS)
library(matrixStats)

rm(list=ls())

bcg<- read.dta13("BCG_merged.dta",
                 convert.dates = TRUE, ## avoid warnings
                 convert.underscore = FALSE,## avoid warnings
                 nonint.factors = FALSE) ## avoid warnings
m1 <- glm(cac_ons35_b~mgini_intx+  max_rdisc +
            maxhighx + maxlowx+
            downatall + powershare + 
            sip2l + lpopl + lgdpcapl + ethfrac+
            l_cac_inc35_b,
          subset=cci_ons_b==1,
          data=bcg, 
          family=binomial("probit"),
          x=T,y=T)
V1 <-  vcovCL(m1, cluster=bcg$cowcode[bcg$cci_ons_b==1])
round(coeftest(m1, V1),3)
## marginal
dhat <- dnorm(m1$linear.predictors)
AME1 <- m1$coef["max_rdisc"] *mean(dhat)

## SE of marginal
phat <- m1$fitted.values
D <- colMeans(m1$coef["max_rdisc"]*dhat*m1$x * (-phat + 2*dhat)/phat)
D["max_rdisc"] <- D["max_rdisc"]+ mean(dhat)
var.ame1 <-  D %*% V1 %*% D
se.ame1 <- sqrt(var.ame1)
c(AME1, se.ame1)


## Selection
m2 <- gjrm(list(cci_ons_b ~ mgini_intx + max_rdisc + 
                  maxhighx + maxlowx+
              downatall + powershare + sip2l + 
                lpopl + lgdpcapl+ ethfrac + l_cci_inc_b,
            cac_ons35_b ~ mgini_intx + 
              max_rdisc + maxhighx + maxlowx+
              downatall + powershare + sip2l + 
              lgdpcapl + ethfrac + l_cac_inc35_b),
       model="BSS",
       margins=c("probit", "probit"),
       data=bcg
  )

obs.sample <- as.numeric(row.names(m2$X1))


m2 <- adjCov(m2,
             id = bcg[obs.sample,]$cowcode)
summary(m2)

## Set up for the marginal
## NOTE: m2$X2 is X but only for the selected
##       We want all of it, that's in m2$X2s
l <- ncol(m2$X1)
k <- ncol(m2$X2)
beta <- m2$coefficients[(l+1):(l+k)]
gamma <- m2$coefficients[1:l]
rho <- m2$theta
bhat <- beta["max_rdisc"]
ghat <- gamma["max_rdisc"]
XB <- drop(m2$X2s %*% beta)
dhatXB <- dnorm(XB)
phatXB <- pnorm(XB)
ZG <- drop(m2$X1 %*% gamma)
dhatZG <- dnorm(ZG)
phatZG <- pnorm(ZG)
mu1 <- (ZG-XB*rho)/sqrt(1-rho^2)
mu2 <- (XB-ZG*rho)/sqrt(1-rho^2)

part1 <- bhat*dhatXB * pnorm(mu1)
part2 <- ghat*dhatZG * pnorm(mu2)
part3 <- ghat*pbivnorm(XB, ZG, rho)*dhatZG
AME2 <- mean((part1+part2)/phatZG - part3/(phatZG^2))


## We'll use the parametric bootstrap for the 
## standard error
Tmat <- mvrnorm(1000, m2$coefficients, m2$Vb)
beta.sim <- Tmat[,(l+1):(l+k)]
gamma.sim <- Tmat[,1:l]
rho.sim <- tanh(Tmat[,(l+k+1)])
bhat.sim <- beta.sim[,"max_rdisc"]
ghat.sim <- gamma.sim[,"max_rdisc"]
XB.sim <- m2$X2s %*% t(beta.sim)
dhatXB.sim <- dnorm(XB.sim)
phatXB.sim <- pnorm(XB.sim)
ZG.sim <- m2$X1 %*% t(gamma.sim)
dhatZG.sim <- dnorm(ZG.sim)
phatZG.sim <- pnorm(ZG.sim)
mu1.sim <- (ZG.sim-XB.sim*rho.sim)/sqrt(1-rho.sim^2)
mu2.sim <- (XB.sim-ZG.sim*rho.sim)/sqrt(1-rho.sim^2)

## sadly pbivnorm does not handle matrices 
## so we will just loop over 
phat.sim <- sapply(1:1000,
                   \(i){
                   pbivnorm(XB.sim[,i], ZG.sim[,i], rho.sim[i])
                   })


part1 <- bhat.sim*dhatXB.sim * pnorm(mu1.sim)
part2 <- ghat.sim*dhatZG.sim * pnorm(mu2.sim)
part3 <- ghat.sim*phat.sim*dhatZG.sim
AME.sim <- colMeans((part1+part2)/phatZG - part3/(phatZG^2))
se.AME2 <- sd(AME.sim)
c(AME2, se.AME2)

## p-value
2*pnorm(AME2/se.AME2, lower=FALSE)

## plot the predicted probability of escalation 
## given a crisis
summary(m2$X2s[,"max_rdisc"])
quantile(m2$X2s[,"max_rdisc"])
## variable of interest is between 0 and 1 
## mostly 0
Xstar <- apply(m2$X2s, 2, 
               \(x){
                 if(all(x %in% c(0,1))){
                   return(rep(median(x),50))
                 }else{
                   return(rep(mean(x),50))
                 }
               })
Xstar[,"max_rdisc"] <- seq(0,1, length=50)
Zstar <- apply(m2$X1, 2, 
               \(x){
                 if(all(x %in% c(0,1))){
                   return(rep(median(x),50))
                 }else{
                   return(rep(mean(x),50))
                 }
               })
Zstar[,"max_rdisc"] <- seq(0,1, length=50)

Pr.y.given <- pbivnorm(drop(Xstar %*% beta), 
                       drop(Zstar %*% gamma),
                       m2$theta)/pnorm(drop(Zstar%*%gamma))

## Confidence interval
XB.sim <- Xstar %*% t(beta.sim)
ZG.sim <- Zstar %*% t(gamma.sim)
phat.sim <- sapply(1:1000,
                   \(i){
                     pbivnorm(XB.sim[,i], 
                              ZG.sim[,i], 
                              rho.sim[i])
                   })
Pr.y.given.sim <- phat.sim/pnorm(ZG.sim)
## 90% CI
CI.sim <- rowQuantiles(Pr.y.given.sim, probs=c(0.05,0.95))

plot.df <- data.frame(Pr=Pr.y.given,
                      lo=CI.sim[,1],
                      hi=CI.sim[,2],
                      LDG=seq(0,1, length=50))
ggplot(plot.df)+
  geom_ribbon(aes(x=LDG, ymin=lo, ymax=hi), alpha=.4)+ ## CIs
  geom_line(aes(x=LDG, Pr))+
  theme_bw(14)+
  ylab("Probability of escalation")+
  geom_rug(aes(x=max_rdisc),data=as.data.frame(m2$X2s))

## Big effects but large uncertainty
