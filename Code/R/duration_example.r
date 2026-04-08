library(readstata13)
library(dplyr)
library(survival) ## functions for duration models
library(MASS)
library(ggplot2)
library(matrixStats)

rm(list=ls())

ms <- read.dta13("datasets/mattesSavun_duration.dta",
                 nonint.factors = TRUE)
ms <- ms %>%
  mutate(y.obs = Surv(peacedur, event=peacefail)) ## flag the censored obs
summary(ms$y.obs)
head(ms$y.obs)


## survreg handles exponential, Weibull, and some other parametric versions
## always in terms of E[y|X] (i.e., always AFT)
exp.reg <- survreg(  y.obs ~ polpssum + milpssum + econpssum + terrpssum + 
                       costinc + guarantee + issue + costs + duration, 
                    data = ms, robust=TRUE,
                    dist = "exponential", x = TRUE)
summary(exp.reg)



wei.reg <- survreg( y.obs ~ polpssum + milpssum + econpssum + terrpssum + 
                      costinc + guarantee + issue + costs + duration, 
                   data = ms, robust=TRUE,
                   dist = "weibull", x = TRUE)


summary(wei.reg)

## scale here is sigma from above, so hat(alpha) is
a.hat <- 1/wei.reg$scale #pretty close to 1, so not that different from exponential
a.hat

## if we wanted to convert to PH
exp.ph <- list(est = -1* exp.reg$coef,
               vcov=vcov(exp.reg))
exp.ph$se <- sqrt(diag(exp.ph$vcov))
with(exp.ph, cbind(est, se, est/se, 2*pnorm(abs(est/se),lower=FALSE))) %>% round(4)



## for Weibull, the transformation function is
## g(beta, eta) = (exp(-eta))*beta, exp(-eta)),
## D_beta g(beta, eta) = [I_k (-1/s)]
##                       [    0     ]
## D_eta g(beta, eta) = (beta/s, -1/s)
## D_theta g(beta, eta) = [D_beta,  ]
Dtheta <- cbind(rbind(diag(rep(-1/wei.reg$scale, length(wei.reg$coef))),0), #Dbeta
                (c(wei.reg$coef,-1)/(wei.reg$scale))) #Deta

wei.ph <- list(est = c((-1/wei.reg$scale)*wei.reg$coef,1/wei.reg$scale),
               vcov= Dtheta %*% vcov(wei.reg) %*% t(Dtheta))
wei.ph$se <- sqrt(diag(wei.ph$vcov))
with(wei.ph, cbind(est, se, est/se, 2*pnorm(abs(est/se),lower=FALSE))) %>% round(4)


#### marginals with respect to coalition size ###

###### Exponential ######
100*(exp(exp.reg$coefficients["costinc"])-1) ## increase in expected duration by 75%
100*(exp(exp.ph$est["costinc"])-1) ## decrease in the hazard by 43% !


exp.reg$coefficients["costinc"]*mean(exp(exp.reg$linear.predictors)) #increase by 464 months!
exp.ph$est["costinc"]*mean(exp(exp.reg$x %*% exp.ph$est)) #decrease prob of returning by about 0.005

###### Weibull #####
100*(exp(wei.reg$coefficients["costinc"])-1) ## increase in expected duration by 78%
100*(exp(wei.ph$est["costinc"])-1) ## decrease in the hazard by 43% !


wei.reg$coefficients["costinc"]*mean(exp(wei.reg$linear.predictors)) *gamma(1+1/a.hat) #increase by 524 months!
wei.ph$est["costinc"]*mean(exp(drop(wei.reg$x %*% wei.ph$est[1:ncol(wei.reg$x)]))*ms$dur^(a.hat-1))* a.hat #decrease prob of returning by about 0.005



## these are bit tough to really consider
table(exp.reg$x[,"costinc"])

## let's look at the effect of going from 0 to 1
X1 <- X0 <- exp.reg$x
X1[,"costinc"] <- 1
X0[,"costinc"] <- 0

## Exponential 
mean(exp(X1%*% exp.reg$coef)) -mean(exp(X0%*% exp.reg$coef)) #duration increases by 226 months
mean(exp(X1%*% exp.ph$est)) -mean(exp(X0%*% exp.ph$est)) #hazard decreases by 0.006

mean(exp(X1%*% wei.reg$coef)) -mean(exp(X0%*% wei.reg$coef)) #duration
mean(exp(X1%*% exp.ph$est)) -mean(exp(X0%*% exp.ph$est)) #hazard


## maybe we can plot the survival function or expected duration to make the point clearer
X<- apply(exp.reg$x, 2, 
          \(x){
            if(all(x %in% c(0,1))){
              rep(median(x), 5)
            }else{
              rep(mean(x), 5)
            }
          })

X[,"costinc"] <- 0:4
yhat <- exp(drop(X %*% exp.reg$coef))

Bmat <- mvrnorm(1000, exp.reg$coef, Sigma=vcov(exp.reg))
ysim <-exp(X %*% t(Bmat))
CI <- rowQuantiles(ysim, probs= c(0.025,0.975))

plot.df <- data.frame(Duration=yhat,
                      Provisions=0:4,
                      lo=CI[,1],
                      hi=CI[,2])
ggplot(plot.df)+
  geom_ribbon(aes(x=Provisions, ymin=lo, ymax=hi), alpha=.3)+
  geom_line(aes(y=Duration,
                x=Provisions))+
  ylab("Expected duration (months)")+
  xlab("Cost increasing provisions")+
  theme_bw(14)

## weibull
yhat.wei <- exp(drop(X %*% wei.reg$coef))
Bmat <- mvrnorm(1000, wei.reg$coef, Sigma=vcov(wei.reg)[1:10, 1:10])
ysim <-exp(X %*% t(Bmat))
CI <- rowQuantiles(ysim, probs= c(0.025,0.975))
plot.df <- data.frame(Duration=yhat.wei,
                      Provisions=0:4,
                      lo=CI[,1],
                      hi=CI[,2])
ggplot(plot.df)+
  geom_ribbon(aes(x=Provisions, ymin=lo, ymax=hi), alpha=.3)+
  geom_line(aes(y=Duration,
                x=Provisions))+
  ylab("Expected duration (months)")+
  xlab("Cost increasing provisions")+
  theme_bw(14)

## Survival and hazard for different provision numbers, over time
X<- apply(exp.reg$x, 2, 
          \(x){
            if(all(x %in% c(0,1))){
              rep(median(x),3)
            }else{
              rep(mean(x),3)
            }
          })

X[,"costinc"] <- 0:2
haz <- exp(X %*% exp.ph$est)
haz
haz <- rep(haz, each=100)
time <- rep(seq(1, 480, length=100),3)
S <- exp(-haz*time)

plot.df <- data.frame(Time=time,
                      Value=c(S, haz),
                      stat=rep(c("Survival", "Hazard"), each=300),
                      Hazard=haz,
                      Provisions=factor(rep(0:2,each=100)))

ggplot(plot.df)+
  geom_line(aes(x=Time, y=Value, color=Provisions), linewidth=1)+
  theme_bw(14)+
  facet_wrap(~stat,scales="free_y")+
  ylab("")




# Weibull
lam <- exp(X %*% wei.ph$est[-length(wei.ph$est)])
lam <- rep(lam, each=100)

haz <- (a.hat) * (time^(a.hat-1))*lam
S <- exp(-time^(1/wei.reg$scale) *lam)

plot.df <- data.frame(Time=time,
                      Value=c(S, haz),
                      stat=rep(c("Survival", "Hazard"), each=300),
                      Hazard=haz,
                      Provisions=factor(rep(0:2,each=100)))

ggplot(plot.df)+
  geom_line(aes(x=Time, y=Value, color=Provisions), linewidth=1)+
  theme_bw(14)+
  facet_wrap(~stat,scales="free_y")+
  ylab("")












##### Non parametric #### 
library(biostat3) ## for smoothing the hazard estimates
library(survminer) ## additional helpful functions

cox.reg <- coxph(y.obs ~ polpssum + milpssum + econpssum +
                   terrpssum + costinc + guarantee + issue +
                   costs + duration, 
                 data = ms, 
                 method="efron",
                 robust=TRUE, x = TRUE)
summary(cox.reg)

cox.reg2 <- coxph(y.obs ~ polpssum + milpssum + econpssum +
                    terrpssum + costinc + guarantee + issue +
                    costs + duration, 
                  data = ms, 
                  method="breslow",
                  robust=TRUE, x = TRUE)
cox.reg3 <- coxph(y.obs ~ polpssum + milpssum + econpssum +
                    terrpssum + costinc + guarantee + issue +
                    costs + duration, 
                  data = ms, 
                  method="exact", x = TRUE)

compareCoefs(cox.reg, cox.reg2, cox.reg3)

## a test for PH
cox.zph(cox.reg)
## costs looks concerning, as does that 
## last overall one, but most looks OK

## we can plot it too, these should all be flat
ggcoxdiagnostics(cox.reg,"schoenfeld")

## Fit the baseline  with KM/NA
base <- survfit(cox.reg)
plot(base, ylab="Survival", xlab="Time")

## estimate Survival and hazard 
## for different values of 
## of costinc
X<- apply(cox.reg$x, 2, 
          \(x){
            if(all(x %in% c(0,1))){
              rep(median(x),3)
            }else{
              rep(mean(x),3)
            }
          })

X[,"costinc"] <- 0:2
X.df <- as.data.frame(X)
X.df
h0 <- survfit(cox.reg,newdata=X.df)
plot(h0, col=c("red","blue","darkgreen"),
     ylab="Survival",
     xlab="time")
legend("topright",
       c("0", "1", "2"), 
       col=c("red", "blue", "darkgreen"),
       lty=1)



## check for proportional hazards by 
## looking at -log(-log(S(t))) for two values
## of costinc
plot(h0, col=c("red","blue","darkgreen"),
     ylab="Survival",
     xlab="time", fun="cloglog")
legend("bottomright",
       c("0", "1", "2"), 
       col=c("red", "blue", "darkgreen"),
       lty=1)

## biostat3 gives us smoothing for the hazard
h.hat <- coxphHaz(cox.reg, newdata=X.df)
plot(h.hat, 
     legend.args=list(legend=c("0", "1", "2")))
