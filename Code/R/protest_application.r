## basics
library(readstata13)
library(dplyr)
## econometric
library(MASS)
library(sandwich)
library(lmtest)
library(car)
library(nonnest2)
## tables and figures
library(xtable)
library(ggplot2)

rm(list=ls())


protests <- read.dta13("datasets/remittances_and_protests1.dta")
protests <- protests %>% 
  filter(sample==1) %>%  ## authors' chosen sample
  mutate(sum.y = sum(banks_protest_count), .by=cowcode) %>% ## the infinite cases
  filter(sum.y >0) 

## we'll keep everything comparable across models by limiting ourselves to the 
## countries with protests



### Looking at our main variables####
ggplot(protests) +
  geom_histogram(aes(x=banks_protest_count, y=after_stat(density)),
                 fill="blue", color="black", 
                 bins = nclass.scott(protests$banks_protest_count))+
  xlab("Protests")+
  ylab("Density")+
  theme_bw(12)

summary(protests$remit)
dim(protests)

length(unique(protests$cowcode))
summary(protests$year)

#### Fit the models ####

## model formula we'll use today
f1 <- banks_protest_count ~ remit*dict + 
  l1gdp + #GDP pc capita
  l1pop+ #population
  l1nbr5 + #neighboring protests
  l12gr+ #Economic growth
  l1migr+ #Net migration
  elec3+ #election
  factor(period)+
  factor(cowcode)-1

## Normal model
normal <- lm(f1, data=protests, x=TRUE, y=TRUE)
length(normal$resid) ## nothing lost to missingness 

## log normals
log.normal1 <- lm(update(f1, log(banks_protest_count+1)~.), data=protests,
                  y=TRUE, x=TRUE)
log.normal2 <- lm(update(f1, log(banks_protest_count+0.01)~.), data=protests,
                  y=TRUE, x=TRUE)

## Poisson and NB (glm and MASS::glm.nb)
poi <- glm(f1, data=protests, family=poisson, y=TRUE, x=TRUE)
nbreg <- glm.nb(f1, data=protests,y=TRUE, x=TRUE)


#### housekeeping before comparing estimates ####

##collect the models
mod.list <- list(normal,
                 log.normal1,
                 log.normal2,
                 poi,
                 nbreg)

## remove all the country and year coefficients to make this easier to look at
factors <- grep(names(normal$coefficients),pattern="factor")
Ests <- sapply(mod.list,\(x){x$coefficients[-factors]})
SE <- sapply(mod.list,\(x){sqrt(diag(vcov(x)))[-factors]})



## ordinary Hessian based standard errors
tab1 <- rbind(c(Ests),c(SE)) %>% 
  matrix(., ncol=5) %>% 
  round(., 3)
tab1[seq(2, 18,2),] <- paste0("(", tab1[seq(2, 18,2),], ")")

tab1 <- cbind(c(rbind(c("Remit",
                        "Autocracy", 
                        "GDP pc",
                        "Pop",
                        "Neighboring protest",
                        "Growth",
                        "Net Migration",
                        "Election", 
                        "Remit $\\times$ Autocracy"), "")),
              
              tab1)
colnames(tab1) <- c(" ",
                    "Normal",
                    "Log-normal (+1)", 
                    "Log-normal (+0.01)",
                    "Poisson",
                    "Negative binomial"
)

## model information rows
## NOTE: R calls the r parameter from the negative binomial theta, we'll follow suit,
## but this is the r we discussed above

tab1 <- rbind(tab1,
              c("Observations", round(sapply(mod.list, \(x){length(x$fitted.values)}), 0)),
              c("Log-Likelihood", round(sapply(mod.list, logLik), 2)),
              c("$\\theta$", rep(" ", 4), round(nbreg$theta,2)))
                    

tab1

print(xtable(tab1, align="llccccc"), 
      include.rownames=FALSE,
      booktabs=TRUE,
      sanitize.text.function = \(x){x},
      caption.placement="top")










#### clustered variance matrices###
Vn <- vcovCL(normal, cluster=protests$cowcode)
Vln <- vcovCL(log.normal1, cluster=protests$cowcode)
Vln2 <- vcovCL(log.normal2, cluster=protests$cowcode)
Vp <- vcovCL(poi, cluster=protests$cowcode)
Vnb <- vcovCL(nbreg, cluster=protests$cowcode)

#### clustered variance matrices###
Vn <- vcovCL(normal, cluster=protests$cowcode)
Vln <- vcovCL(log.normal1, cluster=protests$cowcode)
Vln2 <- vcovCL(log.normal2, cluster=protests$cowcode)
Vp <- vcovCL(poi, cluster=protests$cowcode)
Vnb <- vcovCL(nbreg, cluster=protests$cowcode)
SE.list <- list(Vn, Vln,Vln2, Vp, Vnb)
SE <- sapply(SE.list,\(x){sqrt(diag(x))[-factors]})

## repeat the tabling process
tab1 <- rbind(c(Ests),c(SE)) %>% 
  matrix(., ncol=5) %>% 
  round(., 3)
tab1[seq(2, 18,2),] <- paste0("(", tab1[seq(2, 18,2),], ")")

tab1 <- cbind(c(rbind(c("Remit",
                        "Autocracy", 
                        "GDP pc",
                        "Pop",
                        "Neighbor protest",
                        "Growth",
                        "Net Migration",
                        "Election", 
                        "Remit $\\times$ Autocracy"), "")),
              
              tab1)
colnames(tab1) <- c(" ",
                    "Normal",
                    "Log-normal (+1)", 
                    "Log-normal (+0.01)",
                    "Poisson",
                    "Negative binomial"
)

tab1 <- rbind(tab1,
              c("Observations", round(sapply(mod.list, \(x){length(x$fitted.values)}), 0)),
              c("Log-Likelihood", round(sapply(mod.list, logLik), 2)),
              c("$\\theta$", rep(" ", 4), round(nbreg$theta,2)))

print(xtable(tab1, align="llccccc"), 
      include.rownames=FALSE,
      booktabs=TRUE,
      sanitize.text.function = \(x){x},
      caption.placement="top")


### interaction ###
## deltaMethod from the car package is an easy way to 
## combine covariates from a fitted model in either 
## a linear or nonlinear way. We use the back ticks
## on remit:dict because it has that special symbol :
## the ` ` keeps it together
deltaMethod(normal, "remit+`remit:dict`", vcov=Vn)
deltaMethod(log.normal1, "remit+`remit:dict`", vcov=Vln)
deltaMethod(log.normal2, "remit+`remit:dict`", vcov=Vln2)
deltaMethod(poi, "remit+`remit:dict`", vcov=Vp)
deltaMethod(nbreg, "remit+`remit:dict`", vcov=Vnb)

### Average marginal effects ###

## we can look at these in terms of the logged variable D_r E[y|X]
## or in terms of remittances D_s E[y|X] 
## leaving it in log terms is what we're used to 
normalAME <- c(normal$coefficients['remit'],               
               normal$coefficients['remit'] + normal$coefficients["remit:dict"])

## matricies for counter factuals
X0 <- X1 <- poi$x
X1[,"dict"] <- 1
X1[,"remit:dict"] <- X1[,"remit"]
X0[,"dict"] <- 0
X0[,"remit:dict"] <- 0

log.normal1.s2.hat <- mean(log.normal1$residuals^2)
lnormalAME1 <- c(log.normal1$coefficients['remit']*
                   exp(log.normal1.s2.hat/2) * 
                   mean(exp(X0 %*% log.normal1$coef)),
                 (log.normal1$coefficients['remit'] + 
                    log.normal1$coefficients["remit:dict"])*
                   exp(log.normal1.s2.hat/2) * 
                   mean(exp(X1 %*% log.normal1$coef)))




log.normal2.s2.hat <- mean(log.normal2$residuals^2)
lnormalAME2 <-  c(log.normal2$coefficients['remit']*
                    exp(log.normal2.s2.hat/2) * 
                    mean(exp(X0 %*% log.normal2$coef)),
                  (log.normal2$coefficients['remit'] + 
                     log.normal2$coefficients["remit:dict"])*
                    exp(log.normal2.s2.hat/2) * 
                    mean(exp(X1 %*% log.normal2$coef)))

poissonAME <-  c(poi$coefficients['remit']* mean(exp(X0%*% poi$coefficients)),
                 (poi$coefficients['remit'] + 
                    poi$coefficients["remit:dict"])*
                   mean(exp(X1 %*% poi$coefficients)))

nbAME <-  c(nbreg$coefficients['remit']* mean(exp(X0 %*% nbreg$coefficients)),
            (nbreg$coefficients['remit'] + 
               nbreg$coefficients["remit:dict"])*
              mean(exp(X1 %*% nbreg$coefficients)))

ame.tab <- round(cbind(normalAME, lnormalAME1,lnormalAME2, poissonAME, nbAME), 3)
ame.tab <- cbind( c("Demo", "Auto"), ame.tab)
colnames(ame.tab) <- colnames(tab1)

print(xtable(ame.tab, align="llccccc", caption="AME of logged remittances"), 
      include.rownames=FALSE,
      booktabs=TRUE,
      sanitize.text.function = \(x){x},
      caption.placement="top")


## we can look at these in terms of the logged variable D_r E[y|X]
## or in terms of remittances D_s E[y|X] 
## The latter
s <- mean(exp(protests$remit)/1000000) ## remittances in millions

normalAME <- c(normal$coefficients['remit'],               
               normal$coefficients['remit'] + normal$coefficients["remit:dict"])/s
## matricies for counter factuals
X0 <- X1 <- poi$x
X1[,"dict"] <- 1
X1[,"remit:dict"] <- X1[,"remit"]
X0[,"dict"] <- 0
X0[,"remit:dict"] <- 0

log.normal1.s2.hat <- mean(log.normal1$residuals^2)
lnormalAME1 <- colMeans(cbind(log.normal1$coefficients['remit']*
                                exp(log.normal1.s2.hat/2) * 
                                exp(X0 %*% log.normal1$coef),
                              (log.normal1$coefficients['remit'] + 
                                 log.normal1$coefficients["remit:dict"])*
                                exp(log.normal1.s2.hat/2) * 
                                exp(X1 %*% log.normal1$coef))/s)




log.normal2.s2.hat <- mean(log.normal2$residuals^2)
lnormalAME2 <- colMeans(cbind(log.normal2$coefficients['remit']*
                    exp(log.normal2.s2.hat/2) * 
                    exp(X0 %*% log.normal2$coef),
                  (log.normal2$coefficients['remit'] + 
                     log.normal2$coefficients["remit:dict"])*
                    exp(log.normal2.s2.hat/2) * 
                    exp(X1 %*% log.normal2$coef))/s)

poissonAME <-  colMeans(cbind(poi$coefficients['remit']* exp(X0%*% poi$coefficients),
                 (poi$coefficients['remit'] + poi$coefficients["remit:dict"])*
                   exp(X1 %*% poi$coefficients))/s)

nbAME <-  colMeans(cbind(nbreg$coefficients['remit']* exp(X0%*% nbreg$coefficients),
                         (nbreg$coefficients['remit'] + nbreg$coefficients["remit:dict"])*
                           exp(X1 %*% nbreg$coefficients))/s)

ame.tab <- round(cbind(normalAME, lnormalAME1,lnormalAME2, poissonAME, nbAME), 3)
ame.tab <- cbind( c("Demo", "Auto"), ame.tab)
colnames(ame.tab) <- colnames(tab1)

print(xtable(ame.tab, align="llccccc", caption="AME of a 1 million USD increase in remittances"), 
      include.rownames=FALSE,
      booktabs=TRUE,
      sanitize.text.function = \(x){x},
      caption.placement="top")



#### Expected values ####
## What to vary remittances over?
hist(poi$x[,"remit"],freq = FALSE,breaks = "scott")
quantile(protests$remit, probs=c(0,0.025,.05, .95,.975,1))

## Create an "average" profile 
## whatever that actually means
## I've never been to averagastan, have you?
Xstar <- t(replicate(50, colMeans(poi$x))) ##means for most things
Xstar[,"elec3"] <- median(poi$x) ## one binary covariate we'll set to the median

Xstar[,"dict"] <- c(rep(0,25), rep(1,25)) ## setting regime type to change
Xstar[,"remit"] <- seq(8, 18, length=25) ## setting remittances to vary
Xstar[,"remit:dict"] <- Xstar[,"dict"] *Xstar[,"remit"]

ystar <- drop(exp(Xstar %*% poi$coefficients))

## delta method
DyDb <- Xstar* ystar
V.ystar <- DyDb %*% Vp %*% t(DyDb)
se.ystar <- sqrt(diag(V.ystar))


## create a data frame for plotting
plot.df <- data.frame(Remittances=Xstar[,"remit"],
                      Regime=rep(c("Democracy", "Autocracy"), each=25),
                      Protests=ystar,
                      lo=ystar-1.96*se.ystar,
                      hi=ystar+1.96*se.ystar)

## plot
ggplot(plot.df)+
  geom_ribbon(aes(x=Remittances, ymin=lo, ymax=hi), alpha=.4)+ ## CIs
  geom_line(aes(x=Remittances, Protests))+
  theme_bw(18)+
  facet_grid(~Regime)+
  xlab("logged Remittances")+
  geom_rug(aes(x=remit),data=protests%>%
             filter(remit>=8 & remit <= 18))




#### Model comparisons ####
lrtest(nbreg, poi) ## from lmtest package

##let's double check this one
LRstat <- 2*(logLik(nbreg) -logLik(poi))
LRpval <- pchisq(LRstat, df=1, lower.tail = FALSE)
c(LRstat, LRpval)

## Here we do find some evidence that negative binomial may be preferred
## However, remember that the Poisson is robust to more types of misspecification
## so it comes down to what our goal(s) are.


pairs(data.frame(Protests= protests$banks_protest_count,
                 Normal=normal$fitted.values,
                 Log.normal1=  exp(log.normal1$fitted.values+log.normal1.s2.hat/2)-1,
                 Log.normal2=  exp(log.normal2$fitted.values+log.normal2.s2.hat/2),
                 Poisson=  poi$fitted.values,
                 NegBin=  nbreg$fitted.values),
      upper.panel = \(x,y){
        par(usr = c(0, 1, 0, 1))
        text(.5,.5,round(cor(x,y),2), 
             cex=1.5)
        },
      cex.labels =1.285)



#### comparing normal to Poisson ####
## start with AIC/BIC
c(AIC(normal), AIC(poi))
c(BIC(normal), BIC(poi))

## which is better on these metrics?



###### Vuong test #####
## we need the likelihood ratio at each point
normal.ll <- function(theta, y,X){
  s2 <- theta[length(theta)]
  beta <- theta[-length(theta)]
  mu <- X %*%  beta
  return(-log(2*pi*s2)/2-(y-mu)^2 / (2*s2))
}
poisson.ll <- function(beta, y,X){
  mu <- X %*% beta
  lambda <- exp(X %*% beta)
  return(y*(mu)-lambda -lgamma(y+1))
}
s2.hat <-mean(normal$residuals^2) ##MLE for s2
LR1 <- normal.ll(c(normal$coef,s2.hat), y=normal$y, X=normal$x)
LR2 <- poisson.ll(poi$coef, y=poi$y,X=poi$x)
N <- length(normal$residuals)
LR <- LR1-LR2
omega2 <- mean(LR^2) - mean(LR)^2
LR <- sum(LR)
V <- LR/(sqrt(N*omega2))
V

pnorm(V) ## reject the null of equal fit in favor of model 2 (Poisson)


vuongtest(normal, poi) ## from nonnest2 package (no penalty)

## Here we conclude that Poisson is a better model for y than the normal
