library(readstata13)
library(dplyr)
library(MASS)
library(lmtest)
library(car)
library(ggplot2)
library(sandwich)
library(survival) ## for conditional estimator
library(lme4) ## for (C)RE estimator
library(fastGHQuad)
library(xtable)
rm(list=ls())

source("clogit_functions.R")
emw <- read.dta13("temp-micro.dta",
                  nonint.factors = TRUE) # avoid some warning

var.names <- c("remit", "remitXpro", "cellphone",
               "lage", "education", "wealth", "male",
               "employment", "travel")

emw <- emw %>%
  subset(select=c(var.names, "protest01", "dcode")) %>%
  na.omit() %>%
  mutate(ybar = mean(protest01, na.rm=TRUE), 
         .by=dcode) %>% ## match their coding
  mutate(sub.sam = 1-1*(ybar %in% c(0,1)))

## This is the sub set that is actually informative with 
## fixed effects (either MLDV or CML)
cml.data <- emw %>% filter(sub.sam==1)


##### pooled version ####
pooled <- glm(protest01~remit + remitXpro+
               cellphone+ lage+ education+ wealth+ male+ 
               employment+ travel,
              x=TRUE,
             data=emw, family=binomial())
V0 <- vcovCL(pooled, emw$dcode)
coeftest(pooled, V0)

## marginals
X5.lo <- X0.lo <- X5.hi <- X0.hi <- pooled$x
X5.lo[,"remit"] <- X5.hi[,"remit"] <- 5
X0.lo[,"remit"] <- X0.hi[,"remit"] <- 0
X5.lo[,"remitXpro"] <- 5*.18
X5.hi[,"remitXpro"] <- 5*.8
X0.lo[,"remitXpro"] <- X0.hi[,"remitXpro"] <- 0
FD0 <- c(mean(plogis(X5.lo %*% pooled$coef)-plogis(X0.lo %*% pooled$coef)),
         mean(plogis(X5.hi %*% pooled$coef)-plogis(X0.hi %*% pooled$coef)))



## SE
Jbeta <- rbind(colMeans(X5.lo*dlogis(drop(X5.lo %*% pooled$coef))- 
                          X0.lo*dlogis(drop(X0.lo %*% pooled$coef))),
               colMeans(X5.hi*dlogis(drop(X5.lo %*% pooled$coef)) - 
                          X0.lo*dlogis(drop(X0.lo %*% pooled$coef))))


Vame <-Jbeta %*% V0 %*% t(Jbeta)
SeAme0 <- sqrt(diag(Vame))
rbind(FD0,SeAme0) %>% round(.,3)

#### Dummy variables ####
mldv1 <- glm(protest01~remit + remitXpro+
               cellphone+ lage+ education+ wealth+ male+ 
               employment+ travel+factor(dcode)-1,
             x=TRUE,
             data=cml.data, family=binomial())
V1 <- vcovCL(mldv1, cml.data$dcode)
coeftest(mldv1, V1)[grep("factor", 
                         names(mldv1$coef),
                         invert = TRUE),] %>% 
  round(.,3)

## marginals
X5.lo <- X0.lo <- X5.hi <- X0.hi <- mldv1$x
X5.lo[,"remit"] <- X5.hi[,"remit"] <- 5
X0.lo[,"remit"] <- X0.hi[,"remit"] <- 0
X5.lo[,"remitXpro"] <- 5*.18
X5.hi[,"remitXpro"] <- 5*.8
X0.lo[,"remitXpro"] <- X0.hi[,"remitXpro"] <- 0

FD1 <- c(mean(plogis(X5.lo %*% mldv1$coef)-plogis(X0.lo %*% mldv1$coef)),
         mean(plogis(X5.hi %*% mldv1$coef)-plogis(X0.hi %*% mldv1$coef)))

  

## SE
Jbeta <- rbind(colMeans(X5.lo*dlogis(drop(X5.lo %*% mldv1$coef))- 
                    X0.lo*dlogis(drop(X0.lo %*% mldv1$coef))),
               colMeans(X5.hi*dlogis(drop(X5.lo %*% mldv1$coef)) - 
                       X0.lo*dlogis(drop(X0.lo %*% mldv1$coef))))

  
Vame <-Jbeta %*% V1 %*% t(Jbeta)
SeAme1 <- sqrt(diag(Vame))
rbind(FD1,SeAme1) %>% round(.,3)

## does it change with full sample?
mldv2 <- glm(protest01~remit + remitXpro+
               cellphone+ lage+ education+ wealth+ male+ 
               employment+ travel+factor(dcode)-1,
             x=TRUE,
             data=emw, family=binomial())
## certainly slower
V1a <- vcovCL(mldv2, emw$dcode)
coeftest(mldv2, V1a)[grep("factor", 
                         names(mldv2$coef),
                         invert = TRUE),] %>% 
  round(.,3)

## marginals
X5.lo <- X0.lo <- X5.hi <- X0.hi <- mldv2$x
X5.lo[,"remit"] <- X5.hi[,"remit"] <- 5
X0.lo[,"remit"] <- X0.hi[,"remit"] <- 0
X5.lo[,"remitXpro"] <- 5*.18
X5.hi[,"remitXpro"] <- 5*.8
X0.lo[,"remitXpro"] <- X0.hi[,"remitXpro"] <- 0

FD1a <- c(mean(plogis(X5.lo %*% mldv2$coef)-plogis(X0.lo %*% mldv2$coef)),
         mean(plogis(X5.hi %*% mldv2$coef)-plogis(X0.hi %*% mldv2$coef)))

## estimates don't change but the marginals do... That may concern us
## Does it make sense to average in all those unidentified alphas?
## Does it make sense not to?

## SE
Jbeta <- rbind(colMeans(X5.lo*dlogis(drop(X5.lo %*% mldv2$coef))- 
                          X0.lo*dlogis(drop(X0.lo %*% mldv2$coef))),
               colMeans(X5.hi*dlogis(drop(X5.lo %*% mldv2$coef)) - 
                          X0.lo*dlogis(drop(X0.lo %*% mldv2$coef))))


Vame <-Jbeta %*% V1a %*% t(Jbeta)
SeAme1a <- sqrt(diag(Vame))
rbind(FD1a,SeAme1a) %>% round(.,3)


#### Chamberlain's conditional ML (survival)####

cml <- clogit(protest01~remit + remitXpro+
               cellphone+ lage+ education+ wealth+ male+ 
               employment+ travel+strata(dcode),
             data=cml.data,x=TRUE,y=TRUE)
summary(cml)

V2 <- vcovCL.clogit(cml)
coeftest(cml, vcov=V2) %>% round(.,3)

## marginals.... oh wait
## Fun fact, the authors present marginal effects, but 
## calculate them assuming that each $\alpha_i=0$, which is...a choice



##### CRE#####
cre.dat <- emw %>% 
  group_by(dcode) %>%
  mutate(across(all_of(var.names), 
                .fns=mean,
                .names="{.col}.bar")) %>% 
  ungroup()

## canned version (can be flakey)
cre <- glmer(protest01~remit + remitXpro+
        cellphone+ lage+ education+ wealth+ male+ 
        employment+ travel+
        remit.bar + remitXpro.bar+
        cellphone.bar+ lage.bar+
        education.bar+ wealth.bar+ male.bar+
        employment.bar+ travel.bar+
        (1|dcode),
      data=cre.dat,
      family=binomial(),
      nAGQ = 7, ## number of GH nodes
      control=glmerControl(optCtrl=list(maxfun=5e5)))
summary(cre)
  
  
## Compute clustered SE's for it
X <- model.matrix(~remit + remitXpro+
                    cellphone+ lage+ education+ 
                    wealth+ male+ 
                    employment+ travel+
                    remit.bar + remitXpro.bar+
                    cellphone.bar+ lage.bar+
                    education.bar+ wealth.bar+ male.bar+
                    employment.bar+ travel.bar, 
                  data=cre.dat)
id <- cre.dat$dcode
y <- cre.dat$protest01
theta <- c(cre@beta, log(cre@theta))
names(theta) <- c(colnames(X), "ln.sigma")
R <- 7
GH <- gaussHermiteData(R)

J <-   r.relogit(theta,
                 X=X, y=y, 
                 id=id, GH=GH, sum=FALSE)
OPG <- crossprod(J) ## note that the jacobian already sums up within units
bread <- vcov(cre)
cheese <- OPG[-length(cre@beta), -length(cre@beta)]
V3 <- bread %*% cheese %*% bread
tab.cre <- data.frame(Est=cre@beta, SE=sqrt(diag(V3))) %>% 
  mutate(tstat=Est/SE, 
         p=2*pnorm(abs(tstat),lower=FALSE)) 
tab.cre <- tab.cre[grep("bar|Intercept", 
                        row.names(tab.cre),
                        invert=TRUE),]
round(tab.cre,3)
## specification test
psi.names <- colnames(X)[grep(".bar", colnames(X))]
linearHypothesis(cre, psi.names, vcov=V3)

## marginals
X5.lo <- X0.lo <- X5.hi <- X0.hi <- X
X5.lo[,"remit"] <- X5.hi[,"remit"] <- 5
X0.lo[,"remit"] <- X0.hi[,"remit"] <- 0
X5.lo[,"remitXpro"] <- 5*.18
X5.hi[,"remitXpro"] <- 5*.8
X0.lo[,"remitXpro"] <- X0.hi[,"remitXpro"] <- 0
FD3 <- c(mean(plogis(X5.lo %*% cre@beta)-plogis(X0.lo %*% cre@beta)),
          mean(plogis(X5.hi %*% cre@beta)-plogis(X0.hi %*% cre@beta)))

## SE
Jbeta <- rbind(colMeans(X5.lo*dlogis(drop(X5.lo %*% cre@beta))- 
                          X0.lo*dlogis(drop(X0.lo %*% cre@beta))),
               colMeans(X5.hi*dlogis(drop(X5.lo %*% cre@beta)) - 
                          X0.lo*dlogis(drop(X0.lo %*% cre@beta))))


Vame <-Jbeta %*% V3 %*% t(Jbeta)
SeAme3 <- sqrt(diag(Vame))
rbind(FD3,SeAme3) %>% round(.,3)

## optimize it ourselves
cre.a <- optim(theta, l.relogit, gr=r.relogit,
      X=X, y=y, id=id, GH=GH,method="BFGS",
      hessian=TRUE)
J <-   r.relogit(cre.a$par,
                 X=X, y=y, 
                 id=id, GH=GH, sum=FALSE)
OPG <- crossprod(J)
V3a <- solve(cre.a$hessian) %*% OPG %*% solve(cre.a$hessian) 
tab.cre2 <- data.frame(Est=cre.a$par, 
                       SE=sqrt(diag(V3a))) %>% 
  mutate(tstat=Est/SE, 
         p=2*pnorm(abs(tstat),lower=FALSE)) 
tab.cre2 <- tab.cre2[grep("bar|Intercept|sigma", 
                        row.names(tab.cre2),
                        invert=TRUE),]
round(tab.cre2,3)


## marginals
FD3a <- c(mean(plogis(X5.lo %*% cre.a$par[1:ncol(X)])-plogis(X0.lo %*% cre.a$par[1:ncol(X)])),
         mean(plogis(X5.hi %*% cre.a$par[1:ncol(X)])-plogis(X0.hi %*% cre.a$par[1:ncol(X)])))

## SE
Jbeta <- rbind(colMeans(X5.lo*dlogis(drop(X5.lo %*% cre.a$par[1:ncol(X)]))- 
                          X0.lo*dlogis(drop(X0.lo %*% cre.a$par[1:ncol(X)]))),
               colMeans(X5.hi*dlogis(drop(X5.lo %*% cre.a$par[1:ncol(X)])) - 
                          X0.lo*dlogis(drop(X0.lo %*% cre.a$par[1:ncol(X)]))))


Vame <-Jbeta %*% V3a[1:ncol(X), 1:ncol(X)] %*% t(Jbeta)
SeAme3a <- sqrt(diag(Vame))
rbind(FD3a,SeAme3a) %>% round(.,3)

## Table it
m0 <- list(par=pooled$coef, se=sqrt(diag(V0)), ame=FD0, se.ame=SeAme0, Obs=nrow(pooled$x))
m1 <- list(par=mldv1$coef, se=sqrt(diag(V1)), ame=FD1, se.ame=SeAme1, Obs=nrow(mldv1$x))
m1a <- list(par=mldv2$coef, se=sqrt(diag(V1a)), ame=FD1a, se.ame=SeAme1a, Obs=nrow(mldv2$x))
m2 <- list(par=cml$coef, se=sqrt(diag(V2)),  Obs=nrow(cml$x))
m3 <- list(par=cre.a$par, se=sqrt(diag(V3a)), ame=FD3a, se.ame=SeAme3a, Obs=nrow(X))

model.list <- list(m0, m1, m1a,m2, m3)
ests <- lapply(model.list,
               \(x){rbind(x$par[grep("bar|Intercept|sigma|factor", 
                                           names(x$par),
                                           invert=TRUE)],
                          x$se[grep("bar|Intercept|sigma|factor", 
                                    names(x$se),
                                    invert=TRUE)])})
                          
ests <- lapply(ests, round, 2)
ests <- lapply(ests, \(x){x[2,] <- paste0("(", x[2,] ,")");return(x)})
ests <- sapply(ests, c)
estnames <- rbind(c("Remit", 
                    "Remit $\\times$ ProGov",
                    "Cellphone",
                    "Age (log)",
                    "Education",
                    "Wealth",
                    "Male",
                    "Employment",
                    "Travel"),
                  "")
ests <- cbind(c(estnames), ests)
colnames(ests) <- c(" ", "Pooled", "MLDV", "MLDV-Full sample", "CML", "CRE")
ests <- rbind(ests, c("Obs", sapply(model.list, \(x){x$Obs})))
print(ests)
print(xtable(ests, align="llccccc", caption="Replication EMW (2018)"),
      include.rownames=FALSE,
      sanitize.text.function = \(x){x},
      caption.placement="top")


## Marginals
FDout <- lapply(model.list[-4],
               \(x){rbind(x$ame,
                          x$se.ame)})

FDout <- lapply(FDout, round, 3)
FDout <- lapply(FDout, \(x){x[2,] <- paste0("(", x[2,] ,")");return(x)})
FDout <- sapply(FDout, c)
FDout <- cbind(c("Low government support", "", "High government support", ""),
               FDout)
colnames(FDout) <- c(" ","Pooled", "MLDV", "MLDV-Full sample", "CRE")

print(FDout)
print(xtable(FDout, align="llcccc", caption="Effect of remittances by government support"),
      include.rownames=FALSE,
      sanitize.text.function = \(x){x},
      caption.placement="top")

