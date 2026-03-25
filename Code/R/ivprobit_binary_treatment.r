library(ivreg) #for 2SLS
library(GJRM) #for IV probit
library(sandwich) #standard errors
library(lmtest)
library(matrixStats)
library(MASS)
library(xtable) # tables
library(ggplot2) #plots


rm(list=ls())

matt <- read.csv("datasets/mid_data_for_stata.csv")

## ordinary probit
m1 <- glm(start_mid_US ~ turnover+
            lag4startmid  +lag5startmid + lag6startmid+
            lag4logimp_imputed  +
            totalvac +totalvac2 +totalvac3+
            lag4logexp_imputed+  lag4idealdiff_imputed  + 
            lag4polity_imputed +lag4deltapolity_imputed + 
            lag4cap + lag4ally,
          x=TRUE,y=TRUE,
          data=matt, family=binomial("probit"))
coeftest(m1, vcovCL(m1, matt$cow_code))
ame.probit <- m1$coefficients["turnover"]*mean(dnorm(m1$linear.predictors))
ame.probit

## 2sls
m2 <- ivreg(start_mid_US ~ turnover+
              lag4startmid  +lag5startmid + lag6startmid+
              lag4logimp_imputed  +
              totalvac +totalvac2 +totalvac3+
              lag4logexp_imputed+  lag4idealdiff_imputed  + 
              lag4polity_imputed +lag4deltapolity_imputed + 
              lag4cap + lag4ally|
              lag3enter+ lag4startmid  +lag5startmid + lag6startmid+
              lag4logimp_imputed  +
              totalvac +totalvac2 +totalvac3+
              lag4logexp_imputed+  lag4idealdiff_imputed  + 
              lag4polity_imputed +lag4deltapolity_imputed + 
              lag4cap + lag4ally ,
            data=matt)
summary(m2, vcov=\(x){vcovCL(x, matt$cow_code)})

## bivariate probit
m3 <- gjrm(list(start_mid_US ~ turnover+
                  lag4startmid  +lag5startmid + lag6startmid+
                  lag4logimp_imputed  +
                  totalvac +totalvac2 +totalvac3+
                  lag4logexp_imputed+  lag4idealdiff_imputed  + 
                  lag4polity_imputed +lag4deltapolity_imputed + 
                  lag4cap + lag4ally,
                turnover~  lag3enter+ lag4startmid  +
                  lag5startmid + lag6startmid+
                  lag4logimp_imputed  +
                  totalvac +totalvac2 +totalvac3+
                  lag4logexp_imputed+  lag4idealdiff_imputed  + 
                  lag4polity_imputed +lag4deltapolity_imputed + 
                  lag4cap + lag4ally) ,
           data=matt,
           model="B",
           margins=c("probit", "probit")
)
summary(m3)

## it does not use the sandwich package, but does not it's own SE adjustments
m3 <- adjCov(m3, matt$cow_code)
summary(m3)

round(cbind(probit=m1$coefficients, 
            tsls=m2$coefficients,
            biprobit=m3$coefficients[1:15]),3)

X1<- X0 <- m3$X1
X1[,2] <- 1
X0[,2] <-0
biprobit.ate <- mean(pnorm(X1 %*% m3$coefficients[1:15])-
                       pnorm(X0 %*% m3$coefficients[1:15]))
biprobit.ate

Db <- colMeans(X1* dnorm(drop(X1 %*% m3$coefficients[1:15])) - 
                 X0* dnorm(drop(X0 %*%m3$coefficients[1:15]))) 
V.ate.cre <-Db %*% m3$Vb[1:15,1:15] %*% Db
c(biprobit.ate, sqrt(V.ate.cre))
