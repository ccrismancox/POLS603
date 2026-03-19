library(readstata13) #read the data
library(ivreg) #for 2SLS
library(GJRM) #for IV probit
library(sandwich) #standard errors
library(lmtest)
library(matrixStats)
library(MASS)
library(xtable) # tables
library(ggplot2) #plots


rm(list=ls())

fdi <- read.dta13("datasets/FDI_Conflict.dta")
fdi <- subset(fdi, sample==1)

## Start with an ordinary probit to set a baseline
m1 <- glm(onset2v414 ~ cr_rfdicap+
            lmrpe+ llrgdpcap +llpop+ llgrowthrate+ 
            lpolity2 +lloil_gas_cap+ elf85 +
            relfrac +lmtnest+ ncontig+ coldwar+ 
            t1_onset2+ t2_onset2 +t3_onset2,
          family=binomial("probit"),
          data=fdi)
V1 <- vcovCL(m1, fdi$ccode)
coeftest(m1, vcov=V1)

## no real effect of note, very close to 0


## 2SLS 
m2 <- ivreg(onset2v414 ~ cr_rfdicap+
              lmrpe+ llrgdpcap +llpop+ llgrowthrate+ 
              lpolity2 +lloil_gas_cap+ elf85 +
              relfrac +lmtnest+ ncontig+ coldwar+ 
              t1_onset2+ t2_onset2 +t3_onset2|
              lwdist+# instrument
              lmrpe+ llrgdpcap +llpop+ llgrowthrate+ 
              lpolity2 +lloil_gas_cap+ elf85 +
              relfrac +lmtnest+ ncontig+ coldwar+ 
              t1_onset2+ t2_onset2 +t3_onset2,
            data=fdi)
## ivreg includes the option to change the vcov for summary
## we like this because it also updates the first-stage covariance
## matrix accordingly. It needs a function as its input to apply to both 
## stages
summary(m2, vcov=\(x){vcovCL(x, fdi$ccode)})
## in the expected direction, but not significant
## What do we think about the instrument strength?


### control function probit
step1 <- lm(cr_rfdicap ~ lwdist+# instrument
              lmrpe+ llrgdpcap +llpop+ llgrowthrate+ 
              lpolity2 +lloil_gas_cap+ elf85 +
              relfrac +lmtnest+ ncontig+ coldwar+ 
              t1_onset2+ t2_onset2 +t3_onset2,
            x=TRUE, y=TRUE,
            data=fdi)
fdi$u <- step1$residuals

## time out, we do rely a bit on the residuals being normal
## can we test that at all? Yes!
plot(step1, which=2)
ks.test(step1$residuals, pnorm, 0, summary(step1)$sigma)
## maybe not all that normal, but we persist

step2 <- glm(onset2v414 ~ cr_rfdicap+
               lmrpe+ llrgdpcap +llpop+ llgrowthrate+ 
               lpolity2 +lloil_gas_cap+ elf85 +
               relfrac +lmtnest+ ncontig+ coldwar+ 
               t1_onset2+ t2_onset2 +t3_onset2+u,
             family=binomial("probit"),x=TRUE,y=TRUE,
             data=fdi)
Vstep2.cl <- vcovCL(step2, fdi$ccode)
coeftest(step2, vcov=Vstep2.cl)
## things to remember here
## 1. these standard errors are only right if the estimate on u is 0
## it is not
## These estimates are attenuated

## skip the two step correction, but let's do the attenuation correction
correction <- 1/sqrt(step2$coefficients["u"]^2 * summary(step1)$sigma^2 + 1)
theta.u <- step2$coefficients[1:16]*correction
ame <- theta.u["cr_rfdicap"] * mean(dnorm(step2$x[,1:16] %*% theta.u))
ame 


## full information
## gjrm takes:
##   1. a list of 2 formulas
##   2. data
##   3. model="B" (bivariate model)
##   4. margins= the marginal dist of the data 
##      "N" for the normal equation and "probit" for the binary
m3 <- gjrm(list(onset2v414 ~ cr_rfdicap+
                  lmrpe+ llrgdpcap +llpop+ llgrowthrate+ 
                  lpolity2 +lloil_gas_cap+ elf85 +
                  relfrac +lmtnest+ ncontig+ coldwar+ 
                  t1_onset2+ t2_onset2 +t3_onset2,
                cr_rfdicap ~ lwdist+# instrument
                  lmrpe+ llrgdpcap +llpop+ llgrowthrate+ 
                  lpolity2 +lloil_gas_cap+ elf85 +
                  relfrac +lmtnest+ ncontig+ coldwar+ 
                  t1_onset2+ t2_onset2 +t3_onset2),
           data=fdi,
           model="B",
           margins=c("probit", "N")
)
summary(m3)

## it does not use the sandwich package, but does not it's own SE adjustments
m3 <- adjCov(m3, fdi$ccode)
summary(m3)

round(cbind(probit=m1$coefficients, 
            tsls=m2$coefficients,
            cf=theta.u,
            fiml=m3$coefficients[1:16]),3)


##marginal effect
xb <- drop(m3$X1 %*% m3$coefficients[1:16])
ame.fiml <- m3$coefficients["cr_rfdicap"] * mean(dnorm(xb))
ame.fiml

## with standard error


Db <- colMeans(-m3$X1* m3$coefficients["cr_rfdicap"] *dnorm(xb)*xb)
Db[2] <- Db[2] + mean(dnorm(xb))
Vame <- Db %*% m3$Vb[1:16, 1:16] %*% Db
c(ame.fiml, sqrt(Vame), ame.fiml/sqrt(Vame), 
  2*pnorm(abs(ame.fiml/sqrt(Vame)), lower=F))


## What about with a nice plot
newData <- apply(m3$X1, 2, \(x){
  if(all(x %in% c(0,1))){
    return(median(x))
  }else{
    return(mean(x))
  }})
newData <- t(replicate(50, newData))
newData[,2] <- seq(-5,10,length=50)
phat <- pnorm(newData %*% m3$coefficients[1:16])

## clarify approach
Bmat <- mvrnorm(500, m3$coefficients[1:16], m3$Vb[1:16, 1:16])
pmat <- pnorm(newData%*%t(Bmat ) )
CI <- rowQuantiles(pmat, probs= c(0.025, .975))

plot.df <- data.frame(Est=phat, lo=CI[,1], hi=CI[,2], FDI=newData[,2])
ggplot(plot.df)+
  geom_ribbon(aes(x=FDI, ymin=lo, ymax=hi), alpha=.2)+
  geom_line(aes(x=FDI, y=phat))+
  ylab("Pr. of civil conflict")+
  xlab("Cube root of net FDI")