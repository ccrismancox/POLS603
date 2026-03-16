library(brglm2)
library(detectseparation)
library(car)
L.cauchy.logit <- function(b,Z){
  return(-sum(plogis(Z%*%b, log=TRUE))-
           sum(c(dcauchy(b[1],scale=10,log=TRUE),
                 dcauchy(b[-1],scale=2.5,log=TRUE)
           )))
}
r.cauchy.logit <- function(b,Z){
  score <- colSums(Z*(1-plogis(drop(Z%*%b))))
  penalty <- c(-2*b[1]/(b[1]^2+100), 
               -0.32*b[-1]/(0.16*b[-1]^2+1))
  return(-score-penalty)
}
L.logF.logit <-  function(b,Z){
  return(-sum(plogis(Z%*%b, log=TRUE))-
           sum(b[-1]/2-log(exp(b[-1])+1)))
}
r.logF.logit <- function(b,Z){
  score <- colSums(Z*(1-plogis(drop(Z%*%b))))
  penalty <- c(0,(1/2-plogis(b[-1])))
  return(-score-penalty)
}


fulldata <- read.csv('Code/R/datasets/Enemies_Replication.csv')

vars <- c('FSubs', 'FinalGood', 'CompShareZipCoal', 
          'ZipCoalField', 'MarketShare', 'Productivity', 
          'NetPPE', 'MarketCap',
          'HighManufEnergy',
          'electricpowergeneration', 'Support')
data.use <- subset(fulldata, select=vars)
data.use <- na.omit(data.use)

data.use$MNC1Final0 <- 1*(data.use$FSubs == 1 & data.use$FinalGood == 0)
data.use$MNC0Final1 <- 1*(data.use$FSubs == 0 & data.use$FinalGood == 1)
data.use$MNC1Final1 <- 1*(data.use$FSubs == 1 & data.use$FinalGood == 1)

m1 <- glm(Support ~ MNC1Final0 + MNC0Final1 + MNC1Final1 +
            CompShareZipCoal + ZipCoalField + MarketShare + Productivity + 
            NetPPE + MarketCap + HighManufEnergy + electricpowergeneration,
          family=binomial, data=data.use, x=TRUE, y=TRUE)
summary(m1)
detect_separation(x=m1$x, y=m1$y, family=binomial())

m2 <- glm(Support ~ MNC1Final0 + MNC0Final1 + MNC1Final1 +
            CompShareZipCoal + ZipCoalField + MarketShare + Productivity + 
            NetPPE + MarketCap + HighManufEnergy + electricpowergeneration,
          data = data.use, family =binomial, method="brglmFit")
coeftest(m2, vcovHC)


Z <- as.matrix(m1$x) *(2*m1$y-1)

bhat.br2 <- optim(m2$coef,
                  L.cauchy.logit, gr=r.cauchy.logit, 
                  Z=Z,
                  method="BFGS")$par
bhat.br3 <- optim(m2$coef,
                  L.logF.logit, gr=r.logF.logit, 
                  Z=Z,
                  method="BFGS")$par

round(cbind(m1$coef, m2$coef, bhat.br2, bhat.br3), 3)
linearHypothesis(m1, "MNC1Final0+MNC1Final1-MNC0Final1", vcov=vcovHC)


linearHypothesis(m2, "MNC1Final0+MNC1Final1-MNC0Final1", vcov=vcovHC)
linearHypothesis(m2, "MNC1Final0-MNC1Final1", vcov=vcovHC)
