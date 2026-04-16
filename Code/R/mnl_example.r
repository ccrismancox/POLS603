library(readstata13)
library(sandwich)
library(lmtest)
library(car)
library(MASS)
library(mlogit)

rm(list=ls())
kb  <- read.dta13("datasets/TRdataset_JCR2014_replication_final.dta")
table(kb$technologyrebellion)
kb$ln.milper <- log(kb$milper)

## reshape to match mlogit's demands
kb2 <- dfidx(kb,  shape="wide",choice="technologyrebellion")

## note the mlogit formula is 


## y ~ Z | X so for a basic mnl we have 
## y ~ 0 | X

m1 <-mlogit(technologyrebellion~0|poscoldwar91+
              roughterrain+ethnicwar+
              gdpcapita_fl+ln.milper, data=kb2,
            reflevel="Irregular")
summary(m1)
V1 <- vcovCL(m1, kb[names(m1$fitted),]$countryname)

## comparing conventional to irregular
coeftest(m1, V1)[seq(1,12,by=2),] |> round(4)

## comparing SNC to irregular
coeftest(m1, V1)[seq(2,12,by=2),] |> round(4)


## we can also fit the mnl as separate logits with a loss of efficiency
glm(I(technologyrebellion=="Conventional")~ poscoldwar91+
      roughterrain+ethnicwar+
      gdpcapita_fl+ln.milper, data=kb, family=binomial,
    subset=technologyrebellion!="SNC")
glm(I(technologyrebellion=="SNC")~poscoldwar91+
      roughterrain+ethnicwar+
      gdpcapita_fl+ln.milper, data=kb, family=binomial,
    subset=technologyrebellion!="Irregular")


## MN probit
m2 <-mlogit(technologyrebellion~0|poscoldwar91+roughterrain+
              ethnicwar+gdpcapita_fl+ln.milper,
              data=kb2, probit=TRUE,
            reflevel="Irregular")
summary(m2)

compareCoefs(m1,m2)




## Interpretation: First differences for cold war

X <- model.matrix(~poscoldwar91+
                       roughterrain+ethnicwar+
                       gdpcapita_fl+ln.milper,
                     data=kb)
X1 <- X0 <- X 
X1[,"poscoldwar91post-1990"] <- 1
X0[,"poscoldwar91post-1990"] <- 0

expX1.c <- exp(X1 %*% m1$coef[grep("Conventional",names(m1$coef))])
expX1.s <- exp(X1 %*% m1$coef[grep("SNC",names(m1$coef))])
expX0.c <- exp(X0 %*% m1$coef[grep("Conventional",names(m1$coef))])
expX0.s <- exp(X0 %*% m1$coef[grep("SNC",names(m1$coef))])


prob1 <- cbind(1, expX1.c, expX1.s)/ drop( 1+expX1.c+expX1.s)
prob0 <- cbind(1, expX0.c, expX0.s)/ drop( 1+expX0.c+expX0.s)

FD.coldwar <- colMeans(prob1-prob0)

## simulated SE 
Bmat <- mvrnorm(1000, m1$coef, V1)

expX1.c <- exp(X1 %*% t(Bmat[,grep("Conventional",names(m1$coef))]))
expX1.s <- exp(X1 %*% t(Bmat[,grep("SNC",names(m1$coef))]))
expX0.c <- exp(X0 %*% t(Bmat[,grep("Conventional",names(m1$coef))]))
expX0.s <- exp(X0 %*% t(Bmat[,grep("SNC",names(m1$coef))]))

denom1 <- 1+expX1.c+expX1.s
denom0 <- 1+expX0.c+expX0.s

prob1.i <- 1/denom1
prob1.c <- expX1.c/denom1
prob1.s <- expX1.s/denom1

prob0.i <- 1/denom0
prob0.c <- expX0.c/denom0
prob0.s <- expX0.s/denom0

SE <- c(sd(colMeans(prob1.i-prob0.i)),
        sd(colMeans(prob1.c-prob0.c)),
        sd(colMeans(prob1.s-prob0.s)))
FDs <- rbind(FD.coldwar, SE)
colnames(FDs) <- c("AME Pr(Irregular)", "AME Pr(Conventional)", "AME Pr(SNC)")
FDs
## what do we learn here?


## what about a continuous variable. say logged milper
expX.c <- exp(X %*% m1$coef[grep("Conventional",names(m1$coef))])
expX.s <- exp(X %*% m1$coef[grep("SNC",names(m1$coef))])
denom <- 1+expX.c+expX.s
me.milper.i <- (-m1$coefficients["ln.milper:Conventional"]*expX.c- 
                  m1$coefficients["ln.milper:SNC"]*expX.s)/denom^2
me.milper.c <- (m1$coefficients["ln.milper:Conventional"]*expX.c)/denom + 
  me.milper.i*expX.c
me.milper.s <- (m1$coefficients["ln.milper:SNC"]*expX.c)/denom + 
  me.milper.i*expX.s
AMEs <- colMeans(cbind(me.milper.i,me.milper.c,me.milper.s))


#### standard errors
expX.c <- exp(X %*% t(Bmat[,grep("Conventional",names(m1$coef))]))
expX.s <- exp(X %*% t(Bmat[,grep("SNC",names(m1$coef))]))
denom <- 1+expX.c+expX.s
me.milper.i.sim <- (-m1$coefficients["ln.milper:Conventional"]*expX.c- 
                  m1$coefficients["ln.milper:SNC"]*expX.s)/denom^2
me.milper.c.sim <- (m1$coefficients["ln.milper:Conventional"]*expX.c)/denom + 
  me.milper.i.sim*expX.c
me.milper.s.sim <- (m1$coefficients["ln.milper:SNC"]*expX.c)/denom + 
  me.milper.i.sim*expX.s
se <- c(sd(colMeans(me.milper.i.sim)),
        sd(colMeans(me.milper.c.sim)),
        sd(colMeans(me.milper.s.sim)))
AMEs<-rbind(AMEs,se)
colnames(AMEs)<-c("Irregular", "Conventional", "SNC")
AMEs

## plot predicted probs
summary(kb$ln.milper)
truehist(kb$ln.milper)

Xbar <- apply(X,2, \(x){
  if(all(x %in% c(0,1))){
    rep(median(x),25)
  }else{
    rep(mean(x), 25)
  }
})

Xbar[,"ln.milper"] <- seq(0,8, length=25)
expX.c <- exp(Xbar %*% m1$coef[grep("Conventional",names(m1$coef))])
expX.s <- exp(Xbar %*% m1$coef[grep("SNC",names(m1$coef))])
denom <- drop(1+expX.c+expX.s)
phat <- cbind(1, expX.c,expX.s)/denom
#### CI
expX.c <- exp(Xbar %*% t(Bmat[,grep("Conventional",names(m1$coef))]))
expX.s <- exp(Xbar %*% t(Bmat[,grep("SNC",names(m1$coef))]))
denom <- 1+expX.c+expX.s

prob.i <- 1/denom
prob.c <- expX.c/denom
prob.s <- expX.s/denom
CI <- rbind(rowQuantiles(prob.i, probs=c(0.025, 0.975)),
            rowQuantiles(prob.c, probs=c(0.025, 0.975)),
            rowQuantiles(prob.s, probs=c(0.025, 0.975)))
plot.df <- data.frame(Pr=c(phat),
                      lo=CI[,1],
                      hi=CI[,2],
                      milper=seq(0,8, length=25),
                      Outcome=rep(c("Irregular", "Conventional", "SNC"),
                                  each=25))
ggplot(plot.df)+
  geom_ribbon(aes(x=milper, ymin=lo, ymax=hi,fill=Outcome),
              alpha=.3)+
  geom_line(aes(x=milper, y=Pr, color=Outcome))+
  xlab("Logged military size")+
  ylab("Probability")+
  theme_bw(14)+
  theme(legend.position = "bottom")
 
