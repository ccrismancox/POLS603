library(readstata13)
library(MASS)
library(lmtest)
library(sandwich)
library(matrixStats)
library(dplyr)

rm(list=ls())

jklp <- read.dta13("datasets/Kucik_Peritz_JOP2020.dta")


dat <-subset(jklp, comply_b==1|comply_b==0,
             select = c("ds", "comp1", "resp1",
                        "comply_b", "comply_n", "lnthird", 
                        "trade_share_rctp", "gdp_share", "r_greatpower",
                        "times_ruled", "lnclaims", "article_xxii","remedy","clarity",
                        "dispute_combined")) 
dat$id <- as.numeric(as.factor(dat$resp1))
dat <- dat %>% arrange(id)

dat$ln.trade_share_rctp<-log(dat$trade_share_rctp+.05)

table(dat$comply_n) #non compliance, partial, full

# Ordered logit model with baseline controls 
m1 <- polr(factor(comply_n) ~ lnthird+
             ln.trade_share_rctp+
             times_ruled
           +gdp_share+r_greatpower+
             article_xxii+remedy+clarity,
           data = dat, 
           model=TRUE,
           Hess = TRUE,
           method = "logistic")
summary(m1)
V1 <-vcovCL(m1, ~resp1)
coeftest(m1, vcov=V1)


m2 <- polr(factor(comply_n) ~ lnthird+
             ln.trade_share_rctp+
             times_ruled
           +gdp_share+r_greatpower+
             article_xxii+remedy+clarity,
           data = dat, 
           model=TRUE,
           Hess = TRUE,
           method = "probit")
summary(m2)
V2 <-vcovCL(m2, ~resp1)
coeftest(m2, vcov=V2)

## marginal of lnthird
cuts <- rep(c(-Inf, m1$zeta, Inf), each=nrow(m1$model))
cuts <- matrix(cuts, nrow=nrow(m1$model))
X <- as.matrix(m1$model[,-1]) # drop the outcome
xb <- drop(X%*%m1$coefficients)
b.lnthird <- m1$coef["lnthird"]
marginals <- colMeans(b.lnthird*dlogis(cuts[,-ncol(cuts)] -xb))-colMeans(b.lnthird*dlogis(cuts[,-1]-xb))
marginals


Bmat <- mvrnorm(1000, c(m1$coef, m1$zeta),V1 )
tau.sim <- Bmat[,(ncol(X)+1):ncol(Bmat)]
xb.sim <- drop(X%*%t(Bmat[,1:ncol(X)]))
mar.out <- matrix(0, ncol=3, nrow=1000)
for(i in 1:1000){
  cuts <- rep(c(-Inf, tau.sim[i,], Inf), each=nrow(m1$model))
  cuts <- matrix(cuts, nrow=nrow(m1$model))
  
  b.hat <- Bmat[i,"lnthird"]
  mar.out[i,] <- colMeans(b.hat*dlogis(cuts[,-ncol(cuts)] -xb.sim[,i]))-
    colMeans(b.hat*dlogis(cuts[,-1]-xb.sim[,i]))
  
  
}
rbind(marginals,
      colSds(mar.out))


Xbar <- apply(X,2, \(x){
  if(all(x %in% c(0,1))){
    rep(median(x),12)
  }else{
    rep(mean(x), 12)
  }
})
Xbar[,"lnthird"] <- log(seq(1, 24, by=2))
cuts <- rep(c(-Inf, m1$zeta, Inf), each=12)
cuts <- matrix(cuts, nrow=12)
xb <- drop(Xbar%*%m1$coefficients)

phat <- plogis(cuts[,-1] -xb)-plogis(cuts[,-ncol(cuts)]-xb)

plot.df <- data.frame(Prob=c(phat),
                      Third=seq(1, 24, by=2),
                      Outcome=rep(c("Non compliance", "Partial", "Full"), each=12))


ggplot(plot.df)+
  geom_line(aes(x=Third,y=Prob, color=Outcome))
