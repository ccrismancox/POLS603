library(readstata13)
library(sandwich)
library(car)
library(ggplot2)
rm(list=ls())

pm.dat <- read.dta13("datasets/Meng_Paine_APSR_data.dta")
pm.dat <- subset(pm.dat, sample==1)

## 50 countries observed for 25-47 years
length(unique(pm.dat$ccode))
summary(as.numeric(table(pm.dat$year)))


## LPM 
lpm <- lm(fail_main~ rebelregime +
            ln_gdppc+ gdpgrowth+ ln_oilpop +
            ln_pop+ 
            al_ethnic + relfrac +
            col_british+ col_french+ col_port+
            regime_duration_main,
          data=pm.dat)
summary(lpm$fitted.values)
table(lpm$fitted.values <0)/length(lpm$fitted.values)

## logit
l1 <- glm(fail_main~ rebelregime +
            ln_gdppc+ gdpgrowth+ ln_oilpop +
            ln_pop+ 
            al_ethnic + relfrac +
            col_british+ col_french+ col_port+
            regime_duration_main,
          data=pm.dat, 
          family=binomial,x=TRUE)

## probit
p1 <- glm(fail_main~ rebelregime +
            ln_gdppc+ gdpgrowth+ ln_oilpop +
            ln_pop+ 
            al_ethnic + relfrac +
            col_british+ col_french+ col_port+
            regime_duration_main,
          data=pm.dat, 
          family=binomial(link="probit"), x=TRUE) 


Vlpm <- vcovCL(lpm, cluster=pm.dat$ccode)
Vl1 <- vcovCL(l1, cluster=pm.dat$ccode)
Vp1 <- vcovCL(p1, cluster=pm.dat$ccode)


compareCoefs(lpm,l1, p1, 
             vcov.=list(Vlpm,Vl1, Vp1),
             zvals=TRUE,pvals=TRUE)


## Odds ratio (multiplicative effect)
deltaMethod(l1, "exp(rebelregime)", vcov=Vl1)


## % change in the odds
deltaMethod(l1, "(exp(rebelregime)-1)*100", vcov=Vl1)


## marginal effects for the logit and probit
## Since rebel regime is binary, let's do the first difference
## E[E[y|rebelregime=1,X] - E[y|rebelregime=0,X]]
X0 <- X1 <- l1$x
X1[,"rebelregime"] <- 1
X0[,"rebelregime"] <- 0
lpl1 <- drop(X1 %*% l1$coef)
lpl0 <- drop(X0 %*% l1$coef)
lpp1 <- drop(X1 %*% p1$coef)
lpp0 <- drop(X0 %*% p1$coef)
fd.logit <- mean(plogis(lpl1) - plogis(lpl0))
fd.probit <- mean(pnorm(lpp1) - pnorm(lpp0))


Dbeta.logit <- colMeans(X1*dlogis(lpl1)- X0 * dlogis(lpl0))
Dbeta.probit <- colMeans(X1*dnorm(lpp1)- X0 * dnorm(lpp0))
Vfd.logit <- Dbeta.logit %*% Vl1 %*% Dbeta.logit
Vfd.probit <- Dbeta.probit %*% Vp1 %*% Dbeta.probit

## combine and compare
effects <- rbind(c(lpm$coefficients[2], sqrt(diag(Vlpm)[2])),
                 c(fd.logit, sqrt(Vfd.logit)),
                 c(fd.probit, sqrt(Vfd.probit)))
effects <- cbind(effects, effects[,1]/effects[,2])
effects <- cbind(effects, 2*pnorm(abs(effects[,3]),lower=FALSE))
colnames(effects) <- c("Marginal effect", "St. Err.", "z-stat", "p-value")
rownames(effects) <- c("LPM", "Logit", "Probit")
round(effects, 3)

## Let's try a continous variable with an AME
## growth
summary(pm.dat$gdpgrowth)

## one unit isn't a real plausible change, maybe a 1 sd?
s <- sd(pm.dat$gdpgrowth) 
s 
mu.l <- l1$linear.predictors
mu.p <- p1$linear.predictors

ame.logit <- s*mean(l1$coefficients["gdpgrowth"]*dlogis(mu.l))
ame.probit <- s*mean(p1$coefficients["gdpgrowth"]*dnorm(mu.p))

Dbeta.logit <- l1$coefficients["gdpgrowth"] * (1-2*plogis(mu.l))*dlogis(mu.l)*l1$x
Dbeta.logit[,"gdpgrowth"] <- Dbeta.logit[,"gdpgrowth"]+dlogis(mu.l)
Dbeta.logit <- colMeans(Dbeta.logit)*s

Dbeta.probit <- p1$coefficients["gdpgrowth"] * dnorm(mu.p)*mu.p*p1$x
Dbeta.probit[,"gdpgrowth"] <- Dbeta.probit[,"gdpgrowth"]+dnorm(mu.p)
Dbeta.probit <- colMeans(Dbeta.probit)*s

Vame.logit <- Dbeta.logit %*% Vl1 %*% Dbeta.logit
Vame.probit <- Dbeta.probit %*% Vp1 %*% Dbeta.probit

## combine and compare
effects <- rbind(c(s*lpm$coefficients["gdpgrowth"],
                   sqrt(diag(s^2 * Vlpm)["gdpgrowth"])),
                 c(ame.logit, sqrt(Vame.logit)),
                 c(ame.probit, sqrt(Vame.probit)))
effects <- cbind(effects, effects[,1]/effects[,2])
effects <- cbind(effects, 2*pnorm(abs(effects[,3]),lower=FALSE))
colnames(effects) <- c("Marginal effect", "St. Err.", "z-stat", "p-value")
rownames(effects) <- c("LPM", "Logit", "Probit")
round(effects, 3)


## predicted probabilities for growth
hist(pm.dat$gdpgrowth,breaks = "scott", freq=FALSE)
Xbar  <- t(replicate(25,colMeans(l1$x)))
Xbar[,"gdpgrowth"] <- seq(-0.5, 0.5, length=25)



yhat.lpm <- Xbar %*% lpm$coef
yhat.logit <- plogis(Xbar%*%l1$coef)
yhat.probit <- pnorm(Xbar%*%p1$coef)

se.yhat.lpm <- sqrt(diag(Xbar %*% Vlpm %*% t(Xbar)))

dbeta <- (Xbar*dlogis(drop(Xbar%*% l1$coef)))
se.yhat.logit <- sqrt(diag(dbeta %*%  Vl1 %*% t(dbeta)))

dbeta <- (Xbar*dnorm(drop(Xbar%*% p1$coef)))
se.yhat.probit <- sqrt(diag(dbeta %*%  Vp1 %*% t(dbeta)))

plot.df <- data.frame(growth=Xbar[,"gdpgrowth"],
                      Pr=c(yhat.lpm,yhat.logit,yhat.probit),
                      se=c(se.yhat.lpm,se.yhat.logit,se.yhat.probit),
                      Model=factor(rep(c("LPM", "Logit", "Probit"),each=25),
                                   levels=c("LPM", "Logit", "Probit")))
plot.df$lo <-plot.df$Pr-1.96* plot.df$se
plot.df$hi <-plot.df$Pr+1.96* plot.df$se

ggplot(plot.df)+
  geom_ribbon(aes(x=growth,ymax = hi, ymin=lo, fill=Model),
              alpha=.3)+
  geom_line(aes(x=growth, y=Pr, color=Model))+
  facet_wrap(Model~., nrow=1)+
  theme_bw(14)+
  xlab("GDP growth")+
  ylab("Probability of regime failure")+
  theme(legend.position = "bottom")+
  geom_hline(yintercept = 0, alpha=.2, linetype="dotted")+
  geom_rug(aes(x=gdpgrowth),
           data=subset(pm.dat, gdpgrowth>=-.5 & gdpgrowth <= .5))
