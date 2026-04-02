library(matrixStats)
rm(list=ls())
lalonde <- read.csv("datasets/lalonde.csv")

## Baseline
m1 <- glm(re78~treat+age+I(age^2) 
          +educ +I(educ^2)+factor(race)+married
          +log(re74+1)+log(re75+1),
          family=quasipoisson,
          x=TRUE, y=TRUE,
          data=lalonde)

## effects
X <- m1$x
treated <- X[,"treat"]
X[,"treat"] <- 0

ATE1 <- (exp(m1$coef["treat"])-1)*mean(exp(X %*% m1$coef))
ATT1 <- (exp(m1$coef["treat"])-1)*mean(exp(X[treated==1,] %*% m1$coef))
c(ATE1, ATT1)


## propensity scores
pmod <- glm(treat~age+I(age^2) 
            +educ +I(educ^2)+factor(race)+married
            +log(re74+1)+log(re75+1)-1,
            family=binomial,
            x=TRUE, y=TRUE,
            data=lalonde)
summary(pmod)
w.ate <- treated/pmod$fitted.values + (1-treated)/(1-pmod$fitted.values)
w.att <- treated + (1-treated)*pmod$fitted.values/(1-pmod$fitted.values)
#### look for balance###

## no weights
Xnew <- pmod$x
diffs1 <- (colMeans(Xnew[treated==1,])-
             colMeans(Xnew[treated==0,]))/colSds(Xnew)


X.ate <- Xnew*w.ate
w.ate.std <- w.ate/sum(w.ate)
part1 <- (t(t(Xnew)- colSums(Xnew*w.ate.std)) )^2
varX.ate <- colSums(w.ate.std*part1)/(1-sum(w.ate.std^2))
sdX.ate <- sqrt(varX.ate)
diffs2 <- (colSums(X.ate[treated==1,])/sum(w.ate[treated==1])-
            colSums(X.ate[treated==0,])/sum(w.ate[treated==0]))/sdX.ate




X.att <- Xnew*w.att
w.att.std <- w.att/sum(w.att)
part1 <- (t(t(Xnew)- colSums(Xnew*w.att.std)) )^2
varX.att <- colSums(w.att.std*part1)/(1-sum(w.att.std^2))
sdX.att <- sqrt(varX.att)
diffs3 <- (colSums(X.att[treated==1,])/sum(w.att[treated==1])-
             colSums(X.att[treated==0,])/sum(w.att[treated==0]))/sdX.att



plot.df <- data.frame(Variable=names(diffs1),
                      weight=rep(c("None","ATE", "ATT"),
                                 each=length(diffs1)),
                      val=abs(c(diffs1,diffs2,diffs3))
)
plot.df <- plot.df %>%
  mutate(weight=factor(weight, levels=unique(weight)))%>%
  group_by(weight) %>%
  arrange(weight,val) %>%
  mutate(Variable = factor(Variable, levels=unique(Variable)))

ggplot(plot.df)+
  geom_point(aes(x=val,y=factor(Variable), color=weight))+
  stat_summary(aes(x=val,y=factor(Variable), color=weight, group=weight),
               fun.y=sum, geom="line")+
  geom_vline(aes(xintercept = 0))


## ipw estimators
y <- lalonde$re78
p.hat <- pmod$fitted.values

## HT
w1 <- (treated/p.hat) 
w0 <- (1-treated)/(1-p.hat) 
ate.ipw1 <- mean(w1*y)- mean(w0*y)


w1 <- treated 
w0 <-( (1-treated)*(p.hat)/(1-p.hat) )
att.ipw1 <- mean(w1*y)- mean(w0*y)
c(ate.ipw1, att.ipw1)

### Hajek
w1 <- (treated/p.hat) / mean(treated/p.hat)
w0 <- ((1-treated)/(1-p.hat)) / mean((1-treated)/(1-p.hat))
ate.ipw2 <- mean(w1*y)- mean(w0*y)


w1 <- treated / mean(treated)
w0 <-( (1-treated)*(p.hat)/(1-p.hat) )/ 
  mean((1-treated)*(p.hat)/(1-p.hat))
att.ipw2 <- mean(w1*y)- mean(w0*y)
c(ate.ipw2, att.ipw2)

## Also 
m2 <- lm(re78~treat, data=lalonde,
          x=TRUE, y=TRUE,
          weights=w.ate)
summary(m2)
ATE2 <- m2$coef['treat']



m2a <- glm(re78~treat + age + I(age^2) + 
             educ + I(educ^2) + 
             factor(race) + 
             married + 
             log(re74 + 1) + log(re75 + 1),
           x=TRUE,y=TRUE,
           family=quasipoisson,
           data=lalonde,weights=w.ate)
summary(m2a)
X <- m2a$x
X[,"treat"] <- 0
ATE2a <- (exp(m2a$coef["treat"])-1)*mean(exp(X %*% m2a$coef))
ATE2a


## ATT 
m2b <- lm(re78~treat, data=lalonde,
          x=TRUE, y=TRUE,
          weights=w.att)
summary(m2b)
ATT2 <- m2b$coef['treat']

## check on balance 
m2a <- glm(re78~treat + age + I(age^2) + 
             educ + I(educ^2) + 
             factor(race) + 
             married + 
             log(re74 + 1) + log(re75 + 1),
           x=TRUE,y=TRUE,
           family=quasipoisson,
           data=lalonde,weights=w.att)
summary(m2a)
X <- m2a$x
X[,"treat"] <- 0
ATT2a <- (exp(m2a$coef["treat"])-1)*mean(exp(X %*% m2a$coef))
ATT2a


## What if we went more non-parametric in the regression
## side with separate Poissons for each treatment and control
m3.treated <- glm(re78~age + I(age^2) + 
                    educ + I(educ^2) + 
                    factor(race) + 
                    married + 
                    log(re74 + 1) + log(re75 + 1),
                  subset=treat==1,
                  x=TRUE,y=TRUE,
                  family=quasipoisson,
                  data=lalonde)
m3.control <- glm(re78~age + I(age^2) + 
                    educ + I(educ^2) + 
                    factor(race) + 
                    married + 
                    log(re74 + 1) + log(re75 + 1),
                  subset=treat==0,
                  x=TRUE,y=TRUE,
                  family=quasipoisson,
                  data=lalonde)
X <- m1$x[,-2]
cATE <- exp(X%*% m3.treated$coefficients)-
  exp(X%*%m3.control$coefficients)
ATE3 <- mean(cATE)
ATT3 <- mean(cATE[treated==1])
c(ATE3, ATT3)

m3 <- glm(re78~treat*(age + I(age^2) + 
            educ + I(educ^2) + 
            factor(race) + 
            married + 
            log(re74 + 1) + log(re75 + 1)),
          x=TRUE,y=TRUE,
          family=quasipoisson,
          data=lalonde)
Xnew <- with(lalonde,
             cbind.data.frame(treat, age, educ,
                              race,married,re74,re75))
Xnew$treat <- 1
Ey.treat <- predict(m3, newdata=Xnew, type="response")
Xnew$treat <- 0
Ey.control <- predict(m3, 
                      newdata=Xnew, 
                      type="response")
cATE <- Ey.treat-Ey.control
ATE4 <- mean(cATE)
ATT4 <- mean(cATE[treated==1])
c(ATE4, ATT4)

## Double robust?
mu1 <- exp(X%*% m3.treated$coefficients)
mu0 <- exp(X%*% m3.control$coefficients)
double.ate <- mean(mu1-mu0 +
                     ( treated*(y-mu1)/p.hat- 
                         (1-treated)*(y-mu0)/(1-p.hat)))
double.ate

double.att <- sum(y*treated - (y*(1-treated)*p.hat+mu0*(treated-p.hat))/(1-p.hat))/sum(treated)
double.att

## A quick look at some canned versions
library(WeightIt)
library(cobalt)
# Estimate the pscore
pscore_w <- weightit(pmod$formula,
                     data = lalonde, 
                     method = "glm", 
                     estimand = "ATE")
summary(pscore_w$weights-w.ate)

## balance
love.plot(pscore_w,
          drop.distance = TRUE,
          var.order = "unadjusted",
          abs = TRUE,
          line = TRUE,
          estimand = "ATT", method = "weighting",
          thresholds = c(m = .1))

## estimate ATE
ate_ipw_w <- lm_weightit(re78 ~ treat, 
                         data = lalonde, 
                         weightit = pscore_w,
                         x=TRUE,
                         vcov = "asympt")
summary(ate_ipw_w)
X <- ate_ipw_w$x
X[,"treat"] <- 0
ATE5 <- ate_ipw_w$coef["treat"]
ATE5

## ATT 
# Estimate the pscore
pscore_w <- weightit(pmod$formula,
                     data = lalonde, 
                     method = "glm", 
                     estimand = "ATT")
summary(pscore_w$weights-w.att)

## balance
love.plot(pscore_w,
          drop.distance = TRUE,
          var.order = "unadjusted",
          abs = TRUE,
          line = TRUE,
          estimand = "ATT",
          method = "weighting",
          thresholds = c(m = .1))


## estimate ATT
att_ipw_w <- lm_weightit(re78 ~ treat, 
                          data = lalonde, 
                          weightit = pscore_w,
                          x=TRUE,
                          vcov = "asympt")
summary(att_ipw_w)
ATT5 <- att_ipw_w$coef["treat"]
ATT5


out<- cbind(c(ATE1, ate.ipw1,ate.ipw2,ATE2,ATE3,ATE4,ATE5, double.ate),
      c(ATT1, att.ipw1,att.ipw2,ATT2,ATT3,ATT4,ATT5, double.att))
colnames(out) <- c("ATE", "ATT")
row.names(out) <- c("Poisson", "IPW-HT", "IPW-Hajek",
                    "IPW-weights", "Fully interacted",
                    "Fully interacted 2",
                    "Canned IPW weights", "Double robust")
out      
