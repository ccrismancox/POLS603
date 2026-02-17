## basics
library(readstata13)
library(dplyr)
## econometric
library(lmtest)
library(sandwich)
library(MASS)
library(car)
## tables and figures
library(xtable)
library(ggplot2)

rm(list=ls())


protests <- read.dta13("datasets/remittances_and_protests1.dta")
protests <- subset(protests, sample==1)

f1 <- banks_protest_count ~ remit*dict + 
  l1gdp + #GDP pc capita
  l1pop+ #population
  l1nbr5 + #neighboring protests
  l12gr+ #Economic growth
  l1migr+ #Net migration
  elec3+ #election
  factor(cowcode)-1

## Poisson: ML-Dummy variable
poi.mldv <- glm(f1, data=protests,family=poisson, y=TRUE, x=TRUE)

## Mundlak setup
protests <- protests %>% 
  mutate(remit.dict = remit*dict) %>%
  group_by(cowcode)%>%
  mutate(across(c(names(poi.mldv$coef)[grep("factor|:",names(poi.mldv$coef),
                                            invert = TRUE)],
                  "remit.dict"),
                .fns=mean,
                .names="{.col}.bar")) %>%
  ungroup()

## Mundlak formula
f2 <- update(f1, .~.-factor(cowcode)+remit.bar + dict.bar + 
               l1gdp.bar + l1pop.bar + l1nbr5.bar + l12gr.bar + 
               l1migr.bar + elec3.bar + remit.dict.bar+1)
## fit the Mundlak
poi.mundlak <- glm(f2,  data=protests,family=poisson, y=TRUE, x=TRUE)

cbind(poi.mldv$coef[grep("factor",names(poi.mldv$coef),
                         invert = TRUE)],
      poi.mundlak$coef[grep("bar|Intercept",names(poi.mundlak$coef),
                            invert = TRUE)])

## comparing marginal effects
X0 <- X1 <- poi.mldv$x
X0[,"dict"] <- X0[,"remit:dict"] <- 0
X1[,"dict"] <- 1
X1[, "remit:dict"] <- X1[,"remit"]

AME.mldv <- c(poi.mldv$coef["remit"] * mean(exp(X0 %*% poi.mldv$coef)),
              (poi.mldv$coef["remit"]+poi.mldv$coef["remit:dict"]) *
                mean(exp(X1 %*% poi.mldv$coef)))

## drop the all 0s
protests <- protests %>% 
  mutate(sum.y = sum(banks_protest_count), .by=cowcode)

table(protests$sum.y==0)
X0 <- X1 <- poi.mldv$x[protests$sum.y>0,]
X0[,"dict"] <- X0[,"remit:dict"] <- 0
X1[,"dict"] <- 1
X1[, "remit:dict"] <- X1[,"remit"]

AME.mldv2 <- c(poi.mldv$coef["remit"] * mean(exp(X0 %*% poi.mldv$coef)),
               (poi.mldv$coef["remit"]+poi.mldv$coef["remit:dict"]) *
                 mean(exp(X1 %*% poi.mldv$coef)))



## mundlak version
X0 <- X1 <- poi.mundlak$x
X0[,"dict"] <- X0[,"remit:dict"] <- 0
X1[,"dict"] <- 1
X1[, "remit:dict"] <- X1[,"remit"]

AME.mundlak <- c(poi.mundlak$coef["remit"] * mean(exp(X0 %*% poi.mundlak$coef)),
                 (poi.mundlak$coef["remit"]+poi.mundlak$coef["remit:dict"]) *
                   mean(exp(X1 %*% poi.mundlak$coef)))

compMat <- rbind(AME.mldv,AME.mldv2,AME.mundlak)
colnames(compMat) <- c("Demo.", "Auto.")
compMat





############ Instrument version##################






library(ivreg) ##compare to 2sls
protests <- protests %>% 
  mutate(dist =(log(1+(1/(dist_coast)))) ,
         distwremit = log(1+((richremit/1000000)*(dist))))

tsls <- ivreg(log(banks_protest_count+1) ~ remit*dict + 
                l1gdp + #GDP pc capita
                l1pop+ #population
                l1nbr5 + #neighboring protests
                l12gr+ #Economic growth
                l1migr+ #Net migration
                elec3+#election
                factor(cowcode)-1|
                distwremit*dict + 
                l1gdp + #GDP pc capita
                l1pop+ #population
                l1nbr5 + #neighboring protests
                l12gr+ #Economic growth
                l1migr+ #Net migration
                elec3+#election
                factor(cowcode)-1,
              data=protests)
coeftest(tsls, vcovCL(tsls, protests$cowcode))[grep("factor",
                                                    names(tsls$coefficients),
                                                    invert=TRUE),]


## first stages 
m1 <- lm(remit~ distwremit*dict + 
           l1gdp + #GDP pc capita
           l1pop+ #population
           l1nbr5 + #neighboring protests
           l12gr+ #Economic growth
           l1migr+ #Net migration
           elec3+#election
           factor(cowcode)-1,
         x=TRUE,y=TRUE,
         data=protests)
m2 <- lm(I(remit*dict)~ distwremit*dict + 
           l1gdp + #GDP pc capita
           l1pop+ #population
           l1nbr5 + #neighboring protests
           l12gr+ #Economic growth
           l1migr+ #Net migration
           elec3+#election
           factor(cowcode)-1,
         x=TRUE,y=TRUE,
         data=protests)


protests$v1 <- m1$residuals
protests$v2 <- m2$residuals

poi.iv <- glm(banks_protest_count ~ remit*dict + 
                l1gdp + #GDP pc capita
                l1pop+ #population
                l1nbr5 + #neighboring protests
                l12gr+ #Economic growth
                l1migr+ #Net migration
                elec3+ #election
                v1 + v2+
                factor(cowcode)-1,
              data=protests,
              family=poisson,
              x=TRUE,
              y=TRUE)

coeftest(poi.iv, vcovCL(poi.iv, protests$cowcode))[grep("factor",
                                                        names(poi.iv$coefficients),
                                                        invert=TRUE),]



ell <- length(m1$coef)
Vs1 <- rbind( cbind(vcovCL(m1,~cowcode), matrix(0,ell,ell)),
              cbind(matrix(0,ell,ell),vcovCL(m2,~cowcode)))
Vs2 <- vcovCL(poi.iv, ~cowcode)
V2 <- vcov(poi.iv)

V1 <- rbind( cbind(vcov(m1), matrix(0,ell,ell)),
             cbind(matrix(0,ell,ell),vcov(m2)))

## C is the cross deriv D
C0 <- t(m1$x )%*% diag(poi.iv$fitted.values) %*% (poi.iv$x)
C1 <- poi.iv$coefficients["v1"] * C0
C2 <- poi.iv$coefficients["v2"] * C0
C <- rbind(C1, C2)

## R is the OPG between steps 
split.matrix <- function(x, f){
  lapply(split(x = seq_len(nrow(x)), f = f),
         function(ind) x[ind, , drop = FALSE])
  
}

Jpoi <- split.matrix(poi.iv$x*(poi.iv$y-poi.iv$fitted), protests$cowcode)
Jnormal <- cbind(m1$x*(m1$residuals)/mean(m1$residuals^2),
                 m2$x*(m2$residuals)/mean(m2$residuals^2))
J1 <-  split.matrix(Jnormal, protests$cowcode)
R  <- Reduce(`+`,
             lapply(1:length(Jpoi),
                    \(x){
                      return(colSums(J1[[x]]) %*% t(colSums(Jpoi[[x]])))
                    }
             )
)


V2step <- Vs2 + V2 %*% (t(C) %*% Vs1 %*% C- 
                          t(R) %*% V1 %*% C -
                          t(C) %*% V1 %*% R )%*% V2

coeftest(poi.iv, V2step)[grep("factor",
                              names(poi.iv$coefficients),
                              invert=TRUE),]
