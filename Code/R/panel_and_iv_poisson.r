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
           elec3+ #election
           factor(cowcode)-1,
         x=TRUE,y=TRUE,
         data=protests)
m2 <- lm(I(remit*dict)~ distwremit*dict + 
           l1gdp + #GDP pc capita
           l1pop+ #population
           l1nbr5 + #neighboring protests
           l12gr+ #Economic growth
           l1migr+ #Net migration
           elec3+ #election
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



ell <- length(m1$coef) ## size of the instrument matrix 

## First stage covariance matrix is block diagonal
## [V1,1    0 ]
## [0     V1,2]
## We fit these independently so there is no cross correlation.
## We could have tried to fit these as a system with
V1 <- rbind( cbind(vcov(m1),           matrix(0,ell,ell)),
             cbind(matrix(0,ell,ell),  vcov(m2)))
Vcl1 <- rbind( cbind(vcovCL(m1,~cowcode), matrix(0,ell,ell)),
               cbind(matrix(0,ell,ell),   vcovCL(m2,~cowcode)))


## second stage matricies
Vcl2 <- vcovCL(poi.iv, ~cowcode)
V2 <- vcov(poi.iv)

## C is the cross deriv D
lambda.hat <- poi.iv$fitted.values
bonus.term <- c( t(lambda.hat-poi.iv$y ) %*% m1$x)
C0 <- t(poi.iv$x)%*% diag(lambda.hat) %*% m1$x
C1 <- poi.iv$coefficients["v1"] * C0
C2 <- poi.iv$coefficients["v2"] * C0
C1["v1",] <-C1["v1",]  +bonus.term
C2["v2",] <-C2["v2",]  +bonus.term
C <- -cbind(C1, C2)




## R is the OPG between steps.
## Since we're clustering we want to sum within units first
## create function to split our big matricies into smaller matricies
## by units
split.matrix <- function(x, f){
  lapply(split(x = seq_len(nrow(x)), f = f),
         function(ind) x[ind, , drop = FALSE])
  
}

## Poisson Jacobian split by cowcode
Jpoi <- split.matrix(poi.iv$x*(poi.iv$y-poi.iv$fitted), protests$cowcode)

## Two Normal Jacobians split by cowcode
Jnormal <- cbind(m1$x*(m1$residuals)/mean(m1$residuals^2),
                 m2$x*(m2$residuals)/mean(m2$residuals^2))
J1 <-  split.matrix(Jnormal, protests$cowcode)

## apply over cowcodes, return the outer produce of the within-sums
## Then add these up
R  <- Reduce(`+`,
             lapply(1:length(Jpoi),
                    \(x){
                      return(colSums(Jpoi[[x]]) %*% t(colSums(J1[[x]])))
                    }
             )
)


V2step <- Vcl2 + V2 %*% (C %*% Vcl1 %*% t(C)- 
                           R %*% V1 %*% t(C) -
                           C %*% V1 %*% t(R) )%*% V2

coeftest(poi.iv, V2step)[grep("factor",
                              names(poi.iv$coefficients),
                              invert=TRUE),]


## some packages for parallel computing
library(doParallel) ## main work horse
library(doRNG) ## an update that allows us to set a seed in parallel

## number of available cores
detectCores()

#let's use 3/4 of them
ncores <- floor(.75*detectCores())


cl <- makeCluster(ncores) # creates copies of R that can run side-by-side
registerDoParallel(cl) # this connects our cores to the foreach loop we're about to run

## instead of for we now have foreach
## Differences: 
######  1. save the output from the loop up top rather
######     than in pre-initalized matrix
######  2. b=1:B rather than b in 1:B
######  3. It has additional options for how to combine the output,
######     whether to export any packages to our other processes,
######     what to do with errors, and so on. See ?foreach for more
######  4. We give it a "do" option. %do% is the same as for (not parallel)
######     %dopar% is parallel and %dorng% is parallel but it respects the seed.

set.seed(1)
B <- 500
bootMat <- foreach(b=1:B, .combine = "rbind")%dorng%{
  id2 <-  sample(unique(protests$cowcode), replace=TRUE)
  boot.df <- sapply(id2, \(x){subset(protests, cowcode==x)},
                    simplify=FALSE)
  boot.df <- lapply(1:length(boot.df), 
                    \(x){
                      boot.df[[x]]$cowcode <- x
                      return(boot.df[[x]])
                    })
  boot.df <- do.call(rbind, boot.df)
  
  m1 <- lm(remit~ distwremit*dict + 
             l1gdp + #GDP pc capita
             l1pop+ #population
             l1nbr5 + #neighboring protests
             l12gr+ #Economic growth
             l1migr+ #Net migration
             elec3+ #election
             factor(cowcode)-1,
           x=TRUE,y=TRUE,
           data=boot.df)
  m2 <- lm(I(remit*dict)~ distwremit*dict + 
             l1gdp + #GDP pc capita
             l1pop+ #population
             l1nbr5 + #neighboring protests
             l12gr+ #Economic growth
             l1migr+ #Net migration
             elec3+ #election
             factor(cowcode)-1,
           x=TRUE,y=TRUE,
           data=boot.df)
  
  
  boot.df$v1 <- m1$residuals
  boot.df$v2 <- m2$residuals
  
  step2 <- glm(banks_protest_count ~ remit*dict + 
                 l1gdp + #GDP pc capita
                 l1pop+ #population
                 l1nbr5 + #neighboring protests
                 l12gr+ #Economic growth
                 l1migr+ #Net migration
                 elec3+ #election
                 v1 + v2+
                 factor(cowcode)-1,
               data=boot.df,
               family=poisson,
               x=TRUE,
               y=TRUE)
  out <- step2$coefficients[c(1:10,length(step2$coef))]
  if(any(abs(out)>100)){ #sometimes weird things happen. Poisson coefficients should not be more than 100. definite sign of numeric issue
    out <- rep(NA, length(out))
  }
  out
  
}
stopCluster(cl) ## close those other copies so they're not just creeping in the background

dim(bootMat)
dim(na.omit(bootMat))


cbind(Est=poi.iv$coefficients[c(1:10,length(poi.iv$coef))],
      AnalyticSE=sqrt(diag(V2step)[c(1:10,length(poi.iv$coef))]),
      BootstrapSE=sqrt(diag(var(bootMat,na.rm=TRUE))))

