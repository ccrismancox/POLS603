library(readstata13)
library(sandwich)
library(lmtest)
library(ggplot2)
library(splines)
rm(list=ls())

pm.dat <- read.dta13("datasets/Meng_Paine_APSR_data.dta")
pm.dat <- subset(pm.dat, sample==1)
pm.dat$regime_duration_main <- pm.dat$regime_duration_main+1

pm.dat %>% 
  subset(select=c("country",
                  "year",
                  "regime_duration_main",
                  "fail_main"))%>%
  head()


#### Dummies ####
l1.dummies <- glm(fail_main~ rebelregime +
                    ln_gdppc+ gdpgrowth+ ln_oilpop +
                    ln_pop+ 
                    al_ethnic + relfrac +
                    col_british+ col_french+ col_port+
                    factor(regime_duration_main)-1,
                  data=pm.dat, 
                  family=binomial,x=TRUE)
Vd <-  vcovCL(l1.dummies, pm.dat$ccode)
coeftest(l1.dummies, Vd) %>%
  round(3)

## hazard plot
Xbar <- l1.dummies$x
Xbar <- Xbar[,grep("factor",
                   colnames(Xbar), 
                   invert=TRUE)]
Xbar <- apply(Xbar,2, \(x){
  if(all(x %in% c(0,1))){
    rep(median(x),58)
  }else{
    rep(mean(x), 58)
  }
})
Xbar <- cbind(Xbar, diag(58))
yhat <- plogis(Xbar %*% l1.dummies$coef)
Bmat <- mvrnorm(1000, l1.dummies$coef, Vd)
y.sim <- plogis(Xbar %*% t(Bmat))
CI <- rowQuantiles(y.sim, probs=c(0.025, 0.975))

plot.df <- data.frame(Hazard = yhat,
                      lo=CI[,1],
                      hi=CI[,2],
                      time=1:58)
ggplot(plot.df)+
  geom_pointrange(aes(x=time, y=Hazard,
                      ymin=lo,ymax=hi),
                  size=.3)

#### polynomials ####
l1.poly <- glm(fail_main~ rebelregime +
                 ln_gdppc+ gdpgrowth+ ln_oilpop +
                 ln_pop+ 
                 al_ethnic + relfrac +
                 col_british+ col_french+ col_port+
                 regime_duration_main+
                 I(regime_duration_main^2)+
                 I(regime_duration_main^3),
               data=pm.dat, 
               family=binomial,x=TRUE)
Vp <-  vcovCL(l1.poly, pm.dat$ccode)
coeftest(l1.poly, Vp) %>%
  round(3)
linearHypothesis(l1.poly, 
                 c("regime_duration_main",
                   "I(regime_duration_main^2)",
                   "I(regime_duration_main^3)"),
                 vcov=Vp)

## hazard plot
X <- l1.poly$x[,1:11]
Xbar <- apply(X,2, \(x){
  if(all(x %in% c(0,1))){
    rep(median(x),58)
  }else{
    rep(mean(x), 58)
  }
})
Xbar <- cbind(Xbar, (1:58), (1:58)^2, (1:58)^3)
yhat <- plogis(Xbar %*% l1.poly$coef)
Bmat <- mvrnorm(1000, l1.poly$coef, Vp)
y.sim <- plogis(Xbar %*% t(Bmat))
CI <- rowQuantiles(y.sim, probs=c(0.025, 0.975))

plot.df <- data.frame(Hazard = yhat,
                      lo=CI[,1],
                      hi=CI[,2],
                      time=1:58)
ggplot(plot.df)+
  geom_ribbon(aes(x=time, y=Hazard,
                  ymin=lo,ymax=hi),
              alpha=.3)+
  geom_line(aes(x=time, y=Hazard))

#### splines####
spline.mat <- ns(pm.dat$regime_duration_main,df=3)
## df=3 puts knots at 1, 2 interior points, and at the end (58)
attr(spline.mat, "knots")
attr(spline.mat, "Boundary.knots")


## note that the spline values while unintelligible on the surface
## map into the 58 time periods
s.df <- unique(cbind(pm.dat$regime_duration_main,
                     spline.mat))
head(s.df)
dim(s.df)
colnames(s.df) <- c("regime_duration_main",
                    "sp1", "sp2",
                    "sp3")
s.df <- as.data.frame(s.df)
pm.dat <- merge(pm.dat,s.df, 
                by="regime_duration_main",
                all.x=TRUE)


### polynomials ###
l1.splines <- glm(fail_main~ rebelregime +
                    ln_gdppc+ gdpgrowth+ ln_oilpop +
                    ln_pop+ 
                    al_ethnic + relfrac +
                    col_british+ col_french+ col_port+
                    sp1+sp2+sp3,
                  data=pm.dat, 
                  family=binomial,x=TRUE)
Vs <-  vcovCL(l1.splines, pm.dat$ccode)
coeftest(l1.splines, Vs) %>%
  round(3)
linearHypothesis(l1.splines, 
                 c("sp1", "sp2","sp3"),
                 vcov=Vs)

## hazard plot
X <- l1.splines$x[,1:11]
Xbar <- apply(X,2, \(x){
  if(all(x %in% c(0,1))){
    rep(median(x),58)
  }else{
    rep(mean(x), 58)
  }
})
Xbar <- cbind(Xbar, as.matrix(s.df[,2:4]))
yhat <- plogis(Xbar %*% l1.splines$coef)
Bmat <- mvrnorm(1000, l1.splines$coef, Vs)
y.sim <- plogis(Xbar %*% t(Bmat))
CI <- rowQuantiles(y.sim, probs=c(0.025, 0.975))

plot.df <- data.frame(Hazard = yhat,
                      lo=CI[,1],
                      hi=CI[,2],
                      time=1:58)
ggplot(plot.df)+
  geom_ribbon(aes(x=time, y=Hazard,
                  ymin=lo,ymax=hi),
              alpha=.3)+
  geom_line(aes(x=time, y=Hazard))

```
```{r, results='asis'}
## How much do the estimates change?

mod.list <- list(l1.dummies,l1.poly,l1.splines)
V.list<- list(Vd, Vp,Vs)
se.list <- lapply(V.list, \(x){sqrt(diag(x))})

est.list <- lapply(mod.list,
                   \(x){return(x$coef[grep("factor|Intercept|duration|sp",
                                           names(x$coef),
                                           invert=TRUE)])})
se.list <- lapply(se.list,
                  \(x){return(x[grep("factor|Intercept|duration|sp",
                                     names(x),
                                     invert=TRUE)])})
tab <- rbind(do.call(c,est.list),
             do.call(c,se.list))
tab <- round(tab,2)
tab[2,] <- paste0("(", tab[2,], ")")
tab <- matrix(tab, ncol=3)
names <- c("Rebel regime", 
           "GDP pc","Growth", "Oil",
           "Population","Ethnic frac", 
           "Religious frac", "Brit col", 
           "French col", "Port col")
names <- c(rbind(names, " "))
tab <- cbind(names, tab)
colnames(tab) <- c(" ", "Dummies", "Polynomial", "Spline")

print(xtable(tab, align="llccc",
             caption="Effect of rebellion on autocratic regime survival"),
      include.rownames=FALSE,
      sanitize.text.function = \(x){x},
      caption.placement="top")