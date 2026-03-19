rm(list=ls())

source("likelihood_functions.r")

N <- c(50, 100, 200,500)
B <- 250
beta <- c(2, -1, .25)
gamma <- c(2, 0, 0, -1)
tau <- .3


rho <- seq(-2/3, 2/3, by=1/3)
g1 <- rho*sqrt((.16+(pi**2)/3)/(1-rho^2))


conds <- expand.grid(gamma1=g1,N=N)

tic <- proc.time()[3]
results <- list()
for(k in 1:nrow(conds)){
  gamma[2] <- conds[k,1]
  n <- conds[k,2]
  
  output <- matrix(0, B,8)
  for(b in 1:B){
    X <- cbind(1, rnorm(n), rpois(n, lambda=10))
    Z <- cbind(X, rbinom(n, 1, prob=.8))
    u <- rnorm(n, mean=0, sd=25)
    v <- rlogis(n)
    
    d <- Z %*% gamma + v
    cor(d,X[,2])

    y <- exp(d*tau +  X %*% beta) +u
    y[y<=0] <- 0
    dat <- data.frame(y=y,
                      d=d,
                      X1=X[,2],
                      X2=X[,3],
                      Z=Z[,4])

    dat$ln.y <- log(y+any(y==0) )

    Zobs <- with(dat, cbind(1, X2, Z))
    Xobs <- with(dat, cbind(1, d, X2))
    tsls <- solve(t(Zobs) %*% Xobs) %*% t(Zobs) %*%dat$ln.y
    resid.tsls <- drop(dat$ln.y - Xobs %*% tsls)
    Vtsls <- solve(t(Zobs) %*% Xobs) %*% crossprod(Zobs*resid.tsls) %*%  solve(t(Xobs) %*% Zobs)
    ### canned version
    # tsls <- ivreg(ln.y~d+X2|X2+Z, data=dat)
    # step1 <- lm(d~X2+Z, data=dat, x=TRUE)
    # dat$v <- step1$resid
    # poi.ovb <- glm(y~d+X2,data=dat, family=quasipoisson())
    # poi.iv <- glm(y~d+X2+v, data=dat, family=quasipoisson(),x=TRUE)
    # V1 <- vcovHC(step1,type="HC0")
    # V2 <- vcovHC(poi.iv,type="HC0")
    # V1h <- vcov(step1)
    # V2h <- vcov(poi.iv)/summary(poi.iv)$dispersion
    
    #### by hand ####
    step1 <- solve(t(Zobs) %*% Zobs) %*% t(Zobs) %*% dat$d
    dat$v <- drop(dat$d - Zobs %*% step1)
    x0 <- rep(0,3)
    names(x0) <- c("const", "d", "x2")
    poi.ovb <- optim(x0, 
                     l.poisson, gr=r.poisson,
                     X=Xobs,
                     y=dat$y, 
                     method="BFGS",
                     hessian=TRUE)
    Vovb <- solve(poi.ovb$hessian) %*% 
      crossprod(r.poisson(poi.ovb$par, y=dat$y, X=Xobs, sum=FALSE)) %*%
      solve(poi.ovb$hessian)
    
    x0 <- rep(0,4)
    names(x0) <- c("const", "d", "x2", "v")
    X2 <- cbind(Xobs,dat$v)
    poi.iv <- optim(x0, 
                     l.poisson, gr=r.poisson,
                     X=X2,
                     y=dat$y, 
                     method="BFGS",
                     hessian=TRUE)
    V1h <- mean(dat$v^2) *solve(t(Zobs)%*% Zobs)
    V2h <- solve(poi.iv$hessian)
    V1 <- solve(t(Zobs)%*% Zobs) %*% (crossprod(Zobs*dat$v)) %*% solve(t(Zobs)%*% Zobs)
    V2 <- V2h %*% crossprod(r.poisson(poi.iv$par, y=dat$y, X=X2, sum=FALSE)) %*% V2h
    lhat <- drop(exp(X2 %*% poi.iv$par))
    C <- (t(Zobs)%*% diag(lhat) %*% X2) * poi.iv$par["v"]
    C[,4] <- C[,4]+ (lhat - dat$y) %*% Zobs
    C <- -C
    R <- t(r.normal(c(step1, mean(dat$v^2)),X=Zobs, y=dat$d, sum=FALSE)) %*%
      r.poisson(poi.iv$par,X=X2, y=dat$y, sum=FALSE)
    
    V2step <- V2 + V2h %*% (t(C) %*% V1 %*% C- t(R) %*% V1h %*% C - t(C) %*% V1h %*% R )%*%V2h
    
    out <- c(tsls[2],
             poi.ovb$par[2],
             poi.iv$par[2],
             sqrt(diag(Vtsls))[2],
             sqrt(diag(Vovb))[2],
             sqrt(diag(V2step))[2],
             sqrt(diag(V2))[2],
             cor(dat$X1, dat$d))
    
    if(any(abs(out) > 100)){
      output[b,] <- rep(NA, ncol(output))
    }else{
      output[b,] <- out
    }
    
    
  }
  results[[k]] <- output
}
toc <- proc.time()[3]
cat("Simulations complete in ", (toc-tic), " seconds\n")
save(list=c("conds", "results"), file="ps4_simulations.rdata")



### Analysis###
library(matrixStats)
library(ggplot2)
library(reshape2)
rm(list=ls())

rho <- seq(-2/3, 2/3, by=1/3)
tau <- .3
load("ps4_simulations.rdata")

means <- t(sapply(results, colMeans, na.rm=TRUE))
sds <- t(sapply(results, colSds, na.rm=TRUE))

colnames(means) <- c("TSLS", "Poisson.OVB", "Poisson.IV",
                     "TSLS.SE","Poisson.OVB.SE", "Poisson.IV.SE", "Poisson.IV.SE0",
                      "cor.X1.d")

obs.sds <- sds[,1:3]
colnames(obs.sds) <- c("TSLS.sd", "Poisson.OVB.sd", "Poisson.IV.sd")

plot.df <- cbind(conds, means[,1:7],obs.sds)
plot.df <- melt(plot.df,
                id.vars = c("gamma1","N"),
                variable.name = "label")

r <- factor(rho, levels=c(-2/3, -1/3, 0, 2/3, 1/3),
            labels=c("-2/3","-1/3", "0",  "2/3","1/3"))

plot.df$rho <- r
plot.df$Estimator <- c(rep(c("TSLS", "Poisson", "IV-Poisson"),each=nrow(conds)),
                       rep(c("TSLS",  "Poisson", "Murphy-Topel IV-Poisson","Uncorrected IV-Poisson"),each=nrow(conds)),
                       rep(c("TSLS", "Poisson", "IV-Poisson"),each=nrow(conds)))
plot.df$statistic <- c(rep("hat(tau)",3*nrow(conds)),
                       rep("se(hat(tau))",4*nrow(conds)),
                       rep("std",3*nrow(conds)))

bias.df <- plot.df[plot.df$statistic=="hat(tau)",]  
bias.df$value <- abs(bias.df$value-tau)



ggplot(bias.df)+
  geom_line(aes(x=N, y=value, color=Estimator,linetype=Estimator))+
  facet_wrap(rho~., ncol=3)+
  theme_bw(14)+
  ylab("Absolute bias")+
  xlab("Sample size")+
  ggtitle(expression("Bias in"~hat(tau)))+
  theme(legend.position = "bottom")



se.df <- plot.df[plot.df$statistic=="se(hat(tau))",]  
std.df <- plot.df[plot.df$statistic=="std",]
std.df <- std.df[c(1:60, 41:60),]
se.df$value <- (se.df$value-std.df$value)


ggplot(se.df)+
  geom_line(aes(x=N, y=value, color=Estimator, linetype=Estimator))+
  facet_wrap(rho~., ncol=3)+
  theme_bw(14)+
  ylab("Bias")+
  xlab("Sample size")+
  ggtitle(expression("Bias in standard errors"~hat(tau)))+
  theme(legend.position = "bottom")



cover <- function(X){
  colMeans(
    cbind("TSLS"=X[,1]+1.96*  X[,4] > tau & X[,1]-1.96*  X[,4] < tau,
          "Poisson"=X[,2]+1.96*  X[,5] > tau & X[,2]-1.96*  X[,5] < tau,
          "Poisson-IV (corrected)"=X[,3]+1.96*  X[,6] > tau & X[,2]-1.96*  X[,6] < tau,
          "Poisson-IV (uncorrected)"=X[,3]+1.96*  X[,7] > tau & X[,2]-1.96*  X[,7] < tau),
    na.rm=TRUE)
}

coverage <- cbind(conds, rho=r, t(sapply(results, cover)))
coverage <- melt(coverage,id.vars = c("gamma1","N", "rho"),
                            variable.name = "Estimator")

ggplot(coverage)+
  geom_hline(yintercept = .95, alpha=.5)+
  geom_line(aes(x=N, y=value, color=Estimator,linetype=Estimator))+
  facet_wrap(rho~., ncol=3)+
  theme_bw(14)+
  ylab("Coverage")+
  xlab("Sample size")+
  ggtitle(expression("Coverage of"~hat(tau)))+
  theme(legend.position = "bottom")


cbind(means[, "cor.X1.d"],rep(rho, 4))
cor(cbind(means[, "cor.X1.d"],rep(rho, 4)))
par(mfrow=c(1,2))
plot(means[, "cor.X1.d"]~rep(rho, 4),
     ylab="Observed correlation",
     xlab=expression(rho))
plot(means[, "cor.X1.d"]-rep(rho, 4), 
     ylab="Difference in oberserved and expected correlations", 
     xlab="Experiment")
abline(h=0)
