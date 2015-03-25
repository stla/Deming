setwd("~/Github/Deming/scripts")
#library(Deming)
library(rjags)
library(magrittr)
library(data.table)
library(ggplot2)

## define WinBUGS model
WBM <- "gibbs1.txt" # model file
jmodel <- function(){
  for(i in 1:N){
    for(j in 1:J){
      X[i,j] ~ dnorm(theta[i], inv.gammaX2)
      Y[i,j] ~ dnorm(nu[i], inv.gammaY2)
    }
    theta[i] ~ dnorm(0, 0.00001) 
    nu[i] <- alpha+beta*theta[i]
  }
  alpha ~ dnorm(0,0.001)
  beta ~ dnorm(0,0.001)
  #
  inv.gammaX2 ~ dgamma(0.01,0.01)
  inv.gammaY2 ~ dgamma(0.01,0.01)
  gammaX2 <- 1/inv.gammaX2
  gammaY2 <- 1/inv.gammaY2
  lambda <- gammaY2/gammaX2
}
writeLines(paste("model", paste(body(jmodel), collapse="\n"), "}"), WBM)

### data
set.seed(666) 
N <- 25  # number of groups
J <- 4  # number of repetitions within each group
theta0 <- 1:N  # true theta_i
gammaX0 <- sqrt(4); 
gammaY0 <- gammaX0
alpha <- 0
beta <- 1
simulations <- function(){
  X <- Y <- array(NA,dim=c(N,J))
  for(i in 1:N){ 
    for(j in 1:J){
      X[i,j] <- rnorm(1,theta0[i],gammaX0)
      Y[i,j] <- rnorm(1,alpha+beta*theta0[i],gammaY0)  
    }
  }
#   Xs <- stack(data.frame(t(X))) %>% setNames(c("x","group"))
#   Ys <- stack(data.frame(t(Y))) %>% setNames(c("y","group"))
  return(list(X=X,Y=Y))
}

params <- c("alpha", "beta")
nRuns <- 3000
OUT <-  data.table(expand.grid(Run=factor(1:nRuns),method=c("JAGS", "least-squares", "Deming")))[, sapply(params, function(x) as.double(NA), simplify=FALSE), by="Run,method"]
setkey(OUT, Run, method)

for(run in 1:nRuns){
  sims <- simulations()
  X <- sims$X; Y <- sims$Y
  # jags sampler
  data <- list(X=X, Y=Y, N=N, J=J)
  inits1 <- list(alpha=0, beta=1, inv.gammaX2=1, inv.gammaY2=1, theta=rowMeans(X))
  inits2 <- lapply(inits1, function(x) x*rnorm(length(x),1,.01))
  inits <- list(inits1,inits2)
  
  jags <- jags.model(WBM,
                     data = data, 
                     inits=inits, 
                     n.chains = 2,
                     n.adapt = 1000,
                     quiet=TRUE)
  update(jags, 5000, progress.bar="none")
  samples <- coda.samples(jags, c("alpha", "beta", "lambda"), 10000, progress.bar="none")
  
#   summary(samples)$statistics
#   lm(rowMeans(Y)~I(rowMeans(X))) %>% coefficients
#   MethComp::Deming(rowMeans(X), rowMeans(Y))
  
  OUT[.(as.character(run),"JAGS"), identity(params):=data.table(t(summary(samples)$statistics[1:2,1]))]  
  OUT[.(as.character(run),"least-squares"), identity(params):=data.table(t(lm(rowMeans(Y)~I(rowMeans(X))) %>% coefficients))]  
  OUT[.(as.character(run),"Deming"), identity(params):=data.table(t(MethComp::Deming(rowMeans(X), rowMeans(Y))[1:2]))]  
  
}

saveRDS(OUT, "JAGS_COMPARISONS.rds")

ggplot(OUT, aes(x=Run, y=alpha, color=method)) + geom_point()
setkey(OUT, method)

# alpha
par(pty="s")
layout(t(1:2))
plot(OUT[.("JAGS")]$alpha, OUT[.("least-squares")]$alpha, 
     asp=1, pch=19, col="blue",  
     xlab="JAGS", ylab="least-squares",
     main="estimates of alpha")
abline(0,1)
plot(OUT[.("JAGS")]$alpha, OUT[.("Deming")]$alpha, 
     asp=1, pch=19, col="red",  
     xlab="JAGS", ylab="Deming")
abline(0,1)
# beta
layout(t(1:2))
par(pty="s")
plot(OUT[.("JAGS")]$beta, OUT[.("least-squares")]$beta, 
     asp=1, pch=19, col="blue",  
     xlab="JAGS", ylab="least-squares",
     main="estimates of beta")
abline(0,1)
plot(OUT[.("JAGS")]$beta, OUT[.("Deming")]$beta, 
     asp=1, pch=19, col="red",  
     xlab="JAGS", ylab="Deming")
abline(0,1)

# faire plutôt ggpairs et density pour voir meilleur - ggpairs chiant pour abline(0,1)

ggdata <- reshape2::melt(OUT, measure.vars=c("alpha","beta"))
gg1 <- ggplot(OUT, aes(x=alpha, color=method)) + geom_density() +
  geom_vline(xintercept=0) + theme(legend.position="top")
gg2 <- ggplot(OUT, aes(x=beta, color=method)) + geom_density() +
  geom_vline(xintercept=1) + theme(legend.position="top")   
gridExtra::grid.arrange(gg1,gg2,ncol=2)


ggplot(ggdata, aes(x=value, color=method)) + geom_density() +
  geom_vline(xintercept=0) + 
  facet_grid(variable~., scales="free_y", "free_x")

MethComp::Deming(Xs$x, Ys$y) # pas de sens car on peut permuter

deming.ci(rowMeans(X), rowMeans(Y), lambda=1)




# my Gibbs sampler
OUT <- deming_gibbs1(Xs, Ys, nchains=2, nsims=12000, burnin=2000)
mysims <- reshape2::melt(OUT, id.vars=c("Chain","Iteration"), variable.name="parameter")
mysims[!stringr::str_detect(parameter,"_"), list(mean=mean(value), median=median(value), lwr=quantile(value, .025), upr=quantile(value, .975)), by=parameter]
( lmfit <- lm(rowMeans(Y)~I(rowMeans(X))) )
confint(lmfit)
deming.ci(rowMeans(X), rowMeans(Y), lambda=1)
MethComp::Deming(rowMeans(X), rowMeans(Y))
