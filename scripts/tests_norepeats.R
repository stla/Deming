setwd("~/Github/Deming/scripts")
library(ggplot2)

library(Deming)
#library(abind)
library(magrittr)
library(data.table)

set.seed(666) 
N<-25  # number of groups
J<-1  # number of repetitions within each group
theta0 <- 1:N  # true theta_i
gammaX0 <- sqrt(4); 
lambda <- 1
alpha <- 0
beta <- 1
gammaY0 <- sqrt(lambda)*gammaX0
X<-Y<-array(NA,dim=c(N,J))
for(i in 1:N){ 
  for(j in 1:J){
    X[i,j] <- rnorm(1,theta0[i],gammaX0)
    Y[i,j] <- rnorm(1,alpha+beta*theta0[i],gammaY0)  
  }
}
X <- stack(data.frame(t(X))) %>% setNames(c("x","group")) 
Y <- stack(data.frame(t(Y))) %>% setNames(c("y","group"))

deming.ci(X$x,Y$y,lambda=lambda)

OUT <- deming_gibbs2(X, Y, a=1000, b=1000, nchains=2, nsims=6000, burnin=1000)
plot(density(OUT$alpha))
mean(OUT$alpha)
mean(OUT$beta)
quantile(OUT$alpha, c(2.5,97.5)/100)
quantile(OUT$beta, c(2.5,97.5)/100)
confint(lm(Y$y~I(X$x))) # !!!! idem !!! => j'ai modifié l'étape (alpha,beta) en Deming.. estimates ok mais trop serré

plot(OUT$gamma2X)
plot(OUT$kappa2)