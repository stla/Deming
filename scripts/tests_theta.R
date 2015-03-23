setwd("~/Github/Deming/scripts")
library(ggplot2)

library(Deming)
#library(abind)
library(magrittr)
library(data.table)

set.seed(666) 
N<-25  # number of groups
J<-4  # number of repetitions within each group
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
X$group <- reorder(X$group, 1:nrow(X), FUN=identity)
Y$group <- reorder(Y$group, 1:nrow(Y), FUN=identity)

OUT <- deming_gibbs1(X, Y, nchains=nchains, nsims=6000, burnin=1000)

gammaX0/sqrt(J)
sd(OUT$theta_01)

round(apply(as.matrix(OUT), 2, mean),2)
round(apply(as.matrix(OUT), 2, sd),2)

plot(density(OUT$theta_03)) # problème ordre des theta ! => ordonner dans X et Y 

