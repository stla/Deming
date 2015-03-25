setwd("~/Github/Deming/scripts")
library(magrittr)

N <- 25  # number of groups
J <- 4  # number of repetitions within each group
theta0 <- 1:N  # true theta_i
gammaX0 <- sqrt(4); 
lambda <- 1
gammaY0 <- sqrt(lambda)*gammaX0
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
  return(list(X=X,Y=Y))
}

set.seed(666)
sims <- simulations()
X <- sims$X; Y <- sims$Y

lm(rowMeans(Y)~rowMeans(X)) %>% coefficients
MethComp::Deming(rowMeans(X), rowMeans(Y))[1:2]


alpha <- 0; beta <- 1
for(iter in 1:100){
  theta <- rowMeans(X) + beta/(lambda+beta^2)*(rowMeans(Y)-(alpha+beta*rowMeans(X)))
  ab <- lm(c(t(Y))~I(rep(theta, each=J))) %>% coefficients
  alpha <- ab[1]; beta <- ab[2]
}
ab
MethComp::Deming(rowMeans(X), rowMeans(Y))[1:2]

alpha <- 0; beta <- 1
for(iter in 1:10000){
  theta <- rowMeans(X) + beta/(1/rgamma(1,1)+beta^2)*(rowMeans(Y)-(alpha+beta*rowMeans(X)))
  ab <- lm(c(t(Y))~I(rep(theta, each=J))) %>% coefficients
  alpha <- ab[1]; beta <- ab[2]
}
ab
MethComp::Deming(rowMeans(X), rowMeans(Y))[1:2]


