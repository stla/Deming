library(Deming)
library(abind)
library(magrittr)


set.seed(1975) #set.seed(2009)
N<-25  # number of measurments (pour chaque appareil)
J<-4  # repetitions
theta0 <- 1:N  # true theta_i
gammaX0 <- sqrt(4); 
Lambdas <- c(0.5, 1, 2, 4); l <- length(Lambdas)

nsims <- 200
SIMULATIONS <- SUMMARIES <- vector(mode = "list", length = nsims)


for(num in 1:nsims){  
  SUMMARIES[[num]] <- vector(mode = "list", length = l)
  for(ll in 1:l){ 
    lambda <- Lambdas[ll]
    gammaY0 <- sqrt(lambda)*gammaX0
    X<-Y<-array(NA,dim=c(N,J))
    for(i in 1:N){ # simulated data
      for(j in 1:J){
        X[i,j] <- rnorm(1,theta0[i],gammaX0)
        Y[i,j] <- rnorm(1,theta0[i],gammaY0)	
      }
    }
    X <- stack(data.frame(t(X))) %>% setNames(c("x","group"))
    Y <- stack(data.frame(t(Y))) %>% setNames(c("y","group"))
    params <- c("alpha", "beta", "gamma2X", "gamma2Y")
    OUT <- deming_gibbs1(X, Y, nsims=10000, burnin=1000, params=params)
    SIMS <- lapply(setNames(params,params), function(param) do.call(abind, lapply(OUT, function(out) out[[param]])))
    alpha.sims <- as.numeric(SIMS[["alpha"]])
    beta.sims <- as.numeric(SIMS[["beta"]])
    means <- sapply(SIMS[params], mean) 
    quantiles <- sapply(SIMS[params], quantile, probs = c(0.025, 0.975)) 
    posterior.summaries <- round(rbind(means,quantiles),2) 
    SIMULATIONS[[num]][[ll]] <- list(alpha=alpha.sims,beta=beta.sims)
    SUMMARIES[[num]][[ll]] <- posterior.summaries
  }
}


save(SIMULATIONS, SUMMARIES, file="SIMULATIONS_COVERAGE.Rdata")

# faire cover sup et inf 

Covers <- array(NA,dim=c(2,l))
rownames(Covers) <- c("alpha","beta")
colnames(Covers) <- Lambdas
for(ll in 1:l){
  cover <- 0
  for(num in 1:nsims){
    bounds <- SUMMARIES[[num]][[ll]][2:3,"alpha"]
    cover <- cover + as.numeric(bounds[1]<0 & 0<bounds[2])
  }
  Covers["alpha",ll] <- cover/nsims
  cover <- 0
  for(num in 1:nsims){
    bounds <- SUMMARIES[[num]][[ll]][2:3,"beta"]
    cover <- cover + as.numeric(bounds[1]<1 & 1<bounds[2])
  }
  Covers["beta",ll] <- cover/nsims
}

Covers


Alphas <- array(NA,dim=c(nsims,l)) -> Betas
colnames(Alphas) <- colnames(Betas) <- Lambdas
for(ll in 1:l){
  for(num in 1:nsims){
    Alphas[num,ll] <- SUMMARIES[[num]][[ll]]["means","alpha"]
    Betas[num,ll] <- SUMMARIES[[num]][[ll]]["means","beta"]
  }
}
