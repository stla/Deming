library(Deming)
library(abind)
library(magrittr)


set.seed(666) 
N<-25  # number of groups
J<-4  # number of repetitions within each group
theta0 <- 1:N  # true theta_i
gammaX0 <- sqrt(4); 
Lambdas <- c(0.5, 1, 2, 4)

nsims <- 1000
SIMULATIONS <- SUMMARIES <- vector(mode = "list", length = nsims)


for(num in 1:nsims){  
  SUMMARIES[[num]] <- vector(mode = "list", length = length(Lambdas))
  for(l in 1:length(Lambdas)){ 
    # simulated data
    lambda <- Lambdas[l]
    gammaY0 <- sqrt(lambda)*gammaX0
    X<-Y<-array(NA,dim=c(N,J))
    for(i in 1:N){ 
      for(j in 1:J){
        X[i,j] <- rnorm(1,theta0[i],gammaX0)
        Y[i,j] <- rnorm(1,theta0[i],gammaY0)	
      }
    }
    X <- stack(data.frame(t(X))) %>% setNames(c("x","group"))
    Y <- stack(data.frame(t(Y))) %>% setNames(c("y","group"))
    # monitored parameters
    params <- c("alpha", "beta", "gamma2X", "kappa2")
    # Gibbs sampler
    OUT <- deming_gibbs2(X, Y, nsims=4000, burnin=500, params=params)
    # outputs
    SIMS <- lapply(setNames(params,params), 
                   function(param) do.call(abind, lapply(OUT, 
                                                         function(out) out[[param]])))
    alpha.sims <- as.numeric(SIMS[["alpha"]])
    beta.sims <- as.numeric(SIMS[["beta"]])
    means <- sapply(SIMS[params], mean) 
    quantiles <- sapply(SIMS[params], quantile, probs = c(0.025, 0.975)) 
    posterior.summaries <- rbind(means,quantiles)
    SIMULATIONS[[num]][[l]] <- list(alpha=alpha.sims,beta=beta.sims)
    SUMMARIES[[num]][[l]] <- posterior.summaries
  }
}

# faire aussi un truc qui transforme la matrice theta - a ou plutôt faire ça dans gibbs() ?

head(as.data.table(t(OUT[[1]]$theta)))

theta <- t(OUT[[1]]$theta)
Theta <- setNames(as.data.table(theta), sprintf(paste0("theta_%0", floor(log10(ncol(theta)))+1, "d"), 1:ncol(theta)))

x <- list(a=11, b=12, c=13)
Map(function(x, i) paste(i, x), x, seq_along(x))

iterations <- 5:10 # tronquer au préalable ou pas ?.. 
xx <- rbindlist(
  Map(function(out, i) as.data.table(out)[, "Chain":=rep(i,.N)], OUT, seq_along(OUT))
  )[, "Iteration":=seq_len(.N), by=Chain]

save(SIMULATIONS, SUMMARIES, file="SIMULATIONS_COVERAGE2.Rdata")

# faire cover sup et inf 

Covers <- array(NA,dim=c(2,l))
rownames(Covers) <- c("alpha","beta")
colnames(Covers) <- Lambdas
for(l in 1:length(Lambdas)){
  cover <- 0
  for(num in 1:nsims){
    bounds <- SUMMARIES[[num]][[l]][2:3,"alpha"]
    cover <- cover + as.numeric(bounds[1]<0 & 0<bounds[2])
  }
  Covers["alpha",l] <- cover/nsims
  cover <- 0
  for(num in 1:nsims){
    bounds <- SUMMARIES[[num]][[l]][2:3,"beta"]
    cover <- cover + as.numeric(bounds[1]<1 & 1<bounds[2])
  }
  Covers["beta",l] <- cover/nsims
}

Covers


Alphas <- array(NA,dim=c(nsims,l)) -> Betas
colnames(Alphas) <- colnames(Betas) <- Lambdas
for(l in 1:length(Lambdas)){
  for(num in 1:nsims){
    Alphas[num,l] <- SUMMARIES[[num]][[l]]["means","alpha"]
    Betas[num,l] <- SUMMARIES[[num]][[l]]["means","beta"]
  }
}
