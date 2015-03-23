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
Lambdas <- c(0.5, 1, 2, 4)
alpha <- 0
beta <- 1

nRuns <- 1000
SIMULATIONS <- SUMMARIES <- vector(mode = "list", length = nRuns)

# iterations in Gibbs sampler
nchains <- 2
nsims <- 6000
burnin <- 1000
iters <- nsims-burnin
# monitored parameters
params <- c("alpha","beta","gamma2X","kappa2")

BIGOUT <- data.table(expand.grid(Run=factor(1:nRuns),lambda=Lambdas))[, sapply(params, function(x) as.double(rep(NA,nchains*iters)), simplify=FALSE), by="Run,lambda"]
setkey(BIGOUT, Run, lambda)
# BIGOUT[sim==1 & lambda ==1, identity(params):=dd] 


for(run in 1:nRuns){   
  SUMMARIES[[run]] <- vector(mode = "list", length = length(Lambdas))
  for(l in 1:length(Lambdas)){ 
    # simulated data
    gammaY0 <- sqrt(Lambdas[l])*gammaX0
    X<-Y<-array(NA,dim=c(N,J))
    for(i in 1:N){ 
      for(j in 1:J){
        X[i,j] <- rnorm(1,theta0[i],gammaX0)
        Y[i,j] <- rnorm(1,alpha+beta*theta0[i],gammaY0)	
      }
    }
    X <- stack(data.frame(t(X))) %>% setNames(c("x","group"))
    Y <- stack(data.frame(t(Y))) %>% setNames(c("y","group"))
    # Gibbs sampler
    OUT <- deming_gibbs2(X, Y, nchains=nchains, nsims=nsims, burnin=burnin, params=params)
    #BIGOUT[.(as.character(run),Lambdas[l]), names(OUT):=OUT]  
    BIGOUT[.(as.character(run),Lambdas[l]), identity(params):=OUT[,params,with=FALSE]]  
#     # outputs
#     SIMS <- lapply(setNames(params,params), 
#                    function(param) do.call(abind, lapply(OUT, 
#                                                          function(out) out[[param]])))
#     alpha.sims <- as.numeric(SIMS[["alpha"]])
#     beta.sims <- as.numeric(SIMS[["beta"]])
#     means <- sapply(SIMS[params], mean) 
#     quantiles <- sapply(SIMS[params], quantile, probs = c(0.025, 0.975)) 
#     posterior.summaries <- rbind(means,quantiles)
#     SIMULATIONS[[run]][[l]] <- list(alpha=alpha.sims,beta=beta.sims)
#     SUMMARIES[[run]][[l]] <- posterior.summaries
  }
}

saveRDS(BIGOUT, file="SIMULATIONS_gibbs2.rds")

long <- reshape2::melt(BIGOUT, id.vars=c("Run","lambda"), measure.vars=params, variable.name="parameter")
rm(BIGOUT)

summaries <- long[, list(mean=mean(value), lwr=quantile(value, 2.5/100), upr=quantile(value, 97.5/100)), 
     by="Run,lambda,parameter"][, length:=(upr-lwr)]
saveRDS(summaries, file="SUMMARIES_gibbs2.rds")

summaries <- readRDS("SUMMARIES_gibbs2.rds")
trueparams <- data.table(lambda=Lambdas)[, list(alpha=alpha, beta=beta, gamma2X=gammaX0^2, kappa2=lambda), by=lambda] %>% 
  reshape2::melt(id.vars="lambda", variable.name="parameter", factorsAsStrings = TRUE)
setkey(trueparams, parameter)

summaries <- merge(summaries, trueparams, by=c("lambda","parameter"))[, "lambda":=factor(lambda)]

intervals <- reshape2::melt(summaries, id.vars=c("Run","lambda","parameter"), measure.vars=c("lwr","upr"), 
                  variable.name="bound")
setkey(intervals, parameter)

ggplot(intervals[parameter!="kappa2"], aes(x=Run,y=value,color=bound)) + geom_point() + 
  facet_grid(parameter~lambda, scales="free_y") + 
  geom_hline(mapping=aes(yintercept=value), data=trueparams[parameter!="kappa2"], color="red", linetype="dashed")

covers <- summaries[, list(
  lwr=mean(lwr<value),
  upr=mean(upr>value),
  twosided=mean(lwr<value & upr>value)
  ), by="lambda,parameter"] %>% reshape2::melt(id.vars=c("lambda","parameter"), 
                                               variable.name="coverage", factorsAsStrings = FALSE) %>%
  merge(data.table(coverage=factor(c("lwr","upr","twosided")), nominal=c(.975,.975,.95)), 
        by="coverage") %>% 
  .[,"coverage":=relevel(coverage,"twosided")]


ggplot(covers, aes(x=lambda, y=value)) + geom_point() + 
  facet_grid(coverage~parameter, scales="free_y") + geom_hline(aes(yintercept=nominal), color="red", linetype="dashed")

ggplot(summaries, aes(x=lwr,y=upr, color=lambda)) + geom_point() +
  facet_grid(.~parameter)

means <-  recast(data.frame(summaries), lambda+Run~parameter, 
                 measure.var=c("mean"), id.var=c("lambda","Run","parameter"))
ggplot(means, aes(x=alpha,y=beta,color=lambda)) + geom_point()

# faire aussi un truc qui transforme la matrice theta - a ou plutôt faire ça dans gibbs() ?
# 
# head(as.data.table(t(OUT[[1]]$theta)))
# 
# theta <- t(OUT[[1]]$theta)
# Theta <- setNames(as.data.table(theta), sprintf(paste0("theta_%0", floor(log10(ncol(theta)))+1, "d"), 1:ncol(theta)))
# 
# x <- list(a=11, b=12, c=13)
# Map(function(x, i) paste(i, x), x, seq_along(x))
# 
# iterations <- 5:10 # tronquer au préalable ou pas ?.. 
# xx <- rbindlist(
#   Map(function(out, i) as.data.table(out)[, "Chain":=rep(i,.N)], OUT, seq_along(OUT))
#   )[, "Iteration":=seq_len(.N), by=Chain]
# 
# library(data.table)
# f <- function(stack=FALSE){
#   OUT <- list(A=list(a=1:2,b=3:4), B=list(a=5:6,b=7:8))
#   if(stack) OUT <- rbindlist(
#     Map(function(out, i) as.data.table(out)[, "Chain":=rep(i,.N)], OUT, seq_along(OUT))
#   )
#   return(OUT)
# }
# 
# 
# save(SIMULATIONS, SUMMARIES, file="SIMULATIONS_COVERAGE2.Rdata")
# 
# # faire cover sup et inf 
# 
# Covers <- array(NA,dim=c(2,l))
# rownames(Covers) <- c("alpha","beta")
# colnames(Covers) <- Lambdas
# for(l in 1:length(Lambdas)){
#   cover <- 0
#   for(num in 1:nsims){
#     bounds <- SUMMARIES[[num]][[l]][2:3,"alpha"]
#     cover <- cover + as.numeric(bounds[1]<0 & 0<bounds[2])
#   }
#   Covers["alpha",l] <- cover/nsims
#   cover <- 0
#   for(num in 1:nsims){
#     bounds <- SUMMARIES[[num]][[l]][2:3,"beta"]
#     cover <- cover + as.numeric(bounds[1]<1 & 1<bounds[2])
#   }
#   Covers["beta",l] <- cover/nsims
# }
# 
# Covers
# 
# 
# Alphas <- array(NA,dim=c(nsims,l)) -> Betas
# colnames(Alphas) <- colnames(Betas) <- Lambdas
# for(l in 1:length(Lambdas)){
#   for(num in 1:nsims){
#     Alphas[num,l] <- SUMMARIES[[num]][[l]]["means","alpha"]
#     Betas[num,l] <- SUMMARIES[[num]][[l]]["means","beta"]
#   }
# }
