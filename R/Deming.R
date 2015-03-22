#' Burnin and thin
#' 
#' @export
Extract <- function(sims,burnin,thin){ ## extract from simulations
  if(is.null(dim(sims))){
    l <- length(sims)
    sims <- sims[c(1:l)%%thin==0 & c(1:l)>burnin]
  }else{
    l <- dim(sims)[2]
    sims <- sims[,(c(1:l)%%thin==0 & c(1:l)>burnin)]
  }
}


#' Gibbs sampler 1 for Deming regression
#' 
#' Gibbs sampler with independent priors on variances
#' 
#'@import dplyr
#'@import data.table
#'@export 
deming_gibbs1 <- function(X, Y, nsims=5000, nchains=2, burnin=1000, thin=1, params="all", stack=TRUE, m=rep(0,nlevels(X$group)), tau2=1e4, aX=.1, bX=.001, aY=.1, bY=.001, alpha0=0, beta0=0, B0=diag(.0001, 2)){
  ### checks ###
  allparams <- c("alpha", "beta", "gamma2X", "gamma2Y", "theta")
  if(identical(params,"all")){
    params <- allparams
  }else{
    if(any(!is.element(params, allparams))) stop("non valid parameters")
  }
  if(B0[1,2]!=B0[2,1] || det(B0)<=0) stop("B0 is not symmetric positive")
  ##### data 
  X <- arrange(X, group); Y <- arrange(Y, group)
  x <- X$x; y <- Y$y
  if(nlevels(X$group) != nlevels(Y$group)) stop("")
  N <- nlevels(X$group)
  ####### storing simulations  ###########
  OUT <- vector("list",length=nchains)
  # 
  shapeX <- nrow(X)/2+aX
  shapeY <- nrow(Y)/2+aY
  sizesX <- X %>% group_by(group) %>% summarise(sizes=n()) %>% .$sizes
  sizesY <- Y %>% group_by(group) %>% summarise(sizes=n()) %>% .$sizes
  sumX <- X %>% group_by(group) %>% summarise(sum=sum(x)) %>% .$sum
  sumY <- Y %>% group_by(group) %>% summarise(sum=sum(y)) %>% .$sum
  iterations <- subset(seq_len(nsims), seq_len(nsims)%%thin==0 & seq_len(nsims)>burnin)
  ###### run Gibbs' sampler
  for(chain in 1:nchains){ 
    gamma2X.sims <- gamma2Y.sims <- alpha.sims <- beta.sims <- rep(NA,nsims)
    theta.sims   <-  array(NA,dim=c(N,nsims))
    ### initial values - generates random initial values ###
    # theta
    theta <- Xmeans <- X %>% group_by(group) %>% summarise(mean=mean(x)) %>% .$mean
    # alpha and beta
    Ymeans <- Y %>% group_by(group) %>% summarise(mean=mean(y)) %>% .$mean
    ab <- coef(lm(Ymeans*rnorm(N,1,.01) ~ I(Xmeans*rnorm(N,1,.01))))
    alpha <- ab[1]; beta <- ab[2]
    # variances 
    gamma2X <- X %>% group_by(group) %>% mutate(means=mean(x), resid=(x-means)^2) %>% 
      .$resid %>% sum %>% divide_by(length(x)-N) %>% multiply_by(rnorm(1,1,.01))
    gamma2Y <- Y %>% group_by(group) %>% mutate(means=mean(y), resid=(y-means)^2) %>% 
      .$resid %>% sum %>% divide_by(length(y)-N) %>% multiply_by(rnorm(1,1,.01))
    ### simulations ##
    rgammaX <- rgamma(nsims, shapeX, 1)
    rgammaY <- rgamma(nsims, shapeY, 1)
    rmnormAB <- matrix(rnorm(2*nsims), nrow=2)
    rmnormTheta <- matrix(rnorm(N*nsims), nrow=N)
    for(sim in 1:nsims){
      thetaX <- rep(theta, times=sizesX)
      thetaY <- rep(theta, times=sizesY)
      # draw gamma2X
      gamma2X.sims[sim] <- gamma2X <- (crossprod((x-thetaX))[1,]/2+bX)/rgammaX[sim]
      # draw gamma2Y
      gamma2Y.sims[sim] <- gamma2Y <- (crossprod((y-(alpha+beta*thetaY)))[1,]/2+bY)/rgammaY[sim]
      # draw alpha beta
      XX <- cbind(rep(1,length(y)), thetaY)
      # à améliorer : il y a une formule de cholesky pour les 2x2 - voir wiki Wishart
      M <- gamma2Y*B0+crossprod(XX)
      sqrtM <- sqrt(M)
      rho <- M[1,2]/prod(diag(sqrtM))
      cholBn <- cbind(c(sqrtM[1,1], 0), c(rho*sqrtM[2,2], sqrt(1-rho^2)*sqrtM[2,2]))
      inv.Bn <- chol2inv(cholBn)
      #       inv.Bn <-  chol2inv(chol(gamma2Y*B0+crossprod(XX)))
      Mean <- inv.Bn%*%(gamma2Y*B0%*%c(alpha0,beta0)+crossprod(XX,y))
      #S <- chol(gamma2Y*inv.Bn)
      S <- sqrt(gamma2Y)*t(backsolve(cholBn, diag(2)))
      M <- Mean + crossprod(S,rmnormAB[,sim])
      alpha.sims[sim] <- alpha <- M[1]
      beta.sims[sim] <- beta <- M[2]
      # draw theta_i
      mmean <- (gamma2Y*m + beta^2*tau2*(sumY-sizesY*alpha)/beta)/(gamma2Y + beta^2*sizesY*tau2)
      vvar <- (gamma2Y*tau2)/(gamma2Y + beta^2*sizesY*tau2) 
      mmean <- (gamma2X*mmean + vvar*sumX)/(gamma2X + sizesX*vvar)
      vvar <- (gamma2X*vvar)/(gamma2X + sizesX*vvar) 
      theta.sims[,sim] <- theta <- mmean + sqrt(vvar)*rmnormTheta[,sim]  
    }
    alpha.sims <- alpha.sims[iterations]
    beta.sims <- beta.sims[iterations]
    gamma2X.sims <- gamma2X.sims[iterations]
    gamma2Y.sims <- gamma2Y.sims[iterations]
    theta.sims <- theta.sims[,iterations]
    if(stack && "theta"%in%params) theta.sims <- setNames(as.data.table(t(theta.sims)), 
                                     sprintf(paste0("theta_%0", floor(log10(N)+1), "d"), 1:N))
    OUT[[chain]] <- eval(parse(text=sprintf("list(%s)", paste0(params, sprintf("=%s.sims", params), collapse=","))))
  } # end simulated chain
  #OUT <- lapply(OUT, function(out) lapply(out, function(sims) Extract(sims,burnin,thin)))
  if(stack){
    if("theta"%in%params){
      Theta <- rbindlist(
        Map(function(out, i) out[["theta"]], OUT, seq_along(OUT))
      )
    }
    OUT <- rbindlist(
      Map(function(out, i) as.data.table(out[subset(params, params!="theta")])[, "Chain":=rep(i,.N)], OUT, seq_along(OUT))
    )[, "Iteration":=iterations, by=Chain]
    if("theta"%in%params) OUT <- cbind(OUT,Theta)
  }
  return(OUT)
}



#' Gibbs sampler 2 for Deming regression
#' 
#' Gibbs sampler with independent priors on one variance and the variances ratio
#' 
#'@import dplyr
#'@export 
deming_gibbs2 <- function(X, Y, nsims=5000, nchains=2, burnin=1000, thin=1, params="all", m=rep(0,nlevels(X$group)), tau2=1e4, aX=.1, bX=.001, a=.1, b=.001, alpha0=0, beta0=0, B0=diag(.0001, 2)){
  ### checks ###
  allparams <- c("alpha", "beta", "gamma2X", "kappa2", "theta")
  if(identical(params,"all")){
    params <- allparams
  }else{
    if(any(!is.element(params, allparams))) stop("non valid parameters")
  }
  if(B0[1,2]!=B0[2,1] || det(B0)<=0) stop("B0 is not symmetric positive")
  ##### data 
  X <- arrange(X, group); Y <- arrange(Y, group)
  x <- X$x; y <- Y$y
  if(nlevels(X$group) != nlevels(Y$group)) stop("")
  N <- nlevels(X$group)
  ####### storing simulations  ###########
  gamma2X.sims <- kappa2.sims <- alpha.sims <- beta.sims <- rep(NA,nsims)
  theta.sims   <-  array(NA,dim=c(N,nsims))
  OUT <- vector("list",length=nchains)
  # 
  shapeX <- nrow(X)/2+nrow(Y)/2+aX
  shapeYX <- nrow(Y)/2+a
  sizesX <- X %>% group_by(group) %>% summarise(sizes=n()) %>% .$sizes
  sizesY <- Y %>% group_by(group) %>% summarise(sizes=n()) %>% .$sizes
  sumX <- X %>% group_by(group) %>% summarise(sum=sum(x)) %>% .$sum
  sumY <- Y %>% group_by(group) %>% summarise(sum=sum(y)) %>% .$sum
  
  ###### run Gibbs' sampler
  for(chain in 1:nchains){ 
    ### initial values - generates random initial values ###
    # theta
    theta <- Xmeans <- X %>% group_by(group) %>% summarise(mean=mean(x)) %>% .$mean
    # alpha and beta
    Ymeans <- Y %>% group_by(group) %>% summarise(mean=mean(y)) %>% .$mean
    ab <- coef(lm(Ymeans*rnorm(N,1,.01) ~ I(Xmeans*rnorm(N,1,.01))))
    alpha <- ab[1]; beta <- ab[2]
    # variances 
    gamma2X <- X %>% group_by(group) %>% mutate(means=mean(x), resid=(x-means)^2) %>% 
      .$resid %>% sum %>% divide_by(length(x)-N) %>% multiply_by(rnorm(1,1,.01))
    kappa2 <- Y %>% group_by(group) %>% mutate(means=mean(y), resid=(y-means)^2) %>% 
      .$resid %>% sum %>% divide_by(length(y)-N) %>% multiply_by(rnorm(1,1,.01)) %>%
      divide_by(gamma2X)
    ### simulations ##
    rgammaX <- rgamma(nsims, shapeX, 1)
    rgammaYX <- rgamma(nsims, shapeYX, 1)
    rmnormAB <- matrix(rnorm(2*nsims), nrow=2)
    rmnormTheta <- matrix(rnorm(N*nsims), nrow=N)
    for(sim in 1:nsims){
      thetaX <- rep(theta, times=sizesX)
      thetaY <- rep(theta, times=sizesY)
      # draw gamma2X
      gamma2X.sims[sim] <- gamma2X <- 
        ((crossprod(x-thetaX) + crossprod(y-(alpha+beta*thetaY))/kappa2)[1,]/2+bX)/rgammaX[sim]
      # draw kappa2
      kappa2.sims[sim] <- kappa2 <- (crossprod(y-(alpha+beta*thetaY))[1,]/gamma2X/2+b)/rgammaYX[sim]
      gamma2Y <- kappa2*gamma2X
      # draw alpha beta
      XX <- cbind(rep(1,length(y)), thetaY)
      inv.Bn <-  chol2inv(chol(gamma2Y*B0+crossprod(XX)))
      Mean <- inv.Bn%*%(gamma2Y*B0%*%c(alpha0,beta0)+crossprod(XX,y))
      S <- chol(gamma2Y*inv.Bn)
      M <- Mean + crossprod(S,rmnormAB[,sim])
      alpha.sims[sim] <- alpha <- M[1]
      beta.sims[sim] <- beta <- M[2]
      # draw theta_i
#       mmean <- (gamma2Y/beta^2) %>% { (.*m + tau2*(sumY-sizesY*alpha)/beta)/(. + sizesY*tau2) }
#       vvar <-  (gamma2Y/beta^2) %>% { (.*tau2)/(. + sizesY*tau2) }
#       mmean <- gamma2X %>% { (.*mmean + vvar*sumX)/(. + sizesX*vvar) }
#       vvar <- gamma2X %>% { (.*vvar)/(. + sizesX*vvar) }
      mmean <- (gamma2Y*m + beta^2*tau2*(sumY-sizesY*alpha)/beta)/(gamma2Y + beta^2*sizesY*tau2)
      vvar <- (gamma2Y*tau2)/(gamma2Y + beta^2*sizesY*tau2) 
      mmean <- (gamma2X*mmean + vvar*sumX)/(gamma2X + sizesX*vvar)
      vvar <- (gamma2X*vvar)/(gamma2X + sizesX*vvar) 
      theta.sims[,sim] <- theta <- mmean + sqrt(vvar)*rmnormTheta[,sim]  
    }    
    OUT[[chain]] <- eval(parse(text=sprintf("list(%s)", paste0(params, sprintf("=%s.sims", params), collapse=","))))
  } # end simulated chain
  OUT <- lapply(OUT, function(out) lapply(out, function(sims) Extract(sims,burnin,thin)))
  return(OUT)
}


