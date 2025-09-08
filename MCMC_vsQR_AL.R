#' Bayesian variable selection in quantile regression with random effects - ALD
#'
#' This function estimates a Bayesian variable selection in quantile 
#' regression with random effects, where the response variable is assumed to 
#' follow an asymmetric Laplace distribution.
#'
#' @param y n-dimensional vector of responses.
#' @param x Matrix - nx(p+1) - of predictors (include the intercept).
#' @param idobs n-dimensional vector that indicates groups (categorical numerical).
#' @param tau Quantile of interest (value between 0 and 1).
#' @param n_mcmc Number of iterations.
#' @param burnin_mcmc Number of initial iterations to be discarded.
#' @param thin_mcmc Thinning parameter.
#' 
#' @return A list with the chains of all parameters of interest.

## Packages
require(mvtnorm); require(GIGrvg)

## MCMC
vsQR_re_al <- function(y,x,idobs,tau,n_mcmc,burnin_mcmc,thin_mcmc){
  n <- length(y)
  numcov <- ncol(x)
  n_idobs <- length(unique(idobs))
  resultado <- list()
  # Create auxiliary objects
  m_alpha <- matrix(0,n,n_idobs)
  for(l in 1:n){
    m_alpha[l,as.numeric(base$meso[l])] <- 1 
  }
  beta  <- matrix(0, n_mcmc, numcov)
  sigma <- matrix(NA, n_mcmc, 1)
  alpha <- matrix(NA, n_mcmc, n_idobs)
  psi2  <- matrix(NA, n_mcmc, 1)
  malpha<- matrix(NA, n_mcmc, 1)
  delta <- matrix(NA, n_mcmc, 1)
  PI    <- matrix(NA, n_mcmc,numcov)
  v     <- matrix(NA,n_mcmc,n)
  PROB  <- matrix(NA, n_mcmc,numcov)
  # Set the initial values
  beta[1,]   <- rnorm(numcov,0,1)
  sigma[1,1] <- 1
  v[1,]      <- rgamma(n,2,sigma[1,1])
  alpha[1,]  <- rep(0,n_idobs)
  psi2[1,1]  <- 1
  malpha[1,1]<- 0.5
  delta[1,]  <- c(.5)
  PI[1,]<- rep(1,numcov)
  # Auxiliary constants
  tau2  <- 2 / (tau * (1 - tau))
  theta <- (1 - 2 * tau) / (tau * (1 - tau))
  # Prior (beta)
  b <- rep(0,numcov)
  B <- diag(rep(100,numcov))
  # MCMC
  for(k in 2:n_mcmc){
    alpha_aux  <- m_alpha%*%alpha[k-1,]
    #Variable Selection
    delta[k,1] <- atualizarDELTA(1,1,PI[k-1,])
    PI.aux     <- atualizarPI(numcov,PI[k-1,],y,alpha_aux,x,beta[k-1,],v[k-1,],theta,tau2,sigma[k-1,1],delta[k,1],B)
    PI[k,]     <- PI.aux[1,]
    PROB[k,]   <- PI.aux[2,]
    n_select   <- which(PI[k,]==1)
    beta[k,n_select] <- atualizarBETA(b[n_select],B[n_select,n_select],alpha_aux,y,x[,n_select],sigma[k-1,1],tau2,theta,v[k-1,])
    ######################
    v[k,]      <- atualizarV(y,alpha_aux,x,beta[k,],tau2,theta,sigma[k-1,1],n)
    sigma[k,1] <- atualizarSIGMA(0.01,0.01,alpha_aux,y,x,beta[k,],tau2,theta,v[k,],n)
    alpha[k,]  <- atualizarALPHA(y,x,beta[k,],theta,tau2,v[k,],sigma[k,1],malpha[k-1,1],psi2[k-1,1],t(m_alpha),n_idobs)
    malpha[k,1]<- atualizarMALPHA(0,1,alpha[k,],psi2[k-1,1],n_idobs)
    psi2[k,1]  <- atualizarPSI2(0.01,0.01,alpha[k,],malpha[k,1],n_idobs)
  }
  resultado[['beta']]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['sigma']] <- sigma[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['alpha']] <- alpha[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['pi']]    <- PI[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['qsi']]   <- PROB[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['malpha']]<- malpha[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['psi2']]  <- psi2[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  
  return(resultado)
}


#######################################################################
## Auxiliary functions: Sampling from full conditional distributions ##
#######################################################################

### Full conditional for beta
atualizarBETA<-function(b,B,alpha,y,x,sigma,tau2,theta,v){
  B.inv  <- chol2inv(chol(B))
  covar  <- chol2inv(chol(B.inv+(1/(tau2*sigma)*(t(x/v)%*%x))))
  media  <- covar%*%((B.inv%*%b)+(1/(tau2*sigma))*(t(x/v)%*%(y-alpha-theta*v)))
  beta   <- rmvnorm(1,media,covar)
  return(beta)
}

### Full conditional for sigma
atualizarSIGMA<-function(c,C,alpha,y,x,beta,tau2,theta,v,N){
  alpha1 <- c + 1.5*N
  beta1  <- C + sum(v) + (t((y-alpha-x%*%beta-theta*v)/v)%*%(y-alpha-x%*%beta-theta*v))/(2*tau2)
  sigma  <- 1/rgamma(1, alpha1, beta1)
  return(sigma)
}

### Full conditional for the latent variable
atualizarV<-function(y,alpha,x,beta,tau2,theta,sigma,N){
  p1 <- 0.5
  p2 <- (y-alpha-x%*%beta)^2/(tau2*sigma)
  p3 <- 2/sigma + theta^2/(tau2*sigma)
  v  <- NULL
  for(i in 1:N){
    v[i]  <- rgig(1, chi=p2[i], psi=p3, lambda=p1)
  }
  return(v)
}

### Full conditional for the random effects
atualizarALPHA<-function(y,x,beta,A,B,v,sigma,alpha_mean,psi2,m_aux,N){
  z     <- (y-x%*%beta-A*v)/v
  sd    <- sqrt((psi2*B*sigma)/(B*sigma+psi2*(m_aux%*%(1/v))))
  media <- (B*sigma*alpha_mean + psi2*(m_aux%*%z[,1])) / (B*sigma+psi2*(m_aux%*%(1/v)))
  alpha <- rnorm(N,media,sd)
  return(alpha)
}

### Full conditional for the hierarchical prior
atualizarPSI2<-function(c2,C2,alpha,alpha_mean,N){
  alpha1 <- c2 + 0.5*N
  beta1  <- C2 + 0.5*sum((alpha-alpha_mean)^2) 
  psi2   <- 1/rgamma(1, alpha1, beta1)
  return(psi2)
}
atualizarMALPHA<-function(b2,B2,alpha,psi2,N){
  sd     <- sqrt((psi2*B2)/(N*B2+psi2))
  media  <- (B2*sum(alpha) + psi2*b2) / (N*B2+psi2)
  malpha <- rnorm(1, media, sd)
  return(malpha)
}

### Variable Selection
atualizarDELTA<-function(a1,b1,vetorPI){
  delta<-rbeta(1, a1 + sum(vetorPI), b1 + sum(1 - vetorPI))
  return(delta)
}
atualizarPI<-function(numcov,PI,y,alpha,covariaveis,beta,v,theta,tau2,sigma,delta,Bjj){
  chance <- NULL
  for (j in sample(1:numcov)) {
    
    # get the betas for which beta_j is zero
    PI_aux<- PI
    PI_aux[j] <- 0
    beta_j_aux <- c(beta)*PI_aux
    
    # Compute the z variables and the conditional variance
    w  <- y-alpha-x%*%beta_j_aux-theta*v
    xj <- as.vector(covariaveis[,j])
    
    # Compute the chance parameter of the conditional posterior of pi_j (Bernoulli)
    l0 <- log(1-delta)
    l1 <- log(delta) - 0.5*log(1+Bjj[j,j]*sum(xj^2/v)) +
      (Bjj[j,j]*sum(w*xj/v)^2)/(2*tau2*sigma*(1+Bjj[j,j]*sum(xj^2/v)))  
    
    chance[j] <-ifelse(is.infinite(exp(l1)) == TRUE,1,exp(l1) / (exp(l1) + exp(l0)))
    # sample pi_j from a Bernoulli
    PI_aux[j] <- rbinom(1, 1, chance[j])
  }
  return(rbind(PI_aux[1:numcov],chance[1:numcov]))
}
