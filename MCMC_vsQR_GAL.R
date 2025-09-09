#' Bayesian variable selection in quantile regression with random effects - GAL
#'
#' This function estimates a Bayesian variable selection in quantile 
#' regression with random effects, where the response variable is assumed to 
#' follow a generalized asymmetric Laplace distribution.
#'
#' @param y n-dimensional vector of responses.
#' @param x Matrix - nx(p+1) - of predictors (include the intercept).
#' @param idobs n-dimensional vector that indicates groups (categorical numerical, 
#' from 1 to the number of groups).
#' @param tau Quantile of interest (value between 0 and 1).
#' @param n_mcmc Number of iterations.
#' @param burnin_mcmc Number of initial iterations to be discarded.
#' @param thin_mcmc Thinning parameter.
#' 
#' @return A list with the chains of all parameters of interest.

## Packages
require(mvtnorm); require(GIGrvg); require(msm)

## MCMC
vsQR_re_gal <- function(y,x,idobs,tau,n_mcmc,burnin_mcmc,thin_mcmc){
  n <- length(y)
  numcov <- ncol(x)
  n_idobs <- length(unique(idobs))
  resultado <- list()
  # Create auxiliary objects
  m_alpha <- matrix(0,n,n_idobs)
  for(l in 1:n){
    m_alpha[l,as.numeric(idobs[l])] <- 1 
  }
  beta  <- matrix(0, n_mcmc, numcov)
  sigma <- matrix(NA, n_mcmc, 1)
  gama  <- matrix(NA, n_mcmc, 1)
  theta <- matrix(NA, n_mcmc, 1)
  delta <- matrix(NA, n_mcmc, 1)
  PI    <- matrix(NA, n_mcmc,numcov)
  alpha <- matrix(NA, n_mcmc, n_idobs)
  psi2  <- matrix(NA, n_mcmc, 1)
  malpha<- matrix(NA, n_mcmc, 1)
  v     <- matrix(NA,n_mcmc,n)
  s     <- matrix(NA,n_mcmc,n)
  PROB  <- matrix(NA, n_mcmc,numcov)
  # Set the initial values
  beta[1,]   <- rnorm(numcov,0,1)
  sigma[1,1] <- 1
  v[1,]      <- rgamma(n,2,sigma[1,1])
  s[1,]      <- rtnorm(n,0,1,lower=0,upper=Inf)
  alpha[1,]  <- rnorm(n_idobs,0,1)
  psi2[1,1]  <- 1
  malpha[1,1]<- 0.5
  delta[1,]<-c(.5)
  PI[1,]<- rep(1,numcov)
  gama[1,1]  <- 0
  p     <- (gama[1,1]<0) + (tau-(gama[1,1]<0))/(2*pnorm(-abs(gama[1,1]),0,1)*exp(0.5*gama[1,1]^2))
  B     <- 2 / (p * (1 - p))
  A     <- (1 - 2 * p) / (p * (1 - p))
  C     <- 1/((gama[1,1]>0)-p)
  # Set bounds for gamma
  f <- function(gama1){
    2*pnorm(-abs(gama1), 0, 1)*exp(0.5*gama1^2) - tau
  } 
  UB    <- uniroot(f, c(0,20))$root
  f <- function(gama1){
    2*pnorm(-abs(gama1), 0, 1)*exp(0.5*gama1^2) - (1-tau)
  } 
  LB    <- uniroot(f, c(-20,0))$root
  theta[1,1] <- LogitGama(gama[1,1],LB,UB)
  # Prior (beta)
  mu0 <- rep(0,numcov)
  sigma0 <- diag(rep(100,numcov))
  # MCMC
  for(k in 2:n_mcmc){
    alpha_aux  <- m_alpha%*%alpha[k-1,]
    #Variable Selection
    delta[k,1] <- atualizarDELTA(1,1,PI[k-1,])
    PI.aux     <- atualizarPI(numcov,PI[k-1,],y,x,beta[k-1,],A,B,v[k-1,],sigma[k-1,1],C,delta[k,1],sigma0,gama[k-1,1],s[k-1,],alpha_aux)
    PI[k,]     <- PI.aux[1,]
    PROB[k,]   <- PI.aux[2,]
    n_select   <- which(PI[k,]==1)
    beta[k,n_select] <- atualizarBETA(mu0[n_select],sigma0[n_select,n_select],x[,n_select],sigma[k-1,1],A,B,v[k-1,],y,C,gama[k-1,1],s[k-1,],alpha_aux)
    ######################
    v[k,]      <- atualizarV(y,x,beta[k,],A,B,sigma[k-1,1],n,C,gama[k-1,1],s[k-1,],alpha_aux)
    s[k,]      <- atualizarS(y,x,beta[k,],A,B,v[k,],sigma[k-1,1],n,C,gama[k-1,1],alpha_aux)
    sigma[k,1] <- atualizarSIGMA(0.01,0.01,x,beta[k,],A,B,v[k,],y,n,C,gama[k-1,1],s[k,],alpha_aux)
    theta[k,1] <- atualizarGAMA(theta[k-1,1],y,x,beta[k,],v[k,],sigma[k,1],s[k,],n,LB,UB,alpha_aux)
    gama[k,1]  <- theta_to_gama(theta[k,1],LB,UB)
    p     <- (gama[k,1]<0) + (tau-(gama[k,1]<0))/(2*pnorm(-abs(gama[k,1]),0,1)*exp(0.5*gama[k,1]^2))
    A     <- (1 - 2 * p) / (p * (1 - p))
    B     <- 2 / (p * (1 - p))
    C     <- 1/((gama[k,1]>0)-p)
    alpha[k,]  <- atualizarALPHA(y,x,beta[k,],A,B,v[k,],sigma[k,1],C,gama[k-1,1],s[k,],malpha[k-1,1],psi2[k-1,1],t(m_alpha),n_idobs)
    malpha[k,1]<- atualizarMALPHA(0,1,alpha[k,],psi2[k-1,1],n_idobs)
    psi2[k,1]  <- atualizarPSI2(0.01,0.01,alpha[k,],malpha[k,1],n_idobs)
  }
  resultado[['beta']]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['sigma']] <- sigma[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['alpha']] <- alpha[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['gamma']] <- gama[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['pi']]    <- PI[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['qsi']]   <- PROB[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['malpha']]<- malpha[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['psi2']]  <- psi2[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  
  return(resultado)
}


#######################################################################
## Auxiliary functions: Sampling from full conditional distributions ##
#######################################################################

# Full conditional for beta
atualizarBETA<-function(mu0,sigma0,x,sigma,A,B,v,y,C,gama,s,alpha){
  B.inv  <- chol2inv(chol(sigma0))
  covar  <- chol2inv(chol(B.inv+(t(x)%*%diag(1/v)%*%x)/(B*sigma)))
  media  <- covar%*%((B.inv%*%mu0)+(t(x)%*%diag(1/v)%*%(y-alpha-sigma*C*abs(gama)*s-A*v))/(B*sigma))
  beta   <- rmvnorm(1,media,covar)
  return(beta)
}

# Full conditional for sigma
atualizarSIGMA<-function(c1,C1,x,beta,A,B,v,y,N,C,gama,s,alpha){
  lambda <- -(c1 + 1.5*N)
  chi    <- 2*C1 + 2*sum(v) + sum((y-alpha-x%*%beta-A*v)^2/v)/B
  psi    <- sum((C*abs(gama)*s)^2/(B*v))
  sigma  <- rgig(1, chi=chi, psi=psi, lambda=lambda)
  return(sigma)
}

# Full conditional for latent variables
atualizarV<-function(y,x,beta,A,B,sigma,N,C,gama,s,alpha){
  p1 <- 0.5
  p2 <- (y-alpha-x%*%beta-sigma*C*abs(gama)*s)^2/(B*sigma)
  p3 <- 2/sigma + A^2/(B*sigma)
  v  <- NULL
  for(i in 1:N){
    v[i]  <- rgig(1, chi=p2[i], psi=p3, lambda=p1)
  }
  return(v)
}
atualizarS<-function(y,x,beta,A,B,v,sigma,N,C,gama,alpha){
  sd    <- sqrt(1/(1+(C*gama)^2*sigma/(B*v)))
  media <- sd^2*C*abs(gama)*(y-alpha-x%*%beta-A*v)/(B*v)
  s     <- rtnorm(N,media,sd,lower=0,upper=Inf)
  return(s)
}

# Full conditional for gamma
condicionalGAMA<-function(theta,y,x,beta,v,sigma,s,N,L,U,alpha){
  gama <- theta_to_gama(theta,L,U)
  p <- (gama<0) + (tau-(gama<0))/(2*pnorm(-abs(gama),0,1)*exp(0.5*gama^2))
  A <- (1-2*p)/(p*(1-p))
  B <- 2/(p*(1-p)); C<-1/((gama>0)-p)
  C <- 1/((gama>0)-p)
  
  logJ_proposal <- theta-2*log(1+exp(theta))
  logprior      <- dbeta((gama-L)/(U-L),2,2,log=T) 
  logJ_prior    <- -log(U-L)	
  
  funcao <- -N/2*log(B) - sum((y-alpha-x%*%beta-sigma*C*abs(gama)*s-A*v)^2/v)/(2*sigma*B) + logJ_proposal + logprior + logJ_prior
  return(funcao)
}
# Metropolis-Hasting for gamma 
atualizarGAMA<-function(theta,y,x,beta,v,sigma,s,N,L,U,alpha){
  valoratual    <- theta
  valorproposto <- rnorm(1,valoratual,0.1)
  candidato     <- condicionalGAMA(valorproposto,y,x,beta,v,sigma,s,N,L,U,alpha)-condicionalGAMA(valoratual,y,x,beta,v,sigma,s,N,L,U,alpha)
  chanceaceitar <- min(1,candidato)
  contador<-NULL
  if (runif(1)<chanceaceitar){
    THETAfinal<-valorproposto
  } 
  else{
    THETAfinal<-valoratual 
  }
  return(THETAfinal)
}

# Full conditional for the random effects
atualizarALPHA<-function(y,x,beta,A,B,v,sigma,C,gama,s,alpha_mean,psi2,m_aux,N){
  z     <- (y-x%*%beta-sigma*C*abs(gama)*s-A*v)/v
  sd    <- sqrt((psi2*B*sigma)/(B*sigma+psi2*(m_aux%*%(1/v))))
  media <- (B*sigma*alpha_mean + psi2*(m_aux%*%z[,1])) / (B*sigma+psi2*(m_aux%*%(1/v)))
  alpha <- rnorm(N,media,sd)
  return(alpha)
}

# Full conditional for the hierarchical prior
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

# Variable Selection
atualizarDELTA<-function(a1,b1,vetorPI){
  delta<-rbeta(1, a1 + sum(vetorPI), b1 + sum(1 - vetorPI))
  return(delta)
}
atualizarPI<-function(numcov,PI,y,x,beta,A,B,v,sigma,C,delta,Bjj,gama,s,alpha){
  chance <- NULL
  for (j in sample(1:numcov)) {
    
    # get the betas for which beta_j is zero
    PI_aux<- PI 
    PI_aux[j] <- 0
    beta_j_aux <- c(beta)*PI_aux
    
    # Compute the z variables and the conditional variance
    w  <- y-alpha-x%*%beta_j_aux-sigma*C*abs(gama)*s-A*v
    xj <- as.vector(x[,j])
    
    # Compute the chance parameter of the conditional posterior of pi_j (Bernoulli)
    l0 <- log(1-delta)
    l1 <- log(delta) - 0.5*log(1+Bjj[j,j]*sum(xj^2/v)) +
      (Bjj[j,j]*sum(w*xj/v)^2)/(2*B*sigma*(1+Bjj[j,j]*sum(xj^2/v)))  
    
    chance[j] <-ifelse(is.infinite(exp(l1)) == TRUE,1,exp(l1) / (exp(l1) + exp(l0)))
    # sample pi_j from a Bernoulli
    PI_aux[j] <- rbinom(1, 1, chance[j])
  }
  return(rbind(PI_aux[1:numcov],chance[1:numcov]))
}


##################################################
## Auxiliary functions: Supplementary functions ##
##################################################

LogitGama<-function(gama,L,U){
  LogitGama <- log((gama-L)/(U-gama))
  return(LogitGama)
}

theta_to_gama<-function(theta,L,U){
  gama  <- (U*exp(theta)+L)/(1+exp(theta))
  return(gama)
}


