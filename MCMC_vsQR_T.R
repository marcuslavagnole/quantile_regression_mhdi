#' Bayesian variable selection in quantile regression with random effects - 
#' Student-t
#'
#' This function estimates a Bayesian variable selection in quantile 
#' regression with random effects, where the response variable is assumed to 
#' follow a Student-t distribution.
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
require(mvtnorm); require(GIGrvg); require(msm)

## MCMC
vsQR_re_t <- function(y,x,idobs,tau,n_mcmc,burnin_mcmc,thin_mcmc){
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
  sigma <-matrix(NA, n_mcmc, 1)
  alpha <- matrix(NA, n_mcmc, n_idobs)
  psi2  <- matrix(NA, n_mcmc, 1)
  malpha<- matrix(NA, n_mcmc, 1)
  delta <- matrix(NA, n_mcmc, 1)
  PI    <- matrix(NA, n_mcmc,numcov)
  PROB  <- matrix(NA, n_mcmc,numcov)
  u     <- matrix(NA, n_mcmc, n)
  nu    <- matrix(NA, n_mcmc, 1)
  # Set the initial values
  beta[1,]   <- rnorm(numcov,0,1)
  sigma[1,1] <- 1
  alpha[1,]  <- rep(0,n_idobs)
  psi2[1,1]  <- 1
  malpha[1,1]<- 0.5
  delta[1,]<-c(.5)
  PI[1,]<- rep(1,numcov)
  nu[1,1] <- 4
  u[1,]   <- rgamma(n,nu[1,1]/2,nu[1,1]/2)
  # Create auxiliary objects for the MH
  rmw     <- matrix(NA,n_mcmc,3)
  rmw[1,] <- c(0.8,1,0)
  # Prior (beta)
  b <- rep(0,numcov)
  B <- diag(rep(100,numcov))
  #MCMC
  for(k in 2:n_mcmc){
    alpha_aux  <- m_alpha%*%alpha[k-1,]
    #Variable Selection
    delta[k,1] <- atualizarDELTA(1,1,PI[k-1,])
    PI.aux     <- atualizarPI(numcov,PI[k-1,],y,alpha_aux,x,beta[k-1,],sigma[k-1,1],delta[k,1],B,u[k-1,])
    PI[k,]     <- PI.aux[1,]
    PROB[k,]   <- PI.aux[2,]
    n_select   <- which(PI[k,]==1)
    beta[k,n_select] <- atualizarBETA(b[n_select],B[n_select,n_select],alpha_aux,x[,n_select],sigma[k-1,1],y,u[k-1,])
    ######################
    u[k,]      <- atualizarU(nu[k-1,1],y,alpha_aux,x,beta[k,],sigma[k-1,1],n)
    nu.aux     <- atualizarNU(nu[k-1,1],u[k,],n,rmw[k-1,1],rmw[k-1,2],rmw[k-1,3],k)
    nu[k,1]    <- nu.aux[1] ; rmw[k,] <- nu.aux[2:4]
    sigma[k,1] <- atualizarSIGMA(0.01,0.01,alpha_aux,x,beta[k,],y,n,u[k,])
    alpha[k,]  <- atualizarALPHA(y,x,beta[k,],sigma[k,1],malpha[k-1,1],psi2[k-1,1],t(m_alpha),n_idobs,u[k,])
    malpha[k,1]<- atualizarMALPHA(0,1,alpha[k,],psi2[k-1,1],n_idobs)
    psi2[k,1]  <- atualizarPSI2(0.01,0.01,alpha[k,],malpha[k,1],n_idobs)
  }
  resultado[['beta']]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['sigma']] <- sigma[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['alpha']] <- alpha[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['nu']]    <- sigma[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
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
atualizarBETA<-function(b,B,alpha,x,sigma,y,u){
  B.inv  <- chol2inv(chol(B))
  covar  <- chol2inv(chol(B.inv+(1/(sigma)*(t(u*x)%*%x))))
  media  <- covar%*%((B.inv%*%b)+(1/(sigma))*(t(u*x)%*%(y-alpha)))
  beta   <- rmvnorm(1,media,covar)
  return(beta)
}

# Full conditional for sigma
atualizarSIGMA<-function(c,C,alpha,x,beta,y,N,u){
  alpha1 <- c + 0.5*N 
  beta1  <- C + 0.5*(t(u*(y-alpha-x%*%beta))%*%(y-alpha-x%*%beta))
  sigma  <- 1/rgamma(1, alpha1, beta1)
  return(sigma)
}

# Full conditional for the random effects
atualizarALPHA<-function(y,x,beta,sigma,alpha_mean,psi2,m_aux,N,u){
  z     <- u*(y-x%*%beta)
  sd    <- sqrt((psi2*sigma)/(sigma+psi2*(m_aux%*%u)))
  media <- (sigma*alpha_mean + psi2*(m_aux%*%z[,1])) / (sigma+psi2*(m_aux%*%u))
  alpha <- rnorm(N,media,sd)
  return(alpha)
}

# Full conditional for the latent variable
atualizarU<-function(nu,y,alpha,x,beta,sigma,N){
  alpha <- c(rep(nu/2,N))
  beta  <- nu/2 + (y-alpha-x%*%beta)^2/(2*sigma)
  u     <- rgamma(N,alpha,beta)
  return(u)
}

# Full conditional for the degrees of freedom
condicionalNU<-function(nu,u,N){
  d<-9/(1+sqrt(2))
  priori<- log(nu)-3*log(nu+d)
  verossi<- (N*nu/2)*log(nu/2)-N*log(gamma(nu/2))+(nu/2-1)*sum(log(u))-nu/2*sum(u)
  funcao<-priori+verossi
  return(funcao)
}

# Metropolis-Hasting for the degrees of freedom
atualizarNU<-function(nu,u,N,clap,clap.aux,M0,t){
  valoratual<-nu
  valorproposto<-rtnorm(1,valoratual,sqrt(clap*clap.aux),lower=2,upper=40)
  candidato<-exp(condicionalNU(valorproposto,u,N)-condicionalNU(valoratual,u,N)-dtnorm(valorproposto,valoratual,sqrt(clap*clap.aux),lower=2,upper=40,log=TRUE)+dtnorm(valoratual,valorproposto,sqrt(clap*clap.aux),lower=2,upper=40,log=TRUE))
  chanceaceitar<-min(1,candidato)
  contador<-NULL
  if (runif(1)<chanceaceitar){
    NUfinal<-valorproposto
  } 
  else{
    NUfinal<-valoratual 
  }
  gama1<- 1/t^0.5
  gama2<- 1*gama1
  termometro<- exp(log(clap)+gama2*(chanceaceitar-0.234))
  termometro.aux<- clap.aux+gama1*((NUfinal-M0)^2-clap.aux)
  p.auxiliar<-M0+gama2*(NUfinal-M0)
  return(c(NUfinal,termometro,termometro.aux,p.auxiliar))
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
atualizarPI<-function(numcov,PI,y,alpha,covariaveis,beta,sigma,delta,Bjj,u){
  chance <- NULL
  for (j in sample(1:numcov)) {
    
    # get the betas for which beta_j is zero
    PI_aux<- PI
    PI_aux[j] <- 0
    beta_j_aux <- c(beta)*PI_aux
    
    # Compute the z variables and the conditional variance
    w  <- y-alpha-x%*%beta_j_aux
    xj <- as.vector(covariaveis[,j])
    
    # Compute the chance parameter of the conditional posterior of pi_j (Bernoulli)
    l0 <- log(1-delta)
    l1 <- log(delta) - 0.5*log(1+Bjj[j,j]*sum(xj^2*u)) +
      (Bjj[j,j]*sum(w*xj*u)^2)/(2*sigma*(1+Bjj[j,j]*sum(xj^2*u)))  
    
    chance[j] <-ifelse(is.infinite(exp(l1)) == TRUE,1,exp(l1) / (exp(l1) + exp(l0)))
    # sample pi_j from a Bernoulli
    PI_aux[j] <- rbinom(1, 1, chance[j])
  }
  return(rbind(PI_aux[1:numcov],chance[1:numcov]))
}
