# Full conditional for beta
atualizarBETA<-function(mu0,sigma0,x,sigma,A,B,v,dados,C,gama,s,spatial){
  B.inv  <- chol2inv(chol(sigma0))
  covar  <- chol2inv(chol(B.inv+(t(x)%*%diag(1/v)%*%x)/(B*sigma)))
  media  <- covar%*%((B.inv%*%mu0)+(t(x)%*%diag(1/v)%*%(dados-spatial-sigma*C*abs(gama)*s-A*v))/(B*sigma))
  beta   <- rmvnorm(1,media,covar)
  return(beta)
}

# Full conditional for sigma
atualizarSIGMA<-function(c1,C1,x,beta,A,B,v,dados,N,C,gama,s,spatial){
  lambda <- -(c1 + 1.5*N)
  chi    <- 2*C1 + 2*sum(v) + sum((dados-spatial-x%*%beta-A*v)^2/v)/B
  psi    <- sum((C*abs(gama)*s)^2/(B*v))
  sigma  <- rgig(1, chi=chi, psi=psi, lambda=lambda)
  return(sigma)
}

# Full conditional for latent variables
atualizarV<-function(dados,x,beta,A,B,sigma,N,C,gama,s,spatial){
  p1 <- 0.5
  p2 <- (dados-spatial-x%*%beta-sigma*C*abs(gama)*s)^2/(B*sigma)
  p3 <- 2/sigma + A^2/(B*sigma)
  v  <- NULL
  for(i in 1:N){
    v[i]  <- rgig(1, chi=p2[i], psi=p3, lambda=p1)
  }
  return(v)
}
atualizarS<-function(dados,x,beta,A,B,v,sigma,N,C,gama,spatial){
  sd    <- sqrt(1/(1+(C*gama)^2*sigma/(B*v)))
  media <- sd^2*C*abs(gama)*(dados-spatial-x%*%beta-A*v)/(B*v)
  s     <- rtnorm(N,media,sd,lower=0,upper=Inf)
  return(s)
}

# Full conditional for gamma
condicionalGAMA<-function(theta,dados,x,beta,v,sigma,s,N,L,U,spatial){
  gama <- theta_to_gama(theta,L,U)
  p <- (gama<0) + (p0-(gama<0))/(2*pnorm(-abs(gama),0,1)*exp(0.5*gama^2))
  A <- (1-2*p)/(p*(1-p))
  B <- 2/(p*(1-p)); C<-1/((gama>0)-p)
  C <- 1/((gama>0)-p)
  
  logJ_proposal <- theta-2*log(1+exp(theta))
  logprior      <- dbeta((gama-L)/(U-L),2,2,log=T) 
  logJ_prior    <- -log(U-L)	
  
  funcao <- -N/2*log(B) - sum((dados-spatial-x%*%beta-sigma*C*abs(gama)*s-A*v)^2/v)/(2*sigma*B) + logJ_proposal + logprior + logJ_prior
  return(funcao)
}
# Metropolis-Hasting for gamma 
atualizarGAMA<-function(theta,dados,x,beta,v,sigma,s,N,L,U,spatial){
  valoratual    <- theta
  valorproposto <- rnorm(1,valoratual,0.1)
  candidato     <- condicionalGAMA(valorproposto,dados,x,beta,v,sigma,s,N,L,U,spatial)-condicionalGAMA(valoratual,dados,x,beta,v,sigma,s,N,L,U,spatial)
  chanceaceitar <- min(1,candidato)
  contador<-NULL
  if (runif(1)<chanceaceitar){
    THETAfinal<-valorproposto
    contador<-1
  } 
  else{
    THETAfinal<-valoratual 
    contador<-0
  }
  return(c(THETAfinal,contador))
}

# Full conditional for the random effects
atualizarALPHA<-function(dados,x,beta,A,B,v,sigma,C,gama,s,alpha_mean,psi2,m_aux,N){
  z     <- (dados-x%*%beta-sigma*C*abs(gama)*s-A*v)/v
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
atualizarPI<-function(n_covariaveis,PI,dados,x,beta,A,B,v,sigma,C,delta,Bjj,gama,s,spatial){
  chance <- NULL
  for (j in sample(1:n_covariaveis)) {
    
    # get the betas for which beta_j is zero
    PI_aux<- PI 
    PI_aux[j] <- 0
    beta_j_aux <- c(beta)*PI_aux
    
    # Compute the z variables and the conditional variance
    w  <- dados-spatial-x%*%beta_j_aux-sigma*C*abs(gama)*s-A*v
    xj <- as.vector(x[,j])
    
    # Compute the chance parameter of the conditional posterior of pi_j (Bernoulli)
    l0 <- log(1-delta)
    l1 <- log(delta) - 0.5*log(1+Bjj[j,j]*sum(xj^2/v)) +
      (Bjj[j,j]*sum(w*xj/v)^2)/(2*B*sigma*(1+Bjj[j,j]*sum(xj^2/v)))  
    
    chance[j] <-ifelse(is.infinite(exp(l1)) == TRUE,1,exp(l1) / (exp(l1) + exp(l0)))
    # sample pi_j from a Bernoulli
    PI_aux[j] <- rbinom(1, 1, chance[j])
  }
  return(rbind(PI_aux[1:n_covariaveis],chance[1:n_covariaveis]))
}
