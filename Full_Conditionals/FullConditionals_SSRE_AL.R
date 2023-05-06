# Full conditional for beta
atualizarBETA<-function(b,B,spatial,x,sigma,tau2,theta,v,dados){
  B.inv  <- chol2inv(chol(B))
  covar  <- chol2inv(chol(B.inv+(1/(tau2*sigma)*(t(x/v)%*%x))))
  media  <- covar%*%((B.inv%*%b)+(1/(tau2*sigma))*(t(x/v)%*%(dados-spatial-theta*v)))
  beta   <- rmvnorm(1,media,covar)
  return(beta)
}
# Full conditional for sigma
atualizarSIGMA<-function(c,C,spatial,x,beta,tau2,theta,v,dados,N){
  alpha1 <- c + 1.5*N
  beta1  <- C + sum(v) + (t((dados-spatial-x%*%beta-theta*v)/v)%*%(dados-spatial-x%*%beta-theta*v))/(2*tau2)
  sigma  <- 1/rgamma(1, alpha1, beta1)
  return(sigma)
}
# Full conditional for the latent variable
atualizarV<-function(dados,spatial,x,beta,tau2,theta,sigma,N){
  p1 <- 0.5
  p2 <- (dados-spatial-x%*%beta)^2/(tau2*sigma)
  p3 <- 2/sigma + theta^2/(tau2*sigma)
  v  <- NULL
  for(i in 1:N){
    v[i]  <- rgig(1, chi=p2[i], psi=p3, lambda=p1)
  }
  return(v)
}
# Full conditional for the random effects
atualizarALPHA<-function(dados,x,beta,A,B,v,sigma,alpha_mean,psi2,m_aux,N){
  z     <- (dados-x%*%beta-A*v)/v
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
atualizarPI<-function(n_covariaveis,PI,dados,alpha,covariaveis,beta,v,theta,tau2,sigma,delta,Bjj){
  chance <- NULL
  for (j in sample(1:n_covariaveis)) {
    
    # get the betas for which beta_j is zero
    PI_aux<- PI
    PI_aux[j] <- 0
    beta_j_aux <- c(beta)*PI_aux
    
    # compute the z variables and the conditional variance
    w  <- dados-alpha-x%*%beta_j_aux-theta*v
    xj <- as.vector(covariaveis[,j])
    
    # compute chance parameter of the conditional posterior of pi_j (Bernoulli)
    l0 <- log(1-delta)
    l1 <- log(delta) - 0.5*log(1+Bjj[j,j]*sum(xj^2/v)) +
      (Bjj[j,j]*sum(w*xj/v)^2)/(2*tau2*sigma*(1+Bjj[j,j]*sum(xj^2/v)))  
    
    chance[j] <-ifelse(is.infinite(exp(l1)) == TRUE,1,exp(l1) / (exp(l1) + exp(l0)))
    # sample pi_j from a Bernoulli
    PI_aux[j] <- rbinom(1, 1, chance[j])
  }
  return(rbind(PI_aux[1:n_covariaveis],chance[1:n_covariaveis]))
}