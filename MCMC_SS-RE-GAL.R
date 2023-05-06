# Load libraries
library(msm)
library(mvtnorm)
library(GIGrvg)

# Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Full_Conditionals/FullConditionals_SSRE_GAL.R")
source("Auxiliary_Functions/Auxiliary_Functions_GAL.R")
load("data.RData")

# Define dependent and indepent variables
n      = dim(base)[1]
n_meso = 5
x      = scale(base[,c(-1,-dim(base)[2])])
dados  = log(base$IDHM_renda/(1-base$IDHM_renda))
numcov = ncol(x)

# Set hyperparameters values
mu0    <- rep(0,numcov)
sigma0 <- diag(rep(100,numcov))
c1 <- 0.01
C1 <- 0.01

# Define number of iterations
NN = 300000

# Create auxiliary objects
m_alpha <- matrix(0,n,n_meso)
for(l in 1:n){
  m_alpha[l,as.numeric(base$meso[l])] <- 1 
}

beta  <- matrix(0, NN, numcov)
sigma <- matrix(NA, NN, 1)
gama  <- matrix(NA, NN, 1)
theta <- matrix(NA, NN, 1)
delta <- matrix(NA, NN, 1)
PI    <- matrix(NA, NN,numcov)
alpha <- matrix(NA, NN, n_meso)
psi2  <- matrix(NA, NN, 1)
malpha<- matrix(NA, NN, 1)
v     <- matrix(NA,NN,n)
s     <- matrix(NA,NN,n)
PROB  <- matrix(NA, NN,numcov)

# Set the initial values
beta[1,]   <- rnorm(numcov,0,1)
sigma[1,1] <- 1
v[1,]      <- rgamma(n,2,sigma[1,1])
s[1,]      <- rtnorm(n,0,1,lower=0,upper=Inf)
alpha[1,]  <- rnorm(n_meso,0,1)
psi2[1,1]  <- 1
malpha[1,1]<- 0.5
delta[1,]<-c(.5)
PI[1,]<- rep(1,numcov)
gama[1,1]  <- 0

p0    <- 0.10
p     <- (gama[1,1]<0) + (p0-(gama[1,1]<0))/(2*pnorm(-abs(gama[1,1]),0,1)*exp(0.5*gama[1,1]^2))
B     <- 2 / (p * (1 - p))
A     <- (1 - 2 * p) / (p * (1 - p))
C     <- 1/((gama[1,1]>0)-p)

# Set bounds for gamma
f <- function(gama1){
  2*pnorm(-abs(gama1), 0, 1)*exp(0.5*gama1^2) - p0
} 
UB    <- uniroot(f, c(0,20))$root

f <- function(gama1){
  2*pnorm(-abs(gama1), 0, 1)*exp(0.5*gama1^2) - (1-p0)
} 
LB    <- uniroot(f, c(-20,0))$root

theta[1,1] <- LogitGama(gama[1,1],LB,UB)

# Create auxiliary objects for the MH
contador<-NULL
contador[1]<-0

# MCMC
for(k in 2:NN){
  alpha_aux  <- m_alpha%*%alpha[k-1,]
  #Variable Selection
  delta[k,1] <- atualizarDELTA(1,1,PI[k-1,])
  PI.aux     <- atualizarPI(numcov,PI[k-1,],dados,x,beta[k-1,],A,B,v[k-1,],sigma[k-1,1],C,delta[k,1],sigma0,gama[k-1,1],s[k-1,],alpha_aux)
  PI[k,]     <- PI.aux[1,]
  PROB[k,]   <- PI.aux[2,]
  n_select   <- which(PI[k,]==1)
  beta[k,n_select] <- atualizarBETA(mu0[n_select],sigma0[n_select,n_select],x[,n_select],sigma[k-1,1],A,B,v[k-1,],dados,C,gama[k-1,1],s[k-1,],alpha_aux)
  ######################
  
  v[k,]      <- atualizarV(dados,x,beta[k,],A,B,sigma[k-1,1],n,C,gama[k-1,1],s[k-1,],alpha_aux)
  s[k,]      <- atualizarS(dados,x,beta[k,],A,B,v[k,],sigma[k-1,1],n,C,gama[k-1,1],alpha_aux)
  
  sigma[k,1] <- atualizarSIGMA(c1,C1,x,beta[k,],A,B,v[k,],dados,n,C,gama[k-1,1],s[k,],alpha_aux)
  
  gama.aux   <- atualizarGAMA(theta[k-1,1],dados,x,beta[k,],v[k,],sigma[k,1],s[k,],n,LB,UB,alpha_aux)
  theta[k,1] <- gama.aux[1] ; contador[k]<-gama.aux[2]
  gama[k,1]  <- theta_to_gama(theta[k,1],LB,UB)
  
  p     <- (gama[k,1]<0) + (p0-(gama[k,1]<0))/(2*pnorm(-abs(gama[k,1]),0,1)*exp(0.5*gama[k,1]^2))
  A     <- (1 - 2 * p) / (p * (1 - p))
  B     <- 2 / (p * (1 - p))
  C     <- 1/((gama[k,1]>0)-p)

  alpha[k,]  <- atualizarALPHA(dados,x,beta[k,],A,B,v[k,],sigma[k,1],C,gama[k-1,1],s[k,],malpha[k-1,1],psi2[k-1,1],t(m_alpha),n_meso)
  malpha[k,1]<- atualizarMALPHA(0,1,alpha[k,],psi2[k-1,1],n_meso)
  psi2[k,1]  <- atualizarPSI2(0.01,0.01,alpha[k,],malpha[k,1],n_meso)
}
