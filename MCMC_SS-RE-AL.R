# Load libraries
library(msm)
library(mvtnorm)
library(GIGrvg)

# Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Full_Conditionals/FullConditionals_SSRE_AL.R")
load("data.RData")

# Define dependent and indepent variables
n      = dim(base)[1]
n_meso = 5
x      = scale(base[,c(-1,-dim(base)[2])])
dados  = log(base$IDHM_renda/(1-base$IDHM_renda))
numcov = ncol(x)

# Set hyperparameter values
b <- rep(0,numcov)
B <- diag(rep(100,numcov))
c <- 0.01
C <- 0.01

# Define the number of iterations
NN = 300000

# Create auxiliary objects
m_alpha <- matrix(0,n,n_meso)
for(l in 1:n){
  m_alpha[l,as.numeric(base$meso[l])] <- 1 
}

beta  <- matrix(0, NN, numcov)
sigma <-matrix(NA, NN, 1)

alpha <- matrix(NA, NN, n_meso)
psi2  <- matrix(NA, NN, 1)
malpha<- matrix(NA, NN, 1)

delta <- matrix(NA, NN, 1)
PI    <- matrix(NA, NN,numcov)

v     <- matrix(NA,NN,n)
PROB  <- matrix(NA, NN,numcov)

# Set the initial values
beta[1,]   <- rnorm(numcov,0,1)
sigma[1,1] <- 1
v[1,]      <- rgamma(n,2,sigma[1,1])
alpha[1,]  <- rep(0,n_meso)
psi2[1,1]  <- 1
malpha[1,1]<- 0.5
delta[1,]  <- c(.5)
PI[1,]<- rep(1,numcov)
q     <- 0.10
tau2  <- 2 / (q * (1 - q))
theta <- (1 - 2 * q) / (q * (1 - q))

# MCMC
for(k in 2:NN){
  alpha_aux  <- m_alpha%*%alpha[k-1,]
  #Variable Selection
  delta[k,1] <- atualizarDELTA(1,1,PI[k-1,])
  PI.aux     <- atualizarPI(numcov,PI[k-1,],dados,alpha_aux,x,beta[k-1,],v[k-1,],theta,tau2,sigma[k-1,1],delta[k,1],B)
  PI[k,]     <- PI.aux[1,]
  PROB[k,]   <- PI.aux[2,]
  n_select   <- which(PI[k,]==1)
  beta[k,n_select] <- atualizarBETA(b[n_select],B[n_select,n_select],alpha_aux,x[,n_select],sigma[k-1,1],tau2,theta,v[k-1,],dados)
  ######################
  
  v[k,]      <- atualizarV(dados,alpha_aux,x,beta[k,],tau2,theta,sigma[k-1,1],n)
  sigma[k,1] <- atualizarSIGMA(c,C,alpha_aux,x,beta[k,],tau2,theta,v[k,],dados,n)
  
  alpha[k,]  <- atualizarALPHA(dados,x,beta[k,],theta,tau2,v[k,],sigma[k,1],malpha[k-1,1],psi2[k-1,1],t(m_alpha),n_meso)
  malpha[k,1]<- atualizarMALPHA(0,1,alpha[k,],psi2[k-1,1],n_meso)
  psi2[k,1]  <- atualizarPSI2(0.01,0.01,alpha[k,],malpha[k,1],n_meso)
}
