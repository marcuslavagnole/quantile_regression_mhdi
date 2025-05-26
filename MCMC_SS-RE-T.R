# Load libraries
library(msm)
library(mvtnorm)
library(GIGrvg)

# Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Full_Conditionals/FullConditionals_SSRE_T.R")
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

PROB  <- matrix(NA, NN,numcov)

u     <- matrix(NA, NN, n)
nu    <- matrix(NA, NN, 1)

# Set the initial values
beta[1,]   <- rnorm(numcov,0,1)
sigma[1,1] <- 1
alpha[1,]  <- rep(0,n_meso)
psi2[1,1]  <- 1
malpha[1,1]<- 0.5
delta[1,]<-c(.5)
PI[1,]<- rep(1,numcov)
nu[1,1] <- 4
u[1,]   <- rgamma(n,nu[1,1]/2,nu[1,1]/2)

# Create auxiliary objects for the MH
rmw     <- matrix(NA,NN,3)
rmw[1,] <- c(0.8,1,0)
contador<- 0

#MCMC
for(k in 2:NN){
  alpha_aux  <- m_alpha%*%alpha[k-1,]
  #Variable Selection
  delta[k,1] <- atualizarDELTA(1,1,PI[k-1,])
  PI.aux     <- atualizarPI(numcov,PI[k-1,],dados,alpha_aux,x,beta[k-1,],sigma[k-1,1],delta[k,1],B,u[k-1,])
  PI[k,]     <- PI.aux[1,]
  PROB[k,]   <- PI.aux[2,]
  n_select   <- which(PI[k,]==1)
  beta[k,n_select] <- atualizarBETA(b[n_select],B[n_select,n_select],alpha_aux,x[,n_select],sigma[k-1,1],dados,u[k-1,])
  ######################
  u[k,]      <- atualizarU(nu[k-1,1],dados,alpha_aux,x,beta[k,],sigma[k-1,1],n)
  nu.aux     <- atualizarNU(nu[k-1,1],u[k,],n,rmw[k-1,1],rmw[k-1,2],rmw[k-1,3],k)
  nu[k,1]    <- nu.aux[1] ; contador <- contador + nu.aux[2] ; rmw[k,] <- nu.aux[3:5]
  sigma[k,1] <- atualizarSIGMA(c,C,alpha_aux,x,beta[k,],dados,n,u[k,])
  
  alpha[k,]  <- atualizarALPHA(dados,x,beta[k,],sigma[k,1],malpha[k-1,1],psi2[k-1,1],t(m_alpha),n_meso,u[k,])
  malpha[k,1]<- atualizarMALPHA(0,1,alpha[k,],psi2[k-1,1],n_meso)
  psi2[k,1]  <- atualizarPSI2(0.01,0.01,alpha[k,],malpha[k,1],n_meso)
}
