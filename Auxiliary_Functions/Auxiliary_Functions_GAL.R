LogitGama<-function(gama,L,U){
  LogitGama <- log((gama-L)/(U-gama))
  return(LogitGama)
} 
theta_to_gama<-function(theta,L,U){
  gama  <- (U*exp(theta)+L)/(1+exp(theta))
  return(gama)
}
