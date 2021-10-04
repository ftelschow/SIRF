## This short script demonstrates that the Model C overcoverage

library(SCBfun)
library(RFT)
library(SampleFields)
library(tidyverse)


#set.seed(666)
# number of samples
nSamp = 800
# sampling locations
x        = seq( 0, 1, length.out = 175 )
coords   = list()
coords$x = x
# level of control
alpha = 0.05

model = "ModelC"

if(model ==  "ModelA"){
  # true mean curve
  mu    <- Vectorize( function(x){sin(4 * pi * x) * exp(-3 * x)} )
  # true standard deviation
  sigma <- Vectorize( function(x){((1 - x - 0.4)^2 + 1) / 6} )
  # noise model
  noise = RandomNormalSum  
}else if(model ==  "ModelB"){
  # true mean curve
  mu    <- Vectorize( function(x){sin(4 * pi * x) * exp(-3 * x)} )
  # true standard deviation
  sigma <- Vectorize( function(x){(1.5 - x)/2} )
  # noise model
  covf <- function(x, y) SampleFields::covf.nonst.matern(x,
                                                         y,
                                                         params = c( 1, 1 / 4, 0.4 ))
  noise = function(N, x){ ArbCovProcess( N, x, covf = covf ) }
}else{
  # true mean curve
  mu    <- Vectorize(function(x){sin(4 * pi * x) * exp(-3 * x)} )
  # true standard deviation
  sigma <- Vectorize(function(x){(1.5 - x)} )
  # noise model
  noise = DegrasNonGaussProcess 
}

Msim  = 5e2
biasY = matrix( NA, length(x), Msim )

for( m in 1:Msim ){
  # generate data from the sample
  Y = SignalPlusNoise( N = nSamp,
                       x = x,
                       mu = mu,
                       sigma = sigma,
                       noise = noise )
  
  # get the normalized functional delta residuals
  scbY = scb_moments(  Y$values,
                       level          = .95,
                       transformation = "skewness",
                       moments        = NULL,
                       method  = list(name = "tGKF", field = "t"),
                       bias    = TRUE,
                       se.est  = TRUE,
                       coords  = NULL,
                       mask    = NULL )
  biasY[, m] = scbY$bias
}

plot( x, rowMeans(biasY), main = paste( "N =", nSamp) )
plot( x, matrixStats::rowVars(biasY), main = paste( "N =", nSamp) )

matplot(x, biasY[, 1:20], main = paste( "N =", nSamp) )
