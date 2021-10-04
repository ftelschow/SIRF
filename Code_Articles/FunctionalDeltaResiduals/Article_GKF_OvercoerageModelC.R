## This short script demonstrates that the Model C overcoverage

library(SCBfun)
library(RFT)
library(SampleFields)
library(tidyverse)


#set.seed(666)
# number of samples
nSamp = 200
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

# generate data from the sample
Y = SignalPlusNoise( N = nSamp,
                     x = x,
                     mu = mu,
                     sigma = sigma,
                     noise = noise )

# get the normalized functional delta residuals
res = SCBfun::DeltaMomentResiduals( Y$values, "cohensd" )
res = res$delta.res / as.vector(res$delta.sd)

res = ( Y$values - rowMeans(Y$values) ) / sqrt(matrixStats::rowVars(Y$values))

# get bootstrap residuals
res.boot = MultiplierBootstrap( res )$samples

# estimate the LKCs from the delta residuals
L = c( 1, LKC_integral( res, coords = list( x = x ) ) )

# get the EEC threshold
thresh = GKFthreshold( alpha = alpha/2, L, type = "t", nSamp - 1 )$threshold

# obtain the EEC of the excursion set / the number of connected components
number.components = apply( res, 2, function(x) EulerChar( abs(t(x)), thresh ) )

# Number of connected components of the bootstrapped process
number.components.boots = apply( res.boot, 2, function(x) EulerChar( abs(t(x)), thresh ) )


hist1 <- ggplot( tibble( components = number.components[number.components!=0] ),
                 aes( x = components ) ) + 
  geom_histogram( binwidth = 0.5)
hist2 <- ggplot( tibble( components = number.components.boots[number.components.boots!=0] ),
                 aes( x = components ) ) + 
  geom_histogram( binwidth = 0.5)

multiplot( hist1, hist2, cols = 2 )

print( mean( number.components ) )
print( mean( number.components.boots ) )
