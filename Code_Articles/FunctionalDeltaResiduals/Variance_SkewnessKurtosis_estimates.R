#-------------------------------------------------------------------------------
#
#            Variance estimates of delta residuals
#
#-------------------------------------------------------------------------------
# Preapre workspace
#-------------------------------------------------------------------------------
# Clear Workspace
rm( list = ls() )

# Load packages
library( SIRF )
library( tidyverse )
library( SampleFields )

Vibrant = c( "#0077BB", "#33BBEE", "#009988",
             "#EE7733", "#CC3311", "#EE3377",
             "#BBBBBB" )

path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"

path_pics = "/home/fabian/Seafile/Projects/2019_DeltaResiduals/Pics/"


#-------------------------------------------------------------------------------
# General Simulation setup
#-------------------------------------------------------------------------------
# General Simulation parameters
Msim    = 5e3 # number of simulations
Nvec    = c(50, 100, 200, 400, 800) # sample sizes considered
x       = seq(0, 1, length.out = 175) # locations the process is evaluated at


# Simulation function
sim_varEsts <- function(Msim, method, Nvec, x, mu, noise, sigma){
  sim_var = array( NaN, dim = c( length(Nvec), length(x), Msim ) )
  for( nn in 1:length(Nvec) ){
    for( m in 1:Msim ){
      # generate data
      N = Nvec[nn]
      Y1 <- SignalPlusNoise( N = N,
                             x = x,
                             mu = mu,
                             noise = noise,
                             sigma = sigma )
      # get the deltaresidual
      sim_var[nn, , m] <- (DeltaMomentResiduals( Y = Y1$values,
                                                 transformation = method,
                                                 moments = NULL )$delta.sd)^2
    }
  }
  return(sim_var)
}

# Function computing the true Gaussian variance of Cohen's d, skewness and kurtosis
trueVar_Gauss <- function(Nvec, method, x, mu = NULL, sigma = NULL){
  if(method == "cohensd"){
    nu  = Nvec - 1
    d  = function(v)  sqrt(v) * mu(x) / sigma(x)

    trueVar = t(vapply(nu, function(v) v * (1 + d(v)^2) / (v - 2) -d(v)^2 *
                                           (1 - 3 / (4 * v -1 ))^(-2),
                       FUN.VALUE = rep(1, length(x))
    ))
  }else if(method == "skewness"){
    trueVar = matrix(6 * Nvec * (Nvec - 2) / (Nvec + 1) / (Nvec + 3), length(Nvec), length(x))
  }else{
    trueVar = matrix(24 * Nvec^2 * (Nvec - 1)^2 / (Nvec - 3) /
                       (Nvec - 2) / (Nvec + 3) / (Nvec + 5), length(Nvec), length(x))
  }
  rownames(trueVar) <- Nvec

  return(trueVar)
}

# Plot function
plot_var_ests <- function(method, x, estVar, trueVar, outname = NULL){
  relBias = ( estVar - trueVar ) / trueVar
  NN      = dim(estVar)[1]

  if(is.null(outname)){
    par(mfrow = c(1,2))
    plot( NULL,
          xlim = range(x),
          ylim = c(min( 0, min(c(relBias))), max( 0, max(c(relBias)))),
          xlab = "location",
          ylab = "(est.var - true.var) / true.var",
          main = paste("Relative Bias (", method,")", sep = "") )
    abline( h = 0, lty = 2, lwd = 2 )
    matlines( x = x, y = t(relBias), lty = 1, col = 1:NN )
    legend( "bottomright",   # Coordinates (x also accepts keywords)
            legend = rownames(estVar), # Vector with the name of each group
            col = 1:NN, # Color of lines or symbols
            lty = rep(1, NN)
    )
    plot( NULL,
          xlim = range(x),
          ylim = range(c(c(estVar), c(trueVar))),
          ylab = paste("variance (", method,")", sep = ""),
          main = "Estimated vs true variance" )
    matlines( x, t(estVar), lty = 1, col = 1:NN )
    matlines( x, t(trueVar), lty = 2, col = 1:NN )
  }else{
    png( filename = paste( filename, "_", method, ".png" , sep = "" ),
         width = 1000,
         height = 500
    )
    par(mfrow = c(1,2))
    plot( NULL,
          xlim = range(x),
          ylim = c(min( 0, min(c(relBias))), max( 0, max(c(relBias)))),
          xlab = "location",
          ylab = "(est.var - true.var) / true.var",
          main = paste("Relative Bias (", method,")", sep = "") )
    abline( h = 0, lty = 2, lwd = 2 )
    matlines( x = x, y = t(relBias), lty = 1, col = 1:NN )
    legend( "bottomright",   # Coordinates (x also accepts keywords)
            legend = rownames(estVar), # Vector with the name of each group
            col = 1:NN, # Color of lines or symbols
            lty = rep(1, NN)
    )
    plot( NULL,
          xlim = range(x),
          ylim = range(c(c(estVar), c(trueVar))),
          ylab = paste("variance (", method,")", sep = ""),
          main = "Estimated vs true variance" )
    matlines( x, t(estVar), lty = 1, col = 1:NN )
    matlines( x, t(trueVar), lty = 2, col = 1:NN )
    dev.off()
  }
}

methodVec = c("cohensd", "skewness", "kurtosis")
methodVec = c("skewness (normality)")

#-------------------------------------------------------------------------------
# Simulation model A
#-------------------------------------------------------------------------------
# Functional data model
mu_modelA    = function(x){sin(4 * pi * x) * exp(-3 * x)}
sigma_modelA = function(x){((1 - x - 0.4)^2 + 1) / 6}
noise_modelA = RandomNormalSum

# Filename
filename = paste(path_pics, "DeltaResidualEstimateVariance" , sep = "")
filename = NULL


# Simulate the variance of cohens d
method = "cohensd"
  sim_var_d <- sim_varEsts(Msim, method = method, Nvec, x,
                           mu = mu_modelA,
                           noise = noise_modelA,
                           sigma = sigma_modelA)

  # Get the mean of the simulated variance and the relative bias
  estVar = apply( sim_var_d, 1:2, mean )
  rownames(estVar) <- Nvec
  trueVar <- trueVar_Gauss(Nvec, method = method, x, mu = mu_modelA, sigma = sigma_modelA)

  # Plot the data
  plot_var_ests(method = method, x, estVar, trueVar, outname = filename)


  # Simulate the variance of skewness
method = "skewness"
  sim_var_s <- sim_varEsts(Msim, method = method, Nvec, x,
                           mu = mu_modelA,
                           noise = noise_modelA,
                           sigma = sigma_modelA)

  # Get the mean of the simulated variance and the relative bias
  estVar = apply( sim_var_s, 1:2, mean )
  rownames(estVar) <- Nvec
  trueVar_s <- trueVar_Gauss(Nvec, method = method, x, mu = mu_modelA, sigma = sigma_modelA)

  # Plot the data
  plot_var_ests(method = method, x, estVar, trueVar, outname = filename)


  # Simulate the variance of skewness (normalization)
method = "skewness (normality)"
  sim_var_sn <- sim_varEsts(Msim, method = method, Nvec, x,
                           mu = mu_modelA,
                           noise = noise_modelA,
                           sigma = sigma_modelA)

  # Get the mean of the simulated variance and the relative bias
  estVar_sn = apply( sim_var_sn, 1:2, mean )
  rownames(estVar_sn) <- Nvec
  trueVar_sn <- rep(1, length.out = length(x))

  # Plot the data
  plot_var_ests(method = method, x, estVar, trueVar, outname = filename)


  save.image(paste(path_wd, "Workspaces/Variance_Simulation2.Rdata", sep = ""))


#-------------------------------------------------------------------------------
# Simulation model A
#-------------------------------------------------------------------------------
# Functional data model
mu_modelA    = function(x){sin(4 * pi * x) * exp(-3 * x)}
sigma_modelA = function(x){((1 - x - 0.4)^2 + 1) / 6}
noise_modelA = RandomNormalSum

# Filename
filename = paste( path_pics, "DeltaResidualEstimateVariance" , sep = "" )
filename = NULL

# Simulate the variance
for(method in methodVec){
  sim_var_d <- sim_varEsts(Msim, method = method, Nvec, x,
                           mu = mu_modelA,
                           noise = noise_modelA,
                           sigma = sigma_modelA)

  # Get the mean of the simulated variance and the relative bias
  estVar = apply( sim_var_d, 1:2, mean )
  rownames(estVar) <- Nvec
  trueVar <- trueVar_Gauss(Nvec, method = method, x, mu = mu_modelA, sigma = sigma_modelA)


  # Plot the data
  plot_var_ests(method = method, x, estVar, trueVar, outname = NULL)
}
