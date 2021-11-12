#-------------------------------------------------------------------------------
#
#       Script simulating the covering rate for SCBs with Delta residuals
#
#-------------------------------------------------------------------------------
# Preapre workspace
#-------------------------------------------------------------------------------
library(SIRF)
require(SampleFields)

Model   = "ModelA" # "ModelB",  "ModelC"
Msim    = 1e2 # number of simulations
Nvec    = c( 50, 100, 200, 400, 800 ) # sample sizes considered
x       = seq( 0, 1, length.out = 150 ) # locations the process is evaluated at
level   = .95  # level of simultaneous control
sim_tag = "test_adding_roots"
path_wd = "~/Seafile/Projects/2019_DeltaResiduals/"
date    = "YEAR_MO_DY"

  # General constants
  path_data <- paste( path_wd, "Workspaces/", sep = "" )

  #-------------------------------------------------------------------------------
  # General Simulation parameters
  #-------------------------------------------------------------------------------
  # Describing the parameters of the compared methods
  GKF    = list( name = "GKF", field = "z" )
  tGKF   = list( name = "tGKF", field = "t" )
  Mult   = list( name = "MultBoot", Mboots = 2.5e3, weights = "gaussian", method = "regular" )
  tMult  = list( name = "MultBoot", Mboots = 2.5e3, weights = "gaussian", method = "t"  )
  rMult  = list( name = "MultBoot", Mboots = 2.5e3, weights = "rademacher", method = "regular"  )
  trMult = list( name = "MultBoot", Mboots = 2.5e3, weights = "rademacher", method = "t"  )

  methvec = list( GKF = GKF, tGKF = tGKF,
                  Mult = Mult, tMult = tMult,
                  rMult = rMult, trMult = trMult )
  methvec = list( GKF = GKF, tGKF = tGKF )

  transformationvec = f ~ ( (mu4 - 4 * mu1 * mu3 + 6 * mu1^2 * mu2 - 3 * mu1^4) / (mu2 - mu1^2)^2 - 3 )
  momentsvec <- c("mu1", "mu2", "mu3", "mu4")

  biasvec = TRUE
  se.est <- TRUE

  #-------------------------------------------------------------------------------
  # Model A
  #-------------------------------------------------------------------------------
  if( Model == "ModelB" ){
    mu_model    = function(x){ ( x-0.3 )^2 }
    sigma_model = function(x){ (sin(3*pi*x) + 1.5) / 6 }
    noise_model = ArbCovProcess
    covf <- function(x, y) SampleFields::covf.nonst.matern(x,
                                                           y,
                                                           params = c( 1, 1 / 4, 0.4 ))

    trueValue = rbind(
      matrix( 0, 5, length(x) )
    )

  }else if( Model == "ModelC" ){
    mu_model    = function(x){  sin(4 * pi * x) * exp(-3 * x) }
    sigma_model = function(x){ (1.5 - x) }
    noise_model = DegrasNonGaussProcess

    # Get the true values for the models considered in this simulation
    trueValue = rbind(
      matrix( 0, 5, length(x) )
    )
  }else{
    mu_model    = function(x){ sin(4 * pi * x) * exp(-3 * x) }
    sigma_model = function(x){ ((1 - x - 0.4)^2 + 1) / 6 }
    noise_model = RandomNormalSum


    # Get the true values for the models considered in this simulation
    trueValue = rbind(
      matrix( 0, 5, length(x) )
    )
  }

  for( l in (1:length(transformationvec)) ){
    bias.l           = biasvec[l]
    transformation.l = f ~ (( mu3 - 3*mu2*mu1 + 2*mu1^3 ) / ( mu2 - mu1^2 )^(3/2))^(1/3)
    moments.l = momentsvec
    trueValue.l      = trueValue[l,]
    se.est.l         = se.est[l]
    cov = covering_scb( Msim   = Msim,
                        N      = Nvec,
                        level  = level,
                        method = methvec,
                        x      = x,
                        mu     = mu_model,
                        sigma  = sigma_model,
                        noise  = noise_model,
                        transformation = transformation.l,
                        moments = moments.l,
                        bias   = bias.l,
                        se.est = se.est.l,
                        trueValue = trueValue.l )
    save( list = c("cov", "Msim", "Nvec", "level", "methvec", "x", "mu_model", "sigma_model",
                   "noise_model", "transformation.l", "bias.l", "se.est.l", "trueValue.l"),
          file = paste( path_data,
                        date, "_",
                        Model, "_",
                        transformation.l,
                        "_bias_", bias.l,
                        "_se.est_", se.est.l,
                        "_", sim_tag,
                        ".rdata",
                        sep = "" ) )
  }


  Y1 <- SignalPlusNoise( N  = 8e2,
                         x  = x,
                         mu = mu_model,
                         noise = noise_model,
                         sigma = sigma_model )

  test <- scb_moments(  Y1$values,
                            level          = .95,
                            transformation = transformation.l,
                            moments        = moments.l,
                            method  = list(name = "tGKF", field = "t"),
                            bias    = TRUE,
                            se.est  = TRUE,
                            coords  = NULL,
                            mask    = NULL )
  plot_scb(test)
