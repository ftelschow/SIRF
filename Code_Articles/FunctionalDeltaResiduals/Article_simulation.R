#-------------------------------------------------------------------------------
#
#       Script simulating the covering rate for SCBs with Delta residuals
#
#-------------------------------------------------------------------------------
# Preapre workspace
#-------------------------------------------------------------------------------
Article_simulation <- function(Model   = "ModelA", # "ModelB",  "ModelC",
                               transformation = "skewness",
                               Msim    = 5e3, # number of simulations
                               Nvec    = c(20, 50, 100, 200, 400, 800), # sample sizes considered
                               x       = seq(0, 1, length.out = 150), # locations the process is evaluated at
                               level   = .95,  # level of simultaneous control
                               sim_tag = "",
                               path_wd = "~/Seafile/Projects/2019_DeltaResiduals/",
                               date    = "YEAR_MO_DY",... ){
  # Load packages
  require(SIRF)
  require(tidyverse)
  require(SampleFields)

  # General constants
  path_data <- paste( path_wd, "Workspaces/", sep = "" )

  #-------------------------------------------------------------------------------
  # General Simulation parameters
  #-------------------------------------------------------------------------------
  # Describing the parameters of the compared methods
  GKF    = list( name = "GKF", field = "z" )
  tGKF   = list( name = "tGKF", field = "t" )
  Mult   = list( name = "MultBoot", Mboots = 5e3, weights = "gaussian", method = "regular" )
  tMult  = list( name = "MultBoot", Mboots = 5e3, weights = "gaussian", method = "t"  )
  rMult  = list( name = "MultBoot", Mboots = 5e3, weights = "rademacher", method = "regular"  )
  trMult = list( name = "MultBoot", Mboots = 5e3, weights = "rademacher", method = "t"  )

  #methvec = list(GKF = GKF, tGKF = tGKF,
  #               Mult = Mult, tMult = tMult,
  #               rMult = rMult, trMult = trMult )
  methvec = list( GKF = GKF )

  if(transformation %in% c("skewness", "kurtosis", "cohensd")){
    se.est <- c("estimate", "exact gaussian")
  }else if(transformation %in% c("skewness (normality)",
                                 "kurtosis (normality)",
                                 "Ksquare (normality)")){
    se.est <- c("estimate", 1)
  }else if(transformation == "Ksquare"){
    se.est <- c("estimate", 2)
  }


  biasvec = c("asymptotic gaussian", "estimate")

  #-------------------------------------------------------------------------------
  # Models
  #-------------------------------------------------------------------------------
  if( Model == "ModelB" ){
    mu_model    = function(x){ ( x - 0.3 )^2 }
    sigma_model = function(x){ (sin(3 * pi * x) + 1.5) / 6 }

    covf <- function(x, y) SampleFields::covf.nonst.matern(x,
                                                           y,
                                                           params = c( 1, 1 / 4, 0.4 ))
    noise_model = function(N, x){ ArbCovProcess( N, x, covf = covf ) }

    if(transformation == "cohensd"){
      mu_model(x) / sigma_model(x) / vapply(x, function(x) sqrt(covf(x,x)), 1)
    }else if(transformation %in% c("skewness",
                                   "skewness (normality)",
                                   "kurtosis",
                                   "kurtosis (normality)",
                                   "Ksquare",
                                   "Ksquare (normality)"
                                   )){
      trueValue <- rep(0, length(x))
    }

  }else if( Model == "ModelC" ){
    mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
    sigma_model = function(x){(1.5 - x)}
    noise_model = DegrasNonGaussProcess

    f <- function(x) sqrt(2) / 6 * sin(pi * x)
    g <- function(x) 2 / 3 * (x - 0.5)#

    if(transformation == "cohensd"){
      mu_model(x) / sigma_model(x)
    }else if(transformation == "skewness"){
      trueValue <- (8 * f(x)^3 + 2 * g(x)^3) / (2 * f(x)^2 + g(x)^2)^(3/2)
    }else if(transformation == "skewness (normality)"){
      trueValue <- Z1_transform((8 * f(x)^3 + 2 * g(x)^3) / (2 * f(x)^2 + g(x)^2)^(3/2), Inf)
    }else if(transformation == "kurtosis"){
      trueValue = (60 * f(x)^4 + 12 * f(x)^2 * g(x)^2  + 9 * g(x)^4) /
                      (2 * f(x)^2 + g(x)^2)^2
    }else if(transformation == "kurtosis (normality)"){
      trueValue = Z2_transform((60 * f(x)^4 + 12 * f(x)^2 * g(x)^2  + 9 * g(x)^4) /
                      (2 * f(x)^2 + g(x)^2)^2, Inf)
    }else if(transformation == "Ksquare"){
      trueValue = Z2_transform((60 * f(x)^4 + 12 * f(x)^2 * g(x)^2  + 9 * g(x)^4) /
                                 (2 * f(x)^2 + g(x)^2)^2, Inf)^2 +
        Z1_transform((8 * f(x)^3 + 2 * g(x)^3) / (2 * f(x)^2 + g(x)^2)^(3/2), Inf)^2
    }else if(transformation == "Ksquare (normality)"){
      trueValue =  WF_transform(Z2_transform((60 * f(x)^4 + 12 * f(x)^2 * g(x)^2  + 9 * g(x)^4) /
                                  (2 * f(x)^2 + g(x)^2)^2, Inf)^2 +
        Z1_transform((8 * f(x)^3 + 2 * g(x)^3) / (2 * f(x)^2 + g(x)^2)^(3/2), Inf)^2)
    }

  }else{
    mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
    sigma_model = function(x){((1 - x - 0.4)^2 + 1) / 6}
    noise_model = RandomNormalSum

    if(transformation == "cohensd"){
      mu_model(x) / sigma_model(x)
    }else if(transformation %in% c("skewness",
                                      "skewness (normality)",
                                      "kurtosis",
                                      "kurtosis (normality)",
                                      "Ksquare",
                                      "Ksquare (normality)"
    )){
      trueValue <- rep(0, length(x))
    }
  }
    count = 0
    cov <- list()
    for(bias.est.l in biasvec){
      for(se.est.l in se.est){
        count = count + 1
        cov[[count]] = covering_scb( Msim   = Msim,
                            N      = Nvec,
                            level  = level,
                            method = methvec,
                            x      = x,
                            mu     = mu_model,
                            sigma  = sigma_model,
                            noise  = noise_model,
                            transformation = transformation,
                            bias.est = bias.est.l,
                            se.est   = se.est.l,
                            trueValue = trueValue )
        save( list = c("cov", "Msim", "Nvec", "level", "methvec", "x", "mu_model", "sigma_model",
                       "noise_model", "transformation", "bias.est.l", "se.est.l", "trueValue"),
              file = paste( path_data,
                            date, "_",
                            Model, "_",
                            gsub(" ", "_", transformation),
                            "_bias_", gsub(" ", "_", bias.est.l),
                            "_se.est_", gsub(" ", "_", se.est.l),
                            "_", sim_tag,
                            ".rdata",
                            sep = "" ) )
      }
    }
  return(cov)
}
