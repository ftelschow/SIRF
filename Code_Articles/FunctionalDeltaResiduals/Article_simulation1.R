#-------------------------------------------------------------------------------
#
#       Script simulating the covering rate for SCBs with Delta residuals
#
#-------------------------------------------------------------------------------
# Preapre workspace
#-------------------------------------------------------------------------------

Article_simulation1 <- function(Model   = "ModelA", # "ModelB",  "ModelC"
                                Msim    = 5e3, # number of simulations
                                Nvec    = c(50, 100, 200, 400, 800), # sample sizes considered
                                x       = seq(0, 1, length.out = 150), # locations the process is evaluated at
                                level   = .95,  # level of simultaneous control
                                sim_tag = "",
                                path_wd = "~/Seafile/Projects/2019_DeltaResiduals/",
                                date    = "YEAR_MO_DY" ){
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

  methvec = list( GKF = GKF, tGKF = tGKF,
                  Mult = Mult, tMult = tMult,
                  rMult = rMult, trMult = trMult )
  #methvec = list( GKF = GKF, tGKF = tGKF )

  #-------------------------------------------------------------------------------
  # Model A
  #-------------------------------------------------------------------------------
  if( Model == "ModelB" ){
    mu_model    = function(x){ ( x-0.3 )^2 }
    sigma_model = function(x){ (sin(3*pi*x) + 1.5) / 6 }

    covf <- function(x, y) SampleFields::covf.nonst.matern(x,
                                                           y,
                                                           params = c( 1, 1 / 4, 0.4 ))
    noise_model = function(N, x){ ArbCovProcess( N, x, covf = covf ) }

    transformationvec = c( "cohensd",
                           "cohensd",
                           "cohensd",
                           "cohensd",
                           "skewness",
                           "skewness",
                           "skewness" )
    se.est <- c( TRUE,
                 TRUE,
                 FALSE,
                 FALSE,
                 TRUE,
                 TRUE,
                 FALSE )

    biasvec = c( TRUE,
                 FALSE,
                 TRUE,
                 FALSE,
                 TRUE,
                 FALSE,
                 FALSE )

    trueValue = rbind(
      t(matrix(mu_model(x) / sigma_model(x) / vapply(x, function(x) sqrt(covf(x,x)), 1), length(x), sum(transformationvec=="cohensd"))),
      matrix(0, sum(transformationvec == "skewness"), length(x))
    )

  }else if( Model == "ModelC" ){
    mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
    sigma_model = function(x){(1.5 - x)}
    noise_model = DegrasNonGaussProcess

    transformationvec = c( "cohensd",
                           "cohensd",
                           "skewness",
                           "skewness" )
    se.est <- c(TRUE,
                TRUE,
                TRUE,
                TRUE)

    biasvec = c(TRUE,
                FALSE,
                TRUE,
                FALSE)

    f <- function(x) sqrt(2) / 6 * sin(pi * x)
    g <- function(x) 2 / 3 * (x - 0.5)
    skew_thyC <- (8 * f(x)^3 + 2 * g(x)^3) / (2*f(x)^2 + g(x)^2)^(3/2)

    trueValue = rbind(
      t(matrix(mu_model(x) / sigma_model(x), length(x), sum(transformationvec=="cohensd"))),
      t(matrix(skew_thyC, length(x), sum(transformationvec=="cohensd")))
    )
    rm( skew_true, skewnessMC, f )
  }else{
    mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
    sigma_model = function(x){((1 - x - 0.4)^2 + 1) / 6}
    noise_model = RandomNormalSum

    # Get the true values for the models considered in this simulation
    transformationvec = c( "cohensd",
                           "cohensd",
                           "cohensd",
                           "cohensd",
                           "skewness",
                           "skewness",
                           "skewness" )
    se.est <- c( TRUE,
                 TRUE,
                 FALSE,
                 FALSE,
                 TRUE,
                 TRUE,
                 FALSE )

    biasvec = c( TRUE,
                 FALSE,
                 TRUE,
                 FALSE,
                 TRUE,
                 FALSE,
                 FALSE )

    trueValue = rbind(
      t(matrix(mu_model(x) / sigma_model(x), length(x), sum(transformationvec=="cohensd"))),
      matrix(0, sum(transformationvec == "skewness"), length(x))
    )
  }



  for( l in (1:length(transformationvec)) ){
    bias.l           = biasvec[l]
    transformation.l = transformationvec[l]
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
}
