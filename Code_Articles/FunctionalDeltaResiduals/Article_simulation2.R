#-------------------------------------------------------------------------------
#
#       Script simulating the covering rate for SCBs with Delta residuals
#
#-------------------------------------------------------------------------------
# Preapre workspace
#-------------------------------------------------------------------------------

Article_simulation2 <- function(Model   = "ModelA", # "ModelB",  "ModelC"
                                Msim    = 5e3, # number of simulations
                                Nvec    = c(50, 100, 200, 400, 800), # sample sizes considered
                                x       = seq(0, 1, length.out = 175), # locations the process is evaluated at
                                level   = .95,  # level of simultaneous control
                                obs     = 0.05,
                                sim_tag = "",
                                path_wd = "~/Seafile/Projects/2019_DeltaResiduals/",
                                date    = "YEAR_MO_DY",... ){
  # General constants
  path_data <- paste(path_wd, "Schreibtisch/", sep = "")#"Workspaces/", sep = "")

  #-------------------------------------------------------------------------------
  # General Simulation parameters
  #-------------------------------------------------------------------------------
  # Describing the parameters of the compared methods
  GKF    = list(name = "GKF", field = "z")
  tGKF   = list(name = "GKF", field = "t")
  Mult   = list(name = "MultBoot", Mboots = 5e3, weights = "gaussian", method = "regular")
  tMult  = list(name = "MultBoot", Mboots = 5e3, weights = "gaussian", method = "t")
  rMult  = list(name = "MultBoot", Mboots = 5e3, weights = "rademacher", method = "regular")
  trMult = list(name = "MultBoot", Mboots = 5e3, weights = "rademacher", method = "t")

  methvec = list(tGKF = tGKF, Mult = Mult)
  #methvec = list( GKF = GKF, tGKF = tGKF )
  xeval = seq(0, 1, length.out = 175)
  #-------------------------------------------------------------------------------
  # Model A
  #-------------------------------------------------------------------------------
  if( Model == "ModelB" ){
    mu_model    = function(x){ (x - 0.3)^2 }
    sigma_model = function(x){ (sin(3 * pi * x) + 1.5) / 6 }

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
    se.est <- c( "estimate",
                 "estimate",
                 "exact gaussian",
                 "exact gaussian",
                 "estimate",
                 "estimate",
                 "exact gaussian" )

    biasvec = c( "estimate",
                 "asymptotic gaussian",
                 "estimate",
                 "asymptotic gaussian",
                 "estimate",
                 "asymptotic gaussian",
                 "asymptotic gaussian" )

    trueValue = rbind(
      t(matrix(mu_model(xeval) / sigma_model(xeval) / vapply(xeval, function(x) sqrt(covf(x,x)), 1), length(x), sum(transformationvec=="cohensd"))),
      matrix(0, sum(transformationvec == "skewness"), length(xeval))
    )

  }else if( Model == "ModelC" ){
    mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
    sigma_model = function(x){(1.5 - x)}
    noise_model = DegrasNonGaussProcess

    transformationvec = c( "cohensd",
                           "cohensd",
                           "skewness",
                           "skewness" )
    se.est <- c("estimate",
                "estimate",
                "estimate",
                "estimate")

    biasvec = c("estimate",
                "asymptotic gaussian",
                "estimate",
                "asymptotic gaussian")

    f <- function(x) sqrt(2) / 6 * sin(pi * x)
    g <- function(x) 2 / 3 * (x - 0.5)
    skew_thyC <- (8 * f(xeval)^3 + 2 * g(xeval)^3) / (2 * f(xeval)^2 + g(xeval)^2)^(3/2)

    trueValue = rbind(
      t(matrix(mu_model(xeval) / sigma_model(xeval), length(xeval), sum(transformationvec=="cohensd"))),
      t(matrix(skew_thyC, length(xeval), sum(transformationvec=="skewness")))
    )
    rm( g, f )
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
    se.est <- c( "estimate",
                 "estimate",
                 "exact gaussian",
                 "exact gaussian",
                 "estimate",
                 "estimate",
                 "exact gaussian" )

    biasvec = c( "estimate",
                 "asymptotic gaussian",
                 "estimate",
                 "asymptotic gaussian",
                 "estimate",
                 "asymptotic gaussian",
                 "asymptotic gaussian")

    trueValue = rbind(
      t(matrix(mu_model(xeval) / sigma_model(xeval), length(xeval), sum(transformationvec == "cohensd"))),
      matrix(0, sum(transformationvec == "skewness"), length(xeval))
    )
  }



  for(l in 1:length(transformationvec)){
    bias.est.l       = biasvec[l]
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
                        bias.est = bias.est.l,
                        se.est   = se.est.l,
                        smoothing = list(),
                        obs.noise = obs,
                        trueValue = trueValue.l,... )
    save( list = c("cov", "Msim", "Nvec", "level", "methvec", "x", "mu_model", "sigma_model",
                   "noise_model", "transformation.l", "bias.est.l", "se.est.l", "trueValue.l"),
          file = paste( path_data,
                        date, "_",
                        Model, "_",
                        transformation.l,
                        "_bias_", gsub(" ", "_", bias.est.l),
                        "_se_", gsub(" ", "_", se.est.l),
                        "_", sim_tag,
                        ".rdata",
                        sep = "" ) )
  }
}
