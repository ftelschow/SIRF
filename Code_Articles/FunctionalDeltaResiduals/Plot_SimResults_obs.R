#

# Load packages
library(SIRF)
library( SampleFields )

#
path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"
path_pics <- "/home/fabian/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/Pics/"
# Set the working path
setwd(path_wd)

# Get the simulation file from source
source(paste(path_wd, "Plot_Functions.R", sep = ""))

# Define the simulation parameters

Modelvec = "ModelA"#c("ModelA", "ModelB", "ModelC")

date   = "2021_11_19"
simMax = 20

for(model in Modelvec){
  if( model == "ModelC" ){
    transformationvec = c("cohensd",
                          "cohensd",
                          "skewness",
                          "skewness")

    se.est <- c(TRUE,
                TRUE,
                TRUE,
                TRUE)

    biasvec = c(TRUE,
                FALSE,
                TRUE,
                FALSE)
  }else{
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
  }

for(obs in c(0.05, 0.1)){for(nx in c(50, 100, 175)){
  for(l in 1:length(transformationvec)){
    transformation.l = transformationvec[l]
    bias.est.l = biasvec[l]
    se.est.l = se.est[l]

    for(n in 1:simMax){
      # Simulate the covering rate
      load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                          date, "_",
                          model, "_",
                          transformation.l,
                          "_bias_", gsub(" ", "_", bias.est.l),
                          "_se_", gsub(" ", "_", se.est.l),
                          "_nx_", nx,
                          "_obs_", obs,
                          "_", n,
                          ".rdata",
                          sep = "" ) )
      if(n == 1){
        covAll  <- cov
        MsimAll <- Msim
      }else{
        covAll  <- covAll + cov
        MsimAll <- MsimAll + Msim
      }
  }
  covAll <- covAll / simMax
  colnames(covAll) <- c("tGKF", "gMult")

  # Plot the simulation results
  pngname <- paste( path_pics,
                    model, "_",
                    transformation.l,
                    "_bias_", gsub(" ", "_", bias.est.l),
                    "_se_", gsub(" ", "_", se.est.l),
                    "_nx_", nx,
                    "_obs_", obs,
                    ".png",
                    sep = "" )
  covRates = covAll
  Msim     = MsimAll
  lvl      = level
  title = "Simultaneous Confidence Bands"
  png( pngname, width = 550, height = 450 )
  target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
  xLab <- rownames(covRates)
  yLab <- colnames(covRates)
  covs <- as_tibble(covRates, rownames = "N") %>%
    melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
    as_tibble() %>% mutate_if(is.character, as.numeric)

  # Plot the Covering Rates by Method
  print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
    geom_point() + geom_line() +
    xlab( "Sample Size [N]" ) + ylab( "Covering Rate" ) +
    ggtitle( title ) +
    geom_hline( yintercept = target[2] ) +
    geom_hline( yintercept = target[1], linetype = "dashed" ) +
    geom_hline( yintercept = target[3], linetype = "dashed" ) )
  dev.off()

  # Save the concatinated data
  save( list = c("covRates", "Msim", "lvl"), file = paste( paste(path_wd, "Workspaces/",sep = ""),
                      date, "_",
                      model, "_",
                      transformation.l,
                      "_bias_", gsub(" ", "_", bias.est.l),
                      "_se_", gsub(" ", "_", se.est.l),
                      "_nx_", nx,
                      "_obs_", obs,
                      "_full",
                      ".rdata",
                      sep = "" ) )
}}}}

