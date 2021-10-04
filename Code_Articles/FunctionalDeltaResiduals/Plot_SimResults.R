# 

# Load packages
library( SCBfun )
library( SampleFields )

#
path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SCBfun/DeltaResiduals/"
path_pics <- "/home/fabian/Seafile/Code/Rpackages/SCBfun/DeltaResiduals/Pics/"
# Set the working path
setwd(path_wd)

# Get the simulation file from source
source(paste(path_wd, "Plot_Functions.R", sep = ""))

# Define the simulation parameters

Modelvec = c("ModelA", "ModelB", "ModelC")

date   = "2021_05_31"
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
    transformationvec = c("cohensd",
                          "cohensd",
                          "cohensd",
                          "cohensd",
                          "skewness",
                          "skewness",
                          "skewness")
    
    se.est <- c(TRUE,
                TRUE,
                FALSE,
                FALSE,
                TRUE,
                TRUE,
                FALSE)
    
    biasvec = c(TRUE,
                FALSE,
                TRUE,
                FALSE,
                TRUE,
                FALSE,
                FALSE)    
  }
  
  for(l in 1:length(transformationvec)){
    transformation.l = transformationvec[l]
    bias.l = biasvec[l]
    se.est.l = se.est[l]
  
    for(n in 1:simMax){
      # Simulate the covering rate
      load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                          date, "_",
                          model, "_",
                          transformation.l,
                          "_bias_", bias.l,
                          "_se.est_", se.est.l,
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
  colnames(covAll) <- c("GKF", "tGKF", "gMult", "tgMult", "rMult", "rtMult")
  
  # Plot the simulation results
  pngname <- paste( path_pics,
                    model, "_",
                    transformation.l,
                    "_bias_", bias.l,
                    "_se.est_", se.est.l,
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
                      "_bias_", bias.l,
                      "_se.est_", se.est.l,
                      "_full",
                      ".rdata",
                      sep = "" ) )
}}

