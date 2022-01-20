# Clear workspace
rm(list = ls())

# Required packages
require(RFT)
require(SIRF)
require(SampleFields)
require(tidyverse)
require(gaussplotR)
require(fields)
require(reshape2)

#
path_wd   <- "~/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"
path_data <- "~/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/Workspaces/"
path_pics <- "~/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/Pics/"

# Set the working path
setwd(path_wd)

# Get the simulation file from source
#source(paste(path_wd, "Workspaces/Variance_Simulation.Rdata", sep = ""))

date  = "2022_01_11"
simMax = 20

for(Model in c("ModelA", "ModelB", "ModelC") ){
  for(transformation in c("skewness", "skewness (normality)", "kurtosis", "kurtosis (normality)")){

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

      if( transformation == "skewness" ){
        ylow = 0.75
        yup  = 1
      }else{
        ylow = 0.85
        yup  = 1
      }

      #-------------------------------------------------------------------------------
      sLegend <- 25
      sText   <- 25
      Sylab = 1.5
      Sxlab = 1.5
      sLine = 1.5
      sPch = 4

      theme1 <- theme(legend.position = "none",
                      plot.title   = element_text(face = "bold", size = sText),
                      axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                      axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                      axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                      axis.title.y = element_text(color = "black", size = sText, face = "plain"))

      theme2 <- theme(legend.title = element_blank(),
                      legend.position = c(0.82, 0.3),
                      legend.key.size = unit(1.1, 'cm'),
                      legend.text = element_text(size = sLegend),
                      plot.title   = element_text(face = "bold", size = sText),
                      axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                      axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                      axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                      axis.title.y = element_text(color = "black", size = sText, face = "plain"))


      #-------------------------------------------------------------------------------
      #  Concatinate and plot
      #-------------------------------------------------------------------------------
      for(bias.est.l in biasvec){
        for(se.est.l in se.est){
          sim_Name <- paste(Model, "_",
                            gsub( "[()]", "", gsub(" ", "_", transformation)),
                            "_bias_", gsub(" ", "", bias.est.l),
                            "_seEst_", gsub(" ", "_", se.est.l), sep = "")

          for(sim_tag in 1:simMax){
            # Simulate the covering rate
            load( file = paste( path_data,
                                date, "_",
                                sim_Name,
                                "_", sim_tag,
                                ".rdata",
                                sep = "" ) )
            if(sim_tag == 1){
              covAll      <- cov
              covAll$Msim <- Msim
            }else{
              covAll$rates <- covAll$rates + cov$rates
              covAll$Msim  <- covAll$Msim + Msim
            }
          }
          cov$rates <- covAll$rates / simMax
          cov$Msim  <- covAll$Msim

          # Save the concatinated data
          save( list = c("cov", "Msim", "level"),
                file = paste( path_data,
                            date, "_",
                            sim_Name,
                            "_ALL.rdata",
                            sep = "" ) )

          # Plot the data
          pngname <- paste( path_pics, sim_Name, ".png", sep = "" )
#          png( pngname, width = 550, height = 450 )
          plot_covSim(cov, title = "", legend.position = "none", ylims = c(0.2,1))
          savePlot(filename = pngname)
#          dev.off()
        }
      }
  }
}

