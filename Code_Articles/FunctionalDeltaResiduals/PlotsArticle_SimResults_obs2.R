# Clear workspace
rm( list = ls() )

# Required packages
require(RFT)
require(SampleFields)
require(tidyverse)
require(gaussplotR)
require(fields)
require(reshape2)

#
path_wd   <- "~/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"
path_pics <- "~/Seafile/Projects/2019_DeltaResiduals/Article/Figures/"
# Set the working path
setwd(path_wd)

date  = "2022_01_11"

#-------------------------------------------------------------------------------
plot_covResults <- function(model, transformation, bias.l, se.est.l, legend.on = FALSE){
  sim_Name <- paste(Model, "_",
                    gsub( "[()]", "", gsub(" ", "_", transformation)),
                    "_bias_", gsub(" ", "", bias.est.l),
                    "_seEst_", gsub(" ", "_", se.est.l),
                    "_obs_", 100*obs, sep = "")
  
  if(bias.l == "asymptotic gaussian"){
    title = "No Bias correction,"
  }else{
    title = "Bias correction,"
  }
  
  if(se.est.l == "exact gaussian"){
    title = paste(title, "Gaussian variance")
  }else if(se.est.l == "1"){
    title = paste(title, "Gaussian variance")
  }else{
    title = paste(title, "Estimated variance")
  }
  
  if( transformation %in% c("kurtosis", "kurtosis (normality)") ){
    ylow = 0.65
    yup  = 1
    if(model == "ModelC"){
      ylow = 0.4
      yup  = 1
    }
  }else if( transformation %in% c("skewness", "skewness (normality)") ){
    ylow = 0.8
    yup  = 1
    if(model == "ModelC"){
      ylow = 0.4
      yup  = 1
    }
  }else{
    ylow = 0.85
    yup  = 1
  }
  
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
  
  
  if(legend.on){
    themeP = theme2
  }else{
    themeP = theme1
  }
  
  # Collecting the data
  covRate <- NULL
  for(nx in nx_vec){
    load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                        date, "_",
                        sim_Name,
                        Model, "_",
                        gsub( "[()]", "", gsub(" ", "_", transformation)),
                        "_bias_", gsub(" ", "", bias.est.l),
                        "_seEst_", gsub(" ", "_", se.est.l),
                        "_x_", nx,
                        "_obs_", 100*obs,
                        "_ALL",
                        ".rdata",sep = "") )
    
    covRates_nx <- cbind(covRates, nx)
    colnames(covRates_nx) <- c(colnames(covRates), "T")
    
    covRate <- rbind(covRate, covRates_nx)
  }
  
  covRates <- covRate
  rm(covRate)
  
  # Plot the simulation results
  pngname <- paste( path_pics,
                    sim_Name,
                    ".png",
                    sep = "" )
  Msim <- cov$Msim
  level <- cov$level
  #png( pngname, width = 550, height = 450 )
  target <- sqrt( level * ( 1 - level ) / Msim ) * c(-qnorm(level/2 + 1/2), 0, qnorm(level/2 + 1/2)) + level
  xLab <- rownames(cov$rates)
  yLab <- colnames(cov$rates)
  covs <- as_tibble(cov$rates, rownames = "N") %>%
    melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
    as_tibble() %>% mutate_if(is.character, as.numeric)
  
  X11(width = 1.2*7, height = 7)
  # Plot the Covering Rates by Method
  print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
           xlab( "Sample Size [N]" ) +
           ylab( "Covering Rate" ) +
           ggtitle( title ) +
           geom_hline( yintercept = target[2], size = sLine ) +
           geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
           geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
           geom_point(size = sPch) + geom_line(size = sLine) +
           themeP +
           coord_cartesian(ylim = c(ylow, yup ))
  )
  savePlot(filename = pngname)
  dev.off()
}

#-------------------------------------------------------------------------------
for(model in c("ModelA", "ModelB", "ModelC")){
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
    
    for(bias.l in c("asymptotic gaussian", "estimate")){
      for(se.est.l in se.est){
        # Legend plot or not
        if(bias.l == "asymptotic gaussian" & (se.est.l == "exact gaussian" | se.est.l == "1")){
          legend.on = TRUE
        }else{
          legend.on = FALSE
        }
        # Plot the data
        plot_covResults(model, transformation, bias.l, se.est.l, legend.on)
      }
    }
  }
}