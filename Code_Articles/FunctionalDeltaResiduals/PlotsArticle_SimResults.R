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
#date  = "2022_11_19"
transvec <- c("cohensd", "skewness", "skewness (normality)", "kurtosis", "kurtosis (normality)")
#transvec <-

sLegend <- 25
sText   <- 30
sTitle  <- 35
sLegend <- 25
Sylab = 1.5
Sxlab = 1.5
sLine = 1.5
sPch = 4

wvalue = 1.3*7
hvalue = 7

#-------------------------------------------------------------------------------
plot_covResults <- function(model, transformation, bias.l, se.est.l, legend.on = FALSE){
  sim_Name <- paste(model, "_",
                    gsub( "[()]", "", gsub(" ", "_", transformation)),
                    "_bias_", gsub(" ", "", bias.l),
                    "_seEst_", gsub(" ", "_", se.est.l), sep = "")

  if(bias.l == "asymptotic gaussian"){
    title = "No Bias correction,"
  }else{
    title = "Bias correction,"
  }

  if(transformation == "cohensd"){
    title = "Cohen's d,"
  }else if(transformation == "skewness"){
    title = "Skewness,"
  }else if(transformation == "skewness (normality)"){
    title = "Transf. Skewness,"
  }else if(transformation == "kurtosis"){
    title = "Kurtosis,"
  }else if(transformation == "kurtosis (normality)"){
    title = "Transf. Kurtosis,"
  }

  if(se.est.l == "exact gaussian"){
    title = paste(title, "Gaussian S.E.")
  }else if(se.est.l == "1"){
    title = paste(title, "Gaussian S.E.")
  }else{
    title = paste(title, "Estimated S.E.")
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

  theme1 <- theme(legend.position = "none",
                  plot.title   = element_text(face = "bold", size = sTitle),
                  axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                  axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                  axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                  axis.title.y = element_text(color = "black", size = sText, face = "plain"))

  theme2 <- theme(legend.title = element_blank(),
                  legend.position = c(0.82, 0.3),
                  legend.key.size = unit(1.1, 'cm'),
                  legend.text = element_text(size = sLegend),
                  plot.title   = element_text(face = "bold", size = sTitle),
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
  load( file = paste(path_wd, "Workspaces/",
                     date, "_",
                     sim_Name,
                     "_ALL.rdata",
                     sep = "" ) )

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
  nav <- colnames(cov$rates)
  nav[6] <- "rtMult"
  yLab <- nav
  colnames(cov$rates) <- nav
  rm(nav)


  covs <- as_tibble(cov$rates, rownames = "N") %>%
    melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
    as_tibble() %>% mutate_if(is.character, as.numeric)

  X11(width = wvalue, height = hvalue)
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
  for(transformation in transvec){
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
        if( se.est.l == "exact gaussian"){
          legend.on = TRUE
        }else if(model == "ModelC" & se.est.l == "estimate" & transformation == "skewness"){
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
