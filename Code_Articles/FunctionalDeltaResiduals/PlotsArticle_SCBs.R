#-------------------------------------------------------------------------------
#
# Plot the SCBs
#
#-------------------------------------------------------------------------------
# Clear workspace
rm( list = ls() )
library(SIRF)
library(SampleFields)
library(RFT)
library(tidyverse)
# Path to save the figures
pic_path <- "/home/fabian/Seafile/Projects/2019_DeltaResiduals/Article/Figures/"

set.seed(666)

wvalue = 1.22*7.3
hvalue = 6
sLegend <- 25
sText   <- 30
sTitle  <- 35
sLegend <- 35
Sylab = 1.5
Sxlab = 1.5
sLine = 1.5
sPch = 4

wvalue = 1.3*7
hvalue = 7
sLineMean = 3

theme1 <- theme(legend.position = c(0.55, 0.88),
                legend.title = element_blank(),
                legend.text  = element_text(color = "black", size = sLegend, face = "plain"),
                legend.direction = "horizontal",
                legend.key = element_rect(colour = "transparent", fill = "transparent"),
                legend.background = element_blank(),
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

theme2 <- theme(legend.position = "none",
                legend.title = element_blank(),
                legend.text  = element_text(color = "black", size = sText, face = "plain"),
                legend.direction = "horizontal",
                legend.key = element_rect(colour = "transparent", fill = "transparent"),
                legend.background = element_blank(),
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

#-------------------------------------------------------------------------------
plot_scbExamples <- function(Model, transformation, method = list(name = "GKF", field = "z"), legend.on = T){
  x = seq(0, 1, length.out = 100)
  if(Model == "A"){
    # Model A
    mu_model    = function( x ){ sin( 4 * pi * x ) * exp( -3 * x ) }
    sigma_model = function( x ){ ( ( 1 - x - 0.4 )^2 + 1 ) / 6 }
    noise_model = RandomNormalSum

    if(transformation == "cohensd"){
      trueValue <- mu_model(x) / sigma_model(x)
    }else if(transformation %in% c("skewness",
                                   "skewness (normality)",
                                   "kurtosis",
                                   "kurtosis (normality)",
                                   "Ksquare",
                                   "Ksquare (normality)"
    )){
      trueValue <- rep(0, length(x))
    }
  }else if(Model == "B"){
    # Model B
    mu_model    = function( x ){ ( x-0.3 )^2 }
    sigma_model = function( x ){ (sin(3*pi*x) + 1.5) / 6 }
    covf <- function(x, y) SampleFields::covf.nonst.matern(x,
                                                           y,
                                                           params = c( 1, 1 / 4, 0.4 ))
    noise_model = function(N, x){ ArbCovProcess( N, x, covf = covf ) }

    if(transformation == "cohensd"){
      trueValue <- mu_model(x) / sigma_model(x) / vapply(x, function(x) sqrt(covf(x,x)), 1)
    }else if(transformation %in% c("skewness",
                                   "skewness (normality)",
                                   "kurtosis",
                                   "kurtosis (normality)",
                                   "Ksquare",
                                   "Ksquare (normality)"
    )){
      trueValue <- rep(0, length(x))
    }
  }else if(Model == "C"){
    # Model C
    mu_model    = function( x ){  sin( 4 * pi * x ) * exp( -3 * x ) }
    sigma_model = function( x ){ (1.5-x) }
    noise_model = DegrasNonGaussProcess

    f <- function(x) sqrt(2) / 6 * sin(pi * x)
    g <- function(x) 2 / 3 * (x - 0.5)#
    skewV <- (8 * f(x)^3 + 2 * g(x)^3) / (2 * f(x)^2 + g(x)^2)^(3/2)
    kurtV <- (60 * f(x)^4 + 12 * f(x)^2 * g(x)^2  + 9 * g(x)^4) /
      (2 * f(x)^2 + g(x)^2)^2 - 3

    if(transformation == "cohensd"){
      trueValue <- mu_model(x) / sigma_model(x)
    }else if(transformation == "skewness"){
      trueValue <- skewV
    }else if(transformation == "skewness (normality)"){
      trueValue <- Z1_transform(skewV, Inf)
    }else if(transformation == "kurtosis"){
      trueValue = kurtV
    }else if(transformation == "kurtosis (normality)"){
      trueValue = Z2_transform(kurtV + 3, Inf)
    }else if(transformation == "Ksquare"){
      trueValue = Z2_transform(kurtV, Inf)^2 + Z1_transform(skewV, Inf)^2
    }else if(transformation == "Ksquare (normality)"){
      trueValue =  WF_transform(Z2_transform(kurtV + 3, Inf)^2 + Z1_transform(skewV, Inf)^2)
    }
  }

  # define the titel
  if(transformation == "cohensd"){
    title = paste("Model ", Model, ": Cohen's d", sep = "")
  }else if(transformation == "skewness"){
    title = paste("Model ", Model, ": Skewness", sep = "")
  }else if(transformation == "skewness (normality)"){
    title = paste("Model ", Model, ": Transf. Skewness", sep = "")
  }else if(transformation == "kurtosis"){
    title = paste("Model ", Model, ": Kurtosis", sep = "")
  }else if(transformation == "kurtosis (normality)"){
    title = paste("Model ", Model, ": Transf. Kurtosis", sep = "")
  }


  # Get a data sample
  Y <- SignalPlusNoise( N  = 400,
                        x  = x,
                        mu = mu_model,
                        sigma = sigma_model,
                        noise = noise_model )

  Y1 <- Y
  Y1$values <- Y$values[,1:100]
  scb_ModelA = scb_moments(Y1,
                           level          = .95,
                           transformation = transformation,
                           method  = method )

  Y1 <- Y
  Y1$values <- Y$values[,1:400]
  scb_ModelA2 = scb_moments(Y1,
                            level          = .95,
                            transformation = transformation,
                            method  = method )

  data = cbind(x, scb_ModelA$scb)
  data[,"est"] <- trueValue
  colnames(data)[3] <- "True"
  data = as_tibble(data)
  data2 = cbind(x, scb_ModelA2$scb)
  data[,"est"] <- trueValue
  colnames(data)[3] <- "True"
  data2 = as_tibble(data2)

  pngname <- paste( pic_path, "Model", Model, "_SCBs_", gsub( "[()]", "", gsub(" ", "_", transformation)),"_example.pdf", sep = "" )
  X11(width = wvalue, height = hvalue)
  if(legend.on){
    print(ggplot(data) +
            geom_line(aes(x = x, y = True), size = 1.2) +
            geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = "N = 100"), alpha = 0.25) +
            geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = "N = 400"), data = data2, alpha = 0.25) +
            xlab( "" ) + ylab( "" ) +
            ggtitle(title) +
            theme1)
  }else{
    print(ggplot(data) +
            geom_line(aes(x = x, y = True), size = 1.2) +
            geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = "N = 100"), alpha = 0.25) +
            geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = "N = 400"), data = data2, alpha = 0.25) +
            xlab( "" ) + ylab( "" ) +
            ggtitle(title) +
            theme2)
  }

  savePlot(filename = pngname)
  dev.off()

}

#-------------------------------------------------------------------------------

for(Model in c("A", "B", "C")){
  for(transformation in c("cohensd", "skewness", "skewness (normality)")){
    if(Model == "A"){
      method = list(name = "GKF", field = "z")
    }else{
      method   = list( name = "MultBoot", Mboots = 5e3, weights = "gaussian", method = "regular" )
    }
    if(Model == "B" & transformation == "cohensd"){
      leg = T
    }else{
      leg = F
    }
    plot_scbExamples(Model = Model, transformation = transformation, method = method, legend.on = leg)
  }
}


