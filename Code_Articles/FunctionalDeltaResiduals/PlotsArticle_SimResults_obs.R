# Clear workspace
rm( list = ls() )

# Required packages
require(RFT)
require(SampleFields)
require(tidyverse)
require(gaussplotR)
require(fields)
require(reshape2)
require(latex2exp)

#
path_wd   <- "~/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"
path_pics <- "~/Seafile/Projects/2019_DeltaResiduals/Article/Figures/"
# Set the working path
setwd(path_wd)

# Get the simulation file from source
#source(paste(path_wd, "Workspaces/Variance_Simulation.Rdata", sep = ""))

date  = "2021_11_19"
date  = "2022_01_20"

#-------------------------------------------------------------------------------
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

theme1 <- theme(legend.position = "none",
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

theme2 <- theme(legend.position = c(0.8, 0.35),
                legend.title = element_text(face = "bold", size = sLegend),
                legend.key.size = unit(1.1, 'cm'),
                legend.text  = element_text(size = sLegend),
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

#-------------------------------------------------------------------------------
plot_results <- function(Model, transformation.l, bias.est.l, se.est.l, obs, ylims = c(0.7, 1), legend = F){
  # Plot the simulation results
  pngname <- paste( path_pics,
                    Model, "_",
                    gsub( "[()]", "", gsub(" ", "_", transformation.l)),
                    "_bias_", gsub(" ", "_", bias.est.l),
                    "_se_", gsub(" ", "_", se.est.l),
                    "_obs_", 100*obs,
                    ".png",
                    sep = "" )
  # N vector
  nx_vec = c(50, 100, 175)

  if(transformation.l == "cohensd"){
    title = r'(\textbf{Cohen's d, })'
  }else if(transformation.l == "skewness"){
    title = r'(\textbf{Skewness, })'
  }else if(transformation.l == "skewness (normality)"){
    title = r'(\textbf{Transf. Skewness, })'
  }else if(transformation.l == "kurtosis"){
    title = r'(\textbf{Kurtosis, })'
  }else if(transformation.l == "kurtosis (normality)"){
    title = r'(\textbf{Transf. Kurtosis, })'
  }

  if(se.est.l == "exact gaussian"){
    tt = r'(\textbf{Gaussian S.E.})'
  }else if(se.est.l == "1"){
    tt = r'(\textbf{Gaussian S.E.})'
  }else{
    tt = r'(\textbf{Estimated S.E.})'
  }

  title = TeX(paste( title, tt))
  # Get the data
  covRate <- NULL
  for(nx in nx_vec){
    sim_Name <- paste(Model, "_", transformation.l,
                      "_bias_", gsub(" ", "_", bias.est.l),
                      "_se_", gsub(" ", "_", se.est.l),
                      "_nx_", nx,
                      "_obs_", obs,  sep = "")
    sim_Name <- paste(Model, "_",
                      gsub( "[()]", "", gsub(" ", "_", transformation.l)),
                      "_bias_", gsub(" ", "", bias.est.l),
                      "_seEst_", gsub(" ", "_", se.est.l),
                      "_x_", nx,
                      "_obs_", 100*obs, sep = "")
    load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                        date, "_",
                        sim_Name,
                        "_ALL",
                        ".rdata",sep = "") )

    covRates_nx <- cbind(cov$rates, nx)
    colnames(covRates_nx) <- c(colnames(cov$rates), "T")

    covRate <- rbind(covRate, covRates_nx)
  }
  lvl <- cov$level
  Msim    <- cov$Msim
  covRates <- covRate
  rm(covRate)


  target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
  xLab <- rownames(covRates)
  yLab <- colnames(covRates)
  covs <- as_tibble(covRates, rownames = "N") %>%
    melt(., id.vars = c("N", "T"), variable = "Method", value.name = "CovRate") %>%
    as_tibble() %>% mutate_if(is.character, as.numeric) %>% mutate_at("T", as.factor)

  # Plot the Covering Rates by Method
  X11(width = wvalue, height = hvalue)
  if(legend){
    print(ggplot(covs, aes(N, CovRate, shape = Method, col = T, group = interaction(Method, T))) +
            xlab( "Sample Size [N]" ) +
            ylab( "Covering Rate" ) +
            geom_hline( yintercept = target[2], size = sLine ) +
            geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
            geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
            geom_point(size = sPch) + geom_line(size = sLine) +
            ggtitle( title ) +
            theme2 +
            coord_cartesian(ylim = ylims))
  }else{
    print(ggplot(covs, aes(N, CovRate, shape = Method, col = T, group = interaction(Method, T))) +
            xlab( "Sample Size [N]" ) +
            ylab( "Covering Rate" ) +
            geom_hline( yintercept = target[2], size = sLine ) +
            geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
            geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
            geom_point(size = sPch) + geom_line(size = sLine) +
            ggtitle( title ) +
            theme1 +
            coord_cartesian(ylim = ylims))
  }
    savePlot(filename = pngname)
}


#-------------------------------------------------------------------------------
# No Bias / Estimate plot
bias.est.l = "asymptotic gaussian"
se.est.l   = "estimate"
transformation.l = "skewness"

for(Model in c("ModelA", "ModelB")){
  for(obs in c(0.05, 0.1)){
    plot_results (Model, transformation.l, bias.est.l, se.est.l, obs, ylims = c(0.6, 1), legend = F)
    dev.off()
  }
}

bias.est.l = "asymptotic gaussian"
se.est.l   = "exact gaussian"
transformation.l = "skewness"

for(Model in c("ModelA", "ModelB")){
  for(obs in c(0.05, 0.1)){
    plot_results (Model, transformation.l, bias.est.l, se.est.l, obs, ylims = c(0.6, 1), legend = T)
    dev.off()
  }
}

#-------------------------------------------------------------------------------
# No Bias / Estimate plot
bias.est.l = "asymptotic gaussian"
se.est.l   = "estimate"
transformation.l = "cohensd"

for(Model in c("ModelA", "ModelB")){
  for(obs in c(0.05, 0.1)){
    plot_results (Model, transformation.l, bias.est.l, se.est.l, obs, ylims = c(0.6, 1), legend = F)
    dev.off()
  }
}

#-------------------------------------------------------------------------------
# No Bias / Estimate plot
bias.est.l = "asymptotic gaussian"
se.est.l   = "1"# "estimate"
transformation.l = "skewness (normality)"

for(Model in c("ModelA", "ModelB")){
  for(obs in c(0.05, 0.1)){
    plot_results (Model, transformation.l, bias.est.l, se.est.l, obs, ylims = c(0.6, 1), legend = F)
    dev.off()
  }
}

bias.est.l = "asymptotic gaussian"
se.est.l   = "estimate"
transformation.l = "skewness (normality)"

for(Model in c("ModelA", "ModelB")){
  for(obs in c(0.05, 0.1)){
    plot_results (Model, transformation.l, bias.est.l, se.est.l, obs, ylims = c(0.6, 1), legend = F)
    dev.off()
  }
}
