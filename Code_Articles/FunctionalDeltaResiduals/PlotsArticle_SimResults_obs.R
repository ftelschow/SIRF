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

# Get the simulation file from source
#source(paste(path_wd, "Workspaces/Variance_Simulation.Rdata", sep = ""))

date  = "2021_11_19"
model = "ModelA" #  "ModelA" #  "ModelC" #
transformation.l = "skewness" #   #
nx_vec = c(50, 100, 175)

if( transformation.l == "skewness" ){
  ylow = 0.5
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

theme2 <- theme(legend.title = element_text(face = "bold", size = sText),
                legend.key.size = unit(1.1, 'cm'),
                legend.text  = element_text(size = sLegend),
                plot.title   = element_text(face = "bold", size = sText),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))



#-------------------------------------------------------------------------------
# No Bias / Estimate plot
bias.est.l   = "asymptotic gaussian"
se.est.l = "estimate"

#---- Collecting the data
# obs = 0.1
title = "Observation Noise: sd = 0.1"
obs = 0.1
covRate <- NULL
for(nx in nx_vec){
  load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                      date, "_",
                      model, "_",
                      transformation.l,
                      "_bias_", gsub(" ", "_", bias.est.l),
                      "_se_", gsub(" ", "_", se.est.l),
                      "_nx_", nx,
                      "_obs_", obs,
                      "_full",
                      ".rdata",sep = "") )

  covRates_nx <- cbind(covRates, nx)
  colnames(covRates_nx) <- c(colnames(covRates), "T")

  covRate <- rbind(covRate, covRates_nx)
}

covRates <- covRate
rm(covRate)

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", gsub(" ", "_", bias.est.l),
                  "_se_", gsub(" ", "_", se.est.l),
                  "_obs_", obs,
                  ".png",
                  sep = "" )
png( pngname, width = 550*6.1/5, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = c("N", "T"), variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric) %>% mutate_at("T", as.factor)

# Plot the Covering Rates by Method
print(
  ggplot(covs, aes(N, CovRate, shape = Method, col = T, group = interaction(Method, T))) +
         xlab( "Sample Size [N]" ) +
         ylab( "Covering Rate" ) +
         geom_hline( yintercept = target[2], size = sLine ) +
         geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
         geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
         geom_point(size = sPch) + geom_line(size = sLine) +
         theme2 +
         coord_cartesian(ylim = c(ylow, yup ))
)
dev.off()

#---- Collecting the data
# obs = 0.05
title = "Observation Noise: sd = 0.05"
obs = 0.05
covRate <- NULL
for(nx in nx_vec){
  load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                      date, "_",
                      model, "_",
                      transformation.l,
                      "_bias_", gsub(" ", "_", bias.est.l),
                      "_se_", gsub(" ", "_", se.est.l),
                      "_nx_", nx,
                      "_obs_", obs,
                      "_full",
                      ".rdata",sep = "") )

  covRates_nx <- cbind(covRates, nx)
  colnames(covRates_nx) <- c(colnames(covRates), "T")

  covRate <- rbind(covRate, covRates_nx)
}

covRates <- covRate
rm(covRate)

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", gsub(" ", "_", bias.est.l),
                  "_se_", gsub(" ", "_", se.est.l),
                  "_obs_", obs,
                  ".png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = c("N", "T"), variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric) %>% mutate_at("T", as.factor)

# Plot the Covering Rates by Method
print(
  ggplot(covs, aes(N, CovRate, shape = Method, col = T, group = interaction(Method, T))) +
    xlab( "Sample Size [N]" ) +
    ylab( "Covering Rate" ) +
    geom_hline( yintercept = target[2], size = sLine ) +
    geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
    geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
    geom_point(size = sPch) + geom_line(size = sLine) +
    theme1 +
    coord_cartesian(ylim = c(ylow, yup ))
)
dev.off()

#-------------------------------------------------------------------------------
# No Bias / Estimate plot
bias.est.l   = "estimate"
se.est.l = "estimate"

#---- Collecting the data
# obs = 0.1
title = "Observation Noise: sd = 0.1"
obs = 0.1
covRate <- NULL
for(nx in nx_vec){
  load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                      date, "_",
                      model, "_",
                      transformation.l,
                      "_bias_", gsub(" ", "_", bias.est.l),
                      "_se_", gsub(" ", "_", se.est.l),
                      "_nx_", nx,
                      "_obs_", obs,
                      "_full",
                      ".rdata",sep = "") )

  covRates_nx <- cbind(covRates, nx)
  colnames(covRates_nx) <- c(colnames(covRates), "T")

  covRate <- rbind(covRate, covRates_nx)
}

covRates <- covRate
rm(covRate)

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", gsub(" ", "_", bias.est.l),
                  "_se_", gsub(" ", "_", se.est.l),
                  "_obs_", obs,
                  ".png",
                  sep = "" )
png( pngname, width = 550*6.1/5, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = c("N", "T"), variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric) %>% mutate_at("T", as.factor)

# Plot the Covering Rates by Method
print(
  ggplot(covs, aes(N, CovRate, shape = Method, col = T, group = interaction(Method, nx))) +
    xlab( "Sample Size [N]" ) +
    ylab( "Covering Rate" ) +
    geom_hline( yintercept = target[2], size = sLine ) +
    geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
    geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
    geom_point(size = sPch) + geom_line(size = sLine) +
    theme2 +
    coord_cartesian(ylim = c(ylow, yup ))
)
dev.off()

#---- Collecting the data
# obs = 0.05
title = "Observation Noise: sd = 0.05"
obs = 0.05
covRate <- NULL
for(nx in nx_vec){
  load( file = paste( paste(path_wd, "Workspaces/",sep = ""),
                      date, "_",
                      model, "_",
                      transformation.l,
                      "_bias_", gsub(" ", "_", bias.est.l),
                      "_se_", gsub(" ", "_", se.est.l),
                      "_nx_", nx,
                      "_obs_", obs,
                      "_full",
                      ".rdata",sep = "") )

  covRates_nx <- cbind(covRates, nx)
  colnames(covRates_nx) <- c(colnames(covRates), "nx")

  covRate <- rbind(covRate, covRates_nx)
}

covRates <- covRate
rm(covRate)

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", gsub(" ", "_", bias.est.l),
                  "_se_", gsub(" ", "_", se.est.l),
                  "_obs_", obs,
                  ".png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = c("N", "T"), variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric) %>% mutate_at("T", as.factor)

# Plot the Covering Rates by Method
print(
  ggplot(covs, aes(N, CovRate, shape = Method, col = T, group = interaction(Method, nx))) +
    xlab( "Sample Size [N]" ) +
    ylab( "Covering Rate" ) +
    geom_hline( yintercept = target[2], size = sLine ) +
    geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
    geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
    geom_point(size = sPch) + geom_line(size = sLine) +
    theme1 +
    coord_cartesian(ylim = c(ylow, yup ))
)
dev.off()
