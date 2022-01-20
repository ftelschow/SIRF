#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Load packages
suppressPackageStartupMessages(library(SIRF))
suppressPackageStartupMessages(library(SampleFields))
suppressPackageStartupMessages(library(RFT))
suppressPackageStartupMessages(library(tidyverse))

# Get the working directory
path_wd = args[7]

print(path_wd)

# Set the working path
setwd(path_wd)

# Get the simulation file from source
source(paste(path_wd, "Article_simulation_obs.R", sep = ""))

# needed to make string without empty spaces
if(args[8] == "skewnessN"){
  transform = "skewness (normality)"
}else if(args[8] == "kurtosisN"){
  transform = "kurtosis (normality)"
}else{
  transform = args[8]
}

# Simulate the covering rate
Article_simulation( Model   = args[1], # "ModelA", "ModelB", "ModelC"
                    transformation = args[8],
                     Msim    = as.numeric(args[2]), # number of simulations
                     Nvec    = c( 50, 100, 200, 400, 800 ), # sample sizes considered
                     x       = seq(0, 1, length.out = as.numeric(args[4])), # locations the process is evaluated at
                     level   = .95,  # level of simultaneous control
                     obs     = as.numeric(args[5]),
                     sim_tag = args[3],
                     path_wd = path_wd,
                     date    = args[6] )
