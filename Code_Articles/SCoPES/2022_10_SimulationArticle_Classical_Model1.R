#-------------------------------------------------------------------------------
# Simulations for the article
#-------------------------------------------------------------------------------
# Prepare workspace
rm(list = ls())

setwd("~/Rpackages/SIRF/Code_Articles/SCoPES")

library(SampleFields)
library(tidyverse)
library(SIRF)

source("Auxillary_fcns.R")
today = "2022_11_9_"
#-------------------------------------------------------------------------------
mu_name = "1" # "2" #
SCoPEStype = "classical" # "extraction"

# General Simulation parameters
Msim  = 1e4
alpha = 0.1
C     = c(0,3)

NVec    = c(20, 50, 1e2, 2e2, 5e2, 10e2)
betaVec = c(0.99, 0.2, 0.5, 0.15, 0.1, 0.05)

# Model parameters
if(mu_name == "1"){
  NDelta = c(30, 20, 30, 0, 0)
}else{
  NDelta = c(0, 80, 0, 0, 0)
}

# variables: q estimation
name       = "t" # "gauss" # "mboot" #
mu1est     = "thickening" # NULL #

# Output containers
lN = length(NVec)
bN = length(betaVec)
coverage         <- matrix(NaN, bN, lN)
trueDetect_holm  <- trueDetect_sidak <- trueDetect_scopes <- array(NaN, dim = c(bN, lN, Msim))
falseDetect_holm <- falseDetect_sidak <- falseDetect_scopes <- array(NaN, dim = c(bN, lN, Msim))

# Model specification
muvec = generate_muvec(NDelta)
model <- list(
  x  = 1:length(muvec),
  mu = Vectorize(function(x){ t = muvec; return(t[x]) }),
  noise = SampleFields::IIDProcess,
  sigma = Vectorize(function(x){ return(1) }),
  truesigma = FALSE
)

# Variables for the tubes
if(SCoPEStype == "classical"){
  C = rep(0, sum(NDelta))
}else{
  C = cbind(rep(0, sum(NDelta)), rep(3, sum(NDelta)))
}

for(n in 1:lN){for(b in 1:bN){
  kN = get_SCBquant(betaVec[b], NVec[n], muvec)

  # Generate method list
  method = method_gen(name, SCoPEStype, mu1est, NVec[n], kN )

  # Simulate the coverage etc
  test = sim_SCoPES(Msim = Msim, N = NVec[n], alpha = alpha, C = C,
                    q_method = method, model = model)

  coverage[b,n]               = test$coverage
  trueDetect_scopes[b, n, ]   = test$NDetect[1, ]
  falseDetect_scopes[b, n, ]  = test$NDetect[2, ]
  trueDetect_sidak[b, n, ]    = test$NDetect_sidak[1, ]
  falseDetect_sidak[b, n, ]   = test$NDetect_sidak[2, ]
  trueDetect_holm[b, n, ]     = test$NDetect_holm[1, ]
  falseDetect_holm[b, n, ]    = test$NDetect_holm[2, ]
  }
}

rm(test, n, b, lN, bN, get_SCBquant, generate_muvec, generateData, kN)
save.image(paste(today, SCoPEStype, "_Model", mu_name, ".Rdata", sep = ""))

