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
today = "2022_22_9_"
#-------------------------------------------------------------------------------
mu_name = "1"
SCoPEStype = "classical" # "extraction"

# General Simulation parameters
Msim  = 1e4
alpha = 0.1
C     = c(0,3)

NVec    = c(20, 50, 1e2, 2e2, 5e2, 10e2)
betaVec = c(666, 10, 5, 3, 2, 1-0.01, alpha )

# variables: q estimation
name       = "t" # "gauss" # "mboot" #

# Output containers
lN = length(NVec)
bN = length(betaVec)
coverage       <- matrix(NaN, bN, lN)
trueDetect_BH  <- trueDetect_holm  <- trueDetect_sidak <- trueDetect_scopes <- array(NaN, dim = c(bN, lN, Msim))
falseDetect_BH <- falseDetect_holm <- falseDetect_sidak <- falseDetect_scopes <- array(NaN, dim = c(bN, lN, Msim))

# Model specification
muvec = generate_muvec(mu_name)
model <- list(
  x  = 1:length(muvec),
  mu = Vectorize(function(x){ t = muvec; return(t[x]) }),
  noise = SampleFields::IIDProcess,
  sigma = Vectorize(function(x){ return(1) }),
  truesigma = FALSE
)

# Variables for the tubes
if(SCoPEStype == "classical"){
  C = rep(0, length(muvec))
}else{
  C = cbind(rep(0, length(muvec)), rep(3, length(muvec)))
}

for(n in 1:lN){for(b in 1:bN){
  if(betaVec[b] < 1){
    kN     = get_SCBquant(betaVec[b], NVec[n], muvec)
    mu1est = "thickening"
  }else if(b != bN){
    kN     = log(NVec[n]) / betaVec[b]
    mu1est = "thickening"
  }else{
    kN = 1
    mu1est = NULL
  }

  # Generate method list
  method = method_gen(name, SCoPEStype, mu1est, NVec[n], kN )

  # Simulate the coverage etc
  test = sim_SCoPES(Msim = Msim, N = NVec[n], alpha = alpha, C = C,
                    q_method = method, model = model)

  coverage[b,n]               = test$coverage
  trueDetect_scopes[b, n, ]   = test$NDetect[1, ]
  falseDetect_scopes[b, n, ]  = test$NDetect[2, ]

  if(SCoPEStype == "classical"){
    trueDetect_sidak[b, n, ]  = test$NDetect_sidak[1, ]
    falseDetect_sidak[b, n, ] = test$NDetect_sidak[2, ]
    trueDetect_holm[b, n, ]   = test$NDetect_holm[1, ]
    falseDetect_holm[b, n, ]  = test$NDetect_holm[2, ]
    trueDetect_BH[b, n, ]     = test$NDetect_BH[1, ]
    falseDetect_BH[b, n, ]    = test$NDetect_BH[2, ]
  }
}
}

rm(test, n, b, lN, bN, get_SCBquant, generate_muvec, generateData, kN)
save.image(paste(today, SCoPEStype, "_Model", mu_name, ".Rdata", sep = ""))

