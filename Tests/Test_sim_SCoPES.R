#-------------------------------------------------------------------------------
# This file tests the preimage estimation
#-------------------------------------------------------------------------------
# Prepare workspace
rm(list = ls())

library(SampleFields)
library(tidyverse)
library(SIRF)

source("Tests/Auxillary_fcns.R")
#-------------------------------------------------------------------------------
# Variables General
Msim  = 1e3
# Variables: General
N     = 1e3
Ntrue = 80
alpha = 0.1

mu_name    =   "1" #"2" # "3" #
SCoPEStype = "classical" # "extraction"
mu1est     =  "thickening" # NULL #
kNtype =  "log" # "SCB" #
betaN  = 0.99
k_fac  = 3

B     = c(0, 3)

NVec    = c(20, 50, 1e2, 2e2, 5e2, 10e2)
betaVec = c(0.99, 0.2, 0.5, 0.15, 0.1, 0.05)

# Model parameters
if(mu_name == "1"){
  NDelta = c(30, 20, 30, 0, 0)
  muvec = generate_muvec(NDelta)
}else if(mu_name == "2"){
  NDelta = c(0, Ntrue, 0, 0, 0)
  muvec = generate_muvec(NDelta)
}else{
  muvec = sin((1:100)/2/pi)
}


# variables: q estimation
name       = "t" # "gauss" # "mboot" #
truesigma  = FALSE

kNold = log(N) / k_fac

if(kNtype == "SCB"){
  kN = get_SCBquant(betaN, N, muvec)
}else{
  kN = kNold
}

#-------------------------------------------------------------------------------
# Generate data
y = generateData(N, muvec, B, truesigma, SCoPEStype)
x = y$x; hatmu = y$hatmu; hatsigma = y$hatsigma; R = y$R; model = y$model;
C = y$C; tN = y$tN
#-------------------------------------------------------------------------------
# Generate method list
method = method_gen(name, SCoPEStype, mu1est, N, kN, R )

tmp = PreimageC(cbind(C,C), hatmu, hatsigma, tN, kN = kN, method = "classical")
all(tmp$minus)

#-------------------------------------------------------------------------------
# Main test area
test = sim_SCoPES(Msim = Msim, N = N, alpha = alpha, C = C,
                  q_method = method, model = model)

results <- c(test$coverage)
names(results) <- c("Coverage")
results

if(SCoPEStype == "extraction"){
  index_minus = model$mu(model$x) == C[, 1]
  index_plus = model$mu(model$x) == C[, 2]

  par(mfrow = c(1,2))
  hist(apply(test$mu1C$minus[index_minus,], 2, sum),
       xlab = "Number Detected C-",
       main = paste("True = ", sum(index_minus), sep = "") )

  hist(apply(test$mu1C$plus[index_plus,], 2, sum),
       xlab = "Number Detected C+",
       main = paste( "True = ", sum(index_plus), sep = "") )

}else if(SCoPEStype == "classical"){
  index_minus = model$mu(model$x) == C
  par(mfrow = c(1,1))
  hist(apply(test$mu1C$minus[index_minus,], 2, sum),
       xlab = "Number Detected C",
       main = paste("True = ", sum(index_minus), sep = "") )
}

# Check the testing
mtrue_detect = c(mean(test$NDetect[1,]),
                mean(test$NDetect_holm[1,]),
                mean(test$NDetect_sidak[1,]))
mfalse_detect = c(mean(test$NDetect[2,]),
                 mean(test$NDetect_holm[2,]),
                 mean(test$NDetect_sidak[2,]))
false_inclusion = c(mean(test$NDetect[2,]>0),
                  mean(test$NDetect_holm[2,]>0),
                  mean(test$NDetect_sidak[2,]>0))

summ <- rbind(mtrue_detect, mfalse_detect, false_inclusion)
colnames(summ) <- c("SCoPES", "holm", "sidak")
rownames(summ) <- c("av. true detect", "av. false detect", "SCoPES incl false")

results <- c(test$coverage)
names(results) <- c("Coverage")
results
c(kNold, kN)
summ

#apply(test$mu1C$minus | test$mu1C$plus, 2, sum)
