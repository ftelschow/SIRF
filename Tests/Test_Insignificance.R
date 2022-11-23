  #-------------------------------------------------------------------------------
  # This file tests the preimage estimation
  #-------------------------------------------------------------------------------
  # Prepare workspace
  rm(list = ls())

  library(SampleFields)
  library(tidyverse)
  library(SIRF)
  library(mutoss)
  library(multtest)

  setwd("/home/fabian/Seafile/Code/Rpackages/SIRF/")

  source("Tests/Auxillary_fcns.R")
  #-------------------------------------------------------------------------------
  # Variables: General
  N     = 1e2
  Ntrue = 80
  mcorrect = "sidak" # "holm" #

  mu_name    = "1" #  "3" # "2" # "4" #
  SCoPEStype = "classical" # "extraction"
  mu1est     =  "thickening" # NULL #
  kNtype = "log" # "SCB" #
  betaN  = 0.95
  k_fac  = 2

  # General Simulation parameters
  alpha = 0.1
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
  }else if(mu_name == "3"){
    muvec = sin((1:100)/2/pi)
  }else if(mu_name == "4"){
    NDelta = c(5, 75, 0, 0, 0)
    muvec = generate_muvec(NDelta)
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

  #-------------------------------------------------------------------------------
  # Main test area
  if(SCoPEStype == "classical"){
    I0 = model$mu(model$x) == C
    I1 = model$mu(model$x) != C
  }

  # Get the SCoPES for the method
  scopes <- SCoPES(alpha = alpha, C = C, x = x, hatmu = hatmu,
                   hatsigma = hatsigma, tN = 1 / sqrt(N),
                   method = method,  R = R, mu = model$mu(model$x))

  plot_SCoPES(scopes = scopes, index_C = 1, mu = mu)

  # IVs = IVobs(scopes, method)
  # H0 = any(scopes$hatLC) | any(scopes$hatUC)
  #
  # IV = round(c(H0, IVs$IV0, IVs$IV1, IVs$IV1sim, IVs$IVk, IVs$IVkabs, IVs$IVobs), 3)
  # names(IV) <- c("H0", "IV0", "IV1", "IV1sim", "IVk", "IVkabs", "IVobs")
  # IVloc = round(IVs$IVloc, 3)

  if(SCoPEStype == "classical"){
    # Bonferroni-Holm-correction
    pvals <- apply(y$Y, 1, function(v) t.test(x = v,
                              alternative = "two.sided",
                              conf.level = 1-alpha)$p.value)
   # pvals_adjust = p.adjust(pvals, method = "sidak")
    detect_holm  = FWE_control(alpha, pvals, mcorrect)
    # Benyamini-Hochberg-correction
    detect_BH <- FDR_control(alpha, pvals, "BH")

    t_detect_h <- sum(detect_holm[I1])
    f_detect_h <- sum(detect_holm[I0])

    t_detect_bh <- sum(detect_BH[I1])
    f_detect_bh <- sum(detect_BH[I0])

    test = round(cbind(c(t_detect_bh, f_detect_bh),c(t_detect_h, f_detect_h), c(scopes$NtrueDetect, scopes$NfalseDetect), c(kN, kN)), 3)
    colnames(test) <- c("BH",mcorrect, "SCoPES", "kN")
    rownames(test) <- c("true", "false")

  }else{
    test = round(cbind(c(scopes$NtrueDetect, scopes$NfalseDetect), c(kN, kN)), 3)
    names(test) <- c("true", "false", "kN")
  }

  IV <- IVs_IID(scopes, cdf = pt, cdf_abs = function(q) extraDistr::pht(q = q, nu = N-1),
                rcdf = rt, Msim = 2e4, df = N - 1)

#  IVs = round(c(IV$IV0, IV$IV0kN, IV$IV1,  IV$IV2,  IV$IV3,  IV$IV4, IV$IVkabs, IV$IVk, IV$IVobs), 4)
#  names(IVs) = c("IV0", "IVkn", "IV1", "IV2", "IV3", "IV4", "IVkabs", "IVk", "IVobs")
  IVs = round(rbind(IV$IV, IV$IVo), 4)
  rownames(IVs) = c("IVq", rep("IVobs",dim(IV$IVo)[1]))
  colnames(IVs) = 1:length(IV$IV)

  test
  #IVloc
c(all(!detect_holm == scopes$hatmu1C$minus), sum(!detect_holm), sum(scopes$hatmu1C$minus))

# Get the SCoPES for the method
method2 = method
method2$mu1Cest = list(minus = !detect_BH, plus = !detect_BH)

#scopes2 <- SCoPES(alpha = alpha, C = C, x = x, hatmu = hatmu,
#                 hatsigma = hatsigma, tN = 1 / sqrt(N),
#                 method = method2,  R = R, mu = model$mu(model$x))
#test[1,4] = scopes2$NtrueDetect
test
IVs
