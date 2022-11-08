generateData <- function(N, muvec, B, truesigma, SCoPEStype){
  # Get derived quantities
  x = 1:length(muvec)
  if(SCoPEStype == "extraction"){
    C = cbind(rep(B[1], length(x)), rep(B[2], length(x)))
  }else if(SCoPEStype == "selection"){
    B = seq(B[1], B[2], length.out =  4)
    C = cbind(rep(B[1], length(x)), rep(B[2], length(x)), rep(B[3], length(x)), rep(B[4], length(x)))
  }else{
    C = rep(0, length(x))
  }

  # Model specification
  model <- list(
    x  = x,
    mu = Vectorize(function(x){ t = muvec; return(t[x]) }),
    noise = SampleFields::IIDProcess,
    sigma = Vectorize(function(x){ return(1) }),
    truesigma = truesigma
  )

  # Get the sample for the simulation run
  Y = SampleFields::SignalPlusNoise(N,
                                    x = model$x,
                                    mu = model$mu,
                                    noise = model$noise,
                                    sigma = model$sigma)

  # Estimate the mean and sd
  Y        = Y$values
  hatmu    = rowMeans(Y)
  hatsigma = apply(Y, 1, sd)
  R        = Y - rowMeans(Y)


  return(list(Y = Y, x = x, hatmu = hatmu, hatsigma = hatsigma,
              R = R, model = model, C = C, tN = 1 / sqrt(N)))
}

generate_muvec <- function(NDelta, mu_name = "1"){
    # Variables for the tubes
    if(mu_name == "1"){
      B = c(0, 3)
      muu   = c(-0.3, 0, 0.2, 3, 4)
      muvec = c(rep(muu[1], NDelta[1]),
                rep(muu[2], NDelta[2]),
                rep(muu[3], NDelta[3]),
                rep(muu[4], NDelta[4]),
                rep(muu[5], NDelta[5]))
      muvec = muvec[sample(1:length(muvec), replace = FALSE)]
    }else if(mu_name == "2"){
      B = c(0, 3)
      muvec = c(rep(NDelta[1], NDelta[2]), rep(0,  NDelta[3]))
    }
    return(muvec)
}

get_SCBquant <- function(betaN, N, muvec){
  mm = list(minus = rep(T, length(muvec)),#-nx_close2),
            plus = rep(T, length(muvec)))#-nx_close2))
  qest <- function(q) maxT_p(q, mm, df = N-1) - (1 - betaN)

  return(uniroot(qest, interval = c(-100, 100))$root)
}
