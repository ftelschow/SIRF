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

generate_muvec <- function(mu_name){
      # Model parameters
      if(mu_name == "1"){
        NDelta = c(0, 80, 0, 0, 0)
        muvec = generate_muvec(NDelta)
      }else if(mu_name == "2"){
        NDelta = c(30, 20, 30, 0, 0)
        muvec = generate_muvec(NDelta)
      }else if(mu_name == "3"){
        NDelta = c(5, 75, 0, 0, 0)
        muvec = generate_muvec(NDelta)
      }else if(mu_name == "4"){
        muvec = sin( (1:100) / 2 / pi )
      }

    return(muvec)
}

get_SCBquant <- function(betaN, N, muvec){
  mm = list(minus = rep(T, length(muvec)),#-nx_close2),
            plus = rep(T, length(muvec)))#-nx_close2))
  qest <- function(q) maxT_p(q, mm, df = N-1) - (1 - betaN)

  return(uniroot(qest, interval = c(-100, 100))$root)
}
