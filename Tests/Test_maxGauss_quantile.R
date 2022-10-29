beta = 0.9
Msim = 4e5
df   = 5

muC = list()
muC$minus = FALSE * 1:5
muC$plus  = FALSE * 1:5

mv = c(1,2,3)
pv = NULL  # either NULL or extension of mv

lmv = length(mv)
liv = length(pv) - lmv

muC$minus[mv] = TRUE
muC$plus[pv] = TRUE

qg = maxGauss_quantile(beta, muC)

data1 = matrix(rnorm(lmv * Msim), Msim, lmv)
if(is.null(pv)){
  maxXg = apply(data1, 1, max)
}else{
  data2 = matrix(rnorm(liv * Msim), Msim, liv)
  maxXg = apply(cbind(abs(data1), data2), 1, max)
}

qt = maxT_quantile(beta, muC, df = df)
data1 = matrix(rt(lmv * Msim, df = df), Msim, lmv)
if(is.null(pv)){
  maxXt = apply(data1, 1, max)
}else{
  data2 = matrix(rt(liv * Msim, df = df), Msim, liv)
  maxXt = apply(cbind(abs(data1), data2), 1, max)
}

mean(maxXg < qg)
mean(maxXt < qt)
