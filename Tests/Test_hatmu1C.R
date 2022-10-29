N  = 100
lx = 6

mu = c(0,0,0,4,4,4)

X = mu + matrix(rnorm(N * lx), lx, N)

tN = 1 / sqrt(N)
kN = log(tN^-2) / 1.5
pchisq(kN^2, 1)^length(mu)

# cut-off function sets
C = cbind(rep(0, length(hatmu)), rep(1, length(hatmu)), rep(4, length(hatmu)))

hatmu = rowMeans(X)
hatsigma = apply(X, 1, sd)

test1 = SCoPES(0.05, C, hatmu, hatsigma, tN, method = "extraction", R = X - hatmu)

method <- list(name      = "gauss",
               SCoPEtype = "extraction",
               mu1Cest   = "thickening",
               kN        = log(tN^-2) / 5)
test2 = SCoPES(0.05, C, hatmu, hatsigma, tN, method = method)

method <- list(name      = "t",
               df        = N - 1,
               SCoPEtype = "extraction",
               mu1Cest   = "thickening",
               kN        = log(tN^-2) / 5)
test3 = SCoPES(0.05, C, hatmu, hatsigma, tN, method = method)


# Estimate the preimage of C
test = PreimageC(C, hatmu, hatsigma, tN, kN, method = "extraction")

# Estimate the quantile


hatLc(0, mu, hatsigma, 0, q=3, mu = mu, inclusion = "inequal" )
hatLc(0, hatmu, hatsigma, tN, q=3, mu = mu, inclusion = "equal" )
hatUc(1, mu, hatsigma, 0, q=3, mu = mu, inclusion = "inequal" )
hatUc(0, hatmu, hatsigma, tN, q=3, mu = mu, inclusion = "equal" )

test = SCoPES(C, mu, hatsigma, 0, q=3, mu = mu )

test$hatLC
test$hatUC

PreimageC(C, hatmu, hatsigma, tN, kN, method = "relevance")

hatE_pm(hatmu, C, hatsigma, tN, kN)
