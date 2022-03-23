source("R/02-gradient-mle-functions.R")
load("data/data.rda")

mu0 <- matrix(c(0,0,0), ncol=1)
sig0 <- diag(3)
theta0 <- param_convert(mu0, sig0, t_comp=TRUE)$t

optim(par=theta0, fn=log_like_mvn, gr=gradient, method="BFGS", 
      control = list(fnscale=-1, trace=100, abstol=1e-5))

