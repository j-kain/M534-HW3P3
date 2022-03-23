source("R/02-gradient-mle-functions.R")
load("data/data.rda")

mu0 <- matrix(c(0,0,0), ncol=1)
sig0 <- diag(3)
theta0 <- param_convert(mu0, sig0, t_comp=TRUE)$t
# theta0 <- as.numeric(theta0)
# log_like_mvn(theta0)
# gradient(theta0)


optim(par=theta0, fn=log_like_mvn, gr=gradient, method="BFGS", 
      control = list(fnscale=-1, trace=1, abstol=1e-5))

