source("R/02-gradient-mle-functions.R")

data2 <- read.csv(here("data","censor_data.csv"))


theta <- c(1,1,1)


log_like_weib <- function(theta){
    a = theta[1]
    b0 = theta[2]
    b1 = theta[3]
    t = data2$t
    w = data2$w
    d = data2$d
    like <- sum((w*log(t^a * exp(b0+d*b1))) - t^a * exp(b0+d*b1) + w*log(a/t))
    like
}

gradient <- function(theta){
    a = theta[1]
    b0 = theta[2]
    b1 = theta[3]
    t = data2$t
    w = data2$w
    d = data2$d
    
    da <- sum(w/a + w*log(t) - t^a*exp(b0 + b1*d)*log(t))
    db0 <- sum(w - t^a*exp(b0 + b1*d))
    db1 <- sum(d*w - d*t^a*exp(b0 + b1*d))
    rbind(da,db0,db1)
}

hessian <- function(theta){
    a = theta[1]
    b0 = theta[2]
    b1 = theta[3]
    t = data2$t
    w = data2$w
    d = data2$d
    
    daa <- - w/a^2 - t^a*exp(b0 + b1*d)*log(t)^2
    dab0 <- -t^a*exp(b0 + b1*d)*log(t)
    dab1 <- -d*t^a*exp(b0 + b1*d)*log(t)
    
    db0b0 <- -t^a*exp(b0 + b1*d)
    db0a <- -t^a*exp(b0 + b1*d)*log(t)
    db0b1 <- -d*t^a*exp(b0 + b1*d)

    db1b1 <- -d^2*t^a*exp(b0 + b1*d)
    db1a <- -d*t^a*exp(b0 + b1*d)*log(t)
    db1b0 <- -d*t^a*exp(b0 + b1*d)
}
























