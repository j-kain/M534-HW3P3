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
    
    hess <- matrix(0, nrow=3, ncol=3)
    
    daa <- sum(-w/a^2 - t^a*exp(b0 + b1*d)*log(t)^2)
    dab0 <- sum(-t^a*exp(b0 + b1*d)*log(t))
    dab1 <- sum(-d*t^a*exp(b0 + b1*d)*log(t))
    fa <- cbind(daa,dab0,dab1)
    
    db0a <- sum(-t^a*exp(b0 + b1*d)*log(t) )   
    db0b0 <- sum(-t^a*exp(b0 + b1*d))
    db0b1 <- sum(-d*t^a*exp(b0 + b1*d))
    fb0 <- cbind(dab0, db0b0, db0b1)
    
    db1a <- sum(-d*t^a*exp(b0 + b1*d)*log(t))
    db1b0 <- sum(-d*t^a*exp(b0 + b1*d))    
    db1b1 <- sum(-d^2*t^a*exp(b0 + b1*d))
    fb1 <- cbind(dab1, db0b1, db1b1)
    
    hess[1,] <- fa
    hess[2,] <- fb0
    hess[3,] <- fb1
    
    hess
}

hessian(theta)























