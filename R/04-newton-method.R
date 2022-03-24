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
    cbind(da,db0,db1)
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

stop.criteria <- function(t0, t1, grad, tol.mre=1e-6, tol.grad=1e-5, itr, max.itr){
    mre <- max(abs((t1 - t0)/pmax(1, t1)))
    gradn <- norm(grad,"2")
    stop <- ifelse(itr == max.itr | gradn < tol.grad | mre < tol.mre, TRUE, FALSE)
    list(stop=stop, norm.g=gradn)
}

step.halve <- function(data, mu, sig, t, dir.mu, dir.sig, fun, max.h, itr){
    
    sig1 <- theta[]
    mu1 <- param.convert(theta=t, ms.comp=TRUE)$mu
    
    
    halve = 0
    while(any(eigen(sig1)$values < 0 & halve < max.h)){
        halve = halve + 1
        sig1 <- sig + dir.sig/2^halve
    }
    
    mu1 <- mu + dir.mu/2^halve
    
    if(any(eigen(sig1)$values<0)){stop("Sigma can not be negative definite, try increasing max halves")}
    
    fun1 <- log.like.mvn(data, mu1, sig1)
    
    h = 0
    while(fun1 < fun & h < max.h){
        h = h + 1
        halve = halve + 1
        mu1 <- mu + dir.mu/2^halve
        sig1 <- sig + dir.sig/2^halve
        fun1 <- log.like.mvn(data, mu1, sig1)
    }
    list(mu1 = mu1, sig1 = sig1, halves=hprint, fun1=fun1, llh=llh)
}

#----------------------------------------------------------


#----------------------------------------------------------
newton <- function(theta, tol.grad=1e-5, tol.mre=1e-6, max.itr=50){
    
    it <- 1
    stop <- FALSE
    
    while(!stop){
        obj.fn <- log_like_weib(theta) # set objective function        
        
        grad <- gradient(theta)
        hess <- hessian(theta)
        
        dir <- -solve(hess) * grad
        
        t1 <- theta + dir
        
        r <- step.halve(data, mu, sigma, t1, dm, ds, obj.fn, 20, itr=it)
        
        obj.fn <- r$fun1
        
        t1 <- param.convert(mu=r$mu1, sigma=r$sig1, t.comp=TRUE)$t
        
        stop.cri <- stop.criteria(t0, t1, grad=d, tol.mre=tol.mre, tol.grad=tol.grad, itr=it, max.itr=max.itr)
        stop <- stop.cri$stop
        
        
        it <- it + 1
        
        mu <- r$mu1
        
        sigma <- r$sig1
    }
    list(mu = mu, sigma=sigma)
}




















