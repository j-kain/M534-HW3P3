library(devtools)
library(here)
library(formattable)

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

stop.criteria <- function(t0, t1, grad, tol.mre=1e-6, tol.grad=1e-5, itr, max.itr=100){
    mre <- max(abs((t1 - t0)/pmax(1, t1)))
    gradn <- norm(grad,"2")
    stop <- ifelse(itr == max.itr | gradn < tol.grad | mre < tol.mre, TRUE, FALSE)
    list(stop=stop, norm.g=gradn)
}

#----------------------------------------------------------
newton <- function(theta, tol.grad=1e-5, tol.mre=1e-6, max.itr=50){
    
    it <- 1
    stop <- FALSE
    
    while(!stop){
        obj.fn <- log_like_weib(theta) # set objective function        
        
        grad <- gradient(theta)
        hess <- hessian(theta)
        
        dir <- -solve(hess) %*% grad
        
        t1 <- theta + dir
        fun2 <- log_like_weib(t1)

        halve = 0
        while(fun2 < obj.fn & halve <= 20){
            halve = halve + 1
            t1 <- t0 + dir/2^halve
            fun2 <- log.like.weib(t1)
            print(c(it, halve))
        }
        if (halve >= 20) print('Step-halving failed after 20 halvings')
        
        obj.fn <- fun2
        
        stop.cri <- stop.criteria(theta, t1, grad=grad, itr=it)
        grad_norm <- stop.cri$norm.g
        
        
        stop <- stop.cri$stop
        
        it <- it + 1
        
        theta <- t1
        
        writeLines(paste("\nitr: ", it-1))
        writeLines(paste("halving   ", "log-like     ", "norm   "))
        writeLines(paste(halve, "        ", formattable(obj.fn, digits=4, format="f"), "   ", formattable(grad_norm, digits=1, format="e")))
        
        print ('-----------------------------------------')

    }
    theta
}

theta <- c(1,0,0)
theta_star <- newton(theta)
observed_info <- -solve(hessian(theta_star))
sd_err <- sqrt(diag(observed_info))
cov2cor(V=observed_info)











