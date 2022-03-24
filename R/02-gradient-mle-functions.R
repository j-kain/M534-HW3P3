load("data/datan.rda")

param_convert <- function(mu=NULL, sigma=NULL, theta=NULL, t_comp=FALSE, ms_comp=FALSE){
    #writeLines("\n entered in param_conv function")
    
    # decompose mu and sig into theta vector
    if(t_comp){t <- matrix(c(mu, sigma[upper.tri(sigma, TRUE)]),ncol=1)}
    # decompose theta vector into mu and sig 
    if(ms_comp){
        p <- sqrt(2*length(theta) + 9/4) - 3/2 # formula derived from num of thetas = p + p(p+1)/2
        mu <- matrix(theta[1:p],p,1) # create new mu vector
        
        sig <- matrix(0,p,p) # create new sig matrix
        sig[upper.tri(sig,TRUE)] <- theta[-c(1:p)]
        sig[lower.tri(sig)] <- t(sig)[lower.tri(sig)]
    }
    #writeLines("leave param_conv function")
    list(t=if(t_comp){t}, mu=if(ms_comp){mu}, sig=if(ms_comp){sig})
}

diff_xi_mu <- function(mu, comp_C=FALSE, comp_sxm=FALSE){
    # Compute C defined as sum of (xi-mu)(xi-mu)^T 
    xm <- apply(datan, 1, function(x) x-mu) # subtract mu from rows of x
    C <- xm %*% t(xm)        
    
    # Compute sum of (xi-mu) for use in gradient
    if(comp_sxm){
        sxm <- matrix(rowSums(xm), nrow=ncol(datan),1)
    }
    
    list(C = if(comp_C) C, sxm = if(comp_sxm) sxm)
}

log_like_mvn <- function(theta){

    #writeLines("\n entered in log like function")
    n <- nrow(datan) # num of rows
    p <- ncol(datan) # num of columns
    
    t <- param_convert(theta=theta, ms_comp=TRUE)
    sigma <- t$sig
    mu <- t$mu
                       
    if(any(eigen(sigma)$values <= 0)){
        sigma[,] <- -Inf
    }
    
    C <- diff_xi_mu(mu, comp_C=TRUE)$C # C := sum of (xi-mu)(xi-mu)^T 
    
    log_det_sig <- log(det(sigma)) # log determinant of sigma
    sig_inv <- solve(sigma) # sigma inverse
    
    # Compute log likelihood function for multivariate normal
    # sum(sig.inv * C) = trace(sig.inv %*% C)
    #writeLines("leave log like function")
    log_like <- (-1/2)*(n*p*log(2*pi)+n*log_det_sig + sum(sig_inv * C ))
    log_like
    
}

gradient <- function(theta, t_comp=TRUE, dmu_comp=FALSE, dsig_comp=FALSE){    
    t <- param_convert(theta=theta, ms_comp=TRUE)
    mu <- t$mu
    sigma <- t$sig
    
    if(any(eigen(sigma)$values <= 0)){
        sigma[,] <- -Inf
    }
    
    #writeLines("\n entered in gradient function")
    sig_inv <- solve(sigma)
    diff <- diff_xi_mu(mu, comp_C=TRUE, comp_sxm=TRUE)
    sxm <- diff$sxm
    C <- diff$C
    
    # partial deriv_ formulas for mu indices
    dmu <- sig_inv %*% sxm
    
    # partial deriv_ formulas for sigma indices
    off_d_mat <- -(n * sig_inv - sig_inv %*% C %*% sig_inv) # take off diag elem_
    diag_mat <- (-1/2)*(n*sig_inv - sig_inv %*% C %*% sig_inv) # take diag elem_
    
    # dsig
    dsig <- matrix(0, nrow=p, ncol=p) # initialize pxp matrix for gradient of sigmas
    dsig[upper.tri(dsig)] <- off_d_mat[upper.tri(off_d_mat)] # fill matrix with lower tri of off_d_mat
    dsig[lower.tri(dsig)] <- t(dsig)[lower.tri(dsig)] # fill matrix with upper tri of off_d_mat
    diag(dsig) <- diag(diag_mat) # fill matrix with diag elements diagmat
    
    # theta form
    t <- param_convert(mu=dmu, sigma=dsig, t_comp=TRUE)$t
    t
    #writeLines("leave gradient function")

    # output options for gradient
    #list(t=if(t_comp) t, dmu=if(dmu_comp) dmu, dsig=if(dsig_comp) dsig)
}








