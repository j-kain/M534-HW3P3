
sqrtm <- function (A) {
    # Obtain matrix square root of a matrix A
    a = eigen(A)
    sqm = a$vectors %*% diag(sqrt(a$values)) %*% t(a$vectors)
    sqm = (sqm+t(sqm))/2
}

# Generate data

gen <- function(n,p,mu,sig,seed = 534){
    #---- Generate data from a p-variate normal with mean mu and covariance sigma
    # mu should be a p by 1 vector
    # sigma should be a positive definite p by matrix
    # Seed can be optionally set for the random number generator
    set.seed(seed)
    # generate data from normal mu sigma
    x = matrix(rnorm(n*p),n,p)
    datan = x %*% sqrtm(sig) + matrix(mu,n,p, byrow = TRUE)
    datan
}

# sample data
n <- 200
p <- 3
sigma <- matrix(c(1,.7,.7, .7, 1, .7, .7, .7, 1),p,p)
mu <- matrix(c(-1,1,2),p,1)
datan = gen(n,p,mu,sigma,seed = 2022)

save(datan, file="data/datan.rda")




