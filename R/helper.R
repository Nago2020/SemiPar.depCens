
#' Simulate a toy example for testing
#'
#' @param n sample size
#' @importFrom copula pCopula frankCopula gumbelCopula tau
#' @importFrom stats rbinom runif

genData <- function(n){

beta = c(0.5)
lambd = 0.35
eta = c(0.9,0.4)
lambd = 0.35
X = cbind(rbinom(n,1,0.5))
W = cbind(rep(1,n),rbinom(n,1,0.5))

frank.cop <- copula::frankCopula(param = 5,dim = 2)       # generate dependency structure from Frank
U = copula::rCopula(n,frank.cop)
T1 = (-log(1-U[,1]))/(lambd*exp(X*beta))                  # Survival time
T2 = (-log(1-U[,2]))^(1.1)*exp(W%*%eta)                   # Censoring time
A = runif(n,0,15)                                         # administrative censoring time

Z = pmin(T1,T2,A)
d1 = as.numeric(Z==T1)
d2 = as.numeric(Z==T2)

resData = data.frame("Z" = Z,"d1" = d1, "d2" = d2)
return(list("resData" = resData, "X" = X, "W" = W))
}
