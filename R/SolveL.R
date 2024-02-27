
#' @title Cumulative hazard function of survival time under dependent censoring

#' @description  This function estimates the cumulative hazard function of survival time (T) under dependent censoring (C). The estimation
#' makes use of the estimating equations derived based on martingale ideas.
#'
#' @param theta Estimated parameter values/initial values for finite dimensional parameters
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T
#' @param W Data matrix with covariates related to C. First column of W should be ones
#' @param cop Which copula should be computed to account for dependency between T and C. This argument can take
#'   one of the values from \code{c("Gumbel", "Frank", "Normal")}. The default copula model is "Frank".
#' @param dist The distribution to  be used for the dependent censoring C. Only two distributions are allowed, i.e, Weibull
#' and lognormal distributions. With the value \code{"Weibull"} as the
#'   default.
#' @importFrom copula pCopula frankCopula gumbelCopula
#'
#' @return This function returns an estimated hazard function,  cumulative hazard function and distinct observed survival times;
#'
#' @examples
#' \dontrun{
#' n = 200
#' beta = c(0.5)
#' lambd = 0.35
#' eta = c(0.9,0.4)
#' X = cbind(rbinom(n,1,0.5))
#' W = cbind(rep(1,n),rbinom(n,1,0.5))
#' frank.cop <- copula::frankCopula(param = 5,dim = 2)
#' U = copula::rCopula(n,frank.cop)
#' T1 = (-log(1-U[,1]))/(lambd*exp(X*beta))         # Survival time'
#' T2 = (-log(1-U[,2]))^(1.1)*exp(W%*%eta)          # Censoring time
#' A = runif(n,0,15)                               # administrative censoring time
#' Z = pmin(T1,T2,A)
#' d1 = as.numeric(Z==T1)
#' d2 = as.numeric(Z==T2)
#' resData = data.frame("Z" = Z,"d1" = d1, "d2" = d2)
#' theta = c(0.3,1,0.3,1,2)
#'
#' Estimate cumulative hazard function
#'cumFit <- SolveL(theta, resData,X,W)
#'
#' cumulative hazard function
#'cumhaz = cumFit$cumhaz
#'time = cumFit$times
#'
#'# plot hazard vs time
#'
#'plot(time, cumhaz, type = "l",xlab = "Time",
#'ylab = "Estimated cumulative hazard function")
#'
#'}
#'
#' @export


SolveL = function(theta,resData,X,W,
                  cop = c("Frank", "Gumbel",  "Normal"),
                  dist = c("Weibull", "lognormal")){

  cop <- match.arg(cop)
  dist <- match.arg(dist)


  #  Verify column names of input dataframe resData
  if (!all(c("Z", "d1", "d2") %in% colnames(resData)))
    stop("Z, d1 and d2 arguments must be column names of resData")

  id <- match(c("Z", "d1", "d2"), colnames(resData))
  colnames(resData)[id] <- c("Z", "d1", "d2")

  if (is.vector(X)){
    k = 1
    mu1 = theta[1]*X}
  else if(!is.vector(X)){
    k = dim(X)[2]
    mu1 = X%*%theta[1:k]
  }

  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2

  T1 <- c(0,sort(unique(Z[d1 == 1])))
  m = length(T1)
  k = dim(X)[2]
  beta = theta[1:k]

  L = rep(0,m)
  L[1] = 0
  L[2] = sum((Z==T1[2]))/sum((Z>=T1[2])*exp(mu1))
  for (i in 3:m){
    csum = sum(L[1:(i-1)])
    psi = CompC(theta,T1[i],X,W,csum,cop,dist)
    L[i] <- sum(Z == T1[i])/sum((Z>=T1[i])*exp(psi))
  }
  res <- list(lambda = L,cumhaz = cumsum(L), times = T1)
}





#' @title Cumulative hazard function of survival time under independent censoring
#'
#' @description
#' This function estimates the cumulative hazard function of survival time (T) under the assumption of independent censoring.
#' The estimating equation is derived based on martingale ideas.
#'
#' @param theta Estimated parameter values/initial values for finite dimensional parameters
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T
#'
#'
#' @return This function returns an estimated hazard function,  cumulative hazard function and distinct observed survival times;
#'
#' @examples
#' \dontrun{
#' n = 200
#' beta = c(0.5)
#' lambd = 0.35
#' eta = c(0.9,0.4)
#' X = cbind(rbinom(n,1,0.5))
#' frank.cop <- copula::frankCopula(param = 5,dim = 2)
#' U = copula::rCopula(n,frank.cop)
#' T1 = (-log(1-U[,1]))/(lambd*exp(X*beta))           # Survival time'
#' T2 = (-log(1-U[,2]))^(1.1)*exp(W%*%eta)            # Censoring time
#' A = runif(n,0,15)                                  # administrative censoring time
#' Z = pmin(T1,T2,A)
#' d1 = as.numeric(Z==T1)
#' d2 = as.numeric(Z==T2)
#' resData = data.frame("Z" = Z,"d1" = d1, "d2" = d2)
#' theta = c(0.3,1,0.3,1)
#'
#' Estimate cumulative hazard function
#'
#'cumFit_ind <- SolveLI(theta, resData,X)
#'
#' cumulative hazard function
#' cumhaz = cumFit_ind$cumhaz
#' time = cumFit_ind$times
#'
#'# plot hazard vs time
#'
#'plot(time, cumhaz, type = "l",xlab = "Time",
#' ylab = "Estimated cumulative hazard function")
#'
#'}
#'
#'
#' @export


SolveLI = function(theta,resData,X){

  #  Verify column names of input dataframe resData
  if (!all(c("Z", "d1", "d2") %in% colnames(resData)))
    stop("Z, d1 and d2 arguments must be column names of resData")

  id <- match(c("Z", "d1", "d2"), colnames(resData))
  colnames(resData)[id] <- c("Z", "d1", "d2")

  if (is.vector(X)){
    k = 1
    mu1 = theta[1]*X}
  else if(!is.vector(X)){
    k = dim(X)[2]
    mu1 = X%*%theta[1:k]
  }
  Z = resData$Z
  d1 = resData$d1

  T1 <- c(sort(unique(Z[d1 == 1])))
  m = length(T1)
  k = dim(X)[2]
  beta = theta[1:k]

  L = rep(0,m)
  for (i in 1:m){
    L[i] <- sum(Z == T1[i])/sum((Z>=T1[i])*exp(mu1))
  }
  res <- list(lambda = L,cumhaz = cumsum(L), times = T1)
}



#' @title Compute phi function
#'
#' @description  This function estimates phi function at fixed time point t
#'
#' @inheritParams SolveL
#' @param t A fixed time point
#' @param ld Output of \code{\link{SolveL}} function at a fixed time t
#' @importFrom copula pCopula frankCopula gumbelCopula tau
#' @import pbivnorm


CompC = function(theta,t,X,W,ld,cop,dist){

  if (is.vector(X)){
    k = 1
    beta = theta[1:k]}
  else if(!is.vector(X)){
    k = dim(X)[2]
    beta = theta[1:k]
  }
  l = dim(W)[2]
  eta = theta[(k+1):(l+k)]
  nu = theta[(l+k+1)]
  gm = theta[(l+k+2)]
  lfun = (W%*%eta)


  # Distribution of T
  G1 = 1-exp(-ld*exp(X%*%beta))       # Cox model for T

  # Distribution of C

  if (dist == "Weibull"){               # Weibull margin
    G2 = 1-exp(-exp((log(t)-W%*%eta)/nu))
  }else if (dist == "lognormal")         # Log-normal margin
    G2 = pnorm((log(t)-W%*%eta)/nu)
  else
    stop("Only Weibull and lognormal distributions are supported for C")


  # avoid numerical issues
  G1[is.nan(G1)] = 1e-10
  G2[is.nan(G2)] = 1e-10

  G1[is.na(G1)] <- 1e-10
  G2[is.na(G2)] <- 1e-10
  G1 = pmax(1e-6,G1)
  G2 = pmax(1e-6,G2)
  G1 = pmin(G1,0.999999)
  G2 = pmin(G2,0.999999)

  PD = cbind(G1,G2)                         # Join distributions

  if (cop == "Frank")                      # Frank copula
  {
    frank.cop = frankCopula(gm, dim = 2)
    cp1 = (exp(- gm*G1)*(exp(- gm*G2)-1))/(exp(-gm)-1 + (exp(-gm*G1)-1)*(exp(-gm*G2)-1))
    cp1 = pmax(1e-10,cp1)
    z6 =  1-G1-G2+pCopula(PD,frank.cop)
    z6 =   pmax(z6,1e-10)
  }else if (cop == "Gumbel")                    # Gumbel copula
  {
    gumb.cop = gumbelCopula(gm, dim = 2)
    cp1 = pCopula(PD,gumb.cop)*((-log(G1))^gm+(-log(G2))^gm)^(-1+1/gm)*(-log(G1))^(gm-1)/G1
    cp1[is.na(cp1)] <- 1e-5
    cp1 = pmax(1e-5,cp1)
    z6 =  1-G1-G2+pCopula(PD,gumb.cop)
    z6 =   pmax(z6,1e-10)
  }else if (cop == "Normal")              # normal
  {
    cp1 = pnorm((qnorm(G2)-gm * qnorm(G1))/sqrt(1 - gm^2))
    p1 = qnorm(G1)
    p2 = qnorm(G2)
    c2 = cbind(p1,p2)

    z6 =  1-G1-G2+pbivnorm(c2,rho = gm)       # contribution to the likelihood from A
    z6 =   pmax(z6,1e-20)
  }else
    stop("Entered copula function is not supported")
  B0 = X%*%beta-ld*exp(X%*%beta)
  B1 = log(pmax(1e-10,1-cp1))-log(z6)
  tot = B0+B1
  return(tot)
}


#' @title Search function
#'
#' @description  Function to indicate position of t in observed survival time
#'
#' @param t fixed time t
#' @param T1 distinct observed survival time

SearchIndicate = function(t,T1){
  i = 1;
  m = length(T1);

  while(T1[i+1]<=t){
    i = i+1;
    if (i == m) break;
  }
  return(i)
}


#' @title Long format
#'
#' @description
#'  Change hazard and cumulative hazard to long format
#'
#' @param Z Observed survival time, which is the minimum of T, C and A, where A is the administrative censoring time.
#' @param T1 Distinct observed survival time
#' @param lhat Hazard function estimate
#' @param Lhat Cumulative hazard function estimate

Longfun = function(Z,T1,lhat,Lhat){
  n = length(Z)
  llong = rep(0,n)
  Llong = rep(0,n)

  for(i in 1:n){
    j = SearchIndicate(Z[i],T1)
    llong[i] = lhat[j]
    Llong[i] = Lhat[j]
  }
  long = cbind(llong,Llong)
}


#' @title Distance between vectors
#' @description This function computes distance between two vectors based on L2-norm
#' @param a First vector
#' @param b Second vector

Distance = function(a,b){
  x = b-a
  n = length(x)
  l2norm = sqrt(sum(x^2)/n)
  return(l2norm)
}

