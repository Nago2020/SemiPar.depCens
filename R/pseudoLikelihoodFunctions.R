
#' @title Likelihood function under dependent censoring
#'
#' @description
#' The \code{PseudoL} function  is maximized in order to
#' estimate the finite dimensional model parameters, including the dependency parameter.
#' This function assumes that the cumulative hazard function is known.
#'
#' @inheritParams SolveL
#' @param lhat The estimated hazard function obtained from the output of \code{\link{SolveL}}.
#' @param cumL The estimated cumulative hazard function from the output of \code{\link{SolveL}}.
#' @importFrom copula pCopula frankCopula gumbelCopula tau
#' @import pbivnorm

#' @return maximized log-likelihood value
#'


PseudoL = function(theta,resData,X,W,lhat,cumL,cop,dist){
  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2

  if (is.vector(X)){
    k = 1
    beta = theta[k]
    mu1 = X*beta
  }else if(!is.vector(X)){
    k = dim(X)[2]
    beta = theta[1:k]
    mu1 = X%*%beta}
  l = dim(W)[2]
  eta = theta[(k+1):(l+k)]
  nu = theta[(l+k+1)]
  gm = theta[(l+k+2)]
  y = log(Z)

  # distribution of T -- Cox PH

  g1 = lhat*exp(mu1)*exp(-cumL*exp(mu1))
  G1 = 1-exp(-cumL*exp(mu1))

  # dist of C

  if (dist == "Weibull"){
    g2 = 1/nu*exp((y-W%*%eta)/nu)*exp(-exp((y-W%*%eta)/nu))           # log-time scale
    G2 = 1-exp(-exp((y-W%*%eta)/nu))
  }
  else if (dist == "lognormal"){   # Lognormal
    m = (y-W%*%eta)/nu
    g2 = 1/(sqrt(2*pi)*nu*Z)*exp(-(m^2/2))
    G2 = pnorm(m)
  }else
    stop("The distribution of C is not supported")

  # Joint distribution
  # avoid NA

  g1[is.na(g1)] = 1e-10
  g2[is.na(g2)] = 1e-10
  g1[is.nan(g1)] = 1e-10
  g2[is.nan(g2)] = 1e-10
  G1[is.na(G1)] = 1e-10
  G2[is.na(G2)] = 1e-10
  G1[is.nan(G1)] = 1e-10
  G2[is.nan(G2)] = 1e-10
  g1[!is.finite(g1)] = 0
  g2[!is.finite(g2)] = 0
  G1 = pmax(G1,1e-10)
  G1 = pmin(G1,0.99999999)
  G2 = pmax(G2,1e-10)
  G2 = pmin(G2,0.99999999)

  PD = cbind(G1,G2)

  if (cop == "Frank")                                         # Frank copula
  {
    frank.cop = copula::frankCopula(gm, dim = 2)
    cp1 = (exp(- gm*G1)*(exp(- gm*G2)-1))/(exp(-gm)-1 + (exp(-gm*G1)-1)*(exp(-gm*G2)-1))
    cp2 = (exp(- gm*G2)*(exp(- gm*G1)-1))/(exp(-gm)-1 + (exp(-gm*G1)-1)*(exp(-gm*G2)-1))
    z6 =  1-G1-G2+copula::pCopula(PD,frank.cop)
    z6 =   pmax(z6,1e-10)
  }else if (cop == "Gumbel")                                       # Gumbel copula
  {
    gumb.cop = copula::gumbelCopula(gm, dim = 2)
    cp1 = copula::pCopula(PD,gumb.cop)*((-log(G1))^gm+(-log(G2))^gm)^(-1+1/gm)*(-log(G1))^(gm-1)/G1
    cp2 = copula::pCopula(PD,gumb.cop)*((-log(G1))^gm+(-log(G2))^gm)^(-1+1/gm)*(-log(G2))^(gm-1)/G2
    z6 =  1-G1-G2+copula::pCopula(PD,gumb.cop)
    z6 =   pmax(z6,1e-10)
  }else if (cop == "Normal")     # normal
  {
    cp1 = pnorm((qnorm(G2) - gm * qnorm(G1))/sqrt(1 - gm^2))
    cp2 = pnorm((qnorm(G1) - gm * qnorm(G2))/sqrt(1 - gm^2))
    p1 = qnorm(G1)
    p2 = qnorm(G2)
    c2 = cbind(p1,p2)
    z6 =  1-G1-G2+pbivnorm(c2,rho = gm)          # contribution to the likelihood from A
    z6 =   pmax(z6,1e-10)
  }
  cp1[is.na(cp1)] = 1e-20
  cp2[is.na(cp2)] = 1e-20
  cp1[is.nan(cp1)] = 1e-20
  cp2[is.nan(cp2)] = 1e-20
  cp1[!is.finite(cp1)] = 0
  cp2[!is.finite(cp2)] = 0
  d3 =  (1-d1-d2)             # Censoring
  term1 <- (1 - cp1)*g1       #  Contribution from T
  term2 <- (1 - cp2)*g2       #  Contribution from C

  Logn <- -sum(d1 * log(pmax(term1, 1e-10)) + d2 * log(pmax(term2, 1e-10)) + d3*log(z6))
  return(Logn)
}




#' @title Loglikehood function under independent censoring
#'
#' @inheritParams SolveLI
#' @param W Data matrix with covariates related to C. First column of W should be ones
#' @param lhat The estimated hazard function obtained from the output of \code{\link{SolveLI}}.
#' @param cumL The estimated cumulative hazard function from the output of \code{\link{SolveLI}}.
#' @param dist The distribution to  be used for the dependent censoring C. Only two distributions are allowed, i.e, Weibull
#' and lognormal distributions. With the value \code{"Weibull"} as the
#'   default.
#'  @importFrom stats nlminb pnorm  qnorm sd
#'
#' @return Maximized log-likelihood value
#'


LikCopInd <- function(theta,resData,X,W,lhat,cumL,dist){ # gamma = 0
  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2

  if (is.vector(X)){
    k = 1
    beta = theta[k]
    mu1 = X*beta
  }
  else if(!is.vector(X)){
    k = dim(X)[2]
    beta = theta[1:k]
    mu1 = X%*%beta
  }

  l = dim(W)[2]
  eta = theta[(k+1):(l+k)]
  mu2 = W%*%eta
  nu = theta[(l+k+1)]
  y = log(Z)


  # distribution of T

  g1 = lhat*exp(mu1)*exp(-cumL*exp(mu1))
  G1 = 1-exp(-cumL*exp(mu1))

  # dist of C

  if (dist == "Weibull"){                                           # Weibull
    g2 = 1/nu*exp((y-mu2)/nu)*exp(-exp((y-mu2)/nu))*(1/Z)           # density of C
    G2 = 1-exp(-exp((y-mu2)/nu))                                    # dist of C
  } else if (dist == "lognormal"){                                  # Log-normal
    m = (y-mu2)/nu
    g2 = 1/(sqrt(2*pi)*nu*Z)*exp(-(m^2/2))                          # density of C
    G2 = pnorm(m)                                                   # dist of C
  } else
    stop("The distribution of C is not supported")

  # Joint dist,

  G1 = pmax(1e-10,G1)
  G2 = pmax(1e-10,G2)

  Logn <- -sum(log(pmax(g1[d1==1], 1e-10)))-sum((log(1-G1[(1-d1)==1])))-sum(log(pmax(g2[d2==1], 1e-10)))-sum(log(1-G2[(1-d2)==1]))
  return(Logn)
}

