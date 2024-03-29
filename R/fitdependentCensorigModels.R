
#' @title Fit Dependent Censoring Models
#'
#' @description This function allows to estimate the dependency parameter along all other model parameters. First, estimates the cumulative hazard function, and
#' then at the second stage it estimates other model parameters assuming that the cumulative hazard function is known. The details for
#' implementing the dependent censoring methodology can be found in Deresa and Van Keilegom (2023).
#'
#' @references Deresa and Van Keilegom (2023). Copula based Cox proportional hazards models for dependent censoring, Journal of the American Statistical Association (in press).
#'
#'

#' @param start Initial values for the finite dimensional parameters. If \code{start} is NULL, the initial values will be obtained
#' by fitting a Cox model for survival time T and a Weibull model for dependent censoring C.
#'
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T.
#' @param W Data matrix with covariates related to C. First column of W should be a vector of ones.
#' @param cop Which copula should be computed to account for dependency between T and C. This argument can take
#'   one of the values from \code{c("Gumbel", "Frank", "Normal")}.
#' @param dist The distribution to be used for the censoring time C. Only two distributions are allowed, i.e, Weibull
#' and lognormal distributions. With the value \code{"Weibull"} as the
#' default.
#' @param eps Convergence error. This is set by the user in such away that the desired convergence is met; the default is \code{eps = 1e-3}.
#' @param n.iter Number of iterations; the default is \code{n.iter = 20}. The larger the number of iterations, the longer the computational time.
#' @param bootstrap A boolean indicating whether to compute bootstrap standard errors for making inferences.
#' @param n.boot Number of bootstrap samples to use in the estimation of bootstrap standard errors if \code{bootstrap = TRUE}. The default is n.boot = 50. But, higher
#' values  of \code{n.boot} are recommended for obtaining good estimates of bootstrap standard errors.
#' @importFrom copula pCopula frankCopula gumbelCopula tau
#' @importFrom stats nlminb pnorm  qnorm
#' @importFrom survival coxph survreg Surv
#'
#' @return This function returns a fit of dependent censoring model; parameter estimates, estimate of the cumulative hazard function, bootstrap standard
#' errors for finite-dimensional parameters, the nonparametric cumulative hazard function, etc.
#'
#'
#'
#' @examples
#' \donttest{
#' # Toy data example to illustrate implementation
#' n = 300
#' beta = c(0.5)
#' lambd = 0.35
#' eta = c(0.9,0.4)
#' X = cbind(rbinom(n,1,0.5))
#' W = cbind(rep(1,n),rbinom(n,1,0.5))
#' # generate dependency structure from Frank
#' frank.cop <- copula::frankCopula(param = 5,dim = 2)
#' U = copula::rCopula(n,frank.cop)
#' T1 = (-log(1-U[,1]))/(lambd*exp(X*beta))           # Survival time
#' T2 = (-log(1-U[,2]))^(1.1)*exp(W%*%eta)            # Censoring time
#' A = runif(n,0,15)                                  # administrative censoring time
#' Z = pmin(T1,T2,A)
#' d1 = as.numeric(Z==T1)
#' d2 = as.numeric(Z==T2)
#' resData = data.frame("Z" = Z,"d1" = d1, "d2" = d2)   # should be data frame
#' colnames(W) <- c("ones","cov1")
#' colnames(X) <- "cov.surv"
#'
#' # Fit dependent censoring model
#'
#'fit <- fitDepCens(resData = resData, X = X, W = W, bootstrap = TRUE)
#'
#' # parameter estimates
#'
#' fit$parameterEstimates
#'
#' # summary fit results
#' summary(fit)
#'
#' # plot cumulative hazard vs time
#'
#' plot(fit$observedTime, fit$cumhazardFunction, type = "l",xlab = "Time",
#' ylab = "Estimated cumulative hazard function")
#'}
#'
#' @export


fitDepCens = function(resData,X,W,
                cop = c("Frank","Gumbel", "Normal"),
                dist = c("Weibull", "lognormal"), start = NULL, n.iter = 20, bootstrap = TRUE, n.boot = 50, eps = 1e-3){

  cop <- match.arg(cop)
  dist <- match.arg(dist)

  #  Verify column names of input dataframe resData
  if (!all(c("Z", "d1", "d2") %in% colnames(resData)))
    stop("Z, d1 and d2 arguments must be column names of resData")

  id <- match(c("Z", "d1", "d2"), colnames(resData))
  colnames(resData)[id] <- c("Z", "d1", "d2")

  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2

  if(is.vector(W))
    stop("W should be a matrix of dimension larger than 1")

  if (is.vector(X)) k = 1
  if(!is.vector(X)) k = dim(X)[2]
  l = dim(W)[2]

  if(!is.null(start)){
    if(start[(k+l+1)<=0]|start[(k+l+2)<=0])
      stop("Scale or dependency parameter cannot be negative")
  }
  if(is.null(start)){
    fit1 = coxph(Surv(Z,d1)~X)
    fit2 = survreg(Surv(Z,d2)~W-1)
    }

  if(is.null(start) & (cop == "Frank" | cop == "Gumbel")){                  # obtain initial values by fitting standard models
    start = c(fit1$coefficients,fit2$coefficients,fit2$scale,2)
  }else if (is.null(start) & (cop == "Normal")){
    start = c(fit1$coefficients,fit2$coefficients,fit2$scale,0.2)}

  # lb = lower boundary, ub = upper boundary
  if(cop == "Frank")                    # Frank copula
  {
    lb = c(rep(-Inf,k+l),0,1e-5)
    ub = c(rep(Inf,k+l+2))
  }else if(cop == "Gumbel")              # Gumbel copula
  {lb = c(rep(-Inf,k+l),0,1.0001)
    ub = c(rep(Inf,k+l+2))
  } else if(cop == "Normal")               # Normal
  {lb = c(rep(-Inf,k+l),0,-0.99)
   ub = c(rep(Inf,k+l+1),0.99)
  }else
    stop(paste0("The copula function",cop, "is not supported"))

  if(length(start)!= k+l+2)
    stop("The length of initial values is not equal to the number of parameters need to be estimated")


  res <- SolveL(start,resData,X,W,cop,dist)                          # Estimate the cumulative hazard function
  lsml <- res$lambda
  Llrge <- res$cumhaz
  T1 <- res$times
  longfrm <- Longfun(Z,T1,lsml,Llrge)
  lhat <- longfrm[,1]
  cumL <- longfrm[,2]

  # initial step of estimation

  parhat = nlminb(start = start,PseudoL,resData = resData,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par

  a = start
  b = parhat
  flag = 0

  # then do estimation iteratively until convergence

  while (Distance(b,a)>eps){                                        # doing this while loop until the desired convergence criteria are met
    a = b
    res <- SolveL(a,resData,X,W,cop,dist)
    lsml = res$lambda
    Llrge = res$cumhaz
    T1 = res$times
    longfrm = Longfun(Z,T1,lsml,Llrge)
    lhat = longfrm[,1]
    cumL = longfrm[,2]
    parhat = nlminb(start = a,PseudoL,resData = resData,X = X,W = W,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
    b = parhat
    flag = flag+1;
    if (flag>n.iter)                                  # Stop after iteration 20; this usually gives sufficient convergence results
    {
      flag=0;
      warning("The maximum number of iterations reached before convergence criteria is satisified. Better convergence may be obtained by increasing n.iter")
      warning(paste0("The current convergence error (eps) is ", Distance(b,a)))
      break;
    }
  }
  cumHat <- SolveL(b,resData,X,W,cop,dist)

  # nonparametric bootstrap for making inference

  paramsBootstrap = list("init" = parhat,"resData" = resData, "X" = X, "W" = W,"lhat" = lhat,"cumL" = cumL, "dist" = dist,"k" = k, "lb" = lb, "ub" = ub, "Obs.time" = T1, "cop" = cop, "n.boot" = n.boot, "n.iter" = n.iter, "eps" = eps)

  fitObj <- NULL

  if(bootstrap)                                       # Obtain bootstrap standard error
    fitObj <- do.call(boot.fun, paramsBootstrap)

  # Transform the copula parameter to tau scale
  nn <- length(parhat)
  if (cop == "Frank") {parhat[nn] <- tau(frankCopula(parhat[nn]))
  }else if (cop == "Gumbel"){ parhat[nn] <- tau(gumbelCopula(parhat[nn]))
  }else parhat[nn] <- 2*asin(parhat[nn])/pi
  names(parhat)[nn] <- "tau"
  names(parhat)[nn-1] <- "sigma"

  depObj <- c(list("parameterEstimates" = parhat,"copula" = cop, "censoringDistribution" = dist, "bootstrap" = bootstrap, "dimX" = k, "dimW" = l,"observedTime" = cumHat$times,"hazardFunction" = cumHat$lambda, "cumhazardFunction" = cumHat$cumhaz),fitObj)

  class(depObj) <- append(class(depObj), "depFit")                   # dependent censoring fit object

  return(depObj)
}
