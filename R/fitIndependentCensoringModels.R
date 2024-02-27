
#' @title Fit Independent Censoring Models
#'
#' @description This function allows to estimate all model parameters under the assumption of independent censoring. First, estimates the cumulative hazard function, and
#' then at the second stage it estimates other model parameters assuming that the cumulative hazard is known.
#'
#'

#' @param start Initial values for the finite dimensional parameters. If \code{start} is NULL, the initial values will be obtained
#' by fitting a Cox model for survival time T and a Weibull model for censoring time C.
#'
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T.
#' @param W Data matrix with covariates related to C. First column of W should be ones.
#' @param dist The distribution to be used for the censoring time C. Only two distributions are allowed, i.e, Weibull
#' and lognormal distributions. With the value \code{"Weibull"} as the
#'   default.
#' @param eps Convergence error. This is set by the user in such away that the desired convergence is met; the default is \code{eps = 1e-3}.
#' @param n.iter Number of iterations; the default is \code{n.iter = 20}. The larger the number of iterations, the longer the computational time.
#' @param bootstrap A boolean indicating whether to compute bootstrap standard errors for making inferences.
#' @param n.boot Number of bootstrap samples to use in the estimation of bootstrap standard errors if \code{bootstrap = TRUE}. The default is n.boot = 50. But, higher
#' values  of \code{n.boot} are recommended for obtaining good estimates of bootstrap standard errors.
#' @importFrom stats nlminb pnorm  qnorm sd
#' @importFrom survival coxph survreg
#'
#' @return This function returns a fit of independent censoring model; parameter estimates, estimate of the cumulative hazard function, bootstrap standard
#' errors for finite-dimensional parameters, the nonparametric cumulative hazard function, etc.
#'
#'
#'
#' @examples
#' \dontrun{
#'
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
#' T1 = (-log(1-U[,1]))/(lambd*exp(X*beta))                  # Survival time'
#' T2 = (-log(1-U[,2]))^(1.1)*exp(W%*%eta)                   # Censoring time
#' A = runif(n,0,15)                                         # administrative censoring time
#' Z = pmin(T1,T2,A)
#' d1 = as.numeric(Z==T1)
#' d2 = as.numeric(Z==T2)
#' resData = data.frame("Z" = Z,"d1" = d1, "d2" = d2)     # should be data frame
#'
#' # Fit independent censoring model
#'
#' fitI <- fitIndepCens(resData = resData, X = X, W = W, bootstrap = FALSE)
#'
#' # parameter estimates
#'
#' fitI$parameterEstimates
#'
#' # summary fit results
#' summary(fitI)
#'
#' # plot cumulative hazard vs time
#'
#'  plot(fitI$observedTime, fitI$cumhazardFunction, type = "l",xlab = "Time",
#'  ylab = "Estimated cumulative hazard function")
#'}
#'
#' @export


fitIndepCens = function(resData,X,W,
                      dist = c("Weibull", "lognormal"), start = NULL, n.iter = 20, bootstrap = TRUE, n.boot = 50,  eps = 1e-3){

  dist <- match.arg(dist)

  if (!all(c("Z", "d1", "d2") %in% colnames(resData)))
    stop("Z, d1 and d2 arguments must be column names of resData")

  id <- match(c("Z", "d1", "d2"), colnames(resData))
  colnames(resData)[id] <- c("Z", "d1", "d2")

  Z = resData$Z
  d1 = resData$d1
  d2 = resData$d2

  if (is.vector(X)) k = 1
  else if(!is.vector(X)) k = dim(X)[2]
  l = dim(W)[2]

  if(!is.null(start)){
    if(start[(k+l+1)<=0])
      stop("Scale parameter cannot be negative")
  }

  if(is.null(start)){                                    # obtain initial values by fitting standard models
      fit1 = coxph(Surv(Z,d1)~X)
      fit2 = survreg(Surv(Z,d2)~W-1)
      start = c(fit1$coefficients,fit2$coefficients,fit2$scale)
  }
  if(length(start)!= k+l+1)
    stop("The length of initial values is not equal to the number of parameters need to be estimated")


  res <- SolveLI(start,resData,X)                 # estimate the cumulative hazard function
  lsml <- res$lambda
  Llrge <- res$cumhaz
  T1 <- res$times
  longfrm <- Longfun(Z,T1,lsml,Llrge)
  lhat <- longfrm[,1]
  cumL <- longfrm[,2]


  #lb = lower boundary, ub = upper boundary
  lb = c(rep(-Inf,k+l),0)
  ub = c(rep(Inf,k+l+1))

  # Initial estimate of theta

  parhat  = nlminb(start = start,LikCopInd, resData = resData,X = X,W = W,lhat = lhat,cumL = cumL, dist = dist, lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par

  a = start
  b = parhat
  flag = 0

  while (Distance(b,a)>eps){                             # doing this while loop until the desired convergence criteria are met
    a = b
    res = SolveLI(a,resData,X)
    lsml = res$lambda
    Llrge = res$cumhaz
    T1 = res$times
    longfrm = Longfun(Z,T1,lsml,Llrge)
    lhat = longfrm[,1]
    cumL = longfrm[,2]
    parhat  = nlminb(start = a,LikCopInd,resData = resData,X = X,W = W,lhat = lhat,cumL = cumL, dist = dist, lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par

    b = parhat
    flag = flag+1;
    if (flag>n.iter)                                    # stop after iteration 30; this usually gives sufficient convergence results
    {
      flag=0;
      print("The maximum number of iterations reached before convergence criteria is satisified. Better convergence may be obtained by increasing n.iter")
      print(paste0("The current convergence error (eps) is", Distance(b,a)))
      break;
    }
  }
  cumHat = SolveLI(parhat,resData,X)

  # nonparametric bootstrap for making inference

  paramsBootstrap = list("init" = parhat,"resData" = resData, "X" = X, "W" = W,"lhat" = lhat,"cumL" = cumL, "dist" = dist,"k" = k, "lb" = lb, "ub" = ub, "Obs.time" = T1, "n.boot" = n.boot, "n.iter" = n.iter, "eps" = eps)

  fitObj <- NULL

  if(bootstrap)                                                         # Obtain bootstrap standard error
    fitObj <- do.call(boot.funI, paramsBootstrap)

  nn <- length(parhat)
  names(parhat)[nn] <- "sigma"

  indObj <- c(list("parameterEstimates" = parhat, "censoringDistribution" = dist, "bootstrap" = bootstrap, "dimX" = k, "dimW" = l,"observedTime" = cumHat$times, "hazardFunction" = cumHat$lambda, "cumhazardFunction" = cumHat$cumhaz),fitObj)

  class(indObj) <- append(class(indObj), "indepFit")                                # independent censoring fit object

  return(indObj)
}
