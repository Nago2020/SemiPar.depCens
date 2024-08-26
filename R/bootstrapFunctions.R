
#' @title Nonparametric bootstrap approach for the dependent censoring model

#' @description This function estimates the bootstrap standard errors for the finite-dimensional model parameters and for the non-parametric cumulative
#' hazard function. Parallel computing using foreach has been used to speed up the estimation of standard errors.
#'
#'
#' @param init Initial values for the finite dimensional parameters obtained from the fit of \code{\link{fitDepCens}}
#' @param lhat Initial values for the hazard function obtained from the fit of \code{\link{fitDepCens}} based on the original data.
#' @param cumL Initial values for the cumulative hazard function obtained from the fit of \code{\link{fitDepCens}} based on the original data.
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T
#' @param W Data matrix with covariates related to C. First column of W should be a vector of ones
#' @param k Dimension of X
#' @param lb lower boundary for finite dimensional parameters
#' @param ub Upper boundary for finite dimensional parameters
#' @param cop Which copula should be computed to account for dependency between T and C. This argument can take
#'   one of the values from \code{c("Gumbel", "Frank", "Normal")}.
#' @param dist The distribution to be used for the  dependent censoring time C. Only two distributions are allowed, i.e, Weibull
#' and lognormal distributions. With the value \code{"Weibull"} as the
#'   default.
#' @param Obs.time Observed survival time, which is the minimum of T, C and A, where A is the administrative censoring time.
#' @param eps Convergence error. This is set by the user in such away that the desired convergence is met; the default is \code{eps = 1e-3}
#' @param n.iter Number of iterations; the default is \code{n.iter = 20}. The larger the number of iterations, the longer the computational time.
#' @param n.boot Number of bootstraps to use in the estimation of bootstrap standard errors.
#' @param ncore The number of cores to use for parallel computation is configurable, with the default \code{ncore = 7}.
#' @importFrom stats nlminb pnorm  qnorm sd
#' @importFrom copula pCopula frankCopula gumbelCopula tau
#' @import foreach
#' @import parallel
#' @return Bootstrap standard errors for parameter estimates and for estimated cumulative hazard function.
#'


boot.fun = function(init,resData,X,W,lhat, cumL,dist,k,lb, ub, Obs.time,cop,n.boot, n.iter, ncore, eps){
  B = n.boot                                     # number of bootstrap samples
  n.cores <- ncore

  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)

  boot = foreach(b = 1:B, .packages= c('survival', 'copula', 'pbivnorm'), .export = c("PseudoL", "SolveL","Distance", "CompC", "Longfun", "SearchIndicate")) %dopar% {
    samp1 = sample(length(resData$Z),replace = TRUE)
    resData_b = resData[samp1,]
    Zb = resData_b$Z
    if(k==1){
      Xb = X[samp1]
    }else{
     Xb = X[samp1,]
     }
    Wb = W[samp1,]

    # Initial step of estimation

    parhatb = nlminb(start = init, PseudoL, resData = resData_b,X = Xb,W = Wb,lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par
    aboot = init
    bboot = parhatb
    flag = 0

    while (Distance(bboot,aboot)>eps){                            # doing this while loop until the desired convergence criteria are met
      aboot = bboot
      res = SolveL(aboot,resData_b,Xb,Wb,cop,dist)
      lsml = res$lambda
      Llrge = res$cumhaz
      T1 = res$times
      longfrm = Longfun(Zb,T1,lsml,Llrge)
      lhat = longfrm[,1]
      cumL = longfrm[,2]
      parhatb = nlminb(start = aboot, PseudoL, resData = resData_b ,X = Xb,W = Wb, lhat = lhat,cumL = cumL,cop = cop,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par

      bboot = parhatb
      flag = flag+1;
      if (flag>n.iter)                                                              # Stop after iteration 50; this usually gives sufficient convergence results
      {
        break;
      }
    }
    cumHat = SolveL(bboot,resData_b,Xb,Wb,cop,dist)
    parb = bboot
    nn = length(parb)
    lboot = cumHat$lambda
    Lboot = cumHat$cumhaz
    Tb1 = cumHat$times
    longB = Longfun(Obs.time,Tb1,lboot, Lboot)
    lhatb = longB[,1]
    cumH.boot = longB[,2]

    # Transform the copula parameter to tau scale
    if (cop == "Frank"){ parb[nn] = tau(frankCopula(parb[nn]))
    }else if (cop == "Gumbel"){ parb[nn] = tau(gumbelCopula(parb[nn]))
    }else parb[nn] = 2*asin(parb[nn])/pi
    names(parb)[nn] <- "tau"
    names(parb)[(nn-1)] <- "sigma"

    list("beta.star" = parb, cumH.star = cumH.boot)
  }
  parallel::stopCluster(cl = my.cluster)

  beta.star = t(sapply(boot, function(x) x$beta.star))
  cumH.star =  t(sapply(boot, function(x) x$cumH.star))
  Bootse = apply(beta.star,2,sd)
  cumHse = apply(cumH.star,2,sd)
  bootR = list("parEst.std.error" = Bootse, "cumHazard.std.error" = cumHse)
  return(bootR)
}



#' @title  Nonparametric bootstrap approach for the independent censoring model

#' @description This function estimates the bootstrap standard errors for the finite-dimensional model parameters and for the non-parametric cumulative
#' hazard function under the assumption of independent censoring. Parallel computing using foreach has been used to speed up the computation.
#'
#'
#' @param init Initial values for the finite dimensional parameters obtained from the fit of \code{\link{fitIndepCens}}
#' @param lhat Initial values for the hazard function obtained from the fit of \code{\link{fitIndepCens}} based on the original data
#' @param cumL Initial values for the cumulative hazard function obtained from the fit of \code{\link{fitIndepCens}} based on the original data
#' @param resData Data matrix with three columns;  Z = the observed survival time, d1 = the censoring indicator of T
#' and  d2 =  the censoring indicator of C.
#' @param X Data matrix with covariates related to T
#' @param k Dimension of X
#' @param lb lower boundary for finite dimensional parameters
#' @param ub Upper boundary for finite dimensional parameters
#' @param W Data matrix with covariates related to C. First column of W should be a vector of ones
#' @param dist The distribution to be used for the  dependent censoring time C. Only two distributions are allowed, i.e, Weibull
#' and lognormal distributions. With the value \code{"Weibull"} as the
#'   default
#' @param Obs.time Observed survival time, which is the minimum of T, C and A, where A is the administrative censoring time.
#' @param eps Convergence error. This is set by the user in such away that the desired convergence is met; the default is \code{eps = 1e-3}
#' @param n.iter Number of iterations; the default is \code{n.iter = 20}. The larger the number of iterations, the longer the computational time
#' @param n.boot Number of bootstraps to use in the estimation of bootstrap standard errors.
#' @param ncore The number of cores to use for parallel computation is configurable, with the default \code{ncore = 7}.
#' @importFrom stats nlminb pnorm  qnorm sd
#' @importFrom survival coxph survreg
#' @import foreach
#'
#' @return Bootstrap standard errors for parameter estimates and for estimated cumulative hazard function.
#'


boot.funI = function(init,resData,X,W,lhat, cumL,dist,k,lb,ub, Obs.time,n.boot, n.iter, ncore, eps){
  B = n.boot                                     # number of bootstrap samples
  n.cores <- ncore

  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)

  boot = foreach(b = 1:B, .packages= c('survival'), .export = c("LikCopInd", "SolveLI","Distance", "Longfun", "SearchIndicate")) %dopar% {
    samp1 = sample(length(resData$Z),replace = TRUE)
    resData_b = resData[samp1,]
    Zb = resData_b$Z
    if(k==1) {Xb = X[samp1]
    }else Xb = X[samp1,]
    Wb = W[samp1,]

    # Initial step of estimation
    parhatb  = nlminb(start = init,LikCopInd,resData = resData_b,X = Xb,W = Wb,lhat = lhat,cumL = cumL, dist = dist, lower = lb, upper =  ub, control = list(eval.max=300,iter.max=200))$par
    aboot = init
    bboot = parhatb
    flag = 0

    while (Distance(bboot,aboot)>eps){                                   # doing this while loop until the desired convergence criteria are met
      aboot = bboot
      res = SolveLI(aboot,resData_b,Xb)
      lsml = res$lambda
      Llrge = res$cumhaz
      T1 = res$times
      longfrm = Longfun(Zb,T1,lsml,Llrge)
      lhat = longfrm[,1]
      cumL = longfrm[,2]
      parhatb = nlminb(start = aboot, LikCopInd, resData = resData_b,X = Xb,W = Wb,lhat = lhat,cumL = cumL,dist = dist,lower = lb ,upper =  ub, control = list(eval.max=300,iter.max=200))$par

      bboot = parhatb
      flag = flag+1;
      if (flag>n.iter)
      {
        break;
      }
    }
    cumHat = SolveLI(bboot,resData_b,Xb)
    parb = bboot
    nn <- length(parb)
    names(parb)[nn] <- "sigma"

    lboot = cumHat$lambda
    Lboot = cumHat$cumhaz
    Tb1 = cumHat$times
    longB = Longfun(Obs.time,Tb1,lboot,Lboot)
    lhatb = longB[,1]
    cumH.boot = longB[,2]

    list("beta.star" = parb, cumH.star = cumH.boot)
  }
  parallel::stopCluster(cl = my.cluster)

  beta.star = t(sapply(boot, function(x) x$beta.star))
  cumH.star =  t(sapply(boot, function(x) x$cumH.star))
  Bootse = apply(beta.star,2,sd)
  cumHse = apply(cumH.star,2,sd)
  bootR = list("parEst.std.error" = Bootse, "cumHazard.std.error" = cumHse)
  return(bootR)
}
