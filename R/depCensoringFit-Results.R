

#' Summary of \code{depCensoringFit} object
#'
#' @param object Output of \code{\link{fitDepCens}} function
#' @param ... Further arguments
#' @export
#'
#'

summary.depFit <- function(object, ...) {

  message(rep("-", 100))
  message("Summary of dependent censoring model")



 # Sys.sleep(1)

  k = object$dimX
  l = object$dimW


  if(object$bootstrap){
    stError <- object$parEst.std.error
    stError_ch <- object$cumHazard.std.error
    Z <- object$parameterEstimates/stError
    pvalue <- 2*(1-pnorm(abs(Z)))
    out.res <- cbind(object$parameterEstimates,stError,pvalue)

    message(rep("-", 100))
    cat("\n")
    cat("Survival submodel: Cox proportional hazards model")
    cat("\n \n")
    surv.out <- out.res[1:k,]
    if(dim(surv.out)[1]>=2){
      if(is.null(rownames(surv.out))){
        colName <- cbind(c(rep("X",k)),c(seq(1:k)))
        rownames(surv.out) <- noquote(apply(colName,1,paste, collapse = "_"))
      }else {
        rownames(surv.out) <- gsub("X", "", rownames(surv.out))
      }
      colnames(surv.out) <- c("Estimate", "Boot.SE", "Pvalue")
    }
     surv.out <- round(surv.out, 3)
     print(surv.out)

    cat("\n \n")
    cat("Censoring submodel: ", object$censoringDistribution)
    cat("\n \n")
    cens.out <- out.res[(k+1):(l+k),]

    if(is.null(rownames(cens.out))){
      colName <- cbind(c(rep("W",l)),c(0,seq(1:(l-1))))
      colC <- noquote(apply(colName,1,paste, collapse = "_"))
      cens.out <- rbind(cens.out,out.res[(l+k+1),])
      rownames(cens.out) <- c("Intercept",colC[-1], "sigma")
    }else {
      colC <- rownames(cens.out)
      cens.out <- rbind(cens.out,out.res[(l+k+1),])
      rownames(cens.out) <- c("Intercept",colC[-1], "sigma")
      rownames(cens.out) <- gsub("W", "", rownames(cens.out))
    }
    colnames(cens.out) <- c("Estimate", "Boot.SE", "Pvalue")
    print(round(cens.out, 3))


    message(rep("-", 100))
    cat("Assumed copula model:", object$copula)
    cat("\n \n")

    cat("kendall's tau correlation :\n")
    cat("\n")
    tau = round(out.res[(l+k+2),],3)
    names(tau) <- c("tau","Boot.SE", "Pvalue")
    print(tau)
  }else{
    out.res <- object$parameterEstimates

    message(rep("-", 100))
    cat("Survival submodel: Cox proportional hazards model")
    cat("\n \n")
    cat("Parameter estimates:")
    cat("\n")
    surv.out <- out.res[1:k]
    if(k>=2){
    if(is.null(names(surv.out))){
      colName <- cbind(c(rep("X",k)),c(seq(1:k)))

      names(surv.out)<- noquote(apply(colName,1,paste, collapse = "_"))
    }else{
      names(surv.out) <- gsub("X", "", names(surv.out))
    }
    }else if(k==1){
      if(is.null(names(surv.out))){
        names(surv.out) <- "X1"
      }else  names(surv.out) <- gsub("X","",names(surv.out))
    }
    print(round(surv.out, 3))

    cat("\n \n")
    cat("Censoring submodel: ", object$censoringDistribution)
    cat("\n \n")
    cens.out <- out.res[(k+1):(l+k+1)]
    if(is.null(names(cens.out))){
      colName <- cbind(c(rep("W",l)),c(seq(1:l)))

      colN <- noquote(apply(colName,1,paste, collapse = "_"))
      names(cens.out) <- c("Intercept",colN[-1],"sigma")
    }else{
      colC <- names(cens.out)
      names(cens.out) <- c("Intercept",colC[-1])
      names(cens.out) <- gsub("W","",names(cens.out))
    }
    print(round(cens.out, 3))

    message(rep("-", 100))
    cat("Assumed copula model:", object$copula)
    cat("\n \n")

    cat("kendall's tau correlation :")
    cat("\n")
    tau = round(out.res[(l+k+2)],3)
    names(tau) <- c("tau")
    print(tau)
  }
}


#' Summary of \code{indepCensoringFit} object
#'
#' @param object Output of \code{\link{fitIndepCens}} function
#' @param ... Further arguments
#' @export
#'
#'

summary.indepFit <- function(object, ...) {

  message(rep("-", 100))
  message("Summary of independent censoring model")

   k = object$dimX
   l = object$dimW

  if(object$bootstrap){
    stError <- object$parEst.std.error
    stError_ch <- object$cumHazard.std.error
    Z <- object$parameterEstimates/stError
    pvalue <- 2*( 1-pnorm(abs(Z)))
    out.res <- cbind(object$parameterEstimates,stError,pvalue)

    message(rep("-", 100))
    cat("Survival submodel: Cox proportional hazards model")
    cat("\n \n")
    surv.out <- out.res[1:k,]
    if(dim(surv.out)[1]>=2){
      if(is.null(rownames(surv.out))){
        colName <- cbind(c(rep("X",k)),c(seq(1:k)))
        rownames(surv.out) <- noquote(apply(colName,1,paste, collapse = "_"))
      }else {
        rownames(surv.out) <- gsub("X", "", rownames(surv.out))
      }
      colnames(surv.out) <- c("Estimate", "Boot.SE", "Pvalue")
    }
    print(round(surv.out, 3))

    cat("\n \n")
    cat("Censoring submodel: ", object$censoringDistribution)
    cat("\n \n")
    cens.out <- out.res[(k+1):(l+k),]

    if(is.null(rownames(cens.out))){
      colName <- cbind(c(rep("W",l)),c(0,seq(1:(l-1))))
      colC <- noquote(apply(colName,1,paste, collapse = "_"))
      cens.out <- rbind(cens.out,out.res[(l+k+1),])
      rownames(cens.out) <- c("Intercept",colC[-1], "sigma")
    }else {
      colC <- rownames(cens.out)
      cens.out <- rbind(cens.out,out.res[(l+k+1),])
      rownames(cens.out) <- c("Intercept",colC[-1], "sigma")
      rownames(cens.out) <- gsub("W", "", rownames(cens.out))
    }
    colnames(cens.out) <- c("Estimate", "Boot.SE", "Pvalue")
    print(round(cens.out, 3))
  }else{
    out.res <- object$parameterEstimates

    message(rep("-", 100))
    cat("Survival submodel: Cox proportional hazards model")
    cat("\n\n")
    cat("Parameter estimates:")
    cat("\n")
    surv.out <- out.res[1:k]
    if(k>=2){
      if(is.null(names(surv.out))){
        colName <- cbind(c(rep("X",k)),c(seq(1:k)))

        names(surv.out)<- noquote(apply(colName,1,paste, collapse = "_"))
      }else{
        names(surv.out) <- gsub("X", "", names(surv.out))
      }
    }else if(k==1){
      if(is.null(names(surv.out))){
        names(surv.out) <- "X1"
      }else{
        names(surv.out) <- gsub("X", "", names(surv.out))
      }
    }
    print(round(surv.out, 3))
    cat("\n \n")
   # message(rep("-", 100))
    cat("Censoring submodel: ", object$censoringDistribution)
    cat("\n \n")
    cat("Parameter estimates:")
    cat("\n")
    cens.out <- out.res[(k+1):(l+k+1)]
    if(is.null(names(cens.out))){
      colName <- cbind(c(rep("W",l)),c(seq(1:l)))

      colN <- noquote(apply(colName,1,paste, collapse = "_"))
      names(cens.out) <- c("Intercept",colN[-1],"sigma")
    }else{
      colC <- names(cens.out)
      names(cens.out) <- c("Intercept",colC[-1])
      names(cens.out) <- gsub("W","",names(cens.out))
    }
    print(round(cens.out, 3))
  }
}

