
#' @title Follicular Cell Lymphoma
#'
#' @description Competing risk data set involving follicular cell lymphoma as given in the book of Pintilie (2006).
#' The data consist of 541 patients with early disease stage (I or II) and treated with radiation alone (RT) or with
#' radiation and chemotherapy (CMT). The endpoints of interest are what comes first: relapse of the
#' disease or death in remission.
#'
#' @references Pintilie, M. (2006), Competing Risks: A Practical Perspective, Chichester: Wiley.
#'
#'
#' @name follic
#' @docType data
#' @format A data frame with 541 rows and 7 variables:
#' \itemize{
#'   \item age:	age
#'   \item clinstg: clinical stage: 1=stage I, 2=stage II
#'   \item ch:	chemotherapy
#'   \item rt:	radiotherapy
#'   \item hgb: hemoglobin level in g/l
#'   \item time:	first failure time
#'   \item status:	censoring status; 0=censored, 1=relapse, 2=death
#' }
NULL
