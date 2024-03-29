library(SemiPar.depCens)
library(parallel)
library(copula)
library(survival)
library(pbivnorm)

data <- genData(500)
resData <- data$resData
X = data$X
W = data$W

fit <- fitDepCens(resData = data$resData, X = X, W = W, bootstrap = TRUE)
parE <- fit$parameterEstimates
cop <- fit$copula
dist <- fit$censoringDistribution
haz <- fit$hazardFunction
l = length(parE)





