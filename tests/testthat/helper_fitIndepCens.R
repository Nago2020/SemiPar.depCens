# library(depCens)
# library(parallel)
# library(copula)
# library(survival)
# library(pbivnorm)

data <- genData(500)
resData <- data$resData
X = data$X
W = data$W

fitI <- fitIndepCens(resData = data$resData, X = X, W = W, bootstrap = FALSE)
dist <- fitI$censoringDistribution
haz <- fitI$hazardFunction
parI <- fitI$parameterEstimates
m = length(parI)



