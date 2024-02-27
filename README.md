
# SemiPar.depCens

The goal of *SemiPar.depCens* package is to provide easy to use
functions in *R* for estimation of the dependent censoring methodology
proposed by Deresa and Van Keilegom (2024)
<doi:10.1080/01621459.2022.2161387>. The approach presented in the
latter paper is based on a parametric copula for the relation between
the survival time and the dependent censoring time, and the parameter
defining the copula does not need to be known. Instead, the copula
parameter is estimated jointly with other finite model parameters by
maximizing a Pseudo likelihood function. Available copula functions in
*SemiPar.depCens* package include Frank, Gumbel and Normal copulas. Only
Weibull and Lognormal models are allowed for the censoring model, even
though any parametric model that satisfies certain identifiability
conditions could be used.

## Installation

The development version of *SemiPar.depCens* package is available on
github. To install it in R, use:

``` r
devtools::install_github("Nago2020/SemiPar.depCens")
```

## Example

This is a basic example which shows how to use the package in practice:

``` r
# load packages
library(copula)
library(survival)
library(stats)
library(foreach)
library(pbivnorm)
library(SemiPar.depCens)

## load the data
data("follic")
```

``` r

#Prepare the data in the way that is used by the package

follic = follic[order(follic$time),]                     # order the observed survival time
Z = round(follic$time,digits = 3)
d1 = as.numeric(follic$status==1)                        # censoring indicator for survival time T
d2 = as.numeric(follic$status==2)                        # censoring indicator for dependent censoring C
treat = as.numeric(follic$ch=="Y")                       # treatment indicator
age = (follic$age-mean(follic$age))/sd(follic$age)       # recommended to standardize continuous variables
hgb = (follic$hgb-mean(follic$hgb))/sd(follic$hgb)       # standardized hemoglobin
clinstg = as.numeric(follic$clinstg==1)                  # clinical stage

X = cbind(treat,age,hgb,clinstg)                         # data matrix for T, should be in matrix form
W = cbind(rep(1,length(Z)),treat,age,hgb,clinstg)        # data matrix for C, should be in matrix form
resData = data.frame("Z" = Z,"d1" = d1, "d2" = d2)       # resData should be a data frame
```

### Fit dependent censoring model

The following code fit a default copula, which is Frank, for the
relation between the survival time (T) and dependent censoring time (C).
The default marginal model for C is a Weibull model. Other capabilities
can be explored by typing *?fitDepCens* in the console.

``` r

fitD <- fitDepCens(resData = resData, X = X, W = W, bootstrap = FALSE)    
```

The output for the above code chunk should look as below. Since
bootstrapping = FALSE, it does not make any inference based on p-values;
only parameter estimates are shown.

``` r
summary(fitD)
#> ----------------------------------------------------------------------------------------------------
#> Summary of dependent censoring model
#> ----------------------------------------------------------------------------------------------------
#> Survival submodel: Cox proportional hazards model
#>  
#> Parameter estimates:
#>   treat     age     hgb clinstg 
#>  -0.347   0.352   0.042  -0.646 
#> 
#>  
#> Censoring submodel:  Weibull
#>  
#> Intercept     treat       age       hgb   clinstg     sigma 
#>     2.803     0.125    -0.658    -0.038     0.411     0.618
#> ----------------------------------------------------------------------------------------------------
#> Assumed copula model: Frank
#>  
#> kendall's tau correlation :
#>   tau 
#> 0.336
```

We can do bootstrapping by setting bootstrap = TRUE, but note that the
algorithm may take long time to finish the computations, even after
parallelization is used to speed up the work. The default number of
bootstrap size is 50. Increasing number of bootstrap samples may produce
more precise standard error estimates.

``` r
fitD <- fitDepCens(resData = resData,X = X, W = W, bootstrap = TRUE, n.boot = 50)    
summary(fitD)
#> ----------------------------------------------------------------------------------------------------
#> Summary of dependent censoring model
#> ----------------------------------------------------------------------------------------------------
#> 
#> Survival submodel: Cox proportional hazards model
#>  
#>         Estimate Boot.SE Pvalue
#> treat     -0.347   0.170  0.042
#> age        0.352   0.067  0.000
#> hgb        0.042   0.058  0.463
#> clinstg   -0.646   0.119  0.000
#> 
#>  
#> Censoring submodel:  Weibull
#>  
#>           Estimate Boot.SE Pvalue
#> Intercept    2.803   0.135  0.000
#> treat        0.125   0.135  0.356
#> age         -0.658   0.071  0.000
#> hgb         -0.038   0.051  0.452
#> clinstg      0.411   0.125  0.001
#> sigma        0.618   0.046  0.000
#> ----------------------------------------------------------------------------------------------------
#> Assumed copula model: Frank
#>  
#> kendall's tau correlation :
#> 
#>     tau Boot.SE  Pvalue 
#>   0.336   0.102   0.001
```

### Fit independent censoring model

For independent censoring model, the assumption is that the copula
parameter between T and C is zero. Hence, the model is very simplified
in terms of computational costs. We obtain results very quickly in
comparison to dependent censoring model. The default model for censoring
distribution is Weibull.

``` r

fitI<- fitIndepCens(resData = resData, X = X, W = W, bootstrap = TRUE, n.boot = 50)                       
summary(fitI)
#> ----------------------------------------------------------------------------------------------------
#> Summary of independent censoring model
#> ----------------------------------------------------------------------------------------------------
#> Survival submodel: Cox proportional hazards model
#>  
#>         Estimate Boot.SE Pvalue
#> treat     -0.365   0.174  0.036
#> age        0.329   0.066  0.000
#> hgb        0.041   0.061  0.498
#> clinstg   -0.648   0.123  0.000
#> 
#>  
#> Censoring submodel:  Weibull
#>  
#>           Estimate Boot.SE Pvalue
#> Intercept    3.208   0.143  0.000
#> treat        0.127   0.173  0.465
#> age         -0.703   0.087  0.000
#> hgb         -0.029   0.070  0.677
#> clinstg      0.333   0.160  0.037
#> sigma        0.608   0.054  0.000
```
