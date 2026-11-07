
# Rprobit: Estimation of Multinomial Probit Models in R

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/Rprobit)](https://cran.r-project.org/package=Rprobit)
<!-- badges: end -->

The **Rprobit** package provides a robust set of tools for estimating
multinomial probit (MNP) models for panel data in R. Multinomial probit
models are a flexible approach for modeling discrete choice behavior
when the Independence of Irrelevant Alternatives (IIA) assumption
required by multinomial logit models is too restrictive.

## Features

- **Numerically Optimized Estimation:** Offers maximum approximate
  likelihood estimation and the maximum approximate composite marginal
  likelihood (MaCML) approach of Chandra Bhat using the Solow-Joe
  approximation of the Gaussian CDF. An initialization procedure
  accelerates the numerical search for the likelihood optimum.
- **Correlation Structure:** Allows for general covariance structures in
  the error terms to capture dependencies between alternatives.
- **Unbalanced Panel Data:** The package can handle cross-sectional and
  panel data sets. Each individual may take a different number of
  choices (unbalanced panel setting).
- **Gaussian Mixing:** Normally distributed random coefficients are
  supported. The covariance structure of the mixing distributions can be
  chosen freely.

## Installation

You can install the stable version of **Rprobit** from CRAN:

``` r
install.packages("Rprobit")
```

## Example

This is a basic example that simulates, using randomly drawn parameters,
a panel data set of $N=100$ deciders, each making $T=10$ decisions
between three commuting options `A`, `B`, and `C`, where the choice
depends on price and travel time. The model has the following features:

- It contains randomly distributed alternative specific constants (ASCs)
  for alternatives `B` and `C` with a diagonal variance-covariance
  matrix.
- The error term has a fully general variance matrix (normalized for
  level by restricting the error of the first alternative to zero;
  normalized for scale by fixing the variance corresponding to the
  second alternative to 1).
- The coefficient for price is identical for all alternatives, while the
  travel time is perceived differently for each alternative.
- The model parameters `theta` are drawn randomly from a continuous
  uniform distribution on $[0,1]$. Coefficients for price and time are
  set negative.

Estimation is performed using the full pairwise CML function.

``` r
library(Rprobit)

### model formulation 
formula <- choice ~ price | 1 | time
m_basic <- mod$new(formula, re = "ASC", alt = 3)

### generate a simulated data set from the model 
set.seed(1)
d_sim <- choice_data$new()
theta <- runif(m_basic$length_theta) # DGP parameters
theta[1:4] <- -theta[1:4] # negative price and time coefficients
d_sim$simulate(
  mod = m_basic,
  theta = theta, 
  N = 100,
  Tp = 10,
  sd = 2 # standard deviation of centered normal distribution for covariates
) |> head()
#>   id occasion choice    price_A   price_B    price_C     time_A    time_B
#> 1  1        1      C -1.6409368 -3.325298  1.3928352  1.8665019  1.180957
#> 2  1        2      C  0.9748581  1.621680  0.1070314  3.6121785 -2.826140
#> 3  1        3      C  1.4766494 -3.824692 -2.6205670 -0.1130073  3.220683
#> 4  1        4      C  1.1515627 -2.493507 -4.2461321  3.7718226  3.680885
#> 5  1        5      C -0.6107768  1.996309 -0.4161572  3.1567669  2.736596
#> 6  1        6      B  3.0235623 -1.081745 -0.6255732  1.0045692 -2.511146
#>       time_C
#> 1 -2.4754719
#> 2 -2.4603291
#> 3  1.1955818
#> 4  0.5977288
#> 5 -0.2202787
#> 6 -1.6153500

### set up the model structure 
r <- Rprobit$new(mod = m_basic, choice_data = d_sim)

### estimation
r$fit(cml_type = "FP")
#> Looking for initial values...
#> Starting numerical optimization...
#> Calculating approximated Hessian at estimate...
#> Done.

### result overview
summary(r)
#> 
#> Rprobit model summary
#> 
#> Formula:
#> choice ~ price | 1 | time 
#> 
#> Data:
#>   Deciders: 100 
#>   Choice occasions: 1000 
#>   Choice occasions per decider: 10 
#>   Alternatives: 3 (A, B, C)
#> 
#> Response residuals:
#>                 A         B          C
#> Min    -9.519e-01 -0.950112 -9.800e-01
#> 1Q     -7.167e-02 -0.136775 -4.002e-02
#> Median -3.437e-03 -0.008498  6.119e-08
#> 3Q     -4.266e-09  0.048243  6.429e-02
#> Max     9.814e-01  0.981978  9.362e-01
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> price              -0.30126    0.04077  -7.389 1.48e-13 ***
#> time_A             -0.41507    0.05392  -7.698 1.38e-14 ***
#> time_B             -0.64128    0.07188  -8.921  < 2e-16 ***
#> time_C             -0.99729    0.13390  -7.448 9.46e-14 ***
#> ASC_B               0.17337    0.13920   1.245  0.21296    
#> ASC_C               0.95579    0.12731   7.508 6.02e-14 ***
#> Omega[ASC_B,ASC_B]  0.93990    0.29397   3.197  0.00139 ** 
#> Omega[ASC_C,ASC_B]  0.00000         NA      NA       NA    
#> Omega[ASC_C,ASC_C]  0.38252    0.18484   2.070  0.03850 *  
#> Sigma[A,A]          0.00000         NA      NA       NA    
#> Sigma[B,A]          0.00000         NA      NA       NA    
#> Sigma[C,A]          0.00000         NA      NA       NA    
#> Sigma[B,B]          1.00000         NA      NA       NA    
#> Sigma[C,B]          0.61777    0.07601   8.128 4.37e-16 ***
#> Sigma[C,C]          0.38164    0.09402   4.059 4.93e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random-effect covariance matrix (Omega):
#>        ASC_B  ASC_C
#> ASC_B 0.9399 0.0000
#> ASC_C 0.0000 0.3825
#> 
#> Error covariance matrix (Sigma):
#>   A      B      C
#> A 0 0.0000 0.0000
#> B 0 1.0000 0.6178
#> C 0 0.6178 0.3816
#> 
#> Prediction summary:
#>                                 A      B      C
#> Observed choice frequency  0.2020 0.2890 0.5090
#> Mean predicted probability 0.2036 0.2896 0.5068
#> 
#> Fit statistics:
#>   Log-likelihood: -446.1 
#>   AIC: 912.3 
#>   BIC: 961.3

### comparison between DGP parameters and estimates
rbind(
  "DGP" = theta, 
  "estimate" = r$theta_fit, 
  "deviation" = abs(theta - r$theta_fit)
) |> round(2)
#>            [,1]  [,2]  [,3]  [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> DGP       -0.27 -0.37 -0.57 -0.91 0.20 0.90 0.94 0.66 0.63  0.06
#> estimate  -0.30 -0.42 -0.64 -1.00 0.17 0.96 0.97 0.62 0.62  0.00
#> deviation  0.04  0.04  0.07  0.09 0.03 0.06 0.02 0.04 0.01  0.06
```
