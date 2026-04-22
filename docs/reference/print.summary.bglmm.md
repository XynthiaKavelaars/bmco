# Print Method for summary.bglmm Objects

Print Method for summary.bglmm Objects

## Usage

``` r
# S3 method for class 'summary.bglmm'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A `summary.bglmm` object returned by
  [`summary.bglmm`](https://xynthiakavelaars.github.io/bmco/reference/summary.bglmm.md).

- digits:

  Number of digits to display. Default is `3`.

- ...:

  Additional arguments (not used).

## Value

Invisibly returns `x`.

## See also

[`summary.bglmm`](https://xynthiakavelaars.github.io/bmco/reference/summary.bglmm.md),
[`print.bglmm`](https://xynthiakavelaars.github.io/bmco/reference/print.bglmm.md)

## Examples

``` r
# Uses the pre-computed example object shipped with the package:
print(summary(bglmm_fit))
#> 
#> Bayesian Multilevel Multivariate Logistic Regression
#> ======================================================
#> 
#> Multilevel Structure:  J = 20 clusters    n(placebo) = 150    n(drug) = 150
#> 
#> Group Estimates:
#>    Group mean y1 sd y1      95% CI y1 mean y2 sd y2      95% CI y2
#>  placebo   0.515 0.041 [0.435, 0.595]   0.408 0.037 [0.335, 0.481]
#>     drug   0.617 0.037 [0.540, 0.687]   0.592 0.038 [0.515, 0.666]
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.103, 0.184
#>   Delta SE   (y1, y2): 0.001, 0.001
#>   95% CI delta: y1 [0.002, 0.209]   y2 [0.081, 0.296]
#>   Posterior probability P(drug > placebo): 0.977
#> 
#> Fixed Effects (mean [SD]):
#>                      b11            b10            b01
#> group      1.904 [1.555] -0.031 [1.480]  1.569 [1.588]
#> age        0.025 [0.021] -0.014 [0.019] -0.017 [0.022]
#> group_age -0.012 [0.030]  0.011 [0.029] -0.016 [0.031]
#> 
#> Random Effects (population mean [SD]):
#>                      g11           g10           g01
#> Intercept -1.686 [1.126] 0.791 [0.977] 0.478 [1.109]
#> 
#> Variance Components (posterior mean):
#>   y1=1, y2=1 (b11):
#>           Intercept
#> Intercept     0.324
#>   y1=1, y2=0 (b10):
#>           Intercept
#> Intercept     0.153
#>   y1=0, y2=1 (b01):
#>           Intercept
#> Intercept     0.682
#> 
#> Prior Specification:
#>   Fixed effects -- Normal prior:
#>     Mean:      0, 0, 0 
#>     Variance:  10, 10, 10 
#>   Random effects -- Normal prior:
#>     Mean:      0 
#>     Variance:  10 
#>   Covariance -- Inverse-Wishart: df = 1
#> 
#> Marginalization:
#>   Method: Empirical    (Sub)population: [-Inf, Inf]
#>   Decision rule: All
#> 
#> MCMC Convergence Diagnostics:
#>   Fixed effects -- MPSRF: 1.0218
#>         Parameter   ESS   Rhat
#>      b11_group[1] 400.7 1.0047
#>        b11_age[2] 122.7 1.0044
#>  b11_group_age[3] 523.4 1.0033
#>      b10_group[1] 475.6 1.0037
#>        b10_age[2]  91.7 1.0168
#>  b10_group_age[3] 561.2 1.0029
#>      b01_group[1] 632.2 1.0038
#>        b01_age[2] 167.5 1.0140
#>  b01_group_age[3] 729.9 1.0017
#> 
#>   Random effects -- MPSRF: 1.0247
#>         Parameter   ESS   Rhat
#>  g11_Intercept[1] 109.4 1.0051
#>  g10_Intercept[1]  83.8 1.0202
#>  g01_Intercept[1] 157.4 1.0176
#> 
#>   Variance components -- MPSRF: 1.0037
#>               Parameter    ESS   Rhat
#>  g11_InterceptIntercept 2304.1 1.0005
#>  g10_InterceptIntercept 2261.1 1.0056
#>  g01_InterceptIntercept 2592.1 0.9996
#> 
```
