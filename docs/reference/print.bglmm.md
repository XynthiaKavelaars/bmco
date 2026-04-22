# Print Method for bglmm Objects

Print Method for bglmm Objects

## Usage

``` r
# S3 method for class 'bglmm'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A bglmm object.

- digits:

  Number of digits to display. Default is 3.

- ...:

  Additional arguments (not used).

## Value

Invisibly returns the input object.

## See also

[`summary.bglmm`](https://xynthiakavelaars.github.io/bmco/reference/summary.bglmm.md)

## Examples

``` r
# Uses the pre-computed example object shipped with the package:
print(bglmm_fit)
#> 
#> Bayesian Multilevel Multivariate Logistic Regression
#> ======================================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.515   0.408
#>     drug   0.617   0.592
#>   J = 20 clusters    n(placebo) = 150    n(drug) = 150
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.977
#> Marginalization: Empirical over [-Inf, Inf]
#> MPSRF: fixed = 1.0218    random = 1.0247    variance = 1.0037
#> 
#> Use summary() for full coefficient tables, priors and MCMC diagnostics.
#> 
```
