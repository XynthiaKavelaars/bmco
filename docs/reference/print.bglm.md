# Print Method for bglm Objects

Print Method for bglm Objects

## Usage

``` r
# S3 method for class 'bglm'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A bglm object.

- digits:

  Number of digits to display. Default is 3.

- ...:

  Additional arguments (not used).

## Value

Invisibly returns the input object.

## See also

[`summary.bglm`](https://xynthiakavelaars.github.io/bmco/reference/summary.bglm.md)

## Examples

``` r
# Uses the pre-computed example object shipped with the package:
print(bglm_fit)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.517   0.501
#>     drug   0.660   0.689
#>   n(placebo) = 100    n(drug) = 100
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.979
#> Marginalization: Empirical over [-Inf, Inf]
#> 
#> Use summary() for regression coefficients, priors and MCMC diagnostics.
#> 
```
