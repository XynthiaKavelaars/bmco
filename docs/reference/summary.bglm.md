# Summary Method for bglm Objects

Provides a comprehensive summary of a
[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
analysis, including the regression coefficient table, prior
specification, MCMC diagnostics (effective sample sizes and \\\hat{R}\\
per parameter), and, when the model was fitted with
`return_samples = TRUE`, credible intervals.

## Usage

``` r
# S3 method for class 'bglm'
summary(object, prob = 0.95, ...)
```

## Arguments

- object:

  A `bglm` object returned by
  [`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md).

- prob:

  Numeric. Coverage probability for credible intervals. Default is
  `0.95`.

- ...:

  Additional arguments (not used).

## Value

An object of class `summary.bglm`, a list containing:

- estimates:

  Posterior means and SDs of group probabilities and regression
  coefficients.

- sample_sizes:

  Group sample sizes.

- delta:

  Posterior mean differences, SEs, and posterior probability.

- info:

  Prior specification, test settings, and marginalization details.

- credible_intervals:

  If posterior samples are available: credible interval matrices for
  `theta_a`, `theta_b`, and `delta`.

- effective_n:

  If posterior samples are available: effective sample sizes for
  `theta_a`, `theta_b`, and `delta`.

- mcmc_diags:

  MCMC diagnostics for the regression coefficients (effective sample
  sizes and \\\hat{R}\\).

## See also

[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md),
[`print.bglm`](https://xynthiakavelaars.github.io/bmco/reference/print.bglm.md)

## Examples

``` r
# Uses the pre-computed example object shipped with the package:
summary(bglm_fit)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#> Group Estimates:
#>    Group mean y1 sd y1      95% CI y1 mean y2 sd y2      95% CI y2
#>  placebo   0.517 0.049 [0.421, 0.612]   0.501 0.049 [0.405, 0.598]
#>     drug   0.660 0.047 [0.566, 0.749]   0.689 0.045 [0.597, 0.775]
#> 
#>   n(placebo) = 100    n(drug) = 100
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.143, 0.187
#>   Delta SE   (y1, y2): 0, 0
#>   95% CI delta: y1 [0.008, 0.275]   y2 [0.054, 0.316]
#>   Posterior probability P(drug > placebo): 0.979
#> 
#> Regression Coefficients (mean [SD]):
#>                      b11            b10           b01
#> Intercept  0.233 [1.247]  1.863 [1.184] 0.073 [1.248]
#> group      0.065 [1.734] -0.192 [1.758] 0.782 [1.768]
#> age       -0.003 [0.025] -0.035 [0.024] 0.000 [0.025]
#> group:age  0.034 [0.035]  0.028 [0.036] 0.011 [0.036]
#> 
#> Prior Specification (regression coefficients):
#>   Mean:
#>   Intercept group age group:age
#> 1         0     0   0         0
#>   Variance (diagonal of inverse precision):
#>   Intercept group age group:age
#> 1        10    10  10        10
#> 
#> Marginalization:
#>   Method: Empirical
#>   (Sub)population: [-Inf, Inf]
#>   Decision rule: All
#> 
#> MCMC Diagnostics (regression coefficients):
#>   Multivariate PSRF (MPSRF): 1.0028
#>         Parameter    ESS   Rhat
#>  b11_Intercept[1] 6662.3 1.0000
#>      b11_group[2] 6063.1 1.0007
#>        b11_age[3] 6839.7 1.0000
#>  b11_group_age[4] 4799.1 1.0008
#>  b10_Intercept[1] 6756.0 1.0001
#>      b10_group[2] 5993.6 1.0010
#>        b10_age[3] 6830.4 0.9999
#>  b10_group_age[4] 4613.8 1.0010
#>  b01_Intercept[1] 7002.8 1.0002
#>      b01_group[2] 6039.8 1.0020
#>        b01_age[3] 7194.7 1.0002
#>  b01_group_age[4] 4709.6 1.0023
#> 
```
