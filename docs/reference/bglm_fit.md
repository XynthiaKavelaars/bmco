# Pre-computed bglm Example Fit

A fitted
[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md)
object estimated on
[`bglm_data`](https://xynthiakavelaars.github.io/bmco/reference/bglm_data.md).
Used in package examples and tests so that `print`, `summary`, and
`plot` examples run in well under 5 seconds without re-running the MCMC
sampler.

## Usage

``` r
bglm_fit
```

## Format

An object of class `bglm` as returned by
[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md). See
the **Value** section of
[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md) for
a full description of the list components. Key settings:
`n_burn = 10000`, `n_it = 20000`, `n_chain = 2`,
`return_diagnostics = TRUE`, `return_samples = TRUE`.

## Details

Generated with `set.seed(2024)`. See `data-raw/generate_examples.R` for
the full reproducible script.

## See also

[`bglm`](https://xynthiakavelaars.github.io/bmco/reference/bglm.md),
[`bglm_data`](https://xynthiakavelaars.github.io/bmco/reference/bglm_data.md),
[`bglmm_fit`](https://xynthiakavelaars.github.io/bmco/reference/bglmm_fit.md)

## Examples

``` r
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
