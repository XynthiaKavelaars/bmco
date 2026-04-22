# Summary Method for bmvb Objects

Provides a comprehensive summary of a
[`bmvb`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md)
analysis. When the model was fitted with `return_samples = TRUE`,
credible intervals and effective sample sizes are included.

## Usage

``` r
# S3 method for class 'bmvb'
summary(object, prob = 0.95, ...)
```

## Arguments

- object:

  A `bmvb` object returned by
  [`bmvb`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md).

- prob:

  Numeric. Coverage probability for credible intervals. Default is
  `0.95`.

- ...:

  Additional arguments (not used).

## Value

An object of class `summary.bmvb`, a list containing all fields of
`object` plus:

- credible_intervals:

  If posterior samples are available: a list with `prob` and credible
  interval matrices for `theta_a`, `theta_b`, and `delta`.

- effective_n:

  If posterior samples are available: a list with effective sample sizes
  for `theta_a`, `theta_b`, and `delta`.

## See also

[`bmvb`](https://xynthiakavelaars.github.io/bmco/reference/bmvb.md),
[`print.bmvb`](https://xynthiakavelaars.github.io/bmco/reference/print.bmvb.md)

## Examples

``` r
set.seed(2024)
trial_data <- data.frame(
  treatment = rep(c("placebo", "drug"), each = 50),
  y1 = rbinom(100, 1, rep(c(0.40, 0.60), each = 50)),
  y2 = rbinom(100, 1, rep(c(0.50, 0.70), each = 50))
)
fit <- bmvb(
  data = trial_data, grp = "treatment",
  grp_a = "placebo", grp_b = "drug",
  y_vars = c("y1", "y2"), n_it = 1000,
  return_samples = TRUE
)
summary(fit)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#> Group Estimates:
#>    Group mean y1 sd y1      95% CI y1 mean y2 sd y2      95% CI y2
#>  placebo   0.520 0.067 [0.381, 0.646]   0.481 0.072 [0.343, 0.623]
#>     drug   0.616 0.068 [0.482, 0.743]   0.639 0.067 [0.499, 0.758]
#> 
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Treatment Effect:
#>   Delta mean (y1, y2): 0.096, 0.157
#>   Delta SE   (y1, y2): 0.003, 0.003
#>   95% CI delta: y1 [-0.081, 0.282]   y2 [-0.041, 0.349]
#>   Posterior probability P(drug > placebo): 0.795
#> 
#> Test Information:
#>   Decision rule: All
#>   Hypothesis: P(drug > placebo)
#> 
#> Effective Sample Sizes:
#>   theta (placebo): 1000, 1000
#>   theta (drug): 1098, 1000
#>   delta:      1000, 1000
#> 
```
