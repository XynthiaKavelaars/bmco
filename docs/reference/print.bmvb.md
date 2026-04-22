# Print Method for bmvb Objects

Print Method for bmvb Objects

## Usage

``` r
# S3 method for class 'bmvb'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A bmvb object.

- digits:

  Number of digits to display. Default is 3.

- ...:

  Additional arguments (not used).

## Value

Invisibly returns the input object.

## See also

[`summary.bmvb`](https://xynthiakavelaars.github.io/bmco/reference/summary.bmvb.md)

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
  y_vars = c("y1", "y2"), n_it = 1000
)
print(fit)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#>    Group mean y1 mean y2
#>  placebo   0.520   0.481
#>     drug   0.616   0.639
#>   n(placebo) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > placebo) [All rule]: 0.795
#> 
#> Use summary() for credible intervals and ESS.
#> 
```
