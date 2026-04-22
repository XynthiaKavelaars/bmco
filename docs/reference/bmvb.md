# Bayesian Multivariate Bernoulli Test

Perform a Bayesian test for differences between two groups on multiple
binary outcomes using a Multivariate Bernoulli distribution, as
described in Kavelaars et al. (2020) .

## Usage

``` r
bmvb(
  data,
  grp,
  grp_a,
  grp_b,
  y_vars,
  test = c("right_sided", "left_sided"),
  rule = c("All", "Any", "Comp"),
  w = NULL,
  prior_a = 0.5,
  prior_b = 0.5,
  n_it = 10000,
  return_samples = FALSE
)
```

## Arguments

- data:

  Data frame containing the data.

- grp:

  Character string. Name of the grouping variable.

- grp_a:

  Value of grp indicating first group.

- grp_b:

  Value of grp indicating second group.

- y_vars:

  Character vector. Names of outcome variables (currently supports 2
  outcomes).

- test:

  Character. Direction of test: "left_sided" for P(A\>B) or
  "right_sided" for P(B\>A). Default is "right_sided".

- rule:

  Character. Decision rule: "All" (all outcomes favor hypothesis), "Any"
  (at least one outcome favors hypothesis), or "Comp" (weighted
  combination). Default is "All".

- w:

  Numeric vector. Weights for compensatory rule. Only used if rule =
  "Comp". If NULL and rule = "Comp", equal weights are used. Default is
  NULL.

- prior_a:

  Numeric. Prior hyperparameter (Dirichlet) for group A. Default is 0.5
  (Jeffreys' prior)

- prior_b:

  Numeric. Prior hyperparameter (Dirichlet) for group B. Default is 0.5
  (Jeffreys' prior).

- n_it:

  Integer. Number of MCMC iterations. Default is 10000.

- return_samples:

  Logical. Should posterior samples be returned? Default is FALSE.

## Value

An object of class `bmvb`, a list containing:

- estimates:

  A list with posterior means (`mean_a`, `mean_b`) and standard
  deviations (`sd_a`, `sd_b`) of the category probabilities for both
  groups.

- sample_sizes:

  A list with group sample sizes (`n_a`, `n_b`).

- delta:

  A list with posterior mean differences (`mean_delta`), posterior
  standard errors (`se_delta`), posterior probability of the hypothesis
  (`pop`), and, if `rule = "Comp"`, the weighted difference (`w_delta`).

- info:

  A list with test specifications, including the decision rule, test
  direction, group labels, and weights (if applicable).

- samples:

  If `return_samples = TRUE`, a list containing posterior draws of
  `theta_a`, `theta_b`, and `delta`.

## References

Kavelaars X, Mulder J, Kaptein M (2020). “Decision-making with multiple
correlated binary outcomes in clinical trials.” *Statistical Methods in
Medical Research*, **29**(11), 3265–3277.
[doi:10.1177/0962280220922256](https://doi.org/10.1177/0962280220922256)
.

## Examples

``` r
# \donttest{
# Example with simulated data
# Generate data
set.seed(123)
data <- data.frame(
treatment = rep(c("control", "drug"), each = 50),
 outcome1 = rbinom(100, 1, 0.5),
 outcome2 = rbinom(100, 1, 0.5)
)

# Analyze
result <- bmvb(
 data = data,
 grp = "treatment",
 grp_a = "control",
 grp_b = "drug",
 y_vars = c("outcome1", "outcome2"),
 n_it = 10000
)

print(result)
#> 
#> Multivariate Bernoulli Analysis
#> =================================
#> 
#>    Group mean y1 mean y2
#>  control   0.501   0.501
#>     drug   0.442   0.501
#>   n(control) = 50    n(drug) = 50
#> 
#> Posterior probability P(drug > control) [All rule]: 0.129
#> 
#> Use summary() for credible intervals and ESS.
#> 
# }
```
