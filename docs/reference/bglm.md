# Bayesian Generalized Linear Model

Perform a Bayesian test for differences between two (sub)groups on
multiple binary outcomes using multinomial logistic regression as
described in Kavelaars et al. (2024) .

## Usage

``` r
bglm(
  data,
  grp,
  grp_a,
  grp_b,
  x_var,
  y_vars,
  x_method = c("Empirical", "Analytical", "Value"),
  x_def = c(-Inf, Inf),
  test = c("right_sided", "left_sided"),
  rule = c("All", "Any", "Comp"),
  w = NULL,
  b_mu0 = NULL,
  b_sigma0 = NULL,
  n_burn = 10000,
  n_it = 20000,
  n_thin = 1,
  n_chain = 2,
  start = c(0.5, 1),
  return_diagnostics = TRUE,
  return_diagnostic_plots = FALSE,
  return_samples = FALSE
)
```

## Arguments

- data:

  Data frame containing the data.

- grp:

  Character string. Name of the grouping variable (will be treated as
  factor).

- grp_a:

  Value of grp indicating first group.

- grp_b:

  Value of grp indicating second group.

- x_var:

  Character string. Name of covariate variable (currently supports
  single continuous or binary covariate)

- y_vars:

  Character vector. Names of outcome variables (currently supports 2
  outcomes).

- x_method:

  Character. Method for handling covariate: "Analytical" (numerical
  integration), "Empirical" (empirical marginalization), or "Value"
  (specific value). Default is "Empirical".

- x_def:

  Numeric vector. Defines subpopulation: length-2 vector c(lower, upper)
  for "Analytical"/"Empirical", or scalar for "Value". Default is
  c(-Inf, Inf).

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

- b_mu0:

  Vector of prior means of fixed regression coefficients. Default is
  `rep(0, P)`, where P refers to the number of columns in the model
  matrix.

- b_sigma0:

  Prior covariance matrix (PxP) of regression coefficients. Default is
  `diag(1e-2, P)`, where P refers to the number of columns in the model
  matrix.

- n_burn:

  Integer. Number of burn-in iterations. Default is 10000.

- n_it:

  Integer. Number of MCMC iterations. Default is 20000.

- n_thin:

  Integer. Thinning interval. Default is 1.

- n_chain:

  Integer. Number of MCMC chains to be sampled. Default is 2.

- start:

  Numeric vector. Starting values for chains. Should have length
  `n_chain`. Default is c(0.5, 1).

- return_diagnostics:

  Logical. Return MCMC diagnostics? Default is TRUE.

- return_diagnostic_plots:

  Logical. Should MCMC chains for diagnostic plots (traceplots,
  autocorrelation, density) be returned? Default is `FALSE`. If `TRUE`,
  diagnostics are returned by default.

- return_samples:

  Logical. Should posterior samples be returned? Default is FALSE.

## Value

An object of class `bglm`, a list containing:

- estimates:

  A list with posterior means and standard deviations of group
  probabilities (`mean_a`, `mean_b`, `sd_a`, `sd_b`), as well as
  posterior means (`b`) and standard deviations (`b_sd`) of the
  regression coefficients.

- sample_sizes:

  A list with group sample sizes (`n_a`, `n_b`).

- delta:

  A list with posterior mean differences (`mean_delta`), posterior
  standard errors (`se_delta`), posterior probability of the hypothesis
  (`pop`), and, if `rule = "Comp"`, the weighted difference (`w_delta`).

- info:

  A list with prior specifications, test settings, group labels,
  covariate handling method, and subpopulation definition.

- diags:

  If diagnostics are requested, a list with MCMC diagnostic results for
  the regression coefficients.

- samples:

  If `return_samples = TRUE`, a list containing posterior draws of
  `theta_a`, `theta_b`, `delta`, and regression coefficients.

## References

Kavelaars X, Mulder J, Kaptein M (2024). “Bayesian Multivariate Logistic
Regression for Superiority and Inferiority Decision-Making under
Observable Treatment Heterogeneity.” *Multivariate Behavioral Research*,
**59**(4), 859–882.
[doi:10.1080/00273171.2024.2337340](https://doi.org/10.1080/00273171.2024.2337340)
.

## Examples

``` r
# Example with simulated data
# Generate data
set.seed(123)
n <- 100

data <- data.frame(
 group = rep(c("A", "B"), each = n/2),
 x = rnorm(n),
 stringsAsFactors = FALSE
)

p1 <- p2 <- rep(NA, n)

for (i in 1:n) {
 grpB <- ifelse(data$group[i] == "B", 1, 0)

 p1[i] <- plogis(-0.50 + 0.75 * grpB + 0.10 * data$x[i] + 0.20 * grpB * data$x[i])
 p2[i] <- plogis(-0.50 + 0.80 * grpB + 0.05 * data$x[i] + 0.15 * grpB * data$x[i])

 data$y1[i] <- rbinom(1, 1, p1[i])
 data$y2[i] <- rbinom(1, 1, p2[i])
}

# Analyze
result <- bglm(
 data = data,
 grp = "group",
 grp_a = "A",
 grp_b = "B",
 x_var = "x",
 y_vars = c("y1", "y2"),
 x_method = "Empirical",
 x_def = c(-Inf, Inf),
 test = "right_sided",
 rule = "All",
 n_burn = 100, # Too low for proper MCMC sampling
 n_it = 500 # Too low for proper MCMC sampling
)

print(result)
#> 
#> Bayesian Multivariate Logistic Regression
#> ==========================================
#> 
#>  Group mean y1 mean y2
#>      A   0.278   0.508
#>      B   0.637   0.559
#>   n(A) = 50    n(B) = 50
#> 
#> Posterior probability P(B > A) [All rule]: 0.708
#> Marginalization: Empirical over [-Inf, Inf]
#> 
#> Use summary() for regression coefficients, priors and MCMC diagnostics.
#> 
```
