# Bayesian Generalized Linear Mixed Model

Perform a Bayesian test for differences between two (sub)groups on
multiple binary outcomes using multilevel multinomial logistic
regression, as described in Kavelaars et al. (2023) .

## Usage

``` r
bglmm(
  data,
  grp,
  grp_a = NULL,
  grp_b = NULL,
  id_var,
  x_var,
  y_vars,
  x_method = c("Empirical", "Analytical", "Value"),
  x_def = c(-Inf, Inf),
  test = c("right_sided", "left_sided"),
  rule = c("All", "Any", "Comp"),
  w = NULL,
  n_burn = 10000,
  n_it = 50000,
  start = c(0.5, 1),
  fixed = NULL,
  random = NULL,
  b_mu0 = NULL,
  b_sigma0 = NULL,
  g_mu0 = NULL,
  g_sigma0 = NULL,
  nu0 = NULL,
  tau0 = NULL,
  n_chain = 2,
  return_thinned = TRUE,
  n_thin = 10,
  return_diagnostics = TRUE,
  return_diagnostic_plots = FALSE,
  return_samples = FALSE
)
```

## Arguments

- data:

  Data frame containing the data.

- grp:

  Character string. Name of the grouping variable.

- grp_a:

  Value of `grp` indicating first group (will be determined from factor
  levels if NULL).

- grp_b:

  Value of `grp` indicating second group (will be determined from factor
  levels if NULL).

- id_var:

  Character string. Name of cluster/ID variable.

- x_var:

  Character string. Name of covariate variable.

- y_vars:

  Character vector. Names of outcome variables (currently supports 2
  outcomes).

- x_method:

  Character. Method for handling covariate. Default is "Empirical".

- x_def:

  Numeric. Defines subpopulation. Default is `c(-Inf, Inf)`.

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

- n_burn:

  Integer. Number of burn-in iterations. Default is 10000.

- n_it:

  Integer. Number of MCMC iterations. Default is 50000 (takes long
  running time!).

- start:

  Numeric vector. Starting values for chains. Default is `c(0.5, 1)`.

- fixed:

  Character vector. Names of fixed effect variables. Default is c(x_var,
  grp_x_var).

- random:

  Character vector. Names of random effect variables. Default is
  c("Intercept", grp).

- b_mu0:

  Numeric vector. Prior means for fixed effects. Default is
  `rep(0, length(fixed))`.

- b_sigma0:

  Matrix. Prior covariance for fixed effects. Default is
  `diag(0.1, length(fixed))`.

- g_mu0:

  Numeric vector. Prior means for random effects. Default is
  `rep(0, length(random))`.

- g_sigma0:

  Matrix. Prior covariance for random effects. Default is
  `diag(0.1, length(random))`.

- nu0:

  Numeric. Prior degrees of freedom for inverse-Wishart. Default is
  `length(random)`.

- tau0:

  Matrix. Prior scale matrix of dimension
  `length(random) x length(random)` for inverse-Wishart. Default is
  `diag(1e-1, length(random))`.

- n_chain:

  Integer. Number of MCMC chains. Default is 2.

- return_thinned:

  Logical. Return thinned chains? Default is TRUE.

- n_thin:

  Integer. Thinning interval. Default is 10.

- return_diagnostics:

  Logical. Return MCMC diagnostics? Default is TRUE.

- return_diagnostic_plots:

  Logical. Should MCMC chains for diagnostic plots (traceplots,
  autocorrelation, density) be returned? Default is `FALSE`. If `TRUE`,
  diagnostics are returned by default.

- return_samples:

  Logical. Return posterior samples? Default is FALSE.

## Value

An object of class `bglmm`, a list containing:

- estimates:

  A list with posterior means and standard deviations of group
  probabilities (`mean_a`, `mean_b`, `sd_a`, `sd_b`). If estimated,
  posterior means and standard deviations of fixed effects (`b`, `b_sd`)
  and random effects and variance components (`g`, `g_sd`, `tau`,
  `tau_sd`) are included.

- sample_sizes:

  A list with group sample sizes (`n_a`, `n_b`) and the number of
  clusters (`J`).

- delta:

  A list with posterior mean differences (`mean_delta`), posterior
  standard errors (`se_delta`), posterior probability of the hypothesis
  (`pop`), and, if `rule = "Comp"`, the weighted difference (`w_delta`).

- info:

  A list with prior specifications, model structure (fixed and random
  effects), test settings, group labels, covariate handling method, and
  subpopulation definition.

- diags:

  If diagnostics are requested, a list with MCMC diagnostic results for
  fixed effects, random effects, and variance components.

- samples:

  If `return_samples = TRUE`, a list containing posterior draws of group
  probabilities, differences, fixed effects, random effects, and
  variance components (if applicable).

## References

Kavelaars X, Mulder J, Kaptein M (2023). “Bayesian multilevel
multivariate logistic regression for superiority decision-making under
observable treatment heterogeneity.” *BMC Medical Research Methodology*,
**23**(1).
[doi:10.1186/s12874-023-02034-z](https://doi.org/10.1186/s12874-023-02034-z)
.

## Examples

``` r
# \donttest{
# Example with simulated data
# Generate data
set.seed(123)
J <- 20 # No. clusters
nJ <- 15 # Sample size per cluster

# Generate random intercepts
uj_1 <- rnorm(J)
uj_2 <- rnorm(J)
data <- data.frame(
 id = factor(rep(1:J, each = nJ)),
 group = rep(rep(c("A", "B"), each = J/2), each = nJ),
 x = rnorm(J * nJ),
 stringsAsFactors = FALSE
)

p1 <- p2 <- rep(NA, J * nJ)

for (i in 1:(J * nJ)) {
 j <- as.numeric(data$id[i])
 grpB <- ifelse(data$group[i] == "B", 1, 0)

 p1[i] <- plogis(-0.50 + 0.75 * grpB + 0.10 * data$x[i] + 0.20 * grpB * data$x[i] + uj_1[j])
 p2[i] <- plogis(-0.50 + 0.80 * grpB + 0.05 * data$x[i] + 0.15 * grpB * data$x[i] + uj_2[j])

 data$y1[i] <- rbinom(1, 1, p1[i])
 data$y2[i] <- rbinom(1, 1, p2[i])
}

# Analyze
result <- bglmm(
 data = data,
 grp = "group",
 grp_a = "A",
 grp_b = "B",
 id_var = "id",
 x_var = "x",
 y_vars = c("y1", "y2"),
 x_method = "Empirical",
 x_def = c(-Inf, Inf),
 fixed = c("group", "x", "group_x"),
 random = c("Intercept"), # Random intercept model
 test = "right_sided",
 rule = "All",
 n_burn = 100, # Too low for proper MCMC sampling
 n_it = 500 # Too low for proper MCMC sampling
 )
#> Warning: 20 cluster(s) have observations from only one group. This may affect estimation.
#> Warning: Low effective sample size for some parameters. Results may be unreliable.
#> Warning: Low effective sample size for some parameters. Results may be unreliable.
#> Warning: Low effective sample size for some parameters. Results may be unreliable.
#> Warning: MCMC chains may not have converged (MPSRF > 1.1). Consider increasing n_burn or n_it.

print(result) # Warnings due to low number of MCMC iterations (n_burn and n_it)
#> 
#> Bayesian Multilevel Multivariate Logistic Regression
#> ======================================================
#> 
#>  Group mean y1 mean y2
#>      A   0.448   0.341
#>      B   0.631   0.623
#>   J = 20 clusters    n(A) = 150    n(B) = 150
#> 
#> Posterior probability P(B > A) [All rule]: 1.000
#> Marginalization: Empirical over [-Inf, Inf]
#> MPSRF: fixed = 1.1516    random = 1.1371    variance = 1.0614
#> 
#> Use summary() for full coefficient tables, priors and MCMC diagnostics.
#> 
# }
```
