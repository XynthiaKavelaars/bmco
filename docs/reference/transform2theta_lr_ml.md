# Transform Regression Coefficients to Theta for Multilevel Model

Compute theta (success probabilities) from regression coefficients for
multilevel data via various methods.

## Usage

``` r
transform2theta_lr_ml(
  est_pars,
  X,
  measurement_levels,
  population_var,
  grp_var,
  grp_lvl,
  method = c("Empirical", "Analytical", "Value"),
  range = NULL,
  value = NULL,
  fixed,
  random
)
```

## Arguments

- est_pars:

  List of n_it (P x Q) arrays with posterior regression coefficients.

- X:

  List of J (n_j x P) design matrices.

- measurement_levels:

  Vector of measurement levels per predictor ("discrete" or
  "continuous").

- population_var:

  Optional character string with variable name that defines the
  subpopulation. Required when range is not c(-Inf, Inf).

- grp_var:

  Name of variable that defines groups.

- grp_lvl:

  Vector of group level names.

- method:

  "Value" for fixed values, "Empirical" for empirical marginalization,
  "Analytical" for numerical integration. Default is "Empirical".

- range:

  Optional range that defines the population of interest. Default is
  c(-Inf, Inf).

- value:

  Optional scalar that defines the population of interest. Required when
  method = "Value".

- fixed:

  Character vector with names of fixed variables.

- random:

  Character vector with names of random variables.

## Value

A list with:

- m_theta_a:

  List of n_it vectors with multivariate probabilities for group A

- m_theta_b:

  List of n_it vectors with multivariate probabilities for group B
