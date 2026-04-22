# Transform Regression Coefficients to Success Probabilities

Compute theta (success probabilities) from regression coefficients via
various methods.

## Usage

``` r
transform2theta(
  beta_draw_pg,
  X,
  Y = NULL,
  grp_var,
  population_var,
  measurement_level,
  method,
  range,
  value
)
```

## Arguments

- beta_draw_pg:

  List of n_it (P x Q) matrices with posterior regression coefficients.

- X:

  (n x P) matrix with covariate data.

- Y:

  (n x Q) matrix with multinomial response data. Default is NULL.

- measurement_level:

  Character. "discrete" or "continuous".

- method:

  "Value" for vector of fixed values, "Empirical" for empirical
  marginalization, "Analytical" for numerical marginalization.

- range:

  If method = "Analytical" or "Empirical" and if measurement_level is
  continuous: range that defines the population of interest by a vector
  containing a lower and an upper bound.

- value:

  If method = "Value": value that defines the population of interest by
  a scalar value.

## Value

A list with:

- m_theta_a:

  List of n_it vectors of length K with multivariate probabilities for
  group A

- m_theta_b:

  List of n_it vectors of length K with multivariate probabilities for
  group B
