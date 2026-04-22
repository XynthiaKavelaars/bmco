# Integrand for Analytical Marginalization

Integrand function for integration over a range of a covariate in
numerical marginalization. Works for both GLM and GLMM.

## Usage

``` r
integrand(
  x,
  est_rc,
  mu_x,
  sigma_x,
  range_x,
  grp,
  grp_var,
  population_var,
  q,
  fixed = NULL,
  random = NULL
)
```

## Arguments

- x:

  Scalar. Value of covariate x.

- est_rc:

  (P x Q) matrix of regression coefficients.

- mu_x:

  Scalar. Mean of distribution of covariate.

- sigma_x:

  Scalar. Standard deviation of distribution of covariate.

- range_x:

  Numeric vector of lower and upper bound of range to integrate over.

- grp:

  Scalar. Value of group indicator (0 or 1).

- q:

  Scalar. Response category in 1 to Q.

- fixed:

  Optional character vector with names of fixed variables. Default is
  NULL.

- random:

  Optional character vector with names of random variables. Default is
  NULL.

## Value

Scalar. Joint response probability for response category q.
