# Estimate Theta via Analytical Integration

Estimate success probabilities by integrating over a range of covariate
values.

## Usage

``` r
estimate_theta_analytical(
  est_pars,
  X,
  grp,
  range_x,
  grp_var,
  population_var,
  fixed = NULL,
  random = NULL
)
```

## Arguments

- est_pars:

  List of nIt (P x Q) matrices of posterior regression coefficients.

- X:

  Design matrix of covariate data, where the covariate of interest is in
  the third column.

- grp:

  Scalar. Value of group indicator (0 or 1).

- range_x:

  Numeric vector of lower and upper bound of range to integrate over.

- fixed:

  Optional character vector with names of fixed variables. Default is
  NULL.

- random:

  Optional character vector with names of random variables. Default is
  NULL.

## Value

List of nIt vectors of length K with multivariate probabilities.
Currently supported for K=2 only.
